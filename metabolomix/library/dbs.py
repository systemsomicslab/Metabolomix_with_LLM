"""MS-DIAL の `.dbs` / `.lbm2` を読む（LargeListMessagePack + LZ4 block）。

Key 番号と枠組みの正準は `docs/schema/molecule_ms_reference.md`。推測で直さない。

`.lbm2`（および `.msp2`）はファイル全体が `LargeListMessagePack` のチャンク列。
`.dbs` は ZIP で、内部の `<種別>/<DataBaseID>/DataBase` エントリが同じチャンク列
（`.lbm2` とバイト単位で同一）を持ち、`Storage` エントリだけがチャンク分割の無い
単発の ext99+LZ4 フレーム（`MetabolomicsDataBases` 等を持つマップ）を持つ。
両者は展開ヘッダが同一（ext32 + type 99 + 展開後サイズ）なので同じ関数で剥がせる。

deps: record のみ。mcp_core / session_state / tools_* を import しない。
"""
from __future__ import annotations

import io
import struct
import zipfile
from pathlib import Path
from typing import Iterator

import lz4.block
import msgpack

from metabolomix.library.record import ION_MODES, make_record

EXT_TYPE_CODE = 99
CHUNK_HEADER_SIZE = 11          # ext32 ヘッダ 6 + 展開後サイズ 5
ZIP_MAGIC = b"PK\x03\x04"

# MoleculeMsReference の [Key(N)]。docs/schema/molecule_ms_reference.md が正準。
_K_PRECURSOR_MZ = 1
_K_CHROMXS = 2
_K_ION_MODE = 3
_K_SPECTRUM = 4
_K_NAME = 5
_K_FORMULA = 6
_K_ONTOLOGY = 7
_K_SMILES = 8
_K_INCHIKEY = 9
_K_ADDUCT = 10
_K_COMPOUND_CLASS = 14

# MsRefSearchParameterBase の [Key(N)]（配列位置）。
# 7 / 8 は測定スペクトルの足切りで、MS-DIAL が採点前に掛ける正規化の唯一の可変部分
# （Task 1 の調査結果）。これが取れないと、足切りを変えた run で数値が合わない理由を
# 「移植の誤り」と誤分類することになる。
_SEARCH_PARAM_KEYS = {
    0: "mass_range_begin", 1: "mass_range_end", 2: "rt_tolerance",
    5: "ms1_tolerance", 6: "ms2_tolerance",
    7: "relative_amp_cutoff", 8: "absolute_amp_cutoff",
    9: "squared_weighted_dot_cutoff", 10: "squared_simple_dot_cutoff",
    11: "squared_reverse_dot_cutoff", 12: "matched_peaks_percentage_cutoff",
    13: "total_score_cutoff", 14: "minimum_spectrum_match",
    # 15〜18 は bool。総合スコア（`spectral_match.total_score`）に RT / CCS 項を
    # 入れるかがここで決まる。**上流の既定は 4 つとも False** だが実 run では
    # True のことがあるので、既定値を仮定せず必ずここから読む。
    15: "use_time_for_annotation_filtering", 16: "use_time_for_annotation_scoring",
    17: "use_ccs_for_annotation_filtering", 18: "use_ccs_for_annotation_scoring",
}


def iter_decompressed_chunks(data: bytes) -> Iterator[bytes]:
    """チャンク列を展開して順に返す。

    `.dbs` の `Storage` エントリも、チャンク分割が無いだけで同じ 11 バイト
    ヘッダ（ext32 + type 99 + 展開後サイズ）を持つので、この関数でそのまま
    1 個だけ展開できる（呼び出し側が `next()` すればよい）。
    """
    offset = 0
    while offset < len(data):
        if data[offset] != 0xC9:
            raise ValueError(f"ext32 ヘッダではありません: {data[offset]:#04x} (offset={offset})")
        ext_len = struct.unpack(">I", data[offset + 1:offset + 5])[0]
        type_code = struct.unpack("b", data[offset + 5:offset + 6])[0]
        if type_code != EXT_TYPE_CODE:
            raise ValueError(f"想定外の ext type: {type_code}")
        raw_len = struct.unpack(">i", data[offset + 7:offset + 11])[0]
        body_len = ext_len - 5
        body = data[offset + CHUNK_HEADER_SIZE: offset + CHUNK_HEADER_SIZE + body_len]
        offset += CHUNK_HEADER_SIZE + body_len
        yield lz4.block.decompress(body, uncompressed_size=raw_len)


def read_array_count(chunk: bytes) -> int:
    """展開後チャンクの要素数を返す。

    先頭 5 バイトは配列ヘッダ用に予約されるが、`WriteArrayHeader` は要素数が
    小さいと**短形式**で書く。常に 5 バイトの `0xdd` と仮定すると最終チャンクで
    壊れる。要素は常にオフセット 5 から始まる（`DeserializeList` が固定）。
    """
    first = chunk[0]
    if first == 0xDD:
        return struct.unpack(">I", chunk[1:5])[0]
    if first == 0xDC:
        return struct.unpack(">H", chunk[1:3])[0]
    if 0x90 <= first <= 0x9F:
        return first & 0x0F
    raise ValueError(f"配列ヘッダではありません: {first:#04x}")


def _adduct_name(adduct: object) -> str | None:
    """`AdductIon` 入れ子からインデックス 2（表示名）を取り出す。"""
    if adduct is None:
        return None
    if isinstance(adduct, (list, tuple)):
        return adduct[2] if len(adduct) > 2 else None
    return adduct


def _formula_str(formula: object) -> str | None:
    """`Formula` 入れ子からインデックス 0（組成式文字列）を取り出す。"""
    if formula is None:
        return None
    if isinstance(formula, (list, tuple)):
        return formula[0] if formula else None
    return formula


#: `ChromXs` 入れ子の種別コード（`docs/schema/AlignmentSpotProperty.md`）。
#: 1=RT / 2=RI / 3=m/z / 4=drift。位置ではなくこのコードで一致を取る
#: （`metabolomix/pai2/reader.py` の `_convert_to_times` と同じ流儀）。
_CHROMX_TYPE_RT = 1


def _rt_from_chromxs(chromxs: object) -> float | None:
    """`ChromXs`（`[[種別, [値, ...]], ...]`）から種別 1（RT）を探して返す。

    位置決め打ち（`chromxs[0]`）だと、並び順が変わる・別の種別（m/z 等）が
    先頭に来ると黙って別の値を RT として返してしまう。実データでは末尾に
    入れ子でないスカラーが付くこともあるので、要素ごとに list/tuple 判定を
    してから読む（`pai2/reader.py` と同じ防御）。

    値が `-1`（MS-DIAL の「未設定」番兵）のときは実質未設定なので `None` を返す。
    """
    if not isinstance(chromxs, (list, tuple)):
        return None
    for entry in chromxs:
        if not isinstance(entry, (list, tuple)) or len(entry) < 2:
            continue
        if entry[0] != _CHROMX_TYPE_RT:
            continue
        values = entry[1]
        if not isinstance(values, (list, tuple)) or not values:
            continue
        rt = values[0]
        if rt is None or rt < 0:
            return None
        return rt
    return None


def _spectrum_peaks(spectrum: object) -> list[list[float]]:
    """`SpectrumPeak` 列から使う先頭 2 つ（m/z・強度）だけを取り出す。"""
    if not spectrum:
        return []
    return [[peak[0], peak[1]] for peak in spectrum]


def _get(raw: list, key: int):
    return raw[key] if len(raw) > key else None


def _to_record(raw: list, record_index: int, library_id: str | None) -> dict:
    """1 レコードぶんの `MoleculeMsReference` 配列を正規化レコードへ変換する。"""
    return make_record(
        name=_get(raw, _K_NAME),
        precursor_mz=_get(raw, _K_PRECURSOR_MZ),
        ion_mode=ION_MODES.get(_get(raw, _K_ION_MODE)),
        adduct=_adduct_name(_get(raw, _K_ADDUCT)),
        rt=_rt_from_chromxs(_get(raw, _K_CHROMXS)),
        formula=_formula_str(_get(raw, _K_FORMULA)),
        inchikey=_get(raw, _K_INCHIKEY),
        smiles=_get(raw, _K_SMILES),
        compound_class=_get(raw, _K_COMPOUND_CLASS),
        ontology=_get(raw, _K_ONTOLOGY),
        spectrum=_spectrum_peaks(_get(raw, _K_SPECTRUM)),
        library_id=library_id,
        record_index=record_index,
    )


def _iter_records_from_data(data: bytes, library_id: str | None = None) -> Iterator[dict]:
    """チャンク列（`.lbm2` 全体、または `.dbs` の `DataBase` エントリ）を読む。

    チャンクを 1 個ずつ展開・消費し、次のチャンクに進む前に参照を手放す
    （1 チャンクの展開後サイズは数百 MB〜1 GB になり得るため、複数チャンクを
    同時にメモリへ載せない）。
    """
    idx = 0
    for chunk in iter_decompressed_chunks(data):
        count = read_array_count(chunk)
        # 展開後チャンクは数百 MB〜1 GB になり得るため、msgpack の既定 100 MiB
        # 上限（対・信頼できない入力向けの安全弁）に当たる。上限を外すのではなく、
        # このチャンク自身の実サイズ（＝これ以上大きくなりようがない既知の上限）
        # に縛る。`arf/reader.py` の固定 1 GiB 定数と違い、ここは実測サイズが
        # その場で分かっているのでそれを使う——壊れた/切り詰められたチャンクや
        # read_array_count の誤読を、msgpack 層が早期に BufferFull で捕まえられる。
        unpacker = msgpack.Unpacker(raw=False, max_buffer_size=len(chunk))
        unpacker.feed(memoryview(chunk)[5:])
        for _ in range(count):
            raw = unpacker.unpack()
            yield _to_record(raw, idx, library_id)
            idx += 1


#: `.dbs` の中で本サーバが対象とする種別のみ（`docs/schema/molecule_ms_reference.md`
#: 「`.dbs` の ZIP 構造」）。`ProteomicsDB/` `EadLipidomicsDB/` は明示的に対象外——
#: 中身があっても `MoleculeMsReference` の Key 配置とは限らず、拾うとフィールドが
#: ずれたレコードを黙って返すことになる。
_DATABASE_ENTRY_PREFIX = "MetabolomicsDB/"


def _is_database_entry(info: zipfile.ZipInfo) -> bool:
    """`MetabolomicsDB/` 配下の `DataBase` 実体エントリだけを選ぶ。

    `MetabolomicsDB/MS-FINDER/DataBase` のような空エントリが混在するので、
    「`DataBase` で終わる」「`MetabolomicsDB/` 配下」だけでなく
    「サイズが 0 でない」ことも見る。
    """
    return (
        info.filename.startswith(_DATABASE_ENTRY_PREFIX)
        and info.filename.endswith("DataBase")
        and info.file_size > 0
    )


def _iter_records_from_zip(data: bytes) -> Iterator[dict]:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for info in z.infolist():
            if not _is_database_entry(info):
                continue
            parts = info.filename.split("/")
            library_id = parts[-2] if len(parts) >= 2 else None
            entry_data = z.read(info)
            yield from _iter_records_from_data(entry_data, library_id=library_id)


def iter_records(path: str | Path) -> Iterator[dict]:
    """`.dbs`（ZIP）と `.lbm2`（生）の両方から正規化レコードを読む。

    拡張子ではなく中身（ZIP は `PK\\x03\\x04` で始まる）で見分ける。
    """
    data = Path(path).read_bytes()
    if data[:4] == ZIP_MAGIC:
        yield from _iter_records_from_zip(data)
    else:
        yield from _iter_records_from_data(data)


def _search_params_from_parameter(parameter: list) -> dict:
    return {name: parameter[k] for k, name in _SEARCH_PARAM_KEYS.items() if k < len(parameter)}


def read_storage_meta(path: str | Path) -> dict | None:
    """`.dbs` の `Storage` エントリ（使用ライブラリ・元パス・注釈器パラメータ・注釈器の Key とライブラリの対応）を読む。

    `.lbm2`（ZIP でない、または `Storage` を持たない）には無いので `None` を返す。
    ZIP かどうかは先頭 4 バイトだけで判定する——ZIP でなければそこで打ち切り、
    大きな `.lbm2`/`.msp` に対してファイル全体を読む無駄を作らない
    （呼び出し側の store 構築が拡張子に関わらずこれを毎回呼ぶため、ここで
    効くと全呼び出し側が得をする）。
    """
    path = Path(path)
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head != ZIP_MAGIC:
        return None

    with zipfile.ZipFile(path) as z:
        if "Storage" not in z.namelist():
            return None
        raw = z.read("Storage")

    decompressed = next(iter_decompressed_chunks(raw))
    storage = msgpack.unpackb(decompressed, raw=False)

    databases = storage.get("MetabolomicsDataBases") or []
    if not databases:
        return None
    entry = databases[0]

    database = entry.get("DataBase") or []
    library_name = database[0] if len(database) > 0 else None
    source_path = database[3] if len(database) > 3 else None

    search_params: dict = {}
    pairs = entry.get("Pairs") or []
    if pairs:
        annotator_key = pairs[0][1].get("SerializableAnnotatorKey")
        if annotator_key:
            parameter = annotator_key[1].get("Parameter") or []
            search_params = _search_params_from_parameter(parameter)

    annotator_libraries, annotator_search_params = _annotators(databases)
    return {
        "library_name": library_name,
        "source_path": source_path,
        "search_params": search_params,
        "annotator_libraries": annotator_libraries,
        "annotator_search_params": annotator_search_params,
    }


def _annotators(databases: list) -> tuple[dict[str, str], dict[str, dict]]:
    """注釈器の Key（照合結果の AnnotatorID）→ ライブラリ名（`MetabolomicsDB/<名前>/`）。

    AnnotatorID はライブラリ名と同じとは限らない。GUI は `<名前>_<n>`、MS-DIAL
    Console は LBM ならファイルのパス（後の版は `LbmDB: <stem>`）、MSP なら `.msp` の
    パスや設定ファイルで付けた任意の名前を使う。どの規則でも、MS-DIAL 自身が
    保存したこの対応だけは食い違わない。

    2 つ目は Key → その注釈器の検索パラメータ。RT を使ったか等は注釈器ごとに
    違いうるので、先頭の注釈器（`search_params`）だけで代表させない。
    """
    mapping: dict[str, str] = {}
    params: dict[str, dict] = {}
    for database_entry in databases:
        database = database_entry.get("DataBase") or []
        name = database[0] if database else None
        if not isinstance(name, str):
            continue
        for pair in database_entry.get("Pairs") or []:
            body = pair[1] if isinstance(pair, (list, tuple)) and len(pair) > 1 else {}
            annotator = body.get("SerializableAnnotatorKey") if isinstance(body, dict) else None
            if not isinstance(annotator, (list, tuple)) or len(annotator) < 2:
                continue
            key = annotator[1].get("Key") if isinstance(annotator[1], dict) else None
            if isinstance(key, str) and key:
                mapping[key] = name
                params[key] = _search_params_from_parameter(annotator[1].get("Parameter") or [])
    return mapping, params
