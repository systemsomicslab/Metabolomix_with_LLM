"""参照ライブラリの永続 store（sha256 を鍵とする SQLite キャッシュ）。

`.dbs` / `.msp` の `iter_records` はどちらもジェネレータで正規化レコードを吐くが、
`.dbs` の LZ4 block はチャンク単位でしか展開できず全走査が避けられない一方、
必要な情報は元の 1% 以下（`SpectrumPeak` の 13 フィールドのうち使うのは
`[m/z, intensity]` だけ）。そこで初回に 1 度だけコンパクトな SQLite へ変換し、
以後は precursor m/z 窓 × 極性で候補を引く（spec §5.2）。

キャッシュ先は環境変数 `LIPIDMIX_LIBRARY_CACHE_DIR`、無ければ
`mcp_core.DATA_DIR / ".library-cache"`。ファイル名は `<stem>-<sha256[:16]>.sqlite`
なので、元ファイルが変われば別ファイルになり、無効化の判断が要らない（古い
キャッシュは孤立するだけで、参照されなくなる）。

構築は一時ファイルに書いてから rename する（途中で落ちた store を次回の
`open_store` が掴まないため）。sha256 はサイズ・更新時刻と組でキャッシュ置き場の
`DIGEST_INDEX_NAME` に覚え、変わっていなければ計算し直さない（研究室の `.msp` は
1 GB 級で、全読みを `library_load` のたびに払わないため）。

`python -m lipidmix.library.store --ion-mode negative` で MCP の外から事前構築できる
（`main()`。初回構築が MCP クライアントのタイムアウトに当たるとき用）。

deps: record 経由で `.dbs` / `.msp` の iter_records、mcp_core（leaf の DATA_DIR
解決のみ）。session_state / tools_* は import しない。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import tempfile
import zipfile
from pathlib import Path
from typing import Iterable

import msgpack

from lipidmix.core import mcp_core
from lipidmix.core.atomic_io import atomic_write_json, read_text_stable
from lipidmix.library import dbs as dbs_reader
from lipidmix.library import msp as msp_reader

LIBRARY_CACHE_ENV = "LIPIDMIX_LIBRARY_CACHE_DIR"

#: `record.RECORD_FIELDS` が「store のスキーマであり、採点エンジンの入力契約」
#: （spec `docs/superpowers/specs/2026-09-19-msms-spectral-matching-design.md`）
#: なので、13 フィールド全部を持つ。`library_id` は `.dbs` の `MetabolomicsDB/`
#: 配下に複数 DB が同居するとき「どのライブラリのレコードか」を辿る唯一の手掛かり
#: で、`record_index` と合わせて元レコードへ後から遡れることを保証する
#: （controller裁定 2026-09-20: brief の CREATE TABLE は spec とここが食い違って
#: おり、spec を優先して 2 列を足す）。
_SCHEMA = """
CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE record(
    id INTEGER PRIMARY KEY, name TEXT, precursor_mz REAL NOT NULL,
    ion_mode TEXT, adduct TEXT, rt REAL, formula TEXT, inchikey TEXT,
    smiles TEXT, compound_class TEXT, ontology TEXT, spectrum BLOB NOT NULL,
    library_id TEXT, record_index INTEGER);
"""

#: 索引は全件を入れ終えてから張る（1 件ごとに B-tree を更新するより速い。
#: 研究室の pos ライブラリは 1.2 GB ある）。
_INDEX = "CREATE INDEX record_mz ON record(precursor_mz);"

#: 元ファイルの（絶対パス, サイズ, 更新時刻）→ sha256 の対応表。キャッシュ置き場に置く。
#: サイズと更新時刻が前回と同じなら sha256 を計算し直さない——1.2 GB の全読みを
#: `library_load` のたびに払わないため。内容の同一性を保証するのは構築時に計算した
#: sha256 で、この表はその計算を省く近道にすぎない（壊れていたら黙って捨てて計算し直す）。
DIGEST_INDEX_NAME = ".library-digests.json"

_HASH_CHUNK_SIZE = 1 << 20  # 1 MiB

_RECORD_COLUMNS = (
    "name, precursor_mz, ion_mode, adduct, rt, formula, inchikey, "
    "smiles, compound_class, ontology, spectrum, library_id, record_index"
)
_ANNOTATOR_COUNTER = re.compile(r"_\d+$")


def library_id_from_annotator(annotator_id: str | None) -> str | None:
    """MsScanMatchResult の AnnotatorID（`<.dbs のエントリ名>_<n>`）を store の library_id にする。"""
    if annotator_id is None:
        return None
    return _ANNOTATOR_COUNTER.sub("", str(annotator_id))


def _annotators_from_source(source_path_json: str | None) -> tuple[dict[str, str], dict[str, dict]]:
    """対応表を持たない古いキャッシュ用。元の `.dbs` の `Storage` だけを読み直す。

    レコード本体は読まないので、数百 MB の `.dbs` でも作り直すより軽い。元が
    無い・`.dbs` でない（`.msp`/`.lbm2`）ときは空（＝従来の規則だけで引く）。
    """
    try:
        source_path = json.loads(source_path_json) if source_path_json else None
        if not source_path or not Path(source_path).is_file():
            return {}, {}
        meta = dbs_reader.read_storage_meta(source_path) or {}
    except (OSError, ValueError, zipfile.BadZipFile, StopIteration):
        return {}, {}   # 近道を使えないだけ。従来の規則で引く。
    return meta.get("annotator_libraries") or {}, meta.get("annotator_search_params") or {}


def _row_to_record(row) -> dict:
    (name, row_mz, row_ion_mode, adduct, row_rt, formula, inchikey,
     smiles, compound_class, ontology, spectrum_blob, library_id, record_index) = row
    return {
        "name": name, "precursor_mz": row_mz, "ion_mode": row_ion_mode,
        "adduct": adduct, "rt": row_rt, "formula": formula, "inchikey": inchikey,
        "smiles": smiles, "compound_class": compound_class, "ontology": ontology,
        "spectrum": msgpack.unpackb(spectrum_blob, use_list=True),
        "library_id": library_id, "record_index": record_index,
    }


def cache_dir() -> Path:
    """store の置き場所。環境変数優先、無ければ `mcp_core.DATA_DIR` 配下。"""
    override = os.getenv(LIBRARY_CACHE_ENV)
    if override:
        return Path(override)
    return mcp_core.DATA_DIR / ".library-cache"


def source_sha256(path: str | Path) -> str:
    """元ファイルの sha256 を、全体をメモリへ載せずに計算する。"""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(_HASH_CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _store_path_for_digest(path: str | Path, digest: str, cache_dir: Path | None = None) -> Path:
    base = cache_dir if cache_dir is not None else globals()["cache_dir"]()
    stem = Path(path).stem
    return base / f"{stem}-{digest[:16]}.sqlite"


def store_path_for(path: str | Path, cache_dir: Path | None = None) -> Path:
    """`path` に対応する store ファイルのパス。sha256 が鍵。

    呼ぶたびに現在の内容から sha256 を計算し直す（元ファイルが変われば
    別のパスを返す）。大きなファイルに対して繰り返し呼ぶと重複計算になる
    ので、`open_store` 内部では digest を 1 度だけ計算して使い回す。
    """
    digest = source_sha256(path)
    return _store_path_for_digest(path, digest, cache_dir)


def _iter_records_for(path: Path, stats: dict) -> Iterable[dict]:
    """拡張子でリーダを選ぶ。`.msp` はテキスト、それ以外（`.dbs`/`.lbm2`）は
    `dbs.iter_records` に任せる（ZIP かどうかは向こうが中身で判定する）。
    `.msp` のときだけ `stats` に `non_utf8_lines` が入る。"""
    if path.suffix.lower() == ".msp":
        return msp_reader.iter_records(path, stats=stats)
    return dbs_reader.iter_records(path)


def _build(source_path: Path, dest_path: Path, digest: str) -> None:
    """`source_path` を読んで `dest_path` に SQLite store を作る（一時ファイル経由）。"""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(dest_path.parent), suffix=".sqlite.tmp")
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        conn = sqlite3.connect(str(tmp_path))
        try:
            # 書いている最中のファイルは完成まで誰も開かない（rename で差し替える）ので、
            # ジャーナルと同期書き込みは要らない。途中で落ちたら一時ファイルごと捨てる。
            conn.execute("PRAGMA journal_mode = OFF")
            conn.execute("PRAGMA synchronous = OFF")
            conn.executescript(_SCHEMA)
            ion_mode_counts: dict[str, int] = {}
            record_count = 0
            skipped_no_precursor_mz = 0
            reader_stats: dict = {}

            def rows():
                nonlocal record_count, skipped_no_precursor_mz
                for record in _iter_records_for(source_path, reader_stats):
                    if record["precursor_mz"] is None:
                        # `record.precursor_mz REAL NOT NULL`（Important 3 レビュー）:
                        # reader（`.msp` の `msp.iter_records` / `.dbs` の
                        # `_to_record`）は PRECURSORMZ 欠損・パース失敗を契約どおり
                        # None に潰すだけで例外にしない。公開 `.msp`（MassBank 由来
                        # など）には precursor を持たないレコードが普通に混ざるので、
                        # 1 件の欠損でライブラリ全体を使用不能にせず読み飛ばす。
                        skipped_no_precursor_mz += 1
                        continue
                    record_count += 1
                    ion_mode = record["ion_mode"]
                    if ion_mode:
                        ion_mode_counts[ion_mode] = ion_mode_counts.get(ion_mode, 0) + 1
                    yield (
                        record["name"],
                        record["precursor_mz"],
                        record["ion_mode"],
                        record["adduct"],
                        record["rt"],
                        record["formula"],
                        record["inchikey"],
                        record["smiles"],
                        record["compound_class"],
                        record["ontology"],
                        msgpack.packb(record["spectrum"] or []),
                        record["library_id"],
                        record["record_index"],
                    )

            with conn:
                conn.executemany(
                    "INSERT INTO record(name, precursor_mz, ion_mode, adduct, rt, "
                    "formula, inchikey, smiles, compound_class, ontology, spectrum, "
                    "library_id, record_index) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    rows(),
                )
                conn.execute(_INDEX)

                search_params = None
                annotator_libraries: dict = {}
                annotator_search_params: dict = {}
                storage_meta = dbs_reader.read_storage_meta(source_path)
                if storage_meta is not None:
                    search_params = storage_meta.get("search_params") or {}
                    annotator_libraries = storage_meta.get("annotator_libraries") or {}
                    annotator_search_params = storage_meta.get("annotator_search_params") or {}

                meta_rows = [
                    ("record_count", json.dumps(record_count)),
                    ("ion_modes", json.dumps(ion_mode_counts)),
                    ("source_sha256", json.dumps(digest)),
                    ("source_path", json.dumps(str(source_path))),
                    ("search_params", json.dumps(search_params)),
                    ("skipped_no_precursor_mz", json.dumps(skipped_no_precursor_mz)),
                    ("non_utf8_lines", json.dumps(reader_stats.get("non_utf8_lines", 0))),
                    ("annotator_libraries", json.dumps(annotator_libraries)),
                    ("annotator_search_params", json.dumps(annotator_search_params)),
                ]
                conn.executemany("INSERT INTO meta(key, value) VALUES (?, ?)", meta_rows)
        finally:
            conn.close()
        os.replace(str(tmp_path), str(dest_path))
    finally:
        tmp_path.unlink(missing_ok=True)


class LibraryStore:
    """構築済み SQLite を開いて候補検索するだけの薄い読み取り層。"""

    def __init__(self, conn: sqlite3.Connection, path: Path):
        self._conn = conn
        self._path = path
        meta = self._load_meta()
        self.record_count: int = meta["record_count"]
        self._ion_modes: dict[str, int] = meta["ion_modes"]
        self._source_sha256: str = meta["source_sha256"]
        self._search_params = meta["search_params"]
        self._skipped_no_precursor_mz: int = meta["skipped_no_precursor_mz"]
        self._non_utf8_lines: int = meta["non_utf8_lines"]
        self._annotator_libraries: dict[str, str] = meta["annotator_libraries"]
        self._annotator_search_params: dict[str, dict] = meta["annotator_search_params"]

    def rt_used_for(self, annotator_id: str | None) -> bool | None:
        """その注釈器が RT を絞り込みか採点に使ったか。`.dbs` に記録が無ければ None（不明）。"""
        params = self._annotator_search_params.get(str(annotator_id)) if annotator_id else None
        if not params:
            return None
        return bool(params.get("use_time_for_annotation_filtering")
                    or params.get("use_time_for_annotation_scoring"))

    def library_id_for(self, annotator_id: str | None) -> str | None:
        """照合結果の AnnotatorID を、この store の `library_id`（ライブラリ名）にする。

        `.dbs` が保存した注釈器の Key → ライブラリ名の対応を先に引く（Console の
        LBM はパス、MSP は任意名など、AnnotatorID の形は作り手しだいなので推測しない）。
        対応に無い ID だけ、GUI の `<名前>_<n>` 規則（`library_id_from_annotator`）に戻す。
        """
        if annotator_id is None:
            return None
        mapped = self._annotator_libraries.get(str(annotator_id))
        return mapped if mapped is not None else library_id_from_annotator(annotator_id)

    def _load_meta(self) -> dict:
        rows = dict(self._conn.execute("SELECT key, value FROM meta").fetchall())
        if "annotator_libraries" in rows and "annotator_search_params" in rows:
            annotator_libraries = json.loads(rows["annotator_libraries"])
            annotator_search_params = json.loads(rows["annotator_search_params"])
        else:
            # この対応表より前に構築されたキャッシュ。元の `.dbs` から読み直す。
            annotator_libraries, annotator_search_params = _annotators_from_source(
                rows.get("source_path"))
        return {
            "record_count": json.loads(rows["record_count"]),
            "ion_modes": json.loads(rows["ion_modes"]),
            "source_sha256": json.loads(rows["source_sha256"]),
            "search_params": json.loads(rows["search_params"]),
            # .get(): このコード変更より前に構築されたキャッシュ（sha256 キーが
            # 同じまま残っている）には無いキーなので、rebuild なしで開いたときに
            # KeyError にしない。
            "skipped_no_precursor_mz": json.loads(rows.get("skipped_no_precursor_mz", "0")),
            # 同上。以前のリーダは厳密な UTF-8 で読んでいたので、読めた古いキャッシュは 0 で正しい。
            "non_utf8_lines": json.loads(rows.get("non_utf8_lines", "0")),
            "annotator_libraries": annotator_libraries,
            "annotator_search_params": annotator_search_params,
        }

    def summary(self) -> dict:
        return {
            "record_count": self.record_count,
            "ion_modes": dict(self._ion_modes),
            "source_sha256": self._source_sha256,
            "search_params": self._search_params,
            "skipped_no_precursor_mz": self._skipped_no_precursor_mz,
            # 極性不明のレコード。`candidates(ion_mode=...)` はこれを弾かずに残す。
            "records_without_ion_mode": self.record_count - sum(self._ion_modes.values()),
            "non_utf8_lines": self._non_utf8_lines,
        }

    def compound_class_counts(self, top_n: int = 10) -> list[dict]:
        """化合物クラス分布を件数降順で返す（読み取り専用の集計クエリ）。

        `library_load` の要約が使う。`ion_modes` と違いビルド時の meta には
        含めていない（クラス数は上位 N だけ返せば十分で、meta に全件持たせる
        必要が無い）ので、`record` テーブルへ都度クエリする——ただし
        **キャッシュ済み SQLite への SELECT だけ**であり、元ファイル
        （`.dbs`/`.msp`）を読み直すことはない。`compound_class IS NULL` の
        レコード（`.msp` に `COMPOUNDCLASS` が無い等）は除外する。
        """
        rows = self._conn.execute(
            "SELECT compound_class, COUNT(*) AS n FROM record "
            "WHERE compound_class IS NOT NULL "
            "GROUP BY compound_class ORDER BY n DESC, compound_class ASC LIMIT ?",
            (top_n,),
        ).fetchall()
        return [{"name": name, "count": count} for name, count in rows]

    def candidates(
        self,
        precursor_mz: float,
        *,
        mz_tol: float,
        ion_mode: str | None = None,
        rt: float | None = None,
        rt_tol: float | None = None,
    ) -> list[dict]:
        query = f"SELECT {_RECORD_COLUMNS} FROM record WHERE precursor_mz BETWEEN ? AND ?"
        params: list = [precursor_mz - mz_tol, precursor_mz + mz_tol]
        if ion_mode is not None:
            # COLLATE NOCASE: 呼び出し元の大小表記は揃っていない
            # （`.dcl`/pai2 の feature は IonMode Enum の `.name` = "Positive"、
            # store の列は record.ION_MODES = "positive"。Task 11 レビュー
            # Critical 1 —
            # `peak_verification._spectral_match_for_feature` がここを大小区別の
            # BINARY 照合のまま呼び、常に 0 件になっていた。呼び出し側で `lower()`
            # を撒くのではなくここで正規化する——呼び出し側が増えるたびに同じ穴が
            # 開くのを防ぐため）。
            # `ion_mode IS NULL` は「極性不明」として残す。極性ごとに分かれた
            # `.msp`（研究室ライブラリ）は IONMODE 欄を持たないことがあり、
            # NULL を弾くと候補が黙って 0 件になる。
            query += " AND (ion_mode = ? COLLATE NOCASE OR ion_mode IS NULL)"
            params.append(ion_mode)
        if rt is not None:
            if rt_tol is None:
                # Minor 11: `ABS(rt - ?) <= NULL` は NULL（偽）になり `rt IS NULL` の
                # 行しか返らない黙った縮退を起こす。シグネチャ上 rt_tol は省略できる
                # ように見えるので、呼び出し側の取り違えを実行時に検出する。
                raise ValueError(
                    "rt を指定するときは rt_tol も指定してください "
                    "（None のままだと該当行が rt IS NULL の行以外すべて捨てられます）。"
                )
            query += " AND (rt IS NULL OR ABS(rt - ?) <= ?)"
            params.extend([rt, rt_tol])

        rows = self._conn.execute(query, params).fetchall()
        return [_row_to_record(row) for row in rows]

    def record_by_scan_id(self, scan_id: int, *, library_id: str | None,
                          precursor_mz: float, mz_tol: float) -> dict | None:
        """MsScanMatchResult.LibraryID（上流 ScanID = record_index）で参照を 1 件引く。

        `record_index` には索引が無いので、索引のある precursor_mz の窓で先に絞る
        （参照の precursor はスポットの m/z から高々 0.03 Da 程度しか離れない）。
        `library_id` は `.dbs` の複数エントリを区別する。NULL（`.msp` 由来）の行は残す。
        """
        query = (f"SELECT {_RECORD_COLUMNS} FROM record "
                 "WHERE precursor_mz BETWEEN ? AND ? AND record_index = ?")
        params: list = [precursor_mz - mz_tol, precursor_mz + mz_tol, int(scan_id)]
        if library_id is not None:
            query += " AND (library_id = ? OR library_id IS NULL)"
            params.append(library_id)
        row = self._conn.execute(query + " LIMIT 1", params).fetchone()
        return _row_to_record(row) if row else None

    def close(self) -> None:
        self._conn.close()


def _read_digest_index(base: Path) -> dict:
    try:
        data = json.loads(read_text_stable(base / DIGEST_INDEX_NAME))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _digest_for(source_path: Path, base: Path, *, trust_index: bool) -> str:
    """元ファイルの sha256。サイズと更新時刻が記憶と同じなら計算を省く。"""
    stat = source_path.stat()
    key = str(source_path.resolve())
    fingerprint = {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    index = _read_digest_index(base)
    entry = index.get(key)
    if (trust_index and isinstance(entry, dict) and isinstance(entry.get("sha256"), str)
            and all(entry.get(k) == v for k, v in fingerprint.items())):
        return entry["sha256"]
    digest = source_sha256(source_path)
    index[key] = {**fingerprint, "sha256": digest}
    try:
        base.mkdir(parents=True, exist_ok=True)
        atomic_write_json(base / DIGEST_INDEX_NAME, index)
    except OSError:
        pass  # 近道を覚えられないだけ。次回また計算する。
    return digest


def open_store(path: str | Path, *, cache_dir: Path | None = None, rebuild: bool = False) -> LibraryStore:
    """`path` の store を開く。無ければ（または `rebuild=True` なら）構築してから開く。

    `rebuild=True` は元ファイルの破損を疑うときに使うので、記憶した sha256 も信じない。
    """
    source_path = Path(path)
    base = cache_dir if cache_dir is not None else globals()["cache_dir"]()
    digest = _digest_for(source_path, base, trust_index=not rebuild)
    dest_path = _store_path_for_digest(source_path, digest, cache_dir=base)

    if rebuild or not dest_path.exists():
        _build(source_path, dest_path, digest)

    conn = sqlite3.connect(str(dest_path))
    try:
        return LibraryStore(conn, dest_path)
    except Exception:
        # 構築は原子的なので確率は低いが、LibraryStore.__init__（_load_meta）が
        # 例外を投げたときに Connection を漏らさない。
        conn.close()
        raise


def main(argv: list[str] | None = None) -> int:
    """store を MCP の外で事前構築する CLI。

    研究室ライブラリ（pos ≈ 1.2 GB）の初回構築は数分かかり、MCP クライアントの
    タイムアウトに当たりうる。手元で 1 度走らせておけば、以後の `library_load` は
    キャッシュを開くだけで済む。パスの解決規則は `library_load` と同じ
    （`resolve_library_path`）。出力にはファイル名だけを出し、置き場所は出さない。

        python -m lipidmix.library.store --ion-mode negative
        python -m lipidmix.library.store --file <path> --rebuild
    """
    import argparse
    import sys
    import time

    from lipidmix.core.path_resolvers import LibraryPathError, resolve_library_path

    parser = argparse.ArgumentParser(prog="python -m lipidmix.library.store",
                                     description="参照ライブラリの照合用 store を事前構築する。")
    parser.add_argument("--ion-mode", choices=("positive", "negative"),
                        help="どちらの極性の設定（MSDIAL_MSP_POS / MSDIAL_MSP_NEG か "
                             "lipidmix.local.toml の [library] msp_positive / msp_negative）を使うか")
    parser.add_argument("--file", help="ライブラリのパス（極性の設定より優先）")
    parser.add_argument("--rebuild", action="store_true", help="キャッシュを無視して作り直す")
    args = parser.parse_args(argv)

    try:
        resolved = resolve_library_path(args.file, ion_mode=args.ion_mode)
    except LibraryPathError as exc:
        print(f"{exc.code}: {exc.message}", file=sys.stderr)
        return 2
    if not resolved:
        print("LIBRARY_NOT_FOUND: 参照ライブラリが見つかりませんでした。", file=sys.stderr)
        return 2

    started = time.perf_counter()
    library = open_store(resolved, rebuild=args.rebuild)
    try:
        summary = library.summary()
    finally:
        library.close()
    print(json.dumps({
        "file": Path(resolved).name,
        "record_count": summary["record_count"],
        "ion_modes": summary["ion_modes"],
        "records_without_ion_mode": summary["records_without_ion_mode"],
        "skipped_no_precursor_mz": summary["skipped_no_precursor_mz"],
        "non_utf8_lines": summary["non_utf8_lines"],
        "source_sha256": summary["source_sha256"][:16],
        "elapsed_s": round(time.perf_counter() - started, 1),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
