# MoleculeMsReference — MessagePack Key 番号表、および `.dbs` / `.lbm2` のバイナリ枠組み

MS/MS 参照ライブラリ（スペクトル照合用）1 レコードぶんの MessagePack 配列に対する
Key 番号の正準表。`.lbm2`（脂質 in-silico ライブラリの単独ファイル）と `.msp2`
（同じ枠組みのテキストライブラリ変換版）はこのクラスの `List<MoleculeMsReference>` を
そのままファイル全体として持つ。MS-DIAL の**プロジェクト保存**が書く `.dbs`
（`DataBaseStorage`、ZIP）は、内部の `DataBase` エントリとして**同一の枠組み**を
再利用する（詳細は下記）。

**リーダのインデックスを変える前にこの表を見る。推測で直さない。**
このクラスには 1 つの罠がある。

1. **Key 27（`DatabaseUniqueIdentifier`）は宣言順で Key 20 の直後・Key 21 より前に
   書かれている。** Key 番号としては末尾寄り（29 個中 28 番目）だが、ソース上の
   `[Key(N)]` 属性の並びはそこではない。**MessagePack の並び順は常に Key 番号の昇順
   であり、宣言順ではない**——この点は他の Key 表と同じだが、このクラスは
   宣言順と Key 番号の乖離が最も大きい。宣言順で読んで「そこにあるはずの Key」を
   決め打ちしない。

関連:
- [docs/superpowers/specs/2026-09-19-msms-spectral-matching-design.md](../superpowers/specs/2026-09-19-msms-spectral-matching-design.md)
  §2.2, §2.4–2.6（実測値・ファイルサイズ・チャンク数・レコード件数はここが正準。
  本文書には**丸写ししない**）
- [AlignmentSpotProperty.md](AlignmentSpotProperty.md) / [ChromatogramPeakFeature.md](ChromatogramPeakFeature.md)
  （`ChromXs` / `Formula` / `AdductIon` の入れ子デコードは共通なのでそちらを見る）

## 出典

MS-DIAL（MsdialWorkbench）の C# クラス `CompMs.Common.Components.MoleculeMsReference`
（上流の `src/Common/CommonStandard/Components/MoleculeMsReference.cs`）に付いた
`[Key(N)]` 属性から抽出した。

- 上流: <https://github.com/systemsomicslab/MsdialWorkbench>（LGPL-3.0）
- 照合したコミット: `afd5f9522`（2026-09-20、`master` ブランチの HEAD）
- **この表は Key 番号とメンバ名・型の対応だけを写したもので、上流のソースコードは含まない。**
- MS-DIAL 側でクラスが変わったら、上流の同ファイルを開いて `[Key(N)]` を読み直し、
  この表を更新する。実装（Task 3 で書く `metabolomix/library/*`）から逆算して直しては
  いけない。

表は **Key 番号の昇順**で並べている（上流のソース上の宣言順とは一致しない）。

## Key 番号表（0–28、欠番なし）

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `ScanID` | ライブラリ内の連番。`.dbs` の `MoleculeDataBase.Refer()` が `MsScanMatchResult.LibraryID` との突き合わせに使う参照キー（`ScanID == LibraryID` を優先し、ずれていれば線形探索でフォールバック）。 |
| 1 | `double` | `PrecursorMz` |  |
| 2 | `ChromXs` | `ChromXs` | 入れ子。デコードは [AlignmentSpotProperty.md](AlignmentSpotProperty.md#入れ子型のデコード) と同じ。 |
| 3 | `IonMode` | `IonMode` | 0=Positive / 1=Negative / 2=Both（下記「IonMode の列挙」参照。他クラスの `IonMode` と異なり `-1=Unknown` は無い）。 |
| 4 | `List<SpectrumPeak>` | `Spectrum` | **このレコードの実体**（照合対象の MS/MS フラグメント列）。要素は `SpectrumPeak`（下記「入れ子型: SpectrumPeak」）。 |
| 5 | `string` | `Name` | 既定 `""`。 |
| 6 | `Formula` | `Formula` | 入れ子。index0 が組成式文字列（[AlignmentSpotProperty.md](AlignmentSpotProperty.md#入れ子型のデコード) と同じ）。 |
| 7 | `string` | `Ontology` | 既定 `""`。 |
| 8 | `string` | `SMILES` | 既定 `""`。 |
| 9 | `string` | `InChIKey` | 既定 `""`。 |
| 10 | `AdductIon` | `AdductType` | 入れ子。index2 が付加体文字列（`[M-H]-` など）。[AlignmentSpotProperty.md](AlignmentSpotProperty.md#入れ子型のデコード) と同じ。 |
| 11 | `double` | `CollisionCrossSection` |  |
| 12 | `List<IsotopicPeak>` | `IsotopicPeaks` | 要素は `IsotopicPeak`（下記「入れ子型: IsotopicPeak」）。既定は空リスト。 |
| 13 | `double` | `QuantMass` | GC-MS プロジェクト用。 |
| 14 | `string` | `CompoundClass` | 脂質クラス（PC, TG, …）。`Ontology`（Key7）と役割が重なるプロパティで、`OntologyOrCompoundClass`（`IgnoreMember`）が `Ontology` 優先で選ぶ。 |
| 15 | `string` | `Comment` | 既定 `""`。 |
| 16 | `string` | `InstrumentModel` | 既定 `""`。 |
| 17 | `string` | `InstrumentType` | 既定 `""`。 |
| 18 | `string` | `Links` | 既定 `""`。複数データベースへのリンクをセミコロン区切りで持つ想定（ソースコメント）。 |
| 19 | `float` | `CollisionEnergy` |  |
| 20 | `int` | `DatabaseID` | binbase・fastaDB 等の外部 DB 用（ソースコメント）。 |
| 21 | `int` | `Charge` |  |
| 22 | `int` | `MsLevel` |  |
| 23 | `float` | `RetentionTimeTolerance` | 既定 0.05。テキストライブラリ照合専用（ソースコメント）。 |
| 24 | `float` | `MassTolerance` | 既定 0.05。テキストライブラリ照合専用。 |
| 25 | `float` | `MinimumPeakHeight` | 既定 1000。テキストライブラリ照合専用。 |
| 26 | `bool` | `IsTargetMolecule` | 既定 `true`。テキストライブラリ照合専用。 |
| 27 | `string` | `DatabaseUniqueIdentifier` | binbase・fastaDB 等の外部 DB 用（ソースコメント）。**宣言順は Key 20 と Key 21 の間**（上の「罠」参照）。 |
| 28 | `string` | `FragmentationCondition` |  |

## IonMode の列挙

`src/Common/CommonStandard/Enum/CommonEnums.cs`:

```csharp
public enum IonMode { Positive, Negative, Both }
```

`Positive=0` / `Negative=1` / `Both=2`。既定値の割り当てが無い点に注意——
`AlignmentSpotProperty` / `ChromatogramPeakFeature` の `IonMode` で見られる
`-1=Unknown`（未設定を表す値）はここには存在しない。

## 入れ子型: SpectrumPeak

`List<SpectrumPeak>`（Key 4 `Spectrum`）の要素。出典:
`src/Common/CommonStandard/Components/SpectrumPeak.cs`。

MS-DIAL 本体のスコアリングは `SpectrumPeak.Mass` / `.Intensity` の対だけを使う。
それ以外のプロパティのうち `[Key(N)]` が付いていないもの（`Resolution` `Charge`
`IsotopeFrag` `IsotopeParentPeakID` `IsotopeWeightNumber` `IsMatched`）は
**`[IgnoreMember]` でシリアライズ対象から除外されている**——`ChromatogramPeakFeature`
の「Obsolete だが Key は残っている」罠とは逆で、こちらは本当に欠番（ワイヤフォーマット
に存在しない）。

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `double` | `Mass` |  |
| 1 | `double` | `Intensity` |  |
| 2 | `string` | `Comment` |  |
| 6 | `PeakQuality` | `PeakQuality` |  |
| 7 | `int` | `PeakID` |  |
| 11 | `SpectrumComment` | `SpectrumComment` | `[Flags]` enum（フラグメント種別の注釈: `b` `y` `c` `z` `precursor` 等をビットで持つ）。 |
| 12 | `bool` | `IsAbsolutelyRequiredFragmentForAnnotation` |  |
| 13 | `double` | `FragmentationScore` | 既定 0。 |

欠番（3, 4, 5, 8, 9, 10）: `[Key]` 属性が付いていないプロパティは存在しない
（`Resolution` `Charge` `IsotopeFrag` 等は `[IgnoreMember]` なので、そもそも
Key 番号の割り当て対象ではない）。勝手に埋めない。

## 入れ子型: IsotopicPeak

`List<IsotopicPeak>`（Key 12 `IsotopicPeaks`）の要素。出典:
`src/Common/CommonStandard/DataObj/Property/IsotopeProperty.cs`
（同名で `[Key]` の付いていない `src/Common/CommonStandard/FormulaData/IsotopicPeak.cs`
という別クラスも存在するので、参照先を取り違えない）。

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `double` | `RelativeAbundance` |  |
| 1 | `double` | `Mass` |  |
| 2 | `double` | `MassDifferenceFromMonoisotopicIon` |  |
| 3 | `string` | `Comment` | 既定 `""`。 |
| 4 | `double` | `AbsoluteAbundance` |  |

宣言順は 0, 4, 1, 2, 3（Key 4 が Key 1 より先に書かれている）。**MoleculeMsReference
の Key 27 と同じ「宣言順と Key 番号が一致しない」パターン**。

`AlignmentSpotProperty` の Key 53 `IsotopicPeaks` も同じ型。実データ（kidney pos、
2026-10-09 確認）では全スポットに M, M+1, M+2 の 3 要素が入り、Key 1 `Mass` は
単同位体 + 1.004667751·k（Key 2 と同じ差）の計算値で、実測の質量ではない。Key 0 は M を 100 とする
相対強度、Key 4 は絶対強度。キュレーションの MS1 パネル（`curation/isotope.py`）が Key 0・1 を読む。

## `LargeListMessagePack` の枠組み（`.lbm2` / `.msp2` / `.dbs` の `DataBase` 共通）

出典: `src/Common/CommonStandard/MessagePack/LargeListMessagePack.cs`。

```
ExtensionTypeCode = 99
HeaderSize        = 11
OffsetCutoff      = 1073741824   // 1 GiB. これを超えたら次のチャンクへ切り替える
```

ファイルはチャンクの列。各チャンクのバイト列:

```
c9                    MessagePack ext32 マーカ
<4 bytes BE>          ext 長 = LZ4 本体長 + 5
63                    type = 99（ExtensionTypeCode）
d2 <4 bytes BE>       展開後サイズ（int32、force int32 format）
<本体>                LZ4 block（frame ではない。生の LZ4Codec.Encode 出力）
```

先頭 6 バイト（`c9` + 4B + `63`）が拡張フォーマットヘッダ、続く 5 バイト
（`d2` + 4B）が展開後サイズで、合計 `HeaderSize = 11`。

展開すると MessagePack の配列が現れ、要素が `MoleculeMsReference`（または
`MoleculeDataBase.Save` が直接 `LargeListMessagePack.Serialize` する文脈では
同じく `MoleculeMsReference` のコレクション）。

**罠**: 展開後の先頭 5 バイトは配列ヘッダ用に予約されている
（書き込み側は `offset = 5` から要素のシリアライズを始め、その後
`MessagePackBinary.WriteArrayHeader` で先頭に書き戻す）が、`WriteArrayHeader` は
**要素数が小さいと短形式**（fixarray = 1 バイト、または `0xdc` = 2 バイト長）で
書く。読み込み側（`DeserializeList`）は `ReadArrayHeader` で要素数だけ読み、
**その後のオフセットは実際に書かれたヘッダ長を見ずに無条件で `offset = 5` に
した位置から要素を読み始める**。つまり書き込み・読み込みの双方が「ヘッダの
実サイズに関わらず要素はオフセット 5 から」という同じ前提で動いており、
ヘッダとして使われなかった残りのバイトは未使用のまま埋まる。**5 バイト固定の
`0xdd`（配列 32bit 長）ヘッダだと決め打ってパースすると、要素数が少ない
最終チャンクで壊れる。**

## `.dbs` の ZIP 構造

出典: `src/MSDIAL5/MsdialCore/DataObj/DataBaseStorage.cs`
（`Save` / `Load`）と `DataBaseItem.cs`（`Save` / `TryCurrentLoad`）。
ファイル名の導出は `src/MSDIAL5/MsdialCore/DataObj/MsdialDataStorage.cs` の
`GetNewMspFileName` / `GetNewZippedDatabaseFileName` / `GetDataBasesFileName`。

`.dbs` は `System.IO.Compression.ZipArchive` で書かれた ZIP。エントリ構成:

```
MetabolomicsDB/<DataBaseID>/DataBase     ProteomicsDB/<DataBaseID>/DataBase
MetabolomicsDB/<DataBaseID>/Annotators/<AnnotatorID>/...
EadLipidomicsDB/<DataBaseID>/DataBase
Storage
```

（`<DataBaseID>` は `DataBaseItem.DataBaseID`＝`DataBase.Id`。実データでは
ユーザが指定したライブラリファイル由来の名前や `MS-FINDER` が入る。
`ProteomicsDB` / `EadLipidomicsDB` は本サーバの対象外。）

- **`<種別>/<DataBaseID>/DataBase` エントリ**: `MoleculeDataBase.Save()` が
  `LargeListMessagePack.Serialize(stream, Database)` で直接書く。つまり
  **上の「`LargeListMessagePack` の枠組み」と完全に同一のバイト列**であり、
  スタンドアロンの `.lbm2` / `.msp2` ファイルをそのまま ZIP エントリに
  埋め込んだのと同じもの。Task 3 のパーサはこのエントリと `.lbm2` を
  同じ関数で読める。
- **`Storage` エントリ**: `DataBaseStorage` オブジェクト自身
  （使った DB 一覧・元ファイルパス・注釈器パラメータ）を
  `MessagePackDefaultHandler.SaveToStream` で書く。これは
  `LZ4MessagePackSerializer.Serialize`（`MessagePack.LZ4` 名前空間、
  本リポジトリには無い外部 NuGet パッケージ側の実装）を直接呼ぶ**単発の
  LZ4 圧縮**であり、上のチャンク分割された `LargeListMessagePack` とは
  **実装としては別物**。拡張タイプコードは同じ 99 だが、チャンク分割は無い
  （1 回の ext99+LZ4 ブロックで全体を包む）。

  この単発フレームのバイト単位の並びは、**上流ソースからは**（外部 NuGet の
  ため）確認できないが、**実データでの検証は済んでいる**——これは上流ソースの
  確認とは別種の根拠なので区別して書く。実データ（`.dbs` の `Storage` エントリ、
  実測値は spec §2.2 が正準）を、上の「`LargeListMessagePack` の枠組み」節と
  **まったく同じ 11 バイトヘッダ解釈**（`c9` + ext長4B + `63` + `d2` + 展開後
  サイズ4B + LZ4 block）で読んだところ、`lz4.block.decompress` → `msgpack.unpackb`
  で正しく復号できることを確認済み（controller が実測。本文書には実測値そのもの
  は丸写ししない）。**したがって `Storage` は「チャンクが 1 個だけの
  `LargeListMessagePack`」と同じ読み方で読める**——Task 3 のパーサは `Storage` と
  `DataBase` を同じ展開関数（ext99+LZ4 の 1 チャンク分の展開ロジック）に通せる。
  違いはチャンク分割が無いこと（末尾での次チャンク判定が不要）と、展開後の
  MessagePack が配列ではなく `MetabolomicsDataBases` / `ProteomicsDataBases` /
  `EadLipidomicsDatabases` を持つマップであることの 2 点だけ。
  復号すると `MsRefSearchParameterBase` を含む検索パラメータが読める
  （下記）。

ファイル名の導出:

| 関数 | 出力 | 中身 |
|---|---|---|
| `GetNewMspFileName` | `<title>_Loaded.msp2` | 読み込んだ MspDB のシリアライズ（`.msp` 指定時のみ中身あり。`.lbm2` 経路では 0 バイト） |
| `GetNewZippedDatabaseFileName` | `<title>_Loaded.msp2.zip` | 同名の `.zip`（`GetNewMspFileName` の結果に `.zip` を足しただけ） |
| `GetDataBasesFileName` | `<title>_Loaded.msp2.dbs` | `DataBaseStorage`（本文書が扱う ZIP） |

## `.lbm2` / `.msp2`（スタンドアロンのライブラリファイル）

出典: `src/MSDIAL5/MsdialCore/Utility/LibraryHandler.cs`、
`src/Common/CommonStandard/Parser/MspFileParcer.cs`、
`src/Common/CommonStandard/MessagePack/MoleculeMsRefMethods.cs`。

`LibraryHandler.ReadMsLibrary` は拡張子 `.lbm` / `.lbm2` を
`ReadLipidMsLibrary` → `MspFileParser.ReadSerializedLbmLibrary`
（`.lbm2` のときは `ReadSerializedMspObject`）→
`MoleculeMsRefMethods.LoadMspFromFile` →
`MessagePackDefaultHandler.LoadLargerListFromFile<MoleculeMsReference>` →
`LargeListMessagePack.Deserialize<MoleculeMsReference>` の経路で読む。
**ZIP ラップは無い**——`.lbm2` / `.msp2` はファイル全体がそのまま
「`LargeListMessagePack` の枠組み」節のチャンク列であり、`.dbs` 内の
`DataBase` エントリと**バイト単位で同一の形式**。

`ReadSerializedLbmLibrary` は全件を読んだ後に `queryCheck`（イオンモード・
脂質クラス等）でメモリ上フィルタする。**ファイル自体には絞り込みが
掛かっていない**——絞り込み済みなのは `.dbs` の `DataBase` エントリの方
（run 時の `LipidQueryContainer` で選択された DB が保存し直されるため）。

## `MsRefSearchParameterBase` の Key 表（0–19、欠番なし）

出典: `src/Common/CommonStandard/Parameters/MsRefSearchParameterBase.cs`。
`.dbs` の `Storage` エントリを復号した中に、注釈器パラメータの一部として
このオブジェクトが入っている。

**Key 19（`AndromedaScoreCutOff`）は宣言順で Key 12 の直後・Key 13 より前に
書かれている**（`MoleculeMsReference` の Key 27 と同じ「宣言順と Key 番号が
一致しない」パターン）。

**Key 7（`RelativeAmpCutoff`）と Key 8（`AbsoluteAmpCutoff`）は、MS-DIAL が
採点前に測定スペクトル（クエリ側、ライブラリ側ではない）へ掛ける強度の
足切り閾値**である。`RelativeAmpCutoff` は最大強度に対する相対値、
`AbsoluteAmpCutoff` は絶対強度の下限で、両方を下回るフラグメントは
ドット積・一致本数などのスコア計算に入る前に測定スペクトルから間引かれる。
**この 2 つが正規化（本設計でいう「採点前処理」）の唯一の可変部分**であり、
足切りの設定を変えた run 同士でスコアが合わないときに「移植の誤り」と
誤診断しないための手がかりになる。

**Key 16（`IsUseTimeForAnnotationScoring`）は総合スコアの形を変える。** 真なら
`GetTotalScore` に RT 項（`RtSimilarity`）が 1 つ増え、候補の並び順が変わる。
**既定は false だが実 run では true のことがある**（実測: aging mice kidney neg の
`Dataset_2026_09_09_17_28_59_Loaded.msp2.dbs` は Key 15/16 が両方 true、
`RtTolerance=2.0`）。既定値を仮定せず必ず `Storage` から読むこと——
`metabolomix/library/dbs.py` の `_SEARCH_PARAM_KEYS` が 15〜18 を
`use_time_for_annotation_*` / `use_ccs_for_annotation_*` として拾い、
`library_match_feature` の `scoring.use_rt` に反映される。

| Key | 型 | メンバ | 既定 | 備考 |
|---:|---|---|---|---|
| 0 | `float` | `MassRangeBegin` | 0 |  |
| 1 | `float` | `MassRangeEnd` | 2000 |  |
| 2 | `float` | `RtTolerance` | 100.0 |  |
| 3 | `float` | `RiTolerance` | 100.0 |  |
| 4 | `float` | `CcsTolerance` | 20.0 |  |
| 5 | `float` | `Ms1Tolerance` | 0.01 |  |
| 6 | `float` | `Ms2Tolerance` | 0.025 |  |
| 7 | `float` | `RelativeAmpCutoff` | 0 | **測定スペクトルへの相対強度足切り**（採点前処理。上記参照）。 |
| 8 | `float` | `AbsoluteAmpCutoff` | 0 | **測定スペクトルへの絶対強度足切り**（採点前処理。上記参照）。 |
| 9 | `float` | `SquaredWeightedDotProductCutOff` | 0.36（= 0.6²） | 非二乗値は `IgnoreMember` の `WeightedDotProductCutOff`（`sqrt` で導出）。 |
| 10 | `float` | `SquaredSimpleDotProductCutOff` | 0.36（= 0.6²） | 非二乗値は `IgnoreMember` の `SimpleDotProductCutOff`。 |
| 11 | `float` | `SquaredReverseDotProductCutOff` | 0.64（= 0.8²） | 非二乗値は `IgnoreMember` の `ReverseDotProductCutOff`。 |
| 12 | `float` | `MatchedPeaksPercentageCutOff` | 0.25 |  |
| 13 | `float` | `TotalScoreCutoff` | 0.8 |  |
| 14 | `float` | `MinimumSpectrumMatch` | 3 |  |
| 15 | `bool` | `IsUseTimeForAnnotationFiltering` | false | RT 窓で候補を絞るか。 |
| 16 | `bool` | `IsUseTimeForAnnotationScoring` | false | **総合スコアに RT 項を足すか**（下記参照）。 |
| 17 | `bool` | `IsUseCcsForAnnotationFiltering` | false | IM-MS 用。LC-MS の run では false。 |
| 18 | `bool` | `IsUseCcsForAnnotationScoring` | false | 同上。 |
| 19 | `float` | `AndromedaScoreCutOff` | 0.1 | **宣言順は Key 12 と Key 13 の間**（上記参照）。 |

## ドリフト確認の手順

```bash
git -C "C:/Users/yuu18/source/repos/MsdialWorkbench" diff --stat afd5f9522 origin/master -- \
  "*MoleculeMsReference.cs" "*SpectrumPeak.cs" "*IsotopeProperty.cs" \
  "*LargeListMessagePack.cs" "*DataBaseStorage.cs" "*DataBaseItem.cs" \
  "*MsdialDataStorage.cs" "*LibraryHandler.cs" "*MspFileParcer.cs" \
  "*MoleculeMsRefMethods.cs" "*MessagePackHandler.cs" "*MsRefSearchParameterBase.cs" \
  "*CommonEnums.cs"
```

`git fetch` していない場合は `origin/master` の代わりに手元の `master` と
比較する。差分が出たら、該当ファイルを開いて `[Key(N)]` を読み直し、この表を
更新する（推測で直さない）。
