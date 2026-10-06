# AlignmentChromPeakFeature — MessagePack Key 番号表

`.arf`（アライン後・サンプル別ピーク）1 行ぶんの MessagePack 配列に対する Key 番号の正準表。
この行は `AlignmentSpotProperty` の Key8 `AlignedPeakProperties` の要素として入っている
（→ [AlignmentSpotProperty.md](AlignmentSpotProperty.md)）。

**リーダのインデックスを変える前にこの表を見る。推測で直さない。**
関連: [../output_format/arf.md](../output_format/arf.md)（フィールドの意味・ギャップフィルの扱い）

## 出典

MS-DIAL（MsdialWorkbench）の C# クラス `CompMs.MsdialCore.DataObj.AlignmentChromPeakFeature`
（上流の `src/MSDIAL5/MsdialCore/DataObj/AlignmentChromPeakFeature.cs`）に付いた `[Key(N)]` 属性から抽出した。

- 上流: <https://github.com/systemsomicslab/MsdialWorkbench>（LGPL-3.0）
- 照合したコミット: `45a531c`（2026-09-02）
- **この表は Key 番号とメンバ名・型の対応だけを写したもので、上流のソースコードは含まない。**
- MS-DIAL 側でクラスが変わったら、上流の同ファイルを開いて `[Key(N)]` を読み直し、
  この表を更新する。実装（`metabolomix/*/reader.py`）から逆算して直してはいけない。

表は **Key 番号の昇順**で並べている（上流のソース上の宣言順とは一致しない）。

## Key 番号表

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `FileID` | **真の FileID はここ**（配列先頭から文字列を探す実装とは別）。 |
| 1 | `string` | `FileName` | アライン結果作成時に取り込んだ名前。後のファイルプロパティ編集で現在の解析ファイル名と食い違いうる。 |
| 2 | `int` | `MasterPeakID` | 全ピーク通しの連番。**負値（-2）はギャップフィル＝未検出の補間値**。検出ピークと混ぜてはいけない。 |
| 3 | `int` | `PeakID` | 同一次元内の連番。 |
| 4 | `int` | `ParentPeakID` | LC-IM-MS/MS で親ピークを指す。 |
| 6 | `long` | `SeekPointToDCLFile` | `.dcl` 内のシーク位置。デコンボリューション済み MS/MS の所在。 |
| 7 | `int` | `MS1RawSpectrumID` |  |
| 8 | `int` | `MS1RawSpectrumIDatAccumulatedMS1` | LC-IM-MS/MS 用。 |
| 9 | `int` | `MS2RawSpectrumID` | 代表 ID。 |
| 10 | `Dictionary<int, double>` | `MS2RawSpectrumID2CE` | 衝突エネルギーの写像。**非空なら MS/MS 取得済み**の判定に使える。 |
| 11 | `int` | `ChromScanIdLeft` |  |
| 12 | `int` | `ChromScanIdTop` |  |
| 13 | `int` | `ChromScanIdRight` |  |
| 14 | `ChromXs` | `ChromXsLeft` | 入れ子（ChromXs）。ピーク左端。 |
| 15 | `ChromXs` | `ChromXsTop` | 入れ子（ChromXs）。**ピークトップ。RT はここから取る**。 |
| 16 | `ChromXs` | `ChromXsRight` | 入れ子（ChromXs）。ピーク右端。RT の取得では Key15 が無いときの代替にしかならない。 |
| 17 | `double` | `PeakHeightLeft` |  |
| 18 | `double` | `PeakHeightTop` |  |
| 19 | `double` | `PeakHeightRight` |  |
| 20 | `double` | `PeakAreaAboveZero` |  |
| 21 | `double` | `PeakAreaAboveBaseline` |  |
| 22 | `double` | `Mass` | GC-MS では定量質量。 |
| 23 | `IonMode` | `IonMode` | 0=Positive / 1=Negative / 2=Both / -1=Unknown。 |
| 24 | `string` | `Name` |  |
| 25 | `Formula` | `Formula` | 入れ子。index0 が組成式文字列。 |
| 26 | `string` | `Ontology` |  |
| 27 | `string` | `SMILES` |  |
| 28 | `string` | `InChIKey` |  |
| 29 | `double` | `CollisionCrossSection` |  |
| 31 | `Dictionary<int, List<int>>` | `MSRawID2MspIDs` | 閾値を超えた候補の ID 群。 |
| 33 | `List<int>` | `TextDbIDs` | テキスト DB 由来の候補 ID（任意）。 |
| 34 | `Dictionary<int, MsScanMatchResult>` | `MSRawID2MspBasedMatchResult` |  |
| 35 | `MsScanMatchResult` | `TextDbBasedMatchResult` |  |
| 36 | `IonFeatureCharacter` | `PeakCharacter` | 入れ子。付加体・電荷はここに入る。 |
| 37 | `ChromatogramPeakShape` | `PeakShape` | 入れ子（ChromatogramPeakShape）。index0=EstimatedNoise, index1=SignalToNoise。 |
| 38 | `int` | `MS1RawSpectrumIdTop` |  |
| 39 | `int` | `MS1RawSpectrumIdLeft` |  |
| 40 | `int` | `MS1RawSpectrumIdRight` |  |
| 41 | `int` | `MS1AccumulatedMs1RawSpectrumIdTop` | LC-IM-MS/MS 用。 |
| 42 | `int` | `MS1AccumulatedMs1RawSpectrumIdLeft` | LC-IM-MS/MS 用。 |
| 43 | `int` | `MS1AccumulatedMs1RawSpectrumIdRight` | LC-IM-MS/MS 用。 |
| 44 | `double` | `NormalizedPeakHeight` |  |
| 45 | `double` | `NormalizedPeakAreaAboveZero` |  |
| 46 | `double` | `NormalizedPeakAreaAboveBaseline` |  |
| 47 | `MsScanMatchResultContainer` | `MatchResults` | 多行プロパティ。private バックフィールドには Key が振られていない。 |
| 49 | `int` | `MSDecResultIdUsed` | 既定 -1。-1 のときは MasterPeakID を使う運用。 |
| 50 | `string` | `Protein` | プロテオミクス用。 |
| 51 | `int` | `ProteinGroupID` | プロテオミクス用。既定 -1。 |
| 52 | `bool` | `IsManuallyModifiedForQuant` |  |

## 欠番（勝手に埋めない）

| Key | 状態 |
|---:|---|
| 5 | この版のクラス定義には存在しない。 |
| 30 | `MspID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 32 | `TextDbID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 48 | この版のクラス定義には存在しない。 |

## 入れ子型のデコード

MessagePack 上では、下の型は配列として入れ子で載る。**この内訳は上流のクラス定義には
書かれておらず**、本リポジトリのリーダが実データに対して確認したもの。

| 型 | MessagePack 上の形 | 取り出し方 |
|---|---|---|
| `ChromXs` | `[[種別, [値, ...]], ...]` | 種別 1=RT / 2=RI / 3=m/z / 4=ドリフト時間。値は `[1][0]`（`_convert_to_times()` / `_extract_times()`） |
| `Formula` | `[組成式文字列, 質量, ...]` | index0 |
| `AdductIon` | `[質量差, 電荷, 付加体文字列, ...]` | index2 |
| `ChromatogramPeakShape` | `[EstimatedNoise, SignalToNoise, ...]` | index0 / index1 |
| `IonMode` | int | 0=Positive / 1=Negative / 2=Both / -1=Unknown |

## 本リポジトリのリーダが読む Key

`metabolomix/arf/reader.py` の `_convert_to_alignment_feature()`:
Key 0（FileID）, 2（ギャップフィル判定）, 3, 9, 10（MS/MS 取得判定）, 15→16（RT。Key15 が無い場合のみ
Key16 で代替）, 18, 20, 21, 22, 37（index0=EstimatedNoise / index1=S/N）。
