# ChromatogramPeakFeature — MessagePack Key 番号表

`.pai2`（個別測定で検出されたピーク）1 行ぶんの MessagePack 配列に対する Key 番号の正準表。

**リーダのインデックスを変える前にこの表を見る。推測で直さない。**
このクラスには 2 つの罠がある。

1. **Key0〜10 と Key43 は C# 上 Obsolete だが、シリアライズ対象からは外れていない。**
   ピークの位置・高さ・面積・質量はこの番号群にしか無く、`.pai2` リーダは実際にここを読む。
   「非推奨だから使わない」と判断して番号を空けてはいけない。
2. **Key49 は過去に別メンバへ振られていた番号の再利用**で、Key50 は Key49 の private
   バックフィールドに付いた Key（同じ内容が二重に載る）。

関連: [../output_format/pai2.md](../output_format/pai2.md)（フィールドの意味）、
MS/MS 本体は同名の `.dcl`（[../output_format/dcl.md](../output_format/dcl.md)）

## 出典

MS-DIAL（MsdialWorkbench）の C# クラス `CompMs.MsdialCore.DataObj.ChromatogramPeakFeature`
（上流の `src/MSDIAL5/MsdialCore/DataObj/ChromatogramPeakFeature.cs`）に付いた `[Key(N)]` 属性から抽出した。

- 上流: <https://github.com/systemsomicslab/MsdialWorkbench>（LGPL-3.0）
- 照合したコミット: `45a531c`（2026-09-02）
- **この表は Key 番号とメンバ名・型の対応だけを写したもので、上流のソースコードは含まない。**
- MS-DIAL 側でクラスが変わったら、上流の同ファイルを開いて `[Key(N)]` を読み直し、
  この表を更新する。実装（`metabolomix/*/reader.py`）から逆算して直してはいけない。

表は **Key 番号の昇順**で並べている（上流のソース上の宣言順とは一致しない）。

## Key 番号表

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `ChromScanIdLeft` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 1 | `int` | `ChromScanIdTop` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 2 | `int` | `ChromScanIdRight` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 3 | `ChromXs` | `ChromXsLeft` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 4 | `ChromXs` | `ChromXsTop` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 5 | `ChromXs` | `ChromXsRight` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 6 | `double` | `PeakHeightLeft` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 7 | `double` | `PeakHeightTop` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 8 | `double` | `PeakHeightRight` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 9 | `double` | `PeakAreaAboveZero` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 10 | `double` | `PeakAreaAboveBaseline` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 11 | `int` | `MasterPeakID` | 全ピーク通しの連番。pai2 リーダはこれを `id` として返す。 |
| 12 | `int` | `PeakID` | 同一次元内の連番。 |
| 13 | `int` | `ParentPeakID` | LC-IM-MS/MS で親ピークを指す。 |
| 15 | `long` | `SeekPointToDCLFile` | `.dcl` 内のシーク位置。デコンボリューション済み MS/MS の所在。 |
| 16 | `int` | `MS1RawSpectrumIdTop` |  |
| 17 | `int` | `MS1AccumulatedMs1RawSpectrumIdTop` | LC-IM-MS/MS 用。 |
| 18 | `int` | `MS2RawSpectrumID` | 既定 -1。代表 ID。 |
| 19 | `Dictionary<int, double>` | `MS2RawSpectrumID2CE` | 衝突エネルギーの写像。取得有無の判定に使える。 |
| 20 | `int` | `ScanID` | Key16 MS1RawSpectrumIdTop と同義。 |
| 22 | `IonMode` | `IonMode` | 0=Positive / 1=Negative / 2=Both / -1=Unknown。 |
| 24 | `List<SpectrumPeak>` | `Spectrum` | `List<SpectrumPeak>`。**実データではほぼ空**（MS/MS 本体は同名の `.dcl` にある）。 |
| 25 | `string` | `Name` |  |
| 26 | `Formula` | `Formula` | 入れ子。index0 が組成式文字列。 |
| 27 | `string` | `Ontology` |  |
| 28 | `string` | `SMILES` |  |
| 29 | `string` | `InChIKey` |  |
| 30 | `AdductIon` | `AdductType` | 入れ子。index2 が付加体文字列（`[M-H]-` など）。 |
| 31 | `double` | `CollisionCrossSection` |  |
| 33 | `Dictionary<int, List<int>>` | `MSRawID2MspIDs` | 閾値を超えた候補の ID 群。 |
| 35 | `List<int>` | `TextDbIDs` | テキスト DB 由来の候補 ID（任意）。 |
| 36 | `Dictionary<int, MsScanMatchResult>` | `MSRawID2MspBasedMatchResult` |  |
| 37 | `MsScanMatchResult` | `TextDbBasedMatchResult` |  |
| 38 | `string` | `Comment` |  |
| 39 | `IonFeatureCharacter` | `PeakCharacter` | 入れ子。付加体・電荷はここに入る。 |
| 40 | `ChromatogramPeakShape` | `PeakShape` | 入れ子（ChromatogramPeakShape）。index0=EstimatedNoise, index1=SignalToNoise。S/N はここから取る。 |
| 41 | `FeatureFilterStatus` | `FeatureFilterStatus` |  |
| 42 | `List<ChromatogramPeakFeature>` | `DriftChromFeatures` | LC-IM 用。既定 null。 |
| 43 | `double` | `Mass` | ⚠ C# 上は Obsolete（`PeakFeature` へ委譲）。**シリアライズ対象からは外れていない**。 |
| 44 | `int` | `MS1RawSpectrumIdLeft` |  |
| 45 | `int` | `MS1RawSpectrumIdRight` |  |
| 46 | `int` | `MS1AccumulatedMs1RawSpectrumIdLeft` | LC-IM-MS/MS 用。 |
| 47 | `int` | `MS1AccumulatedMs1RawSpectrumIdRight` | LC-IM-MS/MS 用。 |
| 49 | `MsScanMatchResultContainer` | `MatchResults` | 多行プロパティ。**かつて `TextDbIDWhenOrdered` に振られていた番号の再利用**（下の欠番表を参照）。 |
| 50 | `MsScanMatchResultContainer` | `matchResults` | `private` フィールド。 private バックフィールドにも Key が振られており、Key49 と同じ内容が二重に載る。 |
| 51 | `int` | `MSDecResultIdUsed` | 既定 -1。-1 のときは MasterPeakID を使う運用。 |
| 52 | `string` | `Protein` | プロテオミクス用。 |
| 53 | `int` | `ProteinGroupID` | プロテオミクス用。既定 -1。 |

## 欠番（勝手に埋めない）

| Key | 状態 |
|---:|---|
| 14 | この版のクラス定義には存在しない。 |
| 21 | この版のクラス定義には存在しない。 |
| 23 | この版のクラス定義には存在しない。 |
| 32 | `MspID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 34 | `TextDbID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 48 | `MspIDWhenOrdered`（`int`）に振られていたが上流でコメントアウト済み。 |
| 49 | かつて `TextDbIDWhenOrdered` に振られていたが、**現在は別メンバで再利用されている**（上の表を参照）。 |

## 同ファイル内の別クラス: LinkedPeakFeature

Key 番号は **0 から振り直される**。`ChromatogramPeakFeature` の番号と混ぜてはいけない。

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `LinkedPeakID` |  |
| 1 | `PeakLinkFeatureEnum` | `Character` |  |

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

`metabolomix/pai2/reader.py` の `_convert_to_peakfeature()`:
Key 3, 4, 5（ChromXs 左/トップ/右）, 6, 7, 8（高さ）, 9, 10（面積）, 11（id）, 18, 19, 22, 24,
25, 26（index0）, 27, 28, 29, 30（index2）, 31, 38, 40（index1=S/N）, 43（m/z）。

Obsolete 側の番号（3〜10, 43）に依存していることに注意。
