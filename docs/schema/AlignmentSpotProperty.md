# AlignmentSpotProperty — MessagePack Key 番号表

`.arf2`（アライン後のスポット代表）1 行ぶんの MessagePack 配列に対する Key 番号の正準表。
`.arf` のスポット層も同じクラスで、サンプル別ピークは Key8 `AlignedPeakProperties` に
`AlignmentChromPeakFeature`（→ [AlignmentChromPeakFeature.md](AlignmentChromPeakFeature.md)）
の配列として入る。

**リーダのインデックスを変える前にこの表を見る。推測で直さない。**
関連: [../output_format/arf2.md](../output_format/arf2.md)（フィールドの意味）

## 出典

MS-DIAL（MsdialWorkbench）の C# クラス `CompMs.MsdialCore.DataObj.AlignmentSpotProperty`
（上流の `src/MSDIAL5/MsdialCore/DataObj/AlignmentSpotProperty.cs`）に付いた `[Key(N)]` 属性から抽出した。

- 上流: <https://github.com/systemsomicslab/MsdialWorkbench>（LGPL-3.0）
- 照合したコミット: `45a531c`（2026-09-02）
- **この表は Key 番号とメンバ名・型の対応だけを写したもので、上流のソースコードは含まない。**
- MS-DIAL 側でクラスが変わったら、上流の同ファイルを開いて `[Key(N)]` を読み直し、
  この表を更新する。実装（`metabolomix/*/reader.py`）から逆算して直してはいけない。

表は **Key 番号の昇順**で並べている（上流のソース上の宣言順とは一致しない）。

## Key 番号表

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `MasterAlignmentID` | 全次元通しの連番。 |
| 1 | `int` | `AlignmentID` | 同一次元内（RT×m/z など）の連番。`.arf2` の行 ID。 |
| 2 | `int` | `ParentAlignmentID` | LC-IM-MS/MS で親スポットを指す。 |
| 3 | `int` | `RepresentativeFileID` |  |
| 4 | `ChromXs` | `TimesCenter` | 入れ子（ChromXs）。スポットの RT はここから取る。 |
| 5 | `double` | `MassCenter` |  |
| 6 | `double` | `QuantMass` | GC-MS の定量質量。 |
| 7 | `int` | `InternalStandardAlignmentID` | 参照先の MasterAlignmentID が入る。 |
| 8 | `List<AlignmentChromPeakFeature>` | `AlignedPeakProperties` | **サンプル別ピークの配列**。要素は `AlignmentChromPeakFeature`（別表）。`.arf` の実体はここ。 |
| 9 | `List<AlignmentSpotProperty>` | `AlignmentDriftSpotFeatures` | LC-IM 用。既定 null。 |
| 10 | `IonFeatureCharacter` | `PeakCharacter` |  |
| 11 | `IonMode` | `IonMode` | 0=Positive / 1=Negative / 2=Both / -1=Unknown。 |
| 12 | `string` | `Name` | 同定名。未同定は空または Unknown 相当。 |
| 13 | `Formula` | `Formula` | 入れ子。index0 が組成式文字列。 |
| 14 | `string` | `Ontology` | 脂質クラス（PC, TG, …）。 |
| 15 | `string` | `SMILES` |  |
| 16 | `string` | `InChIKey` |  |
| 17 | `double` | `CollisionCrossSection` |  |
| 19 | `Dictionary<int, List<int>>` | `MSRawID2MspIDs` | 閾値を超えた候補の ID 群（MS raw ID → MSP ID リスト）。 |
| 21 | `List<int>` | `TextDbIDs` | テキスト DB 由来の候補 ID（任意）。 |
| 22 | `int` | `IsotopeTrackTextDbID` | 代表となる text ID。 |
| 23 | `Dictionary<int, MsScanMatchResult>` | `MSRawID2MspBasedMatchResult` |  |
| 24 | `MsScanMatchResult` | `TextDbBasedMatchResult` |  |
| 25 | `string` | `Comment` |  |
| 26 | `string` | `AnnotationCode` |  |
| 27 | `string` | `AnnotationCodeCorrDec` |  |
| 28 | `List<int>` | `CorrDecLibraryIDs` | AIF プロジェクト用。 |
| 29 | `float` | `AnovaPvalue` | プロパティでなくフィールド。MS-DIAL 側が計算した値で、本サーバの差次的解析とは別物。 |
| 30 | `float` | `FoldChange` | 同上（フィールド）。MS-DIAL 側の計算値。 |
| 31 | `float` | `HeightAverage` |  |
| 32 | `float` | `HeightMin` |  |
| 33 | `float` | `HeightMax` |  |
| 34 | `float` | `PeakWidthAverage` |  |
| 35 | `float` | `SignalToNoiseAve` |  |
| 36 | `float` | `SignalToNoiseMax` |  |
| 37 | `float` | `SignalToNoiseMin` |  |
| 38 | `float` | `EstimatedNoiseAve` |  |
| 39 | `float` | `EstimatedNoiseMax` |  |
| 40 | `float` | `EstimatedNoiseMin` |  |
| 41 | `ChromXs` | `TimesMin` | 入れ子（ChromXs）。 |
| 42 | `ChromXs` | `TimesMax` | 入れ子（ChromXs）。 |
| 43 | `double` | `MassMin` |  |
| 44 | `double` | `MassMax` |  |
| 45 | `FeatureFilterStatus` | `FeatureFilterStatus` |  |
| 46 | `bool` | `IsManuallyModifiedForQuant` |  |
| 47 | `IonAbundanceUnit` | `IonAbundanceUnit` |  |
| 49 | `float` | `FillParcentage` | **綴りは上流のまま**（Percentage ではなく Parcentage）。検出されたサンプルの割合。 |
| 50 | `float` | `RelativeAmplitudeValue` |  |
| 51 | `float` | `MonoIsotopicPercentage` |  |
| 52 | `List<AlignmentSpotVariableCorrelation>` | `AlignmentSpotVariableCorrelations` | 要素は `AlignmentSpotVariableCorrelation`（別表）。 |
| 53 | `List<IsotopicPeak>` | `IsotopicPeaks` | 代表ファイル由来の同位体ピーク。 |
| 54 | `AdductIon` | `AdductType` | 入れ子。index2 が付加体文字列（`[M-H]-` など）。setter のみ Obsolete。 |
| 56 | `MsScanMatchResultContainer` | `MatchResults` | 多行プロパティ。private バックフィールドには Key が振られていない。 |
| 58 | `bool` | `IsBlankFilteredByPostCurator` | ポストキュレーションの判定結果。 |
| 59 | `int` | `MSDecResultIdUsed` | 既定 -1。-1 のときは MasterAlignmentID 側を使う運用。 |
| 60 | `string` | `Protein` | プロテオミクス用。 |
| 61 | `int` | `ProteinGroupID` | プロテオミクス用。既定 -1。 |
| 62 | `bool` | `IsMzFilteredByPostCurator` | ポストキュレーションの判定結果。 |
| 63 | `bool` | `IsBlankGhostFilteredByPostCurator` | ポストキュレーションの判定結果。 |
| 64 | `bool` | `IsRsdFilteredByPostCurator` | ポストキュレーションの判定結果。 |
| 65 | `bool` | `IsRmdFilteredByPostCurator` | ポストキュレーションの判定結果。 |

## 欠番（勝手に埋めない）

| Key | 状態 |
|---:|---|
| 18 | `MspID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 20 | `TextDbID`（`int`）に振られていたが上流でコメントアウト済み。 |
| 48 | この版のクラス定義には存在しない。 |
| 55 | この版のクラス定義には存在しない。 |
| 57 | この版のクラス定義には存在しない。 |

## 同ファイル内の別クラス: AlignmentSpotVariableCorrelation

Key 番号は **0 から振り直される**。`AlignmentSpotProperty` の番号と混ぜてはいけない。

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `int` | `CorrelateAlignmentID` |  |
| 1 | `float` | `CorrelationScore` |  |

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

`metabolomix/arf2/reader.py` の `extract_arf2_data()`:
Key 0, 1, 4（RT）, 5, 11, 12, 13（index0）, 14, 15, 16, 31, 32, 33, 34, 35, 36, 37, 43, 44, 49, 51,
54（index2）。

`metabolomix/arf/reader.py` はスポット層から Key8 `AlignedPeakProperties` を取り出し、
その各要素を `AlignmentChromPeakFeature` として読む。
