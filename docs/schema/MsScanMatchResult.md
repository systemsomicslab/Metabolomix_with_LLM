# MsScanMatchResult / MsScanMatchResultContainer の Key 番号

- 出典: `C:\Users\yuu18\source\repos\MsdialWorkbench`（`master`、`afd5f9522`、2026-09-08）
  - `src/Common/CommonStandard/DataObj/Result/MsScanMatchResult.cs`
  - `src/MSDIAL5/MsdialCore/DataObj/MsScanMatchResultContainer.cs`
- `.arf2` の `AlignmentSpotProperty` Key 56（`MatchResults`）がこのコンテナ。

## MsScanMatchResultContainer

| Key | 型 | 名前 | 備考 |
|---|---|---|---|
| 0 | `List<MsScanMatchResult>` | `MatchResults` | 候補の一覧。空なら上流は `UnknownResult` を返す |
| 1 | `Dictionary<int, MsScanMatchResult>` | `MSRawID2MspBasedMatchResult` | |
| 2 | `List<MsScanMatchResult>` | `TextDbBasedMatchResults` | |

代表（GUI が表示する 1 件）はシリアライズされない計算値 `Representative`:
非 decoy かつ `Source != Unknown` の候補の
`argmax(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)`
（`ResultOrder`）。`IsManuallyModified = (Source & Manual) != 0`。

## MsScanMatchResult

| Key | 型 | 名前 |
|---|---|---|
| 0 | `string` | `Name` |
| 1 | `string` | `InChIKey` |
| 2 | `float` | `TotalScore` |
| 3 | `float` | `SquaredWeightedDotProduct`（MS/MS 無しは -1） |
| 4 | `float` | `SquaredSimpleDotProduct` |
| 5 | `float` | `SquaredReverseDotProduct` |
| 6 | `float` | `MatchedPeaksCount`（MS/MS 無しは -1） |
| 7 | `float` | `MatchedPeaksPercentage` |
| 8 | `float` | `EssentialFragmentMatchedScore` |
| 9 | `float` | `RtSimilarity` |
| 10 | `float` | `RiSimilarity` |
| 11 | `float` | `CcsSimilarity` |
| 12 | `float` | `IsotopeSimilarity` |
| 13 | `float` | `AcurateMassSimilarity` |
| 14 | `int` | `LibraryID`（参照の ScanID） |
| 15 | `bool` | `IsPrecursorMzMatch` |
| 16 | `bool` | `IsSpectrumMatch` |
| 17 | `bool` | `IsRtMatch` |
| 18 | `bool` | `IsCcsMatch` |
| 19 | `bool` | `IsLipidClassMatch` |
| 20 | `bool` | `IsLipidChainsMatch` |
| 21 | `bool` | `IsLipidPositionMatch` |
| 22 | `bool` | `IsOtherLipidMatch` |
| 23 | `bool` | `IsRiMatch` |
| 24 | `int` | `LibraryIDWhenOrdered` |
| 25 | — | 欠番 |
| 26 | `SourceType`（byte のビット集合） | `Source`: None=0, Unknown=1, FastaDB=2, MspDB=4, TextDB=16, GeneratedLipid=32, Manual=64 |
| 27 | `string` | `AnnotatorID`（`.dbs` では `<エントリ名>_<n>`） |
| 28 | `int` | `SpectrumID` |
| 29 | `float` | `AndromedaScore` |
| 30 | `bool` | `IsDecoy` |
| 31 | `int` | `Priority` |
| 32 | `float` | `PEPScore` |
| 33 | `bool` | `IsReferenceMatched` |
| 34 | `bool` | `IsAnnotationSuggested` |
| 35 | `bool` | `IsLipidDoubleBondPositionMatch` |
| 36 | `double` | `CollisionEnergy` |
| 37 | `float` | `EnhancedDotProduct` |
| 38 | `float` | `SpectralEntropy` |

このリポジトリの `has_msms`（`metabolomix/arf2/match_results.py`）は `SquaredWeightedDotProduct >= 0`
（Key 3 が番兵 -1 でない）で判定しており、上流の `IsSpectrumComparisonPerformed` を簡略化したもの。

## 注釈名の接頭辞（スポット Key 12）

上流 `StandardAnnotationProcess.SetRepresentativeProperty` / `DataAccess.SetMoleculeMsPropertyAsSuggested`:
`IsReferenceMatched` → 接頭辞なし。そうでなく `IsAnnotationSuggested` → MS/MS なしなら `no MS2: `、
ありなら `low score: `。`Unsettled: ` は GUI の手動操作（`SetMoleculeMsPropertyAsUnsettled`）だけ。
