# AnalysisFileBean — MessagePack Key 番号表

`.mddata`（`MsdialDataStorageBase` の Key0 `AnalysisFiles`）の 1 試料ぶんの MessagePack 配列に対する
Key 番号の正準表。試料クラス・測定順・バッチはここから読む。

**リーダのインデックスを変える前にこの表を見る。推測で直さない。**
関連: [AlignmentSpotProperty.md](AlignmentSpotProperty.md)

## 出典

MS-DIAL（MsdialWorkbench）の C# クラス `CompMs.MsdialCore.DataObj.AnalysisFileBean`
（上流の `src/MSDIAL5/MsdialCore/DataObj/AnalysisFile.cs`）に付いた `[Key(N)]` 属性から抽出した。

- 上流: <https://github.com/systemsomicslab/MsdialWorkbench>（LGPL-3.0）
- 照合したコミット: `45a531c`（2026-09-02）
- **この表は Key 番号とメンバ名・型の対応だけを写したもので、上流のソースコードは含まない。**
- MS-DIAL 側でクラスが変わったら、上流の同ファイルを開いて `[Key(N)]` を読み直し、
  この表を更新する。実装（`metabolomix/msdial/classes.py`）から逆算して直してはいけない。

欠番なし（Key 0〜18 が連続）。

## Key 番号表

| Key | 型 | メンバ | 備考 |
|---:|---|---|---|
| 0 | `string` | `AnalysisFilePath` | 測定ファイルのパス。 |
| 1 | `string` | `AnalysisFileName` | 試料名。 |
| 2 | `AnalysisFileType` | `AnalysisFileType` | 列挙は下記。 |
| 3 | `string` | `AnalysisFileClass` | ユーザー定義の Class ID。 |
| 4 | `int` | `AnalysisFileAnalyticalOrder` | 測定順。 |
| 5 | `int` | `AnalysisFileId` | 試料 ID（ARF の FileID と対応）。 |
| 6 | `bool` | `AnalysisFileIncluded` | 解析に含めるか。 |
| 7 | `string` | `DeconvolutionFilePath` | `*.dcl`。 |
| 8 | `List<string>` | `DeconvolutionFilePathList` | `*.dcl`。 |
| 9 | `string` | `PeakAreaBeanInformationFilePath` | `*.pai`。 |
| 10 | `string` | `RiDictionaryFilePath` | |
| 11 | `int` | `AnalysisBatch` | バッチ番号。既定 1。 |
| 12 | `double` | `ResponseVariable` | PLS 用。 |
| 13 | `double` | `DilutionFactor` | 既定 1.0。 |
| 14 | `string` | `AnalysisFileSuperClass` | |
| 15 | `RetentionTimeCorrectionBean` | `RetentionTimeCorrectionBean` | |
| 16 | `ChromatogramPeaksDataSummaryDto` | `ChromPeakFeaturesSummary` | |
| 17 | `string` | `ProteinAssembledResultFilePath` | `*.prf`。 |
| 18 | `AcquisitionType` | `AcquisitionType` | 既定 `None`。 |

## `AnalysisFileType`（Key2）

上流 `src/Common/CommonStandard/Enum/CommonEnums.cs`。

| 値 | 名前 |
|---:|---|
| 0 | `Sample` |
| 1 | `Standard` |
| 2 | `QC` |
| 3 | `Blank` |
