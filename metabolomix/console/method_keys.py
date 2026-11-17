"""MS-DIAL Console（LC-MS）がメソッドファイルから読むキーの照合表と、上書きの値の正規化。

pipeline の要求の `method_overrides`（`{キー: 値}`）を受けるときに使う。Console の
`ConfigParser`（`tests/MSDIAL5/MsdialCoreTestApp/Parser/ConfigParser.cs`）はキーを小文字にして
`switch` で照らし、**知らないキーは `default: return false` で黙って捨てる**。値も
`int.TryParse` などに通らなければ黙って既定のまま残す（例: `Minimum peak height` は int なので
`50000.0` や `1e5` は効かない）。綴りや型を誤ると、数十分かけて元と同じ解析がもう一度でき、
しかも上書きが効いたように見える。だからここで照らして、誤りは要求の段で断る。

表の出どころ: MsdialWorkbench コミット afd5f9522fa2f11990e1ad51f88b22eef00a087c の
`ConfigParser.cs`。`ReadForLcmsParameter` が呼ぶ `ReadCommonParameter` の `case` と、
メソッドファイルを別に読み直す `ReadAlignmentLightMode` / `ReadLbmAnnotatorPriority` /
`ReadDetailedAlignmentProvenance` / `ReadAnnotationCandidateExport` /
`ReadMsp(Text)AnnotatorSettingsFilePath` の `case`。型は各 `case` の本体
（`int.TryParse` → int、`float` / `double.TryParse` → float、`"true" / "false"` → bool、
`== "false"` だけを見るもの → false_only、列挙型 → enum）。`target` は値を受ける先で、
同じ先を持つキーは別名どうし。表記は `ParameterBase.ParameterToString` の綴りがあればそれ、
無ければ小文字の綴りの頭だけ大文字（Console は小文字にして読むので意味は同じ）。
上流を上げたら、この表を `ConfigParser.cs` と照らし直す（推測で足さない）。

`ReadCommonParameter` の `case "sigma window valueLower"` と
`case "replace true zero valueLowers with 1/2 of minimum peak height over all samples"` は、
上流の置換の誤りで大文字を含み、小文字にしたキーと一致しない（Console はこの 2 つを読まない）。
表には載せず、`UNREADABLE_KEYS` で理由を添えて断る。

stdlib だけに依存する（pipeline の worker も import する）。
"""
from __future__ import annotations

import difflib
import math
import re
from dataclasses import dataclass

__all__ = [
    "ALIGNMENT_REFERENCE_KEY",
    "CONSOLE_KEYS",
    "ConsoleKey",
    "MethodOverrideError",
    "RESERVED_KEYS",
    "UNREADABLE_KEYS",
    "UPSTREAM_COMMIT",
    "alias_keys",
    "normalize_overrides",
]

UPSTREAM_COMMIT = "afd5f9522fa2f11990e1ad51f88b22eef00a087c"


@dataclass(frozen=True)
class ConsoleKey:
    name: str      # 書き出すときの綴り
    kind: str      # int / float / bool / false_only / enum / text
    target: str    # Console の中で値を受ける先（同じ先を持つキーは別名どうし）


def _k(name: str, kind: str, target: str) -> ConsoleKey:
    return ConsoleKey(name, kind, target)


_ROWS: tuple[ConsoleKey, ...] = (
    _k('MS1 data type', 'enum', 'MSDataType'),
    _k('MS2 data type', 'enum', 'MS2DataType'),
    _k('Ion mode', 'enum', 'IonMode'),
    _k('Target omics', 'enum', 'TargetOmics'),
    _k('Acquisition type', 'enum', 'ProjectParam.AcquisitionType'),
    _k('Machine category', 'enum', 'ProjectParam.MachineCategory'),
    _k('Solvent type', 'enum', 'LipidQueryContainer.SolventType'),
    _k('Searched lipid class', 'text', 'LipidQueryContainer.LbmQueries'),
    _k('Msp file path', 'text', 'MspFilePath'),
    _k('Lbm file path', 'text', 'LbmFilePath'),
    _k('Text DB file path', 'text', 'TextDBFilePath'),
    _k('Isotope text DB file path', 'text', 'IsotopeTextDBFilePath'),
    _k('Compounds library file path for target detection', 'text', 'CompoundListInTargetModePath'),
    _k('Compounds library file path for RT correction', 'text', 'RetentionTimeCorrectionCommon.StandardLibrary'),
    _k('RT correction peak selection file path', 'text', 'ReferenceFileParam.RtCorrectionPeakSelectionFilePath'),
    _k('Is private version of tada', 'bool', 'IsLabPrivateVersionTada'),
    _k('Is private version', 'bool', 'IsLabPrivate'),
    _k('Retention time begin', 'float', 'RetentionTimeBegin'),
    _k('Retention time end', 'float', 'RetentionTimeEnd'),
    _k('MS1 mass range begin', 'float', 'MassRangeBegin'),
    _k('MS1 mass range end', 'float', 'MassRangeEnd'),
    _k('MS2 mass range begin', 'float', 'Ms2MassRangeBegin'),
    _k('MS2 mass range end', 'float', 'Ms2MassRangeEnd'),
    _k('Accuracy type', 'enum', 'AccuracyType'),
    _k('MS1 tolerance for centroid', 'float', 'CentroidMs1Tolerance'),
    _k('MS2 tolerance for centroid', 'float', 'CentroidMs2Tolerance'),
    _k('Smoothing method', 'enum', 'SmoothingMethod'),
    _k('Smoothing level', 'int', 'SmoothingLevel'),
    _k('Average peak width', 'int', 'AveragePeakWidth'),
    _k('Minimum peak width', 'int', 'MinimumDatapoints'),
    _k('Minimum peak height', 'int', 'MinimumAmplitude'),
    _k('Mass slice width', 'float', 'MassSliceWidth'),
    _k('Mass accuracy', 'float', 'CentroidMs1Tolerance'),
    _k('Max charge number', 'int', 'MaxChargeNumber'),
    _k('Searched adduct ions', 'text', 'SearchedAdductIons'),
    _k('Amplitude cut off', 'float', 'ChromDecBaseParam.AmplitudeCutoff'),
    _k('Relative amplitude cut off', 'float', 'ChromDecBaseParam.RelativeAmplitudeCutoff'),
    _k('Keep isotope range', 'float', 'KeptIsotopeRange'),
    _k('Exclude after precursor', 'false_only', 'RemoveAfterPrecursor'),
    _k('Keep original precursor isotopes', 'false_only', 'KeepOriginalPrecursorIsotopes'),
    _k('Target CE', 'float', 'TargetCE'),
    _k('RT tolerance for MSP-based annotation', 'float', 'MspSearchParam.RtTolerance'),
    _k('RI tolerance for MSP-based annotation', 'float', 'MspSearchParam.RiTolerance'),
    _k('Ri tolerance for identification', 'float', 'MspSearchParam.RiTolerance'),
    _k('Retention index tolerance for identification', 'float', 'MspSearchParam.RiTolerance'),
    _k('CCS tolerance for MSP-based annotation', 'float', 'MspSearchParam.CcsTolerance'),
    _k('Mass range begin for MSP-based annotation', 'float', 'MspSearchParam.MassRangeBegin'),
    _k('Mass range end for MSP-based annotation', 'float', 'MspSearchParam.MassRangeEnd'),
    _k('Relative amplitude cutoff for MSP-based annotation', 'float', 'MspSearchParam.RelativeAmpCutoff'),
    _k('Absolute amplitude cutoff for MSP-based annotation', 'float', 'MspSearchParam.AbsoluteAmpCutoff'),
    _k('Weighted dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.SquaredWeightedDotProductCutOff'),
    _k('Simple dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.SquaredSimpleDotProductCutOff'),
    _k('Reverse dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.SquaredReverseDotProductCutOff'),
    _k('Weighted dot product cutoff', 'float', 'MspSearchParam.WeightedDotProductCutOff'),
    _k('Square root of weighted dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.WeightedDotProductCutOff'),
    _k('Simple dot product cutoff', 'float', 'MspSearchParam.SimpleDotProductCutOff'),
    _k('Square root of simple dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.SimpleDotProductCutOff'),
    _k('Reverse dot product cutoff', 'float', 'MspSearchParam.ReverseDotProductCutOff'),
    _k('Square root of reverse dot product cutoff for msp-based annotation', 'float', 'MspSearchParam.ReverseDotProductCutOff'),
    _k('Matched peaks percentage cutoff', 'float', 'MspSearchParam.MatchedPeaksPercentageCutOff'),
    _k('Matched peaks percentage cutoff for msp-based annotation', 'float', 'MspSearchParam.MatchedPeaksPercentageCutOff'),
    _k('Minimum spectrum match', 'float', 'MspSearchParam.MinimumSpectrumMatch'),
    _k('Minimum spectrum match for msp-based annotation', 'float', 'MspSearchParam.MinimumSpectrumMatch'),
    _k('Total score cutoff for MSP-based annotation', 'float', 'MspSearchParam.TotalScoreCutoff'),
    _k('MS1 tolerance for MSP-based annotation', 'float', 'MspSearchParam.Ms1Tolerance'),
    _k('MS2 tolerance for MSP-based annotation', 'float', 'MspSearchParam.Ms2Tolerance'),
    _k('Use retention information for MSP-based annotation scoring', 'bool', 'MspSearchParam.IsUseTimeForAnnotationScoring'),
    _k('Use retention information for MSP-based annotation filtering', 'bool', 'MspSearchParam.IsUseTimeForAnnotationFiltering'),
    _k('Use CCS for MSP-based annotation scoring', 'bool', 'MspSearchParam.IsUseCcsForAnnotationScoring'),
    _k('Use CCS for MSP-based annotation filtering', 'bool', 'MspSearchParam.IsUseCcsForAnnotationFiltering'),
    _k('Only report top hit for MSP-based annotation', 'bool', 'OnlyReportTopHitInMspSearch'),
    _k('Execute annotation process only for alignment file', 'bool', 'IsIdentificationOnlyPerformedForAlignmentFile'),
    _k('Execute annotation process only for alignment file for msp-based annotation', 'bool', 'IsIdentificationOnlyPerformedForAlignmentFile'),
    _k('Rt tolerance for lbm-based annotation', 'float', 'LbmSearchParam.RtTolerance'),
    _k('Ri tolerance for lbm-based annotation', 'float', 'LbmSearchParam.RiTolerance'),
    _k('Ccs tolerance for lbm-based annotation', 'float', 'LbmSearchParam.CcsTolerance'),
    _k('Mass range begin for lbm-based annotation', 'float', 'LbmSearchParam.MassRangeBegin'),
    _k('Mass range end for lbm-based annotation', 'float', 'LbmSearchParam.MassRangeEnd'),
    _k('Relative amplitude cutoff for lbm-based annotation', 'float', 'LbmSearchParam.RelativeAmpCutoff'),
    _k('Absolute amplitude cutoff for lbm-based annotation', 'float', 'LbmSearchParam.AbsoluteAmpCutoff'),
    _k('Weighted dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.SquaredWeightedDotProductCutOff'),
    _k('Simple dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.SquaredSimpleDotProductCutOff'),
    _k('Reverse dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.SquaredReverseDotProductCutOff'),
    _k('Square root of weighted dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.WeightedDotProductCutOff'),
    _k('Square root of simple dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.SimpleDotProductCutOff'),
    _k('Square root of reverse dot product cutoff for lbm-based annotation', 'float', 'LbmSearchParam.ReverseDotProductCutOff'),
    _k('Matched peaks percentage cutoff for lbm-based annotation', 'float', 'LbmSearchParam.MatchedPeaksPercentageCutOff'),
    _k('Minimum spectrum match for lbm-based annotation', 'float', 'LbmSearchParam.MinimumSpectrumMatch'),
    _k('Total score cutoff for lbm-based annotation', 'float', 'LbmSearchParam.TotalScoreCutoff'),
    _k('Ms1 tolerance for lbm-based annotation', 'float', 'LbmSearchParam.Ms1Tolerance'),
    _k('Ms2 tolerance for lbm-based annotation', 'float', 'LbmSearchParam.Ms2Tolerance'),
    _k('Use retention information for lbm-based annotation scoring', 'bool', 'LbmSearchParam.IsUseTimeForAnnotationScoring'),
    _k('Use retention information for lbm-based annotation filtering', 'bool', 'LbmSearchParam.IsUseTimeForAnnotationFiltering'),
    _k('Use ccs for lbm-based annotation scoring', 'bool', 'LbmSearchParam.IsUseCcsForAnnotationScoring'),
    _k('Use ccs for lbm-based annotation filtering', 'bool', 'MspSearchParam.IsUseCcsForAnnotationFiltering'),
    _k('Execute annotation process only for alignment file for lbm-based annotation', 'bool', 'IsIdentificationOnlyPerformedForAlignmentFile'),
    _k('RT tolerance for Text-based annotation', 'float', 'TextDbSearchParam.RtTolerance'),
    _k('RI tolerance for Text-based annotation', 'float', 'TextDbSearchParam.RiTolerance'),
    _k('CCS tolerance for Text-based annotation', 'float', 'TextDbSearchParam.CcsTolerance'),
    _k('Total score cutoff for Text-based annotation', 'float', 'TextDbSearchParam.TotalScoreCutoff'),
    _k('Accurate ms1 tolerance for Text-based annotation', 'float', 'TextDbSearchParam.Ms1Tolerance'),
    _k('Use retention information for Text-based annotation scoring', 'bool', 'TextDbSearchParam.IsUseTimeForAnnotationScoring'),
    _k('Use retention information for Text-based annotation filtering', 'bool', 'TextDbSearchParam.IsUseTimeForAnnotationFiltering'),
    _k('Use CCS for Text-based annotation scoring', 'bool', 'TextDbSearchParam.IsUseCcsForAnnotationScoring'),
    _k('Use CCS for Text-based annotation filtering', 'bool', 'TextDbSearchParam.IsUseCcsForAnnotationFiltering'),
    _k('Only report top hit for Text-based annotation', 'bool', 'OnlyReportTopHitInTextDBSearch'),
    _k('Alignment reference file ID', 'int', 'AlignmentReferenceFileID'),
    _k('Retention time tolerance for alignment', 'float', 'RetentionTimeAlignmentTolerance'),
    _k('Retention time factor for alignment', 'float', 'RetentionTimeAlignmentFactor'),
    _k('Spectrum similarity tolerance for alignment', 'float', 'SpectrumSimilarityAlignmentTolerance'),
    _k('Spectrum similarity factor for alignment', 'float', 'SpectrumSimilarityAlignmentFactor'),
    _k('MS1 tolerance for alignment', 'float', 'Ms1AlignmentTolerance'),
    _k('MS1 factor for alignment', 'float', 'Ms1AlignmentFactor'),
    _k('Force insert peaks in gap filling', 'bool', 'IsForceInsertForGapFilling'),
    _k('Together with alignment', 'bool', 'TogetherWithAlignment'),
    _k('Peak count filter', 'float', 'PeakCountFilter'),
    _k('N percent detected in one group', 'float', 'NPercentDetectedInOneGroup'),
    _k('Remove feature based on peak height fold-change', 'bool', 'IsRemoveFeatureBasedOnBlankPeakHeightFoldChange'),
    _k('Blank filtering', 'enum', 'BlankFiltering'),
    _k('Sample max / blank average', 'float', 'SampleMaxOverBlankAverage'),
    _k('Sample average / blank average', 'float', 'SampleAverageOverBlankAverage'),
    _k('Keep reference matched metabolites', 'bool', 'IsKeepRefMatchedMetaboliteFeatures'),
    _k('Keep suggested metabolites', 'bool', 'IsKeepSuggestedMetaboliteFeatures'),
    _k('Keep removable features and assigned tag for checking', 'bool', 'IsKeepRemovableFeaturesAndAssignedTagForChecking'),
    _k('Execute RT correction', 'bool', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.ExcuteRtCorrection'),
    _k('RT correction with smoothing for RT diff', 'bool', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.doSmoothing'),
    _k('User setting intercept', 'float', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.UserSettingIntercept'),
    _k('RT diff calc method', 'enum', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.RtDiffCalcMethod'),
    _k('Interpolation method', 'enum', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.InterpolationMethod'),
    _k('Extrapolation method (begin)', 'enum', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.ExtrapolationMethodBegin'),
    _k('Extrapolation method (end)', 'enum', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.ExtrapolationMethodEnd'),
    _k('RT correction peak selection mode', 'enum', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.PeakSelectionMode'),
    _k('RT correction peak selection RT weight', 'float', 'RetentionTimeCorrectionCommon.RetentionTimeCorrectionParam.PeakSelectionRtWeight'),
    _k('Tracking isotope label', 'bool', 'TrackingIsotopeLabels'),
    _k('Set fully labeled reference file', 'bool', 'SetFullyLabeledReferenceFile'),
    _k('Non labeled reference ID', 'int', 'NonLabeledReferenceID'),
    _k('Fully labeled reference ID', 'int', 'FullyLabeledReferenceID'),
    _k('Number of threads', 'int', 'NumThreads'),
    _k('Isotope tracking dictionary ID', 'int', 'IsotopeTrackingDictionary.SelectedID'),
    _k('CorrDec execute', 'false_only', 'CorrDecParam.CanExcute'),
    _k('CorrDec MS2 tolerance', 'float', 'CorrDecParam.MS2Tolerance'),
    _k('CorrDec minimum MS2 peak height', 'int', 'CorrDecParam.MinMS2Intensity'),
    _k('CorrDec minimum number of detected samples', 'int', 'CorrDecParam.MinNumberOfSample'),
    _k('CorrDec exclude highly correlated spots', 'float', 'CorrDecParam.MinCorr_MS1'),
    _k('CorrDec minimum correlation coefficient (MS2)', 'float', 'CorrDecParam.MinCorr_MS2'),
    _k('CorrDec margin 1 (target precursor)', 'float', 'CorrDecParam.CorrDiff_MS1'),
    _k('CorrDec margin 2 (coeluted precursor)', 'float', 'CorrDecParam.CorrDiff_MS2'),
    _k('CorrDec minimum detected rate', 'float', 'CorrDecParam.MinDetectedPercentToVisualize'),
    _k('CorrDec minimum MS2 relative intensity', 'float', 'CorrDecParam.MinMS2RelativeIntensity'),
    _k('CorrDec remove peaks larger than precursor', 'bool', 'CorrDecParam.CorrDecRemoveAfterPrecursor'),
    _k('Alignment light mode', 'bool', 'ReadAlignmentLightMode'),
    _k('Alignment light', 'bool', 'ReadAlignmentLightMode'),
    _k('Console alignment light mode', 'bool', 'ReadAlignmentLightMode'),
    _k('Lbm annotator priority', 'int', 'ReadLbmAnnotatorPriority'),
    _k('Lbm annotation priority', 'int', 'ReadLbmAnnotatorPriority'),
    _k('Detailed alignment provenance', 'bool', 'ReadDetailedAlignmentProvenance'),
    _k('Export detailed alignment provenance', 'bool', 'ReadDetailedAlignmentProvenance'),
    _k('Annotation candidates', 'bool', 'ReadAnnotationCandidateExport'),
    _k('Export annotation candidates', 'bool', 'ReadAnnotationCandidateExport'),
    _k('Msp annotator settings file path', 'text', 'ReadMspAnnotatorSettingsFilePath'),
    _k('Msp annotation settings file path', 'text', 'ReadMspAnnotatorSettingsFilePath'),
    _k('Msp search settings file path', 'text', 'ReadMspAnnotatorSettingsFilePath'),
    _k('Text annotator settings file path', 'text', 'ReadTextAnnotatorSettingsFilePath'),
    _k('Text library annotator settings file path', 'text', 'ReadTextAnnotatorSettingsFilePath'),
    _k('Text db annotator settings file path', 'text', 'ReadTextAnnotatorSettingsFilePath'),
    _k('Text annotation settings file path', 'text', 'ReadTextAnnotatorSettingsFilePath'),
)

#: 小文字のキー → 表の行。
CONSOLE_KEYS: dict[str, ConsoleKey] = {row.name.lower(): row for row in _ROWS}

#: 列挙型のキーが受ける値（上流の enum の綴り。照合は大文字小文字を問わず、書き出しはこの綴り）。
#: `Machine category` と `Solvent type` は Console が値を大文字小文字どおりに照らすので、
#: 綴りをそろえて書くことに意味がある。`Blank filtering` は Console が
#: `samplemaxoverblankave` しか受けない。
_ENUM_VALUES: dict[str, tuple[str, ...]] = {
    "ms1 data type": ("Centroid", "Profile"),
    "ms2 data type": ("Centroid", "Profile"),
    "ion mode": ("Positive", "Negative"),
    "target omics": ("Metabolomics", "Lipidomics"),
    "acquisition type": ("DDA", "SWATH", "AIF"),
    "machine category": ("GCMS", "LCMS", "IMMS", "LCIMMS", "IFMS", "IIMMS", "IDIMS"),
    "solvent type": ("CH3COONH4", "HCOONH4", "NH4HCO3"),
    "accuracy type": ("IsNominal", "IsAccurate"),
    "smoothing method": ("SimpleMovingAverage", "LinearWeightedMovingAverage", "SavitzkyGolayFilter",
                         "BinomialFilter", "LowessFilter", "LoessFilter",
                         "TimeBasedLinearWeightedMovingAverage"),
    "blank filtering": ("SampleMaxOverBlankAve",),
    "rt diff calc method": ("SampleMinusSampleAverage", "SampleMinusReference"),
    "interpolation method": ("Linear",),
    "extrapolation method (begin)": ("UserSetting", "FirstPoint", "LinearExtrapolation"),
    "extrapolation method (end)": ("LastPoint", "LinearExtrapolation"),
    "rt correction peak selection mode": ("HighestIntensity", "ClosestToReferenceRt", "Weighted"),
}

_LIBRARY_FIELDS = ("method_file / lbm_file / msp_file / library_mode（ライブラリと外部ファイルの"
                   "パスは計画が解決して書く）")

#: 要求のほかの欄が受け持つキー（小文字 → 代わりに使う欄）。上書きで変えると計画の解決
#: （極性・ライブラリ・パスの絶対化）と食い違う。
RESERVED_KEYS: dict[str, str] = {
    "ion mode": "polarity（極性は単位の polarity で決める）",
    **{key: _LIBRARY_FIELDS for key in CONSOLE_KEYS if key.endswith("file path")},
}

#: `alignment_reference="method"` のときだけ上書きできるキー（既定の "median_injection" では
#: 計画が注入順の中央の試料を書く）。
ALIGNMENT_REFERENCE_KEY = "alignment reference file id"

#: Console が読まないキー（小文字 → 理由）。
UNREADABLE_KEYS: dict[str, str] = {
    "sigma window value": (
        "Console はこのキーを読みません（上流 ConfigParser の case が 'sigma window valueLower' と"
        "大文字を含み、小文字にしたキーと一致しない）。既定値のまま走ります。"),
    "replace true zero values with 1/2 of minimum peak height over all samples": (
        "Console はこのキーを読みません（上流 ConfigParser の case が '... valueLowers ...' と"
        "大文字を含み、小文字にしたキーと一致しない）。既定値のまま走ります。"),
}

_INT_TEXT = re.compile(r"^[+-]?\d+$")


class MethodOverrideError(ValueError):
    """上書きの誤り。`details` は要求のエラーの details に載せる。"""

    def __init__(self, message: str, **details):
        super().__init__(message)
        self.message = message
        self.details = details


def _canonical_key(raw_key: str) -> str:
    """空白をつぶした小文字（Console の照合と同じく大文字小文字を問わない）。"""
    return " ".join(raw_key.split()).lower()


def _suggest(key: str) -> list[str]:
    """綴りの近いキー（最大 3 件）。"""
    names = difflib.get_close_matches(key, list(CONSOLE_KEYS), n=3, cutoff=0.6)
    return [CONSOLE_KEYS[name].name for name in names]


def _format_float(number: float) -> str:
    if number.is_integer() and abs(number) < 1e15:
        return str(int(number))
    return repr(number)


def _normalize_value(row: ConsoleKey, value) -> str:
    """値を Console が読める文字列にする。読めない値は MethodOverrideError。"""
    key, kind = row.name, row.kind
    if value is None or isinstance(value, (dict, list, tuple)):
        raise MethodOverrideError(f"{key} の値は文字列・数値・真偽値のどれかです: {value!r}",
                                  key=key, value=value)
    if kind in ("bool", "false_only"):
        if isinstance(value, bool):
            return "True" if value else "False"
        if isinstance(value, str) and value.strip().lower() in ("true", "false"):
            return "True" if value.strip().lower() == "true" else "False"
        raise MethodOverrideError(f"{key} は true / false です: {value!r}",
                                  key=key, value=value, kind=kind)
    if isinstance(value, bool):
        raise MethodOverrideError(f"{key} に真偽値は書けません（{kind}）: {value!r}",
                                  key=key, value=value, kind=kind)
    if kind == "int":
        if isinstance(value, int):
            return str(value)
        if isinstance(value, float) and math.isfinite(value) and value.is_integer():
            return str(int(value))
        if isinstance(value, str) and _INT_TEXT.match(value.strip()):
            return str(int(value.strip()))
        raise MethodOverrideError(
            f"{key} は整数です（Console は int.TryParse で読み、通らない値は黙って捨てる。"
            f"'1e5' や '50000.5' も通らない）: {value!r}", key=key, value=value, kind=kind)
    if kind == "float":
        number = math.nan
        if isinstance(value, (int, float)):
            number = float(value)
        elif isinstance(value, str):
            try:
                number = float(value.strip())
            except ValueError:
                pass
        if not math.isfinite(number):
            raise MethodOverrideError(f"{key} は有限の数です: {value!r}",
                                      key=key, value=value, kind=kind)
        return _format_float(number)
    if kind == "enum":
        options = _ENUM_VALUES[key.lower()]
        if isinstance(value, str):
            for option in options:
                if option.lower() == value.strip().lower():
                    return option
        raise MethodOverrideError(f"{key} の値は {list(options)} のどれかです: {value!r}",
                                  key=key, value=value, kind=kind, options=list(options))
    # text
    if isinstance(value, float):
        return _format_float(value)
    if isinstance(value, int):
        return str(value)
    text = value.strip()
    if not text:
        raise MethodOverrideError(f"{key} の値が空です。", key=key, value=value, kind=kind)
    return text


def normalize_overrides(overrides, *, alignment_reference_from_method: bool) -> dict[str, str]:
    """`method_overrides` を検査し、`{表記の綴り: 文字列の値}`（表記の昇順）にして返す。

    断るもの（MethodOverrideError）: 辞書でない・Console の読まないキー（近い綴りを候補に添える）・
    ほかの欄が受け持つキー・同じキーを 2 度（大文字小文字違いも）・同じ受け先を持つ別名を 2 つ・
    型に合わない値・ASCII でない／改行を含む値（Console はメソッドを ASCII で 1 行ずつ読む）。
    `Alignment reference file ID` は `alignment_reference_from_method` のときだけ受ける。
    """
    if not isinstance(overrides, dict):
        raise MethodOverrideError(
            f"method_overrides は {{キー: 値}} のオブジェクトです: {overrides!r}", value=overrides)
    out: dict[str, str] = {}
    seen_targets: dict[str, str] = {}
    for raw_key, value in overrides.items():
        if not isinstance(raw_key, str) or not raw_key.strip():
            raise MethodOverrideError(f"method_overrides のキーが空です: {raw_key!r}", key=raw_key)
        key = _canonical_key(raw_key)
        if key in UNREADABLE_KEYS:
            raise MethodOverrideError(UNREADABLE_KEYS[key], key=raw_key)
        if key in RESERVED_KEYS:
            raise MethodOverrideError(
                f"{raw_key} は method_overrides では変えられません。代わりに {RESERVED_KEYS[key]}。",
                key=raw_key, use=RESERVED_KEYS[key])
        if key == ALIGNMENT_REFERENCE_KEY and not alignment_reference_from_method:
            raise MethodOverrideError(
                f"{raw_key} は alignment_reference='method' のときだけ上書きできます（既定の "
                "'median_injection' では計画が注入順の中央の試料を書く）。",
                key=raw_key, use="alignment_reference")
        row = CONSOLE_KEYS.get(key)
        if row is None:
            raise MethodOverrideError(
                f"MS-DIAL Console が読まないキーです（Console は知らないキーを黙って捨てる）: {raw_key!r}",
                key=raw_key, suggestions=_suggest(key), upstream_commit=UPSTREAM_COMMIT)
        if row.name in out:
            raise MethodOverrideError(f"同じキーが 2 度あります: {row.name}", key=raw_key)
        if row.target in seen_targets:
            raise MethodOverrideError(
                f"{raw_key} と {seen_targets[row.target]} は Console の中で同じ値を決める別名です。"
                "どちらか 1 つにしてください。", key=raw_key, other=seen_targets[row.target])
        text = _normalize_value(row, value)
        if not text.isascii() or "\n" in text or "\r" in text:
            raise MethodOverrideError(
                f"{row.name} の値に ASCII 以外の文字か改行があります（Console はメソッドを ASCII で"
                "1 行ずつ読む）。", key=row.name)
        seen_targets[row.target] = row.name
        out[row.name] = text
    return dict(sorted(out.items()))


def alias_keys(name: str) -> tuple[str, ...]:
    """`name` と同じ受け先を持つ、ほかのキーの表記（`name` 自身は含まない）。"""
    row = CONSOLE_KEYS.get(_canonical_key(name))
    if row is None:
        return ()
    return tuple(other.name for other in _ROWS if other.target == row.target and other is not row)
