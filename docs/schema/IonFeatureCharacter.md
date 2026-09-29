# IonFeatureCharacter（`AlignmentSpotProperty` Key 10）

上流: `MsdialWorkbench` `afd5f9522`（2026-09-08）
`src/MSDIAL5/MsdialCore/DataObj/IonFeatureCharacter.cs`、
`LinkedPeakFeature` は `src/MSDIAL5/MsdialCore/DataObj/ChromatogramPeakFeature.cs`、
`PeakLinkFeatureEnum` は `src/Common/CommonStandard/Enum/CommonEnums.cs`。

| Key | 型 | 名前 | 備考 |
|---|---|---|---|
| 0 | AdductIon | AdductType | 本システムは `.arf2` Key 54 の AdductType を使う |
| 1 | AdductIon | AdductTypeByAmalgamationProgram | |
| 2 | int | Charge | 既定 1 |
| 3 | List<LinkedPeakFeature> | PeakLinks | 下表 |
| 4 | int | IsotopeWeightNumber | 0 = 単同位体。-1 既定 |
| 5 | int | IsotopeParentPeakID | |
| 6 | int | PeakGroupID | |
| 7 | bool | IsLinked | |
| 8 | int | AdductParent | |

欠番なし（Key 0〜8 が連続）。

LinkedPeakFeature: Key 0 `LinkedPeakID`（int）、Key 1 `Character`（PeakLinkFeatureEnum）。

PeakLinkFeatureEnum: 0 SameFeature, 1 Isotope, 2 Adduct, 3 ChromSimilar, 4 FoundInUpperMsMs, 5 CorrelSimilar。

## 実データでの確認（2026-09-29、kidney neg/pos）

- IsotopeWeightNumber はアラインメントの全スポットで 0（neg 2207 / pos 3703 スポット。
  同位体は別スポットとして残らない）。
- LinkedPeakID は **同じ `.arf2` の MasterAlignmentID（Key 0）** を指す。リンク先が存在しない
  ものは 0 件（neg 5970 本 / pos 6360 本）。リンク元とリンク先の RT の差は、|ΔRT| の中央値が
  neg 0.004 分・pos 0.007 分、95 パーセンタイルが neg 0.048 分・pos 0.053 分。
  スポットの並びは MasterAlignmentID = 配列の添字だった（両極性で確認）。
- リンクの種類の内訳（neg / pos）: CorrelSimilar 4738 / 4896、FoundInUpperMsMs 998 / 1206、
  ChromSimilar 206 / 218、Adduct 28 / 40。SameFeature・Isotope は 0 本。
- **リンクは全て相互**（A→B があれば必ず同じ種類の B→A がある。全種類で 100%）。
  FoundInUpperMsMs の向き: リンク先の m/z がリンク元より高いものと低いものが
  ちょうど同数（neg 499 / 499、pos 603 / 603）。つまり**片方向の「断片側が上位 precursor へ張る」
  情報ではなく、無向の対**として書かれている。どちらが上位イオンかはリンクだけからは
  分からないので、m/z の大小で自分で決める（リンク先の m/z が高い側を上位とみなす）。
  リンク先との |Δm/z| は neg で 7.7〜1082、中央値 197。
