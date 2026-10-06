# analysis-job.v3 — 検証済みプロファイル実行のジョブ記録

MS-DIAL Console実行の永続スナップショット`analysis-job.json`の第3版。実装は
[`metabolomix/handoff/schema.py`](../../metabolomix/handoff/schema.py)。検証済みLC–MS
メタボロミクス経路（spec
[2026-09-15-validated-lcms-metabolomics-design.md](../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md)
§6, §6.1）が書く版で、v2の全フィールドをそのまま引き継ぎ、**profile snapshot**を1つ足す。

## 版の決め方（reader / writer）

| 場面 | 規則 |
|---|---|
| 読込 | `AnalysisJob.load`は`SUPPORTED_SCHEMA_VERSIONS`（`analysis-job.v1` / `analysis-job.v2` / `analysis-job.v3`）を受ける。未知のschemaは`ValueError`（`metabolomix.mztab.loading`はこれを`MZTAB_NOT_FOUND`の`DomainError`へ畳む） |
| 書出 | `AnalysisJob.save`は**ジョブの内容**で版を決める（`_schema_version_for`）。profile snapshotを持つジョブだけが`analysis-job.v3`、持たないジョブは従来どおり`analysis-job.v2` |

内容で決めるので、リピドミクスv1経路が書くジョブのバイト列は一切変わらない
（`"profile"`キーはsnapshotがあるときだけ出力される）。読み込んだv1のジョブを
保存し直すとv2になる従来の挙動もそのまま。

## v3で追加されるフィールド

| フィールド | 型 | 内容 |
|---|---|---|
| `profile` | object | `metabolomix.console.profiles.snapshot_profile()`の戻り値そのまま。dataclass側の属性名は`profile_snapshot` |
| `dependencies` | list | 依存ファイルmanifest。**`profile.dependencies`の読取専用view**（`AnalysisJob.dependencies`プロパティ） |
| `environment` | object | 実行環境manifest。**`profile.execution_environment`の読取専用view**（`AnalysisJob.environment`プロパティ） |

`dependencies` / `environment`はJSONの独立したキーではなく、snapshotの中の1箇所を
指すプロパティである。同じ一覧をジョブ側にも複製すると、どちらかだけ更新された記録が
生まれて「hashは合うのに中身が違う」状態を作れてしまうため、**正準は常にsnapshot側の
1箇所**にする。

## profile（snapshot）の主なキー

`metabolomix/console/profiles.py`の`resolve_profile_inputs` → `snapshot_profile`が組み立てる。

| キー | 内容 |
|---|---|
| `profile_id` / `profile_revision` / `profile_content_hash` | プロファイルの同一性（spec §6.1の上流指紋の構成要素） |
| `adapter` | adapter版と`dependency_keys` allowlist |
| `polarity` / `omics` / `acquisition_type` / `measure` | 実行条件の宣言 |
| `method` | method原本の`source_path`と`sha256` |
| `dependencies` | 依存ファイルごとの`method_key` / `kind` / `source_path` / `sha256` / `present` |
| `execution_environment` | 実行体hash・adapter版などの実行環境manifest |
| `raw_files` | 生データ全構成ファイルの指紋 |
| `effective_method_relative_path` / `effective_method_sha256` | `run_dir`配下へ書いた**実効メソッド**（依存の絶対パスで書き換えた実行用コピー）の位置と内容hash |
| `method_overrides` | method_keyごとの原本宣言値 → 実効値の差分 |
| `plan_identity_hash` | `run_dir`に依存しない計画のcanonical hash。出力先だけが違う2回の計画は同じ値になる |

## 完了検証での使われ方

`metabolomix.mztab.loading._verify_source`は、既存の検証（終了証跡の存在・`job_id`一致・
`termination=exited` かつ `exit_code=0`・mzTabのhash一致）を通ったあとに、v3のジョブだけ
**実効メソッドの再照合**を行う（`_verify_effective_method`）。`run_dir`配下の
`effective_method_relative_path`が消えている、または`effective_method_sha256`と
一致しない場合は`verified`を名乗らせず`legacy_unverified`にして警告を付ける。

spec §6.1が上流指紋の確認を「計画時、実行直前、**完了検証時**」の3点で求めるうちの
完了検証時にあたる。原本ではなく実効コピーを見るのは、実際にConsoleへ渡ったのが
そちらだから（spec §5.1「原本hash、書換え後hash、差分を記録する」）。

profile snapshotを持たないジョブ（v1 / v2）はこの検査の対象外で、判定は従来どおり。
依存ファイル本体（LBMは数百MBになりうる）はここでhashし直さない——読み込みのたびに
読み直すのは重く、その照合はpipelineの`validate_outputs`工程の持ち場である。

## 関連

* 実行記録: [pipeline-run.v2](pipeline-run-v2.md)
* 要求契約: [pipeline-request.v2](pipeline-request-v2.md)
* プロファイル契約: [lcms-profile.v1](lcms-profile-v1.md)
