"""Console 実行層の MCP ツール: console_plan, console_run, console_status, job_list。

spec §8 参照。
"""
from __future__ import annotations

from pathlib import Path

from mcp.types import ToolAnnotations

from lipidmix.core import session_state
from lipidmix.core.mcp_core import mcp
from lipidmix.core.mcp_errors import console_error, mztab_error
from lipidmix.core.serialization import json_payload

__all__ = ["console_plan", "console_prepare_input", "console_method_template",
           "console_method_candidates", "console_run", "console_status",
           "console_cleanup", "job_list"]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False),
          structured_output=False)
def console_plan(
    dataset_root: str,
    method_file: str | None = None,
    polarity: str = "positive",
    measure: str = "peak_height",
    omics: str = "lipidomics",
    save_project: bool = True,
    timeout_s: int = 21600,
    lbm_file: str | None = None,
    keep_extension: str | None = None,
) -> str:
    """MS-DIAL Console の実行計画を作成し、analysis-job.json を生成します。

    Parameters
    ----------
    dataset_root:
        生データフォルダのパス。MS-DIAL が読む 12 形式
        （.abf / .ibf / .cdf / .mzml / .wiff / .raw / .d / .wiff2 /
        .qgd / .lcd / .lrp / .imzml）が入っているフォルダを指します。
        Agilent・Bruker の .d と Waters の .raw は**フォルダ**が 1 検体です。
        リポジトリ外のパスを指定してください。
    method_file:
        MS-DIAL Console のパラメータファイル（ASCII テキスト。`key: value` 形式）。
        **省略できます** — 省略時は dataset_root 直下・兄弟フォルダ・過去 run
        （`runs/*/analysis-job.json`）まで探し、`Ion mode` が polarity と
        一致する最新のものを使います。一致が無く別極性の候補があれば
        METHOD_FILE_CHOICE_REQUIRED で止まります（候補の列挙だけなら
        console_method_candidates）。
        **`.mdproject` / `.mddata` は使えません**（ZIP なので Console は中身を
        読めず、全パラメータが既定値のまま実行されます）。
    lbm_file:
        脂質ライブラリ（`.lbm2`）のパス。省略時は「メソッドファイルの宣言 →
        ビルド生成物（MsdialWorkbench をソースからビルドしている場合） →
        環境変数 MSDIAL_LBM → MSDIAL_EXE と同じフォルダ」の順に、MS-DIAL GUI と
        同じ規則で自動解決します。GUI 由来のパラメータは `Lbm file path:` が
        必ず空なので、この自動解決が無いと**警告なしで同定 0 件**になります。
        自動補完（宣言より後ろの 3 つ）は lipidomics だけです。引数と宣言は
        metabolomics でも使います（Console は omics を問わず LBM で脂質を同定する）。
        解決した LBM と、相対で宣言されたライブラリ類のパス（`Msp file path`
        `Text DB file path` 等）は、メソッドファイル基準の絶対パスにして
        `run_dir/effective-method.txt` へ書きます（Console は相対パスを
        自分の作業フォルダ基準で読み、見つからなければ黙って飛ばすため）。
    polarity:
        "positive" または "negative"。
    measure:
        "peak_height"（既定）または "peak_area_above_zero"。
        現行 MS-DIAL Console は peak_height のみ正確に出力します。
        "peak_area_above_zero" を指定すると UNSUPPORTED_AREA_CONSOLE で停止します。
    omics:
        "lipidomics"（既定）または "metabolomics"。
    save_project:
        True（既定）なら MS-DIAL Console に -p を渡し、GUI で開ける
        .mdproject を出力フォルダに生成させます。
    timeout_s:
        MS-DIAL Console のタイムアウト秒数（既定 21600 ＝ 6 時間）。
        実測では 4 サンプルで約 3 分。60 サンプル規模では 1 時間を超え得ます。
    keep_extension:
        アライメントする計測形式（例 "wiff" / "wiff2"）。dataset_root に形式が
        混在しているとき（SCIEX は 1 測定につき .wiff と .wiff2 ができる）に
        使います。省略して混在していれば、選べる形式を `choices` に入れた
        MIXED_RAW_FORMATS を返すので、そこから選んで呼び直してください。
        指定すると、その形式の計測ファイルと随伴ファイルだけを集めたフォルダを
        `<元フォルダ>_<拡張子>` にハードリンクで作り（既にあれば欠けている分だけ
        補う）、そちらで計画します。元フォルダは変更しません。

    成功すると session.current_job_path にジョブパスが設定され、
    console_run でそのまま実行できます。
    """
    execution_error = _execution_options_error(save_project, timeout_s)
    if execution_error:
        return execution_error

    if measure == "peak_area_above_zero":
        return mztab_error(
            "UNSUPPORTED_AREA_CONSOLE",
            "現行 MS-DIAL Console は peak_area_above_zero を正確に出力できません。"
            "peak_height を使用するか、GUI から Area 出力を別途実行してください。",
        )

    if polarity not in ("positive", "negative"):
        return console_error("JOB_NOT_PLANNED", f"polarity は 'positive' または 'negative' です: {polarity!r}")
    if omics not in ("lipidomics", "metabolomics"):
        return console_error("JOB_NOT_PLANNED", f"omics は 'lipidomics' または 'metabolomics' です: {omics!r}")

    root = Path(dataset_root).expanduser()
    if not root.is_dir():
        return console_error("JOB_NOT_PLANNED", f"dataset_root が存在しません: {dataset_root}")

    from lipidmix.console import method_file as method_file_mod

    discovered_from: str | None = None
    if method_file is None:
        candidates, searched = method_file_mod.discover_method_candidates(
            root, polarity=polarity, omics=omics)
        direct = [c for c in candidates if c.usable == "direct"]
        if direct:
            discovered_from = direct[0].path
            mf = Path(discovered_from)
        elif candidates:
            # 別極性を黙って採ると Ion mode と Searched adduct ions が違うまま走り、
            # 別の解析になる。console_method_template を通して caveat を出させる。
            return console_error(
                "METHOD_FILE_CHOICE_REQUIRED",
                f"極性 {polarity} に一致するパラメータファイルはありませんが、"
                f"別極性の候補が {len(candidates)} 件あります。"
                "console_method_template(based_on=<選んだ path>, polarity=...) で"
                "その極性用に変換してから console_plan に渡してください。"
                "検出・アライメント条件は元のまま引き継がれます。",
                {"dataset_root": str(root), "polarity": polarity,
                 "searched": searched,
                 **_capped_candidates(candidates)},
                required_tools=["console_method_template", "console_plan"])
        else:
            return console_error(
                "METHOD_FILE_NOT_GIVEN",
                "method_file が省略され、使えるパラメータファイルも"
                f"見つかりませんでした（極性 {polarity}）: {root}  "
                "MS-DIAL GUI は解析のたびに `<project>_param_<終了時刻>.txt` を"
                "プロジェクトフォルダへ自動保存します。その極性でも別極性でも"
                "GUI 実行が無い場合は、他のデータセットのパラメータを "
                "console_method_template の based_on で明示してください。",
                {"dataset_root": str(root), "polarity": polarity,
                 "searched": searched},
                required_tools=["console_method_template", "console_method_candidates"])
    else:
        mf = Path(method_file).expanduser()
    if not mf.is_file():
        return console_error("METHOD_FILE_NOT_FOUND", f"メソッドファイルが見つかりません: {method_file}")
    if not _looks_like_method_text(mf):
        return console_error(
            "METHOD_FILE_NOT_TEXT",
            f"メソッドファイルが MS-DIAL Console の読める形式ではありません: {method_file}  "
            "Console は ASCII のテキストを `key: value`（例 `Ion mode: Negative`）として"
            "1 行ずつ読みます。.mdproject / .mddata は ZIP なので、渡しても"
            "エラーにならず全パラメータが既定値のまま実行されます。"
            "MS-DIAL GUI の Export > Parameter で出したパラメータファイルを指定してください。",
            {"method_file": str(mf)},
        )

    method_keys = method_file_mod.read_method_keys(mf)
    declared_mode = (method_keys.get("ion mode") or "").strip().lower()
    if declared_mode and declared_mode != polarity:
        return console_error(
            "METHOD_FILE_POLARITY_MISMATCH",
            f"メソッドファイルの `Ion mode: {method_keys.get('ion mode')}` と "
            f"polarity={polarity!r} が食い違っています: {mf}  "
            "Console はメソッドファイル側を使うため、このまま実行すると"
            "宣言と別の極性の結果が analysis-job に記録されます。",
            {"method_file": str(mf), "method_ion_mode": declared_mode, "declared": polarity})

    try:
        from lipidmix.console import runner as console_runner
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        return _exe_error(exc)
    if not console_runner.is_console_exe(exe):
        return console_error(
            "MSDIAL_EXE_NOT_CONSOLE",
            f"MSDIAL_EXE が MS-DIAL Console ではありません: {exe}  "
            "--help にサブコマンド `lcms` が現れませんでした。GUI の MSDIAL.exe を"
            "指している可能性があります（GUI はコマンドラインを解釈せずウィンドウを"
            "開いたままになります）。MsdialWorkbench の Console 実行体"
            "（MSDIALCUI.exe）のパスを設定してください。",
            {"exe": exe},
        )

    from lipidmix.console.job_manager import create_job, raw_input_summary
    formats = raw_input_summary(root)
    if not formats:
        return console_error(
            "MIXED_RAW_FORMATS",
            f"データフォルダに MS-DIAL が読める計測ファイルがありません: {dataset_root}  "
            "対象拡張子: abf / ibf / cdf / mzml / wiff / raw / d / wiff2 / qgd / lcd / lrp / imzml  "
            "（.d と .raw はフォルダ 1 つが 1 検体。それ以外はファイルであることが要ります）",
            {"formats": formats},
        )

    # 混在フォルダでは、選ばれた形式だけのフォルダを兄弟に作ってそちらで計画する。
    # method_file の自動探索は上で元フォルダに対して済ませてある（GUI の
    # パラメータは元フォルダ側に保存されるため）。
    prepared = None
    source_root = root
    if keep_extension:
        ext = keep_extension.lower().lstrip(".")
        if ext not in formats:
            return console_error(
                "MIXED_RAW_FORMATS",
                f"指定された形式 '{ext}' の計測ファイルがデータフォルダにありません: {root}  "
                "実在する形式: " + ", ".join(f"{e}×{n}" for e, n in sorted(formats.items())),
                {"formats": formats, "requested": ext},
            )
        if len(formats) > 1:
            from lipidmix.console.input_prep import prepare_single_format_input
            try:
                prepared = prepare_single_format_input(root, ext, _prepared_input_dir(root, ext))
            except (ValueError, OSError) as exc:
                return console_error("INPUT_PREP_FAILED", str(exc),
                                     {"dataset_root": str(root), "keep_extension": ext})
            root = Path(prepared.out_dir)
            formats = raw_input_summary(root)
    if len(formats) > 1:
        return console_error(
            "MIXED_RAW_FORMATS",
            "データフォルダに MS-DIAL が対象とする拡張子が 2 種類以上あります: "
            + ", ".join(f"{ext}×{n}" for ext, n in sorted(formats.items()))
            + "。どちらをアライメントするかを選び、console_plan を "
            "keep_extension=<choices の値> 付きで呼び直してください。"
            "選んだ形式の計測ファイルと随伴ファイルだけを集めたフォルダを"
            "`<元フォルダ>_<拡張子>` にハードリンクで作り、そこで計画します"
            "（実体コピーなし・元フォルダは無変更）。"
            "SCIEX の出力は 1 測定につき .wiff と .wiff2 が両方できるのが普通です。"
            "混在のまま MS-DIAL Console に渡すと対話プロンプトで異常終了し、"
            "続行できても同じ測定が複数の解析ファイルとして扱われます。",
            {"formats": formats,
             "choices": [{"keep_extension": ext, "count": n}
                         for ext, n in sorted(formats.items())]},
            required_tools=["console_plan"],
        )
    input_count = sum(formats.values())

    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        method_keys, mf, omics=omics, exe_path=exe, lbm_setting=lbm_setting, override=lbm_file)
    if lbm.error_code:
        return console_error(lbm.error_code, lbm.message or "",
                             {"candidates": list(lbm.candidates)} if lbm.candidates else None)

    # 解決した LBM と、相対で宣言されたライブラリ等のパスは絶対パスにして
    # run_dir の実効メソッドファイルに書く。LC-MS の Console は宣言パスを
    # 解決せず cwd（run_dir）基準で読むので、相対のまま渡すと見つからずに
    # 黙って飛ばされる（同定 0 件で完走）。ジョブを作る前に組んで ASCII 検査する。
    overrides = method_file_mod.relative_path_overrides(method_keys, mf)
    if lbm.path:
        overrides[method_file_mod.LBM_KEY] = lbm.path
    encoding_error = _non_ascii_override_error(overrides)
    if encoding_error:
        return encoding_error

    try:
        job, job_path = create_job(
            dataset_root=root,
            method_file=mf,
            polarity=polarity,  # type: ignore[arg-type]
            measure=measure,    # type: ignore[arg-type]
            omics=omics,        # type: ignore[arg-type]
            input_count=input_count,
        )
    except ValueError as exc:
        return console_error("DATASET_ROOT_IN_REPO", str(exc))

    from lipidmix.console.job_manager import save_job
    job.save_project = save_project
    job.timeout_s = timeout_s

    # 実効メソッドを書いてジョブをそちらへ向ける。ユーザーのパラメータファイルは
    # 触らない（GUI が次に開いたときの整合が崩れる）。
    if overrides:
        effective = method_file_mod.write_effective_method_file(
            mf, Path(job.run_dir) / "effective-method.txt", overrides)
        job.method_file = str(effective)
    save_job(job, job_path)
    session_state.session.current_job_path = str(job_path)

    warnings: list[str] = []
    if any(p.is_file() and (p.name.startswith("AlignResult-") or "AlignmentResult" in p.name)
           for p in root.iterdir()):
        warnings.append(
            "既存のアライメント結果がデータフォルダにあります。MS-DIAL Console は"
            "実行のたびに別タイムスタンプの一式を同じフォルダへ追加するため、"
            "複数バッチが混在します。どれが今回の生成物かは console_status の"
            "artifacts で確認してください。")

    return json_payload({
        "status": "planned",
        "job_id": job.job_id,
        "job_path": str(job_path),
        "run_dir": job.run_dir,
        "polarity": polarity,
        "measure": measure,
        "omics": omics,
        "input_count": input_count,
        "save_project": save_project,
        "timeout_s": timeout_s,
        "method_file": job.method_file,
        "method_source": {
            "given": method_file,
            "discovered_from": discovered_from,
            "effective": job.method_file,
        },
        "lbm": {"path": lbm.path, "source": lbm.source},
        "keep_extension": next(iter(formats)),
        **({"prepared_input": {"source": str(source_root), "out_dir": prepared.out_dir,
                               "linked": prepared.linked, "copied": prepared.copied}}
           if prepared is not None else {}),
        "warnings": warnings,
        "next": "console_run を呼び出して実行を開始してください",
    })


def _non_ascii_override_error(overrides: dict[str, str]) -> str | None:
    """実効メソッドへ書く値に ASCII で表せない文字があれば封筒を返す。

    Console はメソッドを `Encoding.ASCII` で読む（`ConfigParser`）。非 ASCII の
    フォルダ名は `?` に化けてライブラリが見つからず、黙って同定 0 件になる。
    pipeline の `METHOD_ENCODING_UNSUPPORTED` と同じコードで止める。
    """
    bad = {key: value for key, value in overrides.items() if not value.isascii()}
    if not bad:
        return None
    return console_error(
        "METHOD_ENCODING_UNSUPPORTED",
        "実効メソッドへ書くパスに ASCII 以外の文字があります。MS-DIAL Console は"
        "メソッドファイルを ASCII で読むため、このパスのライブラリは見つからず"
        "黙って同定 0 件になります。ASCII だけのフォルダへ置き直してください: "
        + ", ".join(f"{key}={value}" for key, value in sorted(bad.items())),
        {"keys": sorted(bad)})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False),
          structured_output=False)
def console_run(job_path: str | None = None, detach: bool = False) -> str:
    """MS-DIAL Console を実行します。

    job_path: analysis-job.json へのパス。省略時は session.current_job_path を使用します。
    事前に console_plan を実行しておく必要があります。

    detach: True にすると**監視ワーカーを親から切り離して起動し、待たずに戻ります**。
        実データ 60 サンプルは約 44 分かかるため、既定の同期実行では 1 ツール
        呼び出しがその間ずっと戻らず、呼び出し元が中断されると生成物ごと失います。
        **長い実行ではこちらを使ってください。**

        返る `pid` は**監視ワーカーのもの**で、MS-DIAL Console のものではありません
        （Console の pid は run_dir の execution-result.json に載ります）。ワーカーが
        終了まで見張り、生成物の収集とジョブの確定まで行うので、`console_status` を
        呼ばなくても結果は確定します。

    どちらの経路でも run_dir に終了証跡（execution-result.json）が残ります。
    終了コードが 0 でも、主 mzTab-M が一意に選べ、定量行列に有限値があり、予定した
    入力が全て assay に対応していなければ completed にはなりません（partial / failed）。

    実行完了後、session.current_job_path は同じジョブを指し続けます。
    結果は console_status または dataset_load で確認してください。
    """
    resolved = _resolve_job_path(job_path)
    if isinstance(resolved, str):
        return resolved  # error envelope

    from lipidmix.console.job_manager import load_job, update_status
    try:
        job = load_job(resolved)
    except FileNotFoundError:
        return console_error("JOB_NOT_FOUND", f"analysis-job.json が見つかりません: {resolved}")
    except ValueError as exc:
        return console_error("JOB_NOT_FOUND", str(exc))

    pipeline_block = _pipeline_owner_block(Path(job.run_dir))
    if pipeline_block is not None:
        return console_error(
            "JOB_OWNED_BY_PIPELINE",
            "このジョブは pipeline が所有しています。単体 console_run では"
            "実行できません。取消・再実行は所有 pipeline 側の操作から行って"
            f"ください（pipeline_root={pipeline_block.get('pipeline_path')}）。",
            {"job_id": job.job_id, **pipeline_block},
        )

    if job.status != "planned":
        return console_error(
            "JOB_NOT_PLANNED",
            f"実行できるのは status=planned のジョブだけです。現在のステータス: {job.status}",
            {"job_id": job.job_id, "status": job.status},
        )

    execution_error = _execution_options_error(job.save_project, job.timeout_s)
    if execution_error:
        return execution_error

    from lipidmix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        update_status(resolved, "failed", error=str(exc))
        return _exe_error(exc)
    try:
        is_console = console_runner.is_console_exe(exe, raise_on_os_error=True)
    except OSError as exc:
        update_status(resolved, "failed", error=str(exc))
        return console_error("MSDIAL_EXE_NOT_FOUND", str(exc))
    if not is_console:
        return console_error(
            "MSDIAL_EXE_NOT_CONSOLE",
            f"MSDIAL_EXE が MS-DIAL Console ではありません: {exe}  "
            "--help にサブコマンド `lcms` が現れませんでした。GUI の MSDIAL.exe を"
            "指している可能性があります（GUI はコマンドラインを解釈せずウィンドウを"
            "開いたままになります）。MsdialWorkbench の Console 実行体"
            "（MSDIALCUI.exe）のパスを設定してください。",
            {"exe": exe},
        )

    run_dir = Path(job.run_dir)
    dataset_root = Path(job.dataset_root)

    # 何を入力として実行するかを、起動より前に固定する。実行後に mzTab の
    # ms_run[N]-location と 1 対 1 で突き合わせるのはこの目録で、無いまま走らせると
    # 「予定した検体が全部入っているか」を誰も検証できない。
    from lipidmix.console.execution import write_supervision_inputs
    from lipidmix.console.job_manager import list_raw_inputs
    try:
        write_supervision_inputs(run_dir, {
            "raw_inventory": [str(p.resolve()) for p in list_raw_inputs(dataset_root)],
            "method_sha256": _file_sha256_or_none(Path(job.method_file)),
            "exe_path": exe,
            "exe_sha256": _file_sha256_or_none(Path(exe)),
        })
    except OSError as exc:
        update_status(resolved, "failed", error=repr(exc))
        return console_error(
            "JOB_POST_RUN_FAILED",
            f"実行前の監視入力を run_dir へ書けませんでした: {exc!r}",
            {"run_dir": str(run_dir)})

    from lipidmix.console import worker as console_worker
    from lipidmix.core.atomic_io import DomainError

    if detach:
        try:
            launched = console_worker.launch_console_worker(resolved)
        except (DomainError, OSError) as exc:
            update_status(resolved, "failed", error=repr(exc))
            return console_error(
                "WORKER_LAUNCH_FAILED",
                f"監視ワーカーを起動できませんでした: {exc}",
                {"job_path": str(resolved),
                 "worker_log": str(console_worker.worker_log_path(run_dir))})
        # 起動できてから running にする。起動前に書くと、失敗したときに
        # 誰も走っていないジョブが running のまま残る。
        #
        # ただし**無条件には書かない**。ワーカーは切り離して起動しており、
        # ここへ来るまでに自分で running を書き、走り切って completed まで
        # 書き終えていることがある（短い実行・偽 Console）。そこへ親が
        # running を上書きすると、終わった実行が永久に running のまま残る
        # ——`console_status` は読取専用なので、あとから誰も直せない。
        # 所有記録も同じ理由で、ワーカーが自分の identity を書いていたら
        # 触らない（親が知っているのはワーカーの pid だけで、そちらのほうが
        # 情報が古い）。
        warnings: list[str] = []
        current = load_job(resolved)
        if current.status == "planned":
            update_status(resolved, "running")
        try:
            if console_worker.read_owner(run_dir) is None:
                # 起動できた事実を先に記録する。ここで失敗しても走り出した
                # ワーカーは取り消せないので、応答では必ず pid を返す。
                console_worker.write_owner(run_dir, {
                    "kind": console_worker.OWNER_KIND_CONSOLE,
                    "pid": launched["pid"],
                    "identity": launched["identity"],
                    "job_path": str(resolved),
                    "status": "running",
                })
        except OSError as exc:
            warnings.append(
                f"監視ワーカー（pid={launched['pid']}）は起動しましたが、"
                f"所有記録を worker.json へ書けませんでした: {exc!r}")
        return json_payload({
            "status": load_job(resolved).status,
            "job_id": job.job_id,
            "job_path": str(resolved),
            "pid": launched["pid"],
            # pid はワーカーのもの。Console の pid は起動後に証跡へ載る。
            "pid_of": "worker",
            "run_dir": job.run_dir,
            "log": str(run_dir / "msdial.log"),
            "worker_log": str(console_worker.worker_log_path(run_dir)),
            "warnings": warnings,
            "next": "console_status で完了を確認してください"
                    "（ワーカーが最後まで監視し、生成物の収集まで行います）",
        })

    try:
        receipt = console_worker.run_job(resolved)
    except DomainError as exc:
        if exc.code == "LOCK_TIMEOUT":
            # 状態は書き換えない。このジョブを持っているのは別のプロセスで、
            # そちらが終端状態を書く。ここで planned へ戻すと、勝者が書いた
            # running / completed を敗者が消してしまう。
            return console_error(
                "JOB_BUSY",
                "このジョブは別のプロセスが実行中です。"
                "console_status で進行を確認してください。",
                {"job_path": str(resolved), "run_dir": job.run_dir,
                 "owner": console_worker.owner_summary(run_dir)},
                required_tools=["console_status"])
        if exc.code == "JOB_ALREADY_FINISHED":
            # ロック待ちの間に別プロセスが走り切っていた。状態は書き換えない
            # ——勝者が書いた終端状態を敗者が壊してはいけない（LOCK_TIMEOUT と
            # 同じ理由）。二度目の MS-DIAL 起動は worker 側で止まっている。
            return console_error(
                "JOB_NOT_PLANNED", exc.message,
                {"job_id": job.job_id, **(exc.details or {})},
                required_tools=["console_status"])
        _record_finalization_failure(resolved, repr(exc))
        return console_error("JOB_POST_RUN_FAILED", str(exc),
                             {"job_path": str(resolved)})
    except Exception as exc:  # 想定外。running に固着させないことが最優先
        _record_finalization_failure(resolved, repr(exc))
        return console_error(
            "JOB_POST_RUN_FAILED",
            f"MS-DIAL Console の実行中に想定外のエラーが発生しました: {exc!r}",
            {"log": str(run_dir / "msdial.log")})

    return _console_run_result(resolved, receipt)


def _prepared_input_dir(src: Path, keep_extension: str) -> Path:
    """単一形式フォルダの既定の置き場所: `<元フォルダ>_<拡張子>` の兄弟。"""
    return src.parent / f"{src.name}_{keep_extension.lower().lstrip('.')}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def console_prepare_input(
    dataset_root: str,
    keep_extension: str = "wiff",
    out_dir: str | None = None,
) -> str:
    """計測フォーマットが 1 種類だけの入力フォルダを作ります（MIXED_RAW_FORMATS の解消）。

    MS-DIAL は `.wiff` と `.wiff2` を別フォーマットとして数えるため、SCIEX の
    生データフォルダはそのままでは対話プロンプトが出て実行できません。
    このツールは指定した拡張子の計測ファイルと**その随伴ファイル**
    （`.wiff.scan` / `.timeseries.data` 等）だけを集めたフォルダを作ります。

    Agilent・Bruker の `.d` と Waters の `.raw` のように**フォルダそのものが
    1 検体**の形式も扱えます。その場合はフォルダ構造を実体として作り直し、
    中のファイルだけをリンクします（フォルダ自体はリンクにしません）。

    実体はコピーせずハードリンクを張ります（同一ボリュームで無い場合のみコピー）。
    **元フォルダは一切変更しません。** MS-DIAL の生成物はこの新しいフォルダ側に
    出るので、生データ本体を汚さずに済みます。

    dataset_root: 生データフォルダ。
    keep_extension: 残す計測拡張子（既定 "wiff"）。MS-DIAL が読む 12 形式
        （abf / ibf / cdf / mzml / wiff / raw / d / wiff2 / qgd / lcd / lrp /
        imzml）から選びます。`MIXED_RAW_FORMATS` の `formats` が実際の内訳です。
    out_dir: 出力先。省略時は `<元フォルダ>_<拡張子>` を兄弟として作ります。
    """
    src = Path(dataset_root).expanduser()
    dest = Path(out_dir).expanduser() if out_dir else _prepared_input_dir(src, keep_extension)

    from lipidmix.console.input_prep import prepare_single_format_input
    try:
        result = prepare_single_format_input(src, keep_extension, dest)
    except (ValueError, OSError) as exc:
        return console_error("INPUT_PREP_FAILED", str(exc),
                             {"dataset_root": str(src), "keep_extension": keep_extension})

    from lipidmix.console.job_manager import raw_input_summary
    return json_payload({
        "status": "prepared",
        "out_dir": result.out_dir,
        "primary": result.primary,
        "companions": result.companions,
        "linked": result.linked,
        "copied": result.copied,
        "skipped_existing": result.skipped_existing,
        "mode": result.mode,
        "formats": raw_input_summary(Path(result.out_dir)),
        "next": f"console_plan(dataset_root={result.out_dir!r}, ...) を実行してください",
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def console_status(job_path: str | None = None, include_artifacts: bool = False) -> str:
    """ジョブの現在のステータスを返します。

    **保存済みの状態を読むだけ**で、完了処理はしません。実行の監視・生成物の
    収集・ジョブの確定は切り離しワーカー（console_run が起こす）が最後まで行うので、
    このツールを呼ばなくても結果は確定します。

    job_path: analysis-job.json へのパス。省略時は session.current_job_path を使用します。
    include_artifacts: 生成物の全文一覧（TSV）を含めます。既定 False。
        生成物は 1 サンプルにつき 5 件出るため 60 サンプルで 300 行を超え、
        全文を返すと後ろの warnings / error が埋没します。既定では件数
        (`artifact_count`) と役割別内訳 (`artifacts_by_role`) だけを返します。
        個々のパスが必要なとき（欠落の特定など）だけ True にしてください。
    """
    resolved = _resolve_job_path(job_path)
    if isinstance(resolved, str):
        return resolved

    from lipidmix.console.job_manager import load_job
    try:
        job = load_job(resolved)
    except (FileNotFoundError, ValueError) as exc:
        return console_error("JOB_NOT_FOUND", str(exc))

    legacy = _legacy_detached_note(job)
    if isinstance(legacy, str):
        return legacy  # error envelope（旧 detach の終了が未解決）

    from lipidmix.core import version
    update = version.update_status()
    return json_payload({
        "server_version": version.server_version(),
        # 最新なら鍵ごと出さない（戻り値は LLM の文脈をそのまま食う）。
        **({"update_available": update} if update else {}),
        "job_id": job.job_id,
        "status": job.status,
        **({"detached": legacy} if legacy else {}),
        **({"execution_receipt": _receipt_summary(job)}
           if _receipt_summary(job) else {}),
        "polarity": job.polarity,
        "measure": job.measure,
        "omics": job.omics,
        "dataset_root": job.dataset_root,
        "run_dir": job.run_dir,
        "mztab_files": [
            {"path": e.path, "polarity": e.polarity, "measure": e.measure, "root": e.root}
            for e in job.primary_mztab_files
        ],
        "artifact_count": len(job.artifacts),
        # 生成物は 1 サンプルにつき複数出る（.pai2 / .dcl / _tags.xml / .mdpeak /
        # .mdmsp）。60 サンプルの実走で 313 件になり、全文 TSV は数万字に達して
        # 後ろの warnings / error を埋没させた。既定は役割別の件数だけにし、
        # 全文は include_artifacts=True のときだけ返す。
        "artifacts_by_role": _artifacts_by_role(job.artifacts),
        **({"artifacts": _artifacts_tsv(job.artifacts)} if include_artifacts else {}),
        "execution": {"save_project": job.save_project, "timeout_s": job.timeout_s},
        "warnings": job.warnings,
        "error": job.error,
        "updated_at": job.updated_at,
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def console_method_template(
    out_path: str,
    polarity: str,
    based_on: str | None = None,
    dataset_root: str | None = None,
    omics: str = "lipidomics",
) -> str:
    """既存のパラメータから、別極性用のメソッドファイルを作ります（最後の手段）。

    **まず console_plan の method_file 省略を試してください。** MS-DIAL GUI は解析の
    たびに `<project>_param_<終了時刻>.txt` を自動保存するので、その極性で一度でも
    GUI 実行があれば、作る必要はありません。

    このツールが要るのは「その極性で一度も処理していない」場合だけです
    （別極性は処理済み、というのが典型）。やることは 2 行の差し替えです:
    `Ion mode` と `Searched adduct ions` をその極性の標準セットにする。
    **検出・アライメント条件は元のまま引き継ぎます**（勝手に変えると別の解析になる）。
    加えて、GUI 由来では必ず空の `Lbm file path` を解決して埋めます。

    out_path: 書き出し先。
    polarity: "positive" / "negative"（作りたい側）。
    based_on: 元にするパラメータファイル。省略時は dataset_root から探します
        （**極性は問いません** — 別極性から作るのがこのツールの用途なので）。
    dataset_root: based_on 省略時の探索先。dataset_root 直下・兄弟フォルダ・
        過去 run まで探します。候補が複数あれば based_on で選ぶよう
        METHOD_FILE_CHOICE_REQUIRED で止まります（土台の選択は解析条件そのもの
        なので、最新を黙って採りません）。
    """
    if polarity not in ("positive", "negative"):
        return console_error("JOB_NOT_PLANNED",
                             f"polarity は 'positive' または 'negative' です: {polarity!r}")

    from lipidmix.console import method_file as method_file_mod

    if based_on:
        src = Path(based_on).expanduser()
    elif dataset_root:
        # 極性で絞らない。「別極性から作る」のがこのツールの用途なので、
        # polarity=None のまま候補を集める。
        candidates, searched = method_file_mod.discover_method_candidates(
            Path(dataset_root).expanduser(), polarity=None, omics=omics)
        if not candidates:
            return console_error(
                "METHOD_FILE_NOT_GIVEN",
                "元にできるパラメータファイル（`*_param_<ts>.txt`）が見つかりません: "
                f"{dataset_root}  MS-DIAL GUI で一度も解析していないフォルダには"
                "存在しません（兄弟フォルダと過去 run も探しました）。"
                "他のデータセットのパラメータを based_on で明示してください。",
                {"dataset_root": dataset_root, "searched": searched})
        if len(candidates) > 1:
            return console_error(
                "METHOD_FILE_CHOICE_REQUIRED",
                f"元にできる候補が {len(candidates)} 件あります。based_on で 1 つ"
                "選んでください。検出・アライメント条件は選んだファイルのものが"
                "そのまま引き継がれるため、どれを土台にするかは解析条件の選択です。",
                {"dataset_root": dataset_root, "searched": searched,
                 **_capped_candidates(candidates)},
                required_tools=["console_method_template"])
        src = Path(candidates[0].path)
    else:
        return console_error("METHOD_FILE_NOT_GIVEN",
                             "based_on か dataset_root のどちらかを指定してください。")

    if not src.is_file():
        return console_error("METHOD_FILE_NOT_FOUND", f"元ファイルが見つかりません: {src}")
    if not _looks_like_method_text(src):
        return console_error(
            "METHOD_FILE_NOT_TEXT",
            f"元ファイルが MS-DIAL Console の読める ASCII テキストではありません: {src}  "
            ".mdproject / .mddata は ZIP で、解析パラメータを含んでいません"
            "（中身は .mddata へのポインタだけです）。",
            {"based_on": str(src)})

    src_keys = method_file_mod.read_method_keys(src)
    # 別フォルダへ書き出すので、相対宣言は原本基準の絶対パスに固定する。
    overrides = {
        **method_file_mod.relative_path_overrides(src_keys, src),
        method_file_mod.ION_MODE_KEY: polarity.capitalize(),
        method_file_mod.ADDUCT_KEY: method_file_mod.STANDARD_ADDUCTS[polarity],
    }

    from lipidmix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        # exe はビルド生成物・exe フォルダからの LBM 推定にだけ使う。未設定でも
        # メソッドの宣言や [msdial] lbm で解決できれば足りるので、ここでは止めない。
        if getattr(exc, "code", "") == "CONFIG_INVALID":
            return _exe_error(exc)
        exe = None
    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        src_keys, src, omics=omics, exe_path=exe, lbm_setting=lbm_setting)
    if lbm.error_code:
        return console_error(lbm.error_code, lbm.message or "",
                             {"candidates": list(lbm.candidates)} if lbm.candidates else None)
    if lbm.path:
        overrides[method_file_mod.LBM_KEY] = lbm.path
    encoding_error = _non_ascii_override_error(overrides)
    if encoding_error:
        return encoding_error

    dest = Path(out_path).expanduser()
    try:
        method_file_mod.write_effective_method_file(src, dest, overrides)
    except OSError as exc:
        return console_error("METHOD_FILE_NOT_FOUND", f"書き出せませんでした: {exc}")

    return json_payload({
        "status": "written",
        "out_path": str(dest),
        "based_on": str(src),
        "polarity": polarity,
        "lbm": {"path": lbm.path, "source": lbm.source},
        "changed_keys": sorted(overrides),
        "caveat": "検出・アライメント条件は元ファイルのまま引き継いでいます。"
                  "その極性に妥当かは実行前に確認してください。",
        "next": f"console_plan(dataset_root=..., method_file={str(dest)!r}, polarity={polarity!r})",
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def console_method_candidates(
    dataset_root: str,
    polarity: str | None = None,
    omics: str | None = "lipidomics",
    search_dirs: list[str] | None = None,
) -> str:
    """使えるメソッドファイルの候補を列挙します（console_plan が失敗する前に呼べます）。

    MS-DIAL GUI は解析のたびに `<project>_param_<終了時刻>.txt` をプロジェクト
    フォルダへ自動保存します。その極性で一度も GUI 実行が無いフォルダには存在
    しないため、**兄弟フォルダ**（POS の隣の NEG 等）と**過去 run** まで探します。

    dataset_root: 生データフォルダ。
    polarity: "positive" / "negative"。省略すると極性で区別せず全件返します。
        指定すると各候補に `usable` が付き、一致しないものは
        `needs_polarity_conversion`（console_method_template を経由させる）。
    omics: 既定 "lipidomics"。None で絞り込みません。
    search_dirs: 追加で探すフォルダ。再帰はしません。

    候補には比較用の `key_params`（検出・アライメント条件のうち結果を変える少数）が
    付きます。候補が 10 件を超えるときは付きません（戻り値が肥大するため）。
    """
    root = Path(dataset_root).expanduser()
    if not root.is_dir():
        return console_error("DATASET_ROOT_NOT_FOUND",
                             f"データフォルダが見つかりません: {dataset_root}",
                             {"dataset_root": str(dataset_root)})

    from lipidmix.console import method_file as method_file_mod

    candidates, searched = method_file_mod.discover_method_candidates(
        root, polarity=polarity, omics=omics, search_dirs=search_dirs)
    return json_payload({
        "dataset_root": str(root),
        "polarity": polarity,
        "omics": omics,
        "searched": searched,
        "n_candidates": len(candidates),
        **_capped_candidates(candidates),
        "next": _candidates_next_hint(polarity),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True),
          structured_output=False)
def console_cleanup(job_path: str | None = None, dry_run: bool = True) -> str:
    """あるジョブが生成したファイルだけを一覧・削除します（既定は一覧のみ）。

    MS-DIAL Console は生データフォルダ側にも生成物を出すため、再実行のたびに
    別タイムスタンプのアライメント一式が同じフォルダへ積まれます。放置すると
    「複数バッチ混在フォルダ」になり、どれが今回の結果か分からなくなります。

    **消す対象は analysis-job.json が記録した生成物だけ**です。タイムスタンプの
    推測では消しません（別バッチの成果物を巻き込むため）。記録が無いジョブは
    NO_JOB_OUTPUT で拒否します。生データ（.wiff 等）には触れません。

    dry_run: True（既定）なら一覧を返すだけ。False で実際に削除します。
        削除後、ジョブの status は `cleaned` になります（生成物を指したまま
        completed で残ると、dataset_load が存在しないファイルを読もうとします）。
    """
    resolved = _resolve_job_path(job_path)
    if isinstance(resolved, str):
        return resolved

    from lipidmix.console.job_manager import load_job, update_status
    try:
        job = load_job(resolved)
    except (FileNotFoundError, ValueError) as exc:
        return console_error("JOB_NOT_FOUND", str(exc))

    # 生成物は run_dir と dataset_root に書き分かれる。記録された出所から戻す
    # 解決は mztab_tools が正準（関数レベル import で循環を避ける）。
    from lipidmix.tools.mztab_tools import _artifact_abs_path

    targets: list[Path] = []
    for entry in job.primary_mztab_files:
        targets.append(_artifact_abs_path(job, getattr(entry, "root", "run_dir"), entry.path))
    for art in job.artifacts:
        targets.append(_artifact_abs_path(job, getattr(art, "root", "run_dir"), art.path))

    if not targets:
        return console_error(
            "NO_JOB_OUTPUT",
            "このジョブは生成物を 1 件も記録していないため、何を消してよいか決められません"
            f"（status={job.status}）。中断された実行はここに該当します。"
            "ファイル名のタイムスタンプで推測すると別バッチの成果物を巻き込むため、"
            "この場合は console_status と実フォルダを見て手で片付けてください。",
            {"job_id": job.job_id, "status": job.status,
             "dataset_root": job.dataset_root, "run_dir": job.run_dir})

    # 実行中のジョブの生成物を消すと、監視ワーカーが書いている最中のファイルを
    # 足元から抜くことになる。所有者が今も生きているかは process identity で見る
    # （pid だけでは pid 再利用を見分けられない）。
    from lipidmix.console.worker import owner_is_active, owner_summary
    owned = owner_is_active(Path(job.run_dir))
    pipeline_block = _pipeline_owner_block(Path(job.run_dir))

    if dry_run:
        warnings: list[str] = []
        if owned:
            warnings.append("このジョブは監視ワーカー（worker）が実行中です。"
                            "実行が終わるまで削除は拒否されます。")
        if pipeline_block is not None:
            warnings.append(
                "このジョブは pipeline が所有しています "
                f"（pipeline_root={pipeline_block.get('pipeline_path')}）。"
                "所有 pipeline が活動中である限り削除は拒否されます。")
        return json_payload({
            "status": "dry_run",
            "dry_run": True,
            "job_id": job.job_id,
            "count": len(targets),
            "files": [str(p) for p in targets],
            "warnings": warnings,
            "next": "実際に削除するには dry_run=False を指定してください",
        })

    if owned:
        return console_error(
            "JOB_BUSY",
            "このジョブは監視ワーカーが実行中です。生成物を消すと、書き込み中の"
            "ファイルを実行中のプロセスから奪うことになります。"
            "console_status で完了を確認してからやり直してください。",
            {"job_id": job.job_id, "status": job.status, "run_dir": job.run_dir,
             "owner": owner_summary(Path(job.run_dir))},
            required_tools=["console_status"])

    if pipeline_block is not None:
        return console_error(
            "JOB_OWNED_BY_PIPELINE",
            "このジョブは pipeline が所有しています。所有 pipeline が活動中、"
            "または状態を確認できない間は削除できません"
            f"（pipeline_root={pipeline_block.get('pipeline_path')}）。"
            "取消・再開は所有 pipeline 側の操作から行ってください。",
            {"job_id": job.job_id, "status": job.status, "run_dir": job.run_dir,
             **pipeline_block})

    deleted = absent = failed = 0
    errors: list[str] = []
    for path in targets:
        if not path.exists():
            absent += 1
            continue
        try:
            path.unlink()
            deleted += 1
        except OSError as exc:
            failed += 1
            errors.append(f"{path}: {exc}")

    update_status(resolved, "cleaned")
    return json_payload({
        "status": "cleaned",
        "dry_run": False,
        "job_id": job.job_id,
        "deleted": deleted,
        "already_absent": absent,
        "failed": failed,
        "errors": errors,
        "note": "生データ（.wiff 等）と analysis-job.json は残しています。",
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def job_list(dataset_root: str) -> str:
    """dataset_root/runs/ 以下のジョブ一覧を新しい順に返します。

    dataset_root: 生データフォルダのパス（console_plan で指定したものと同じ）。
    """
    root = Path(dataset_root).expanduser()
    if not root.is_dir():
        return console_error("JOB_NOT_FOUND", f"dataset_root が存在しません: {dataset_root}")

    from lipidmix.console.job_manager import list_jobs, load_job
    paths = list_jobs(root)
    if not paths:
        return json_payload({"dataset_root": str(root), "jobs": [], "count": 0})

    jobs = []
    for p in paths:
        try:
            j = load_job(p)
            jobs.append({
                "job_id": j.job_id,
                "status": j.status,
                "polarity": j.polarity,
                "measure": j.measure,
                "updated_at": j.updated_at,
                "job_path": str(p),
            })
        except Exception:
            jobs.append({"job_path": str(p), "status": "unreadable"})

    return json_payload({"dataset_root": str(root), "jobs": jobs, "count": len(jobs)})


# ---------- 内部ヘルパ ----------

def _pipeline_owner_block(run_dir: Path) -> dict | None:
    """このrun_dirを所有するpipelineが活動中／判定不能なら拒否理由を返す。

    None なら単体 console_run / console_cleanup を続行してよい。判定は
    lipidmix.pipeline.store（Task 14）へ委譲する——所有権sidecarの読み書き・
    pipeline-run.json の状態解釈はそちらが正準。
    """
    from lipidmix.pipeline import store as pipeline_store
    return pipeline_store.pipeline_owner_block_reason(run_dir)


def _file_sha256_or_none(path: Path) -> str | None:
    """読めれば SHA-256、読めなければ None（証跡側が読めなかった標識を置く）。"""
    from lipidmix.handoff.schema import sha256_file
    try:
        return sha256_file(path)
    except OSError:
        return None


def _console_run_result(job_path: Path, receipt: dict) -> str:
    """終了証跡と保存済みジョブから、同期実行の応答を組む。

    分岐は「どう終わったか」の因果順に並べる。起動できなかった実行を非ゼロ終了と
    呼んだり、収集に失敗した実行を出力ゼロと呼んだりすると、復旧手順が変わる。
    どの経路でも終了証跡（execution-result.json）は run_dir に残っているので、
    封筒には必ずその場所を入れる。
    """
    from lipidmix.console.execution import receipt_path
    from lipidmix.console.job_manager import load_job

    job = load_job(job_path)
    if receipt["job_save"]["status"] != "succeeded":
        # supervise が終端状態を書けなかった。封筒を組む前にここで直す——
        # running のまま残ると以降の console_run が全部 JOB_NOT_PLANNED で
        # 拒否され、誰もこのジョブを直せなくなる。終了の理由（timeout 等）は
        # 保存失敗より前の事実なので、両方を 1 行に残す。
        _record_finalization_failure(
            job_path,
            f"termination={receipt['termination']} exit_code={receipt['exit_code']}: "
            f"ジョブを保存できませんでした: {receipt['job_save'].get('error')}")
        job = load_job(job_path)
    run_dir = Path(job.run_dir)
    validation = receipt.get("validation") or {}
    details = {
        "job_path": str(job_path),
        "status": job.status,
        "termination": receipt["termination"],
        "exit_code": receipt["exit_code"],
        "mztab_files": len(job.primary_mztab_files),
        "other_artifacts": len(job.artifacts),
        "errors": validation.get("errors", []),
        # 完了しなかった実行こそ所見が要る。warnings は数件で、原因の説明を
        # 別の呼び出し（console_status）に取りに行かせない。
        "warnings": job.warnings,
        "log": str(run_dir / "msdial.log"),
        "receipt": str(receipt_path(run_dir)),
    }

    if job.status == "completed":
        return json_payload({
            "status": "completed",
            "job_id": job.job_id,
            "job_path": str(job_path),
            "mztab_files": len(job.primary_mztab_files),
            "other_artifacts": len(job.artifacts),
            "run_dir": job.run_dir,
            "warnings": job.warnings,
            "next": "dataset_load でmzTab-M を読み込み、解析を開始してください",
        })

    if receipt["termination"] == "launch_failed":
        return console_error(
            "MSDIAL_EXE_NOT_FOUND",
            f"MS-DIAL Console を起動できませんでした: {receipt.get('error', '')}",
            {**details, "exe": _configured_exe()})
    if receipt["termination"] == "timeout":
        return console_error(
            "MSDIAL_TIMEOUT",
            f"MS-DIAL Console がタイムアウトしました（{receipt['timeout_s']}s）",
            details)
    if receipt["termination"] == "cancelled":
        return console_error("MSDIAL_CANCELLED", "実行が取り消されました", details)
    if receipt["termination"] == "worker_lost":
        return console_error(
            "EXECUTION_UNRESOLVED",
            "Console の停止も終了コードの回収もできませんでした。終了の事実が"
            "確定していないため、成果物だけで完了とは判定しません。",
            details)
    if receipt["exit_code"] != 0:
        return console_error(
            "MSDIAL_NONZERO_EXIT",
            f"MS-DIAL Console が終了コード {receipt['exit_code']} で終了しました",
            details)
    if receipt["collection"]["status"] != "succeeded":
        return console_error(
            "JOB_POST_RUN_FAILED",
            "MS-DIAL Console の実行後処理（生成物の収集）でエラーが発生しました: "
            f"{receipt['collection'].get('error')}",
            details)
    if receipt["job_save"]["status"] != "succeeded":
        return console_error(
            "JOB_POST_RUN_FAILED",
            "MS-DIAL Console の実行後処理（ジョブの保存）でエラーが発生しました: "
            f"{receipt['job_save'].get('error')}",
            {**details, "recovery": receipt["job_save"].get("recovery")})
    if not job.primary_mztab_files and not job.artifacts:
        return console_error(
            "NO_JOB_OUTPUT",
            "MS-DIAL Console が終了しましたが、出力ファイルが生成されませんでした。"
            f"ログを確認してください: {run_dir / 'msdial.log'}",
            details)
    return console_error(
        "OUTPUT_VALIDATION_FAILED",
        "MS-DIAL Console は終了しましたが、出力が完了条件を満たしませんでした"
        "（details.errors に検証コードが入ります）。",
        details)


def _legacy_detached_note(job) -> dict | str | None:
    """旧 `.detached-state.json` を持つジョブの扱いを決める（状態は書き換えない）。

    旧経路は誰も監視しないまま Console を放流していたので、この state には
    **終了コードも process identity も無い**。プロセスが消えていても分かるのは
    「もう走っていない」ことだけで、成功したかどうかは分からない。ファイルが
    増えたという事実だけで completed へ進めると、途中で落ちた実行が完了に化ける
    ——それが `console-execution.v1` の証跡を導入した理由そのもの。

    Returns
    -------
    None
        旧 state は無い（新しい実行は証跡で判定する）。
    dict
        まだ生きている旧実行。実行中として表示する。
    str
        終了済みだが結果が確定できない旧実行。`EXECUTION_UNRESOLVED` の封筒。
    """
    from lipidmix.console.detached import read_detached_state

    run_dir = Path(job.run_dir)
    state = read_detached_state(run_dir)
    if state is None:
        return None

    from lipidmix.console import runner as console_runner
    pid = state["pid"]
    if console_runner.is_process_running(pid):
        return {"pid": pid, "alive": True, "legacy": True}

    return console_error(
        "EXECUTION_UNRESOLVED",
        "監視されていない旧方式の実行（.detached-state.json）が残っています。"
        "終了コードも実行の同一性も記録されていないため、生成物の有無だけでは"
        "完了と判定できません。msdial.log を確認し、必要なら console_plan から"
        "実行し直してください（新しい実行は終了証跡を残します）。",
        {"job_id": job.job_id, "status": job.status, "pid": pid,
         "run_dir": job.run_dir, "log": str(run_dir / "msdial.log"),
         "detached_state": str(run_dir / ".detached-state.json")},
        required_tools=["console_plan"])


def _receipt_summary(job) -> dict | None:
    """終了証跡があれば、その要点だけを返す（無い・壊れているなら None）。

    `validate_execution_record` を通してから読む。検証していない記録を状態表示に
    流用すると、壊れた証跡がそのまま「実行の事実」として伝播する。
    """
    import json as _json

    from lipidmix.console.execution import receipt_path, validate_execution_record
    from lipidmix.core.atomic_io import DomainError

    try:
        data = _json.loads(receipt_path(Path(job.run_dir)).read_text(encoding="utf-8"))
        record = validate_execution_record(data)
    except (OSError, ValueError, DomainError):
        return None
    return {
        "execution_id": record["execution_id"],
        "termination": record["termination"],
        "exit_code": record["exit_code"],
        "ended_at": record["ended_at"],
    }


def _looks_like_method_text(path: Path) -> bool:
    """MS-DIAL Console の ConfigParser が読める形かを判定する。

    ConfigParser は ASCII のテキストを 1 行ずつ読み、`#` 始まりを飛ばして
    最初の `:` または `=` で key/value に割る。ここでは先頭 8192 バイトだけを
    調べ、テキストであることと key/value 行が 1 つ以上あることを確認する。

    **先頭だけを読む**。この判定が弾く相手（.mdproject / .mddata）は 60 サンプル
    規模で GB 級になり、全体を読むとエラー封筒を返す前にメモリを潰す。
    """
    try:
        with path.open("rb") as fh:
            head = fh.read(8192)
    except OSError:
        return False
    if b"\x00" in head:
        return False
    text = head.decode("ascii", errors="replace")
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        positions = [i for i in (stripped.find(":"), stripped.find("=")) if i > 0]
        if positions:
            return True
    return False


def _artifacts_by_role(artifacts) -> dict[str, int]:
    """生成物を役割別に数える。役割名は handoff の `Artifact.role`。

    「何が何件出たか」は取りこぼしの検出に足りる（60 サンプルなら各役割 60 件）。
    個々のパスは include_artifacts に譲る。
    """
    counts: dict[str, int] = {}
    for artifact in artifacts:
        counts[artifact.role] = counts.get(artifact.role, 0) + 1
    return dict(sorted(counts.items()))


def _artifacts_tsv(artifacts) -> str:
    """生成物一覧を TSV（列名 1 回）で返す。生成物ゼロなら空文字。

    列名だけの行を返すと「1 件ある」と読めるため、空のときは何も返さない
    （件数は artifact_count が持つ）。
    """
    if not artifacts:
        return ""
    lines = ["path\trole\tformat\troot"]
    lines.extend(f"{a.path}\t{a.role}\t{a.format}\t{a.root}" for a in artifacts)
    return "\n".join(lines)


def _candidate_payload(candidate) -> dict:
    """MethodCandidate を UI が読める辞書にする。mtime は ISO 文字列で返す。

    `usable` は polarity を指定せずに探索したとき None になる（「一致するか」を
    判定しようがないため）。"direct" と取り違えられないよう、`key_params` と
    同じくキーごと省く。
    """
    from datetime import datetime
    payload = {
        "path": candidate.path,
        "origin": candidate.origin,
        "ion_mode": candidate.ion_mode,
        "omics": candidate.omics,
        "has_lbm": candidate.has_lbm,
        "mtime": datetime.fromtimestamp(candidate.mtime).isoformat(timespec="seconds"),
    }
    if candidate.usable is not None:
        payload["usable"] = candidate.usable
    if candidate.key_params is not None:
        payload["key_params"] = candidate.key_params
    return payload


def _candidates_next_hint(polarity: str | None) -> str:
    """console_method_candidates の `next` 文面。polarity 省略時は usable が
    候補に付かない（finding 2）ので、判断材料を ion_mode に差し替える。
    """
    if polarity is None:
        return ("polarity を省略したため候補に usable は付きません。"
                "ion_mode を見て、目的の極性と一致するものを "
                "console_plan(method_file=...) に、違うものは "
                "console_method_template(based_on=..., polarity=...) を"
                "通してから console_plan に渡してください。")
    return ("usable=direct なら console_plan(method_file=...)、"
            "usable=needs_polarity_conversion なら console_method_template("
            "based_on=..., polarity=...) を通してから console_plan")


def _capped_candidates(candidates) -> dict:
    """候補一覧を `MAX_REPORTED_CANDIDATES` で切る。

    候補は `discover_method_candidates` が direct 優先・新しい順にソート済みなので、
    先頭から切れば最も有用な候補が残る。切ったときだけ `truncated: true` を立て、
    切っていないときはキー自体を出さない。呼び出し側が別途持つ「総数」
    （`n_candidates` や封筒メッセージの件数）はここでは変えない —
    切った件数と混同させないため。
    """
    from lipidmix.console import method_file as method_file_mod

    emitted = candidates[:method_file_mod.MAX_REPORTED_CANDIDATES]
    section: dict = {"candidates": [_candidate_payload(c) for c in emitted]}
    if len(candidates) > method_file_mod.MAX_REPORTED_CANDIDATES:
        section["truncated"] = True
    return section


def _exe_error(exc: EnvironmentError) -> str:
    """`get_exe_path` の失敗を封筒にする。

    未設定（`MSDIAL_EXE_NOT_FOUND`）には人間がすべき手順を、設定ファイルが読めない
    （`CONFIG_INVALID`）ときはその行・列を載せる。
    """
    code = getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND")
    details = dict(getattr(exc, "details", {}) or {})
    if code == "MSDIAL_EXE_NOT_FOUND":
        details = {**_msdial_exe_setup_help(), **details}
    return console_error(code, str(exc), details or None)


def _lbm_setting():
    """`[msdial] lbm` / MSDIAL_LBM の設定。設定ファイルが読めなければ封筒（str）を返す。"""
    from lipidmix.core import user_config
    try:
        return user_config.get_setting("msdial.lbm")
    except user_config.ConfigInvalidError as exc:
        return console_error(exc.code, exc.message, exc.details())


def _configured_exe() -> str:
    """封筒に載せるための設定済み実行体（決められなければ空文字）。"""
    from lipidmix.console import runner as console_runner
    try:
        return console_runner.get_exe_path()
    except EnvironmentError:
        return ""


def _msdial_exe_setup_help() -> dict:
    """Console の実行体が未設定の封筒に、誰が何をすべきかを機械可読で載せる。

    設定ファイル（`lipidmix.local.toml`）の編集も環境変数の設定も MCP クライアント
    からはできない。LLM に「設定してください」とだけ返すと、設定を試みて失敗するか
    黙って諦める。人間の作業であることを型で示し、コピペできる手順を渡す。
    設定ファイルは呼ばれるたびに読まれるので再起動は要らないが、環境変数は
    起動中のサーバに反映されない。
    """
    from lipidmix.console.runner import msdial_exe_candidates
    from lipidmix.core import user_config

    try:
        candidates = msdial_exe_candidates()
    except OSError:
        candidates = []
    config = user_config.config_file_path()
    if config.is_file():
        config_step = (f"{config} の [msdial] に exe = '<MSDIALCUI.exe のパス>' と書く"
                       "（単一引用符で囲む。サーバの再起動は不要）")
    else:
        config_step = (f"{user_config.EXAMPLE_FILENAME} を {config} として複製し、"
                       "[msdial] exe = '<MSDIALCUI.exe のパス>' と書く"
                       "（単一引用符で囲む。サーバの再起動は不要）")
    return {
        "human_action_required": True,
        "why": "設定ファイルの編集も環境変数の設定も MCP クライアントからはできません。",
        "candidates": candidates[:10],
        "how_to_set": [
            config_step,
            "PowerShell（恒久設定）: [Environment]::SetEnvironmentVariable("
            "'MSDIAL_EXE','<MSDIALCUI.exe のパス>','User')（設定後に MCP サーバの再起動が必要）",
        ],
        "restart_required": False,
        "restart_note": "設定ファイルなら次の呼び出しから効きます。環境変数で設定した場合だけ、"
                        "起動中のプロセスが変更を読み直さないので MCP サーバの再起動が必要です。",
    }


def _execution_options_error(save_project: object, timeout_s: object) -> str | None:
    """新しい実行オプションを、ジョブ作成・実行の両入口で同じ規則で検査する。"""
    if type(save_project) is not bool:
        return console_error(
            "JOB_NOT_PLANNED",
            f"save_project は bool です: {save_project!r}",
        )
    if type(timeout_s) is not int or timeout_s <= 0:
        return console_error(
            "JOB_NOT_PLANNED",
            f"timeout_s は 0 より大きい int です: {timeout_s!r}",
        )
    return None


def _record_finalization_failure(job_path: Path, error: str) -> None:
    """最終化に失敗したジョブを、書き込める場合は failed として残す。"""
    from lipidmix.console.job_manager import update_status

    try:
        update_status(job_path, "failed", error=error)
    except Exception:
        # ロック等で failed を書けない場合も、MCP 境界から例外を漏らさない。
        pass


def _resolve_job_path(job_path: str | None) -> Path | str:
    """引数またはセッションから job_path を解決する。失敗はエラーエンベロープ文字列を返す。"""
    target = job_path or session_state.session.current_job_path
    if not target:
        return console_error(
            "JOB_NOT_FOUND",
            "job_path が指定されておらず、セッションにジョブもありません。"
            "先に console_plan を実行してください。",
            {"required_tools": ["console_plan"]},
        )
    return Path(target).expanduser()
