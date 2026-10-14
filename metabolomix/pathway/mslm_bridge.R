# metabolomix/pathway/mslm_bridge.R
# MSLipidMapper（R パッケージ）を呼ぶ薄いブリッジ。契約は bridge.py の docstring と計画書を参照。
# 引数: [1] renv プロジェクト（リポジトリの r/）、[2] request.json。結果は <out_dir>/bridge-result.json。
args <- commandArgs(trailingOnly = TRUE)
# 最初に renv を有効にする（R 本体には jsonlite も無い。パッケージは r/ の仮想環境だけから読む）。
local({
  owd <- setwd(args[1])
  on.exit(setwd(owd))
  source(file.path("renv", "activate.R"))
})
req <- jsonlite::fromJSON(args[2], simplifyVector = TRUE)
out_dir <- normalizePath(req$out_dir, winslash = "/", mustWork = TRUE)
result_path <- file.path(out_dir, "bridge-result.json")

write_result <- function(x) {
  jsonlite::write_json(x, result_path, auto_unbox = TRUE, null = "null", pretty = FALSE)
}
fail <- function(code, message) {
  write_result(list(status = "error", code = code, message = message))
  quit(status = 2, save = "no")
}

if (!requireNamespace("MSLipidMapper", quietly = TRUE)) {
  fail("MSLM_R_PACKAGES_MISSING",
       "R の仮想環境に MSLipidMapper がありません。リポジトリ直下で `Rscript r/setup.R` を実行してください。")
}
needed <- c("mztab_pathway_measurements", "render_pathway_from_config")
missing_fns <- needed[!needed %in% getNamespaceExports("MSLipidMapper")]
if (length(missing_fns)) {
  fail("MSLM_R_PACKAGES_MISSING",
       paste0("MSLipidMapper が古く、次の関数がありません: ", paste(missing_fns, collapse = ", ")))
}
version <- as.character(utils::packageVersion("MSLipidMapper"))
sha <- utils::packageDescription("MSLipidMapper")$RemoteSha
if (is.null(sha) || is.na(sha)) sha <- NULL
network_path <- function(p) {
  f <- system.file("extdata", "examples", req$networks[[p]], package = "MSLipidMapper")
  if (!nzchar(f)) fail("PATHWAY_BAD_REQUEST", paste0("同梱パスウェイが見つかりません: ", p))
  f
}

outputs <- list(); mapped <- list()
if (identical(req$mode, "check")) {
  write_result(list(status = "ok", code = NULL, message = "", package_version = version, package_sha = sha,
                    outputs = outputs, mapped_nodes = mapped))
  quit(status = 0, save = "no")
}

if (identical(req$mode, "measure")) {
  for (p in req$pathways) {
    out_file <- file.path(out_dir, paste0(p, "-measurements.json"))
    res <- tryCatch(
      MSLipidMapper::mztab_pathway_measurements(
        req$mztab, network_path(p), output_file = out_file, group_col = "class",
        acyl_chains = as.character(unlist(req$acyl_chains)), aggregate = req$aggregation,
        normalization = req$normalization, unit = req$unit,
        dataset_id = paste0(req$dataset_id, "-", p)),
      error = function(e) fail(
        if (grepl("acyl", conditionMessage(e), ignore.case = TRUE)) "PATHWAY_BAD_REQUEST" else "MSLM_RUN_FAILED",
        conditionMessage(e)))
    outputs[[p]] <- out_file
    nm <- res$node_mapping
    mapped[[p]] <- if (is.data.frame(nm) && "matched" %in% names(nm)) sum(nm$matched, na.rm = TRUE) else NA
  }
} else if (identical(req$mode, "render")) {
  cfg <- list(
    input = list(format = "mztab-m", file = req$mztab),
    analysis = list(normalization = "none", aggregation = req$aggregation, group_column = "class"),
    pathway = list(browser = req$browser, timeout_seconds = req$timeout_seconds),
    plot = list(type = req$plot_type, x_order = as.list(req$x_order)),
    output = list(paper_size = "A4", orientation = "landscape"))
  cfg_path <- file.path(out_dir, "render-config.yml")
  yaml::write_yaml(cfg, cfg_path)
  for (p in req$pathways) {
    pdf <- file.path(out_dir, paste0(p, ".pdf"))
    res <- tryCatch(
      MSLipidMapper::render_pathway_from_config(cfg_path, network_path = network_path(p),
                                                output_path = pdf),
      error = function(e) fail("PATHWAY_RENDER_FAILED", conditionMessage(e)))
    outputs[[p]] <- pdf
    mapped[[p]] <- unname(unlist(res$mapped_nodes))[1]
  }
} else {
  fail("PATHWAY_BAD_REQUEST", paste0("未知の mode: ", req$mode))
}
write_result(list(status = "ok", code = NULL, message = "", package_version = version, package_sha = sha,
                  outputs = outputs, mapped_nodes = mapped))
