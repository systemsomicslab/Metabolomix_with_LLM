# limma の参照値（tests/stats_reference/limma_expected.json）を作る。一度だけ実行する。
#
# 使い方（リポジトリのルートで）:
#   Rscript tests/stats_reference/make_limma_expected.R <limma と jsonlite が入ったライブラリ>
# または環境変数 LIMMA_R_LIB にライブラリの場所を置いて引数なしで実行する。
# ライブラリは一時の場所（スクラッチ）に BiocManager::install("limma") したものを指す。
# R 本体のライブラリにも r/ の renv にも入れない。pytest はこの JSON だけを読み、R を要らない。
#
# 入力: limma_input.tsv（特徴 30 x 試料 12。log2 の値。NA は欠測）と limma_design.tsv
# （sample, group, batch, subject）。入力は numpy の default_rng(20261107) で一度だけ作り、
# 小数 6 桁の固定値で保存したもの（6 特徴に 1〜2 個の NA、f30 は分散がほぼ 0）。
#
# 4 つの設計を当てる（どれも計画行列は列フルランク）:
#   cell_means           ~0+group、対比 B-A・C-A（全 30 特徴。欠測で残差の自由度が揃わない）
#   cell_means_complete  ~0+group、対比 B-A・C-A（欠測の無い 24 特徴。自由度が揃う＝fitFDist の経路）
#   batch_covariate      ~0+group+batch、対比 B-A
#   paired_block         B と C の 8 試料、~0+group+subject、対比 C-B
#   missing_cell         cell_means と同じ ~0+group、対比 B-A・C-A。ただし f05 の C を全部 NA、
#                        f12 の C2〜C4 を NA にする（na_cells に記録）。C の係数が推定できない特徴・
#                        1 試料だけから推定される特徴が当てはめと事前分布に入る経路（lm.series の NA 係数）を縛る
# 各設計で lmFit -> contrasts.fit -> eBayes（trend=FALSE, robust=FALSE）を回し、
# topTable(coef=, confint=0.95, number=Inf, sort.by="none") の logFC, CI.L, CI.R, t, P.Value と
# fit の df.total, s2.prior, df.prior, s2.post, F, F.p.value などを書く。

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else Sys.getenv("LIMMA_R_LIB", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "tests/stats_reference"
})

y_all <- as.matrix(read.delim(file.path(here, "limma_input.tsv"), row.names = 1,
                              na.strings = "NA", check.names = FALSE))
info <- read.delim(file.path(here, "limma_design.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(y_all), info$sample))

run_design <- function(samples, features, formula, contrast_list, na_cells = list()) {
  d <- info[match(samples, info$sample), , drop = FALSE]
  d$group <- factor(d$group)
  d$batch <- factor(d$batch)
  d$subject <- factor(d$subject)
  X <- model.matrix(formula, data = d)
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  stopifnot(is.fullrank(X))
  y <- y_all[features, samples, drop = FALSE]
  for (f in names(na_cells)) y[f, na_cells[[f]]] <- NA
  C <- sapply(contrast_list, function(v) v[colnames(X)])
  C <- matrix(C, nrow = ncol(X), dimnames = list(colnames(X), names(contrast_list)))
  fit <- lmFit(y, X)
  cfit <- contrasts.fit(fit, C)
  eb <- eBayes(cfit)
  per_contrast <- list()
  for (k in colnames(C)) {
    tt <- topTable(eb, coef = k, confint = 0.95, number = Inf, sort.by = "none")
    stopifnot(identical(rownames(tt), features))
    per_contrast[[k]] <- list(logFC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R,
                              t = tt$t, P.Value = tt$P.Value)
  }
  list(
    samples = samples,
    features = features,
    column_names = colnames(X),
    design = unname(X),
    contrast_names = I(colnames(C)),
    contrasts = unname(C),
    df_residual = unname(eb$df.residual),
    sigma = unname(eb$sigma),
    stdev_unscaled = unname(eb$stdev.unscaled),
    s2_prior = eb$s2.prior,
    df_prior = eb$df.prior,
    s2_post = unname(eb$s2.post),
    df_total = unname(eb$df.total),
    F = unname(eb$F),
    F_p_value = unname(eb$F.p.value),
    contrast = per_contrast,
    na_cells = lapply(na_cells, I)
  )
}

all_samples <- info$sample
all_features <- rownames(y_all)
complete <- all_features[rowSums(is.na(y_all)) == 0]
bc <- info$sample[info$group %in% c("B", "C")]

out <- list(
  `_meta` = list(
    r_version = R.version.string,
    limma_version = as.character(packageVersion("limma")),
    generator = "tests/stats_reference/make_limma_expected.R",
    note = "eBayes(trend=FALSE, robust=FALSE); topTable(confint=0.95, sort.by='none')"
  ),
  cell_means = run_design(all_samples, all_features, ~0 + group,
                          list(BvsA = c(groupA = -1, groupB = 1, groupC = 0),
                               CvsA = c(groupA = -1, groupB = 0, groupC = 1))),
  cell_means_complete = run_design(all_samples, complete, ~0 + group,
                                   list(BvsA = c(groupA = -1, groupB = 1, groupC = 0),
                                        CvsA = c(groupA = -1, groupB = 0, groupC = 1))),
  batch_covariate = run_design(all_samples, all_features, ~0 + group + batch,
                               list(BvsA = c(groupA = -1, groupB = 1, groupC = 0, batchb2 = 0))),
  paired_block = run_design(bc, all_features, ~0 + group + subject,
                            list(CvsB = c(groupB = -1, groupC = 1, subjectS2 = 0,
                                          subjectS3 = 0, subjectS4 = 0))),
  missing_cell = run_design(all_samples, all_features, ~0 + group,
                            list(BvsA = c(groupA = -1, groupB = 1, groupC = 0),
                                 CvsA = c(groupA = -1, groupB = 0, groupC = 1)),
                            na_cells = list(f05 = c("C1", "C2", "C3", "C4"),
                                            f12 = c("C2", "C3", "C4")))
)

writeLines(toJSON(out, auto_unbox = TRUE, digits = NA, matrix = "rowmajor", pretty = FALSE),
           file.path(here, "limma_expected.json"))
cat("wrote", file.path(here, "limma_expected.json"), "\n")
