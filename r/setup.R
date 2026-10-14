# R の仮想環境（renv）を r/renv.lock から復元する。
# リポジトリ直下で:  Rscript r/setup.R
# R 本体のライブラリには何も入れない（renv 自身も r/ の中に自動で取得される）。
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
here <- normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/")
setwd(here)
source(file.path("renv", "activate.R"))
renv::restore(prompt = FALSE)
cat("MSLipidMapper", as.character(utils::packageVersion("MSLipidMapper")),
    "ready in", here, "\n")
