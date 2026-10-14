# R の仮想環境（renv）を r/renv.lock から復元する。
# リポジトリ直下で:  Rscript r/setup.R
# R 本体のライブラリには何も入れない（renv 自身も r/ の中に自動で取得される）。
# 更新で r/renv.lock が変わったら、もう一度実行する（ブリッジが目印の食い違いで止める）。
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
here <- normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/")
setwd(here)
source(file.path("renv", "activate.R"))
renv::restore(prompt = FALSE)
# 復元に使った lock の md5 を目印に残す（metabolomix/pathway/bridge.py の RESTORE_STAMP）。
# library/ ごと git の追跡外。restore が失敗したらここまで来ないので、目印は成功時だけ書かれる。
stamp <- file.path("renv", "library", "restored-lock.md5")
dir.create(dirname(stamp), recursive = TRUE, showWarnings = FALSE)
writeLines(unname(tools::md5sum("renv.lock")), stamp)
cat("MSLipidMapper", as.character(utils::packageVersion("MSLipidMapper")),
    "ready in", here, "\n")
