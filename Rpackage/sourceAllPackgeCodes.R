
# ==============================================================================
# Compile every cpp files in /src and every R files in R.
# ==============================================================================
tmp = list.files('src')[tools::file_ext( list.files('src')) == "cpp"]
tmp = paste0("src/", tmp[tmp != "RcppExports.cpp"])
dir.create("../tempFiles", showWarnings = FALSE)
for (x in tmp) Rcpp::sourceCpp(x, cacheDir = "../tempFiles", verbose = 1)
tmp = list.files('R')[tools::file_ext( list.files('R')) == "R"]
tmp = paste0("R/", tmp[tmp != "RcppExports.R"])
for (x in tmp) source(x)

