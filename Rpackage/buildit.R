
pkgpath = "/finance_develop/Charlie/BayesianFitDistr/Rpackage"
# if ( .Platform$OS.type != "unix") pkgpath = 


setwd(pkgpath)
doc = readLines(con = "../slides/slides.tex")
removeBegin = which(
  doc == "% ---- first slide is here. Do not delete or change this line (set for making R package)")
removeEnd = which(doc == "% ---- New slides begin for R package")
tmp = "
\\begin{frame}
\\centering\\Huge Graphics in the slide tutorial are from analyses on the New Zealand earthquake data.
\\end{frame}
"
doc = c(doc[1:removeBegin], tmp, doc[removeEnd:length(doc)])
setwd('../tempFiles')
fs = list.files()
which2delete = which(tools::file_path_sans_ext(fs) == "slides-1159912635")
unlink( fs[which2delete], recursive = T )
writeLines(doc, con = "slides-1159912635.tex")
system("/home/i56087/TeXLive/bin/x86_64-linux/pdflatex -halt-on-error slides-1159912635.tex")
system("/home/i56087/TeXLive/bin/x86_64-linux/pdflatex -halt-on-error slides-1159912635.tex")
file.copy("slides-1159912635.pdf", "../Rpackage/vignettes/slides.pdf",
          overwrite = T)
setwd(pkgpath)


system("mv descAndNamespace/DESCRIPTION DESCRIPTION")
system("mv descAndNamespace/NAMESPACE NAMESPACE")


# ==============================================================================
# Remove all the object files in the source code.
# ==============================================================================


# Change all hpp to cpp
if (F)
{
  sfs = list.files("src", full.names = T)
  dir.create("src/obsoleteHppFiles", showWarnings = F)
  file.copy(sfs[tools::file_ext(sfs) == "hpp"], to = "obsoleteHppFiles")
  file.rename( sfs[tools::file_ext(sfs) == "hpp"], 
               paste0(tools::file_path_sans_ext(
                 sfs[tools::file_ext(sfs) == "hpp"]), ".cpp"))
}


# ==============================================================================
# Remove all the object files in the source code and rebuild
# ==============================================================================
setwd('/finance_develop/Charlie/BayesianFitDistr/Rpackage')
remove.packages('NGFMfitDistr')
try(system("sudo rm -r /home/i56087/R/x86_64-redhat-linux-gnu-library/4.3/NGFMfitDistr"))
sfs = list.files("src", full.names = T)
unlink( c(sfs[tools::file_ext(sfs) %in% c("o", "so", "dll", "rds")], 
          "src/RcppExports.cpp") )
Rcpp::compileAttributes('./')
roxygen2::roxygenize()


sfs = list.files("src", full.names = T)
unlink( c(sfs[tools::file_ext(sfs) %in% c("o", "so", "dll", "rds")], 
          "src/RcppExports.cpp") )
Rcpp::compileAttributes('./')
setwd("../")
system("R CMD build Rpackage")
# system("R CMD check NGFMfitDistr_1.0.5.tar.gz")
system("R CMD INSTALL NGFMfitDistr_1.0.5.tar.gz")
setwd("Rpackage")


system("R CMD INSTALL --build NGFMfitDistr_1.0.5.tar.gz")

# For windows
system("R CMD INSTALL --build NGFMfitDistr_1.0.5.tar.gz")












