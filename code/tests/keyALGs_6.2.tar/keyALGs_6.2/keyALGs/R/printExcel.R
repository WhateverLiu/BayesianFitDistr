

#'
#' Print Excel sheet.
#' 
#' Print a rectangle block of cells in a sheet in an Excel book.
#' 
#' @param excelBookFullPath Full path to the Excel book. Relative path will not
#' work since Excel does not recognize your current Python working directory. 
#' 
#' @param sheetName A character string. Name of the sheet in the Excel book.
#' 
#' @param printArea A character string in the format of "A1:B2" to identify the
#' top-left and bottom-right cells.
#' 
#' @param pdfFullPath Full path to the printed pdf.
#' 
#' @return Nothing.
#' 
#' @details The function uses \code{\link[reticulate]{source_python}} to 
#' source "excelPrint.py" in /inst/python.
#' 
#'
printExcel <- function(excelBookFullPath, sheetNames, 
                       printArea = 'A1:B2', pdfFullPath = "sample.pdf")
{
  pypath = paste0(find.package('keyALGs'), "/python/excelPrint.py")
  reticulate::source_python(pypath)
  printExcelSheet(excelBookFullPath, sheetNames, printArea, pdfFullPath)
}












