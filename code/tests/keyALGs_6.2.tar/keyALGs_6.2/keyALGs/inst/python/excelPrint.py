

def printExcelSheet(excelBookFullPath, sheetNames, printArea = 'A1:B2', 
                    pdfFullPath = "excelPrint.pdf"):
  from win32com import client
  import os
  if os.path.exists(pdfFullPath): os.remove(pdfFullPath)
  o = client.Dispatch("Excel.Application")
  o.Visible = False
  wb = o.Workbooks.Open(excelBookFullPath, ReadOnly = True)
  allSheetNames = dict([(x.Name, i + 1) for i, x in enumerate(wb.Worksheets)])
  sheetInds = [allSheetNames[x] for x in sheetNames]
  wb.WorkSheets(sheetInds).Select()
  wb.ActiveSheet.PageSetup.PrintArea = printArea
  wb.ActiveSheet.PageSetup.LeftMargin = 0
  wb.ActiveSheet.PageSetup.RightMargin = 0
  wb.ActiveSheet.PageSetup.TopMargin = 0
  wb.ActiveSheet.PageSetup.BottomMargin = 0
  wb.ActiveSheet.PageSetup.HeaderMargin = 0
  wb.ActiveSheet.PageSetup.FooterMargin = 0
  # wb.ActiveSheet.PageSetup.Zoom = False
  # wb.ActiveSheet.PageSetup.FitToPagesWide = 1
  # wb.ActiveSheet.PageSetup.FitToPagesTall = False
  # wb.ActiveSheet.PageSetup.Orientation = 2
  wb.ActiveSheet.ExportAsFixedFormat(
    Type = 0, Filename = pdfFullPath, Quality = 1, 
    IncludeDocProperties = 0, IgnorePrintAreas = 0)
  wb.Close(True)
  

if True:
  printExcelSheet(
    "C:/Users/i56087/Desktop/tsunamiML/Charlie/figure/figure001.xlsx",
    ["Sheet17", "Sheet18"], "A1:HU25", 
    "C:/Users/i56087/Desktop/tsunamiML/Charlie/figure/explainUnet.pdf")

  
if False:
  printExcelSheet(
    "C:/Users/i56087/Desktop/grossLoss/grosslossExample20210216/figure/figureVSexplain003.xlsx",
    ["VS06", "VS07", "VS08"], "A1:AM22", 
    "C:/Users/i56087/Desktop/grossLoss/grosslossExample20210216/figure/tmp.pdf")





