


# Rcpp::sourceCpp("src/cppfuns.cpp", verbose = T)
# Copy files to linux directory. File names preserved.
win2lnxZipCopy = function(winFiles, lnxDir, winLnxPathAdd = "//RSGrid", linuxGrid = "rscgrid1finance1", passwdPath = "data/passwd.Rdata", removeZipOnWin = TRUE, removeZipOnLnx = TRUE)
{
  winFileName = paste0("winFile", sample(2e9L, 1), ".zip")
  tmp = zip(winFileName, winFiles)
  lnxDirWinview = lnxDir
  lnxDir = gsub(winLnxPathAdd, "", lnxDir)
  cmd = paste0("powershell Copy-Item '", winFileName, "' -Destination '", lnxDirWinview, "'")
  system(cmd)
  if(removeZipOnWin) unlink(winFileName, recursive = T)
  load(passwdPath)
  session = ssh::ssh_connect(linuxGrid, passwd = passwd)
  rm(passwd); gc()


  if(removeZipOnLnx) qsub = paste0("cd ", lnxDir, "; unzip -o ", winFileName, "; unlink ", winFileName)
  else qsub = paste0("cd ", lnxDir, "; unzip -o ", winFileName)
  ssh::ssh_exec_internal(session, command = qsub)
  ssh::ssh_disconnect(session)
}




lnx2winZipCopy = function(lnxFiles, winDir, winLnxPathAdd = "//RSGrid", linuxGrid = "rscgrid1finance1", passwdPath = "data/passwd.Rdata", removeZipOnWin = TRUE, removeZipOnLnx = TRUE)
{


  tmp = strsplit(lnxFiles[1], split = "/")[[1]]
  lnxDirWinView = paste0(tmp[-length(tmp)], collapse = "/")
  lnxDir = gsub(winLnxPathAdd, "", lnxDirWinView)
  lnxFiles = gsub(winLnxPathAdd, "", lnxFiles)
  lnxFiles = gsub(paste0(lnxDir, "/"), "", lnxFiles)


  load(passwdPath)
  session = ssh::ssh_connect(linuxGrid, passwd = passwd)
  rm(passwd); gc()
  lnxFilezipName = paste0("lnxFile", sample(2e9L, 1), ".zip")
  qsub = paste0("cd ", lnxDir, "; zip -r ", lnxFilezipName, " ", paste0(lnxFiles, collapse = " "))
  ssh::ssh_exec_internal(session, command = qsub)


  lnxFilezipPath = paste0(winLnxPathAdd, lnxDir, "/", lnxFilezipName)
  cmd = paste0("powershell Copy-Item '", lnxFilezipPath, "' -Destination '", winDir, "'")
  system(cmd)


  if(removeZipOnLnx) qsub = paste0("cd ", lnxDir, "; unlink ", lnxFilezipName)
  ssh::ssh_exec_internal(session, command = qsub)
  ssh::ssh_disconnect(session)


  zipFileOnWin = paste0(winDir, "/", lnxFilezipName)
  unzip(zipFileOnWin, exdir = winDir, junkpaths = F)
  if(removeZipOnWin) unlink(zipFileOnWin)
}



















