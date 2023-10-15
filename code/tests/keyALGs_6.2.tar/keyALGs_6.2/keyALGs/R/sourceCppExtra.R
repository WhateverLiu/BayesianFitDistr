

# MUST have installed Rcpp.
# Codes wraped by ============= are modifications.
sourceCppExtra = function (file = "", code = NULL, env = globalenv(), embeddedR = TRUE,
                           rebuild = FALSE, cacheDir = getOption("rcpp.cache.dir", tempdir()),
                           cleanupCacheDir = FALSE, showOutput = verbose,
                           verbose = TRUE, dryRun = FALSE, windowsDebugDLL = FALSE,
                           echo = TRUE, externObjLibFiles = NULL, retrieveSLL = FALSE,
                           dllsNeeded = NULL)
{
  # ============================================================================
  if (!is.null(externObjLibFiles))
    externObjLibFiles = normalizePath(externObjLibFiles, winslash = "/")
  if (!is.null(dllsNeeded))
    dllsNeeded = normalizePath(dllsNeeded, winslash = "/")
  # ============================================================================
  
  
  cacheDir <- path.expand(cacheDir)
  cacheDir <- Rcpp:::.sourceCppPlatformCacheDir(cacheDir)
  cacheDir <- normalizePath(cacheDir)
  if (!missing(code)) {
    rWorkingDir <- getwd()
    file <- tempfile(fileext = ".cpp", tmpdir = cacheDir)
    con <- file(file, open = "w")
    writeLines(code, con)
    close(con)
  }
  else {
    rWorkingDir <- normalizePath(dirname(file))
  }
  file <- normalizePath(file, winslash = "/")
  if (!tools::file_ext(file) %in% c("cc", "cpp")) {
    stop("The filename '", basename(file), "' does not have an ",
         "extension of .cc or .cpp so cannot be compiled.")
  }
  if (.Platform$OS.type == "windows") {
    if (grepl(" ", basename(file), fixed = TRUE)) {
      stop("The filename '", basename(file), "' contains spaces. This ",
           "is not permitted.")
    }
  }
  else {
    if (windowsDebugDLL) {
      if (verbose) {
        message("The 'windowsDebugDLL' toggle is ignored on ",
                "non-Windows platforms.")
      }
      windowsDebugDLL <- FALSE
    }
  }
  context <- .Call("sourceCppContext", PACKAGE = "Rcpp",
                   file, code, rebuild, cacheDir, .Platform)
  if (context$buildRequired || rebuild) {
    if (verbose)
      Rcpp:::.printVerboseOutput(context)
    succeeded <- FALSE
    output <- NULL
    depends <- Rcpp:::.getSourceCppDependencies(context$depends,
                                                file)
    Rcpp:::.validatePackages(depends, context$cppSourceFilename)
    envRestore <- Rcpp:::.setupBuildEnvironment(depends, context$plugins,
                                                file)
    cwd <- getwd()
    setwd(context$buildDirectory)
    fromCode <- !missing(code)
    if (!Rcpp:::.callBuildHook(context$cppSourcePath, fromCode,
                               showOutput)) {
      Rcpp:::.restoreEnvironment(envRestore)
      setwd(cwd)
      return(invisible(NULL))
    }
    on.exit({
      if (!succeeded) Rcpp:::.showBuildFailureDiagnostics()
      Rcpp:::.callBuildCompleteHook(succeeded, output)
      setwd(cwd)
      Rcpp:::.restoreEnvironment(envRestore)
    })
    if (file.exists(context$previousDynlibPath)) {
      try(silent = TRUE, dyn.unload(context$previousDynlibPath))
      file.remove(context$previousDynlibPath)
    }
    r <- paste(R.home("bin"), "R", sep = .Platform$file.sep)
    lib <- context$dynlibFilename
    deps <- context$cppDependencySourcePaths
    src <- context$cppSourceFilename
    args <- c(r, "CMD", "SHLIB", if (windowsDebugDLL) "-d",
              if (rebuild) "--preclean", if (dryRun) "--dry-run",
              "-o", shQuote(lib), if (length(deps))
                paste(shQuote(deps), collapse = " "), shQuote(src))


    # ==========================================================================
    if (!is.null(externObjLibFiles))
    {
      file.copy(from = externObjLibFiles, to = getwd(), overwrite = TRUE)
      externObjLibFiles = paste0('"', paste0(
        basename(externObjLibFiles), collapse = '" "'), '"')
      args = args[args != "--preclean"]
      args = c(args, externObjLibFiles)
    }
    
    
    if (!is.null(dllsNeeded))
    {
      for (x in dllsNeeded) dyn.load(x)
      tmp = sapply(dllsNeeded, function(x)
      {
        fname = basename(x)
        fname = tools::file_path_sans_ext(fname)
        dname = dirname(x)
        paste0("-L", dname, " -l", fname)
      })
      args = args[args != "--preclean"]
      args = c(args, tmp)
    }
    # ==========================================================================


    cmd <- paste(args, collapse = " ")


    if (showOutput) cat(cmd, "\n")

    result <- suppressWarnings(system(cmd, intern = !showOutput))


    # ==========================================================================
    if(retrieveSLL)
    {
      allFilesInBuild = list.files(getwd(), full.names = T)
      allFilesInBuild = allFilesInBuild[tools::file_ext(allFilesInBuild) == "o"]
      file.copy(allFilesInBuild, to = rWorkingDir, overwrite = TRUE)
    }
    # ==========================================================================


    if (!showOutput) {
      output <- result
      attributes(output) <- NULL
      status <- attr(result, "status")
      if (!is.null(status)) {
        cat(result, sep = "\n")
        succeeded <- FALSE
        stop("Error ", status, " occurred building shared library.")
      }
      else if (!file.exists(context$dynlibFilename)) {
        cat(result, sep = "\n")
        succeeded <- FALSE
        stop("Error occurred building shared library.")
      }
      else {
        succeeded <- TRUE
      }
    }
    else if (!identical(as.character(result), "0")) {
      succeeded <- FALSE
      stop("Error ", result, " occurred building shared library.")
    }
    else {
      succeeded <- TRUE
    }
  }
  else {
    cwd <- getwd()
    on.exit({
      setwd(cwd)
    })
    if (verbose)
      cat("\nNo rebuild required (use rebuild = TRUE to ",
          "force a rebuild)\n\n", sep = "")
  }
  if (dryRun)
    return(invisible(NULL))
  if (length(context$exportedFunctions) > 0 || length(context$modules) >
      0) {
    exports <- c(context$exportedFunctions, context$modules)
    removeObjs <- exports[exports %in% ls(envir = env, all.names = T)]
    remove(list = removeObjs, envir = env)
    scriptPath <- file.path(context$buildDirectory, context$rSourceFilename)
    source(scriptPath, local = env)
  }
  else if (getOption("rcpp.warnNoExports", default = TRUE)) {
    warning("No Rcpp::export attributes or RCPP_MODULE declarations ",
            "found in source")
  }
  if (embeddedR && (length(context$embeddedR) > 0)) {
    srcConn <- textConnection(context$embeddedR)
    setwd(rWorkingDir)
    source(file = srcConn, local = env, echo = echo)
  }
  if (cleanupCacheDir)
    Rcpp:::cleanupSourceCppCache(cacheDir, context$cppSourcePath,
                                 context$buildDirectory)
  invisible(list(functions = context$exportedFunctions, modules = context$modules,
                 cppSourcePath = context$cppSourcePath, buildDirectory = context$buildDirectory))
}




# Test it
if(F)
{

  dir.create("C:/Users/i56087/Desktop/tmp/helpVictor/tmp", showWarnings = F)


  source("R/sourceCppExtra.R")
  # test.o will be copied back to externLib.
  sourceCppExtra("C:/Users/i56087/Desktop/tmp/helpVictor/btConv.cpp", rebuild = T, verbose = T, retrieveSLL = TRUE, cacheDir = "C:/Users/i56087/Desktop/tmp/helpVictor/tmp")
  # There should already have been a .h file included in testLinkLib.cpp. This .h file
  # declare functions declared in test.cpp
  sourceCppExtra("testLinkLib.cpp", verbose = T, rebuild = T, externObjLibFiles = "externLib/test.o")
  x = runif(10); y = runif(10)
  range(timesTwo(x, y) - (x + y) * x)




}










