hicdc2hic_custom <- function (gi_list, hicfile, mode = "normcounts", chrs = NULL, 
                              gen_ver = "hg19", memory = 8) 
{
  options(scipen = 9999)
  if (.Platform$OS.type == "windows" & Sys.getenv("R_ARCH") == 
      "/i386") {
    gc(reset = TRUE, full = TRUE)
    utils::memory.limit(size = 4095)
  }
  HiCDCPlus::gi_list_validate(gi_list)
  binsize <- HiCDCPlus::gi_list_binsize_detect(gi_list)
  if (is.null(chrs)) 
    chrs <- names(gi_list)
  tmpfile <- paste0(base::tempfile(), ".txt")
  HiCDCPlus::gi_list_write(gi_list, tmpfile, columns = "minimal_plus_score", 
                score = mode)
  hicdc2hicoutput <- path.expand(hicfile)
  hicdc2hicoutputdir <- gsub("/[^/]+$", "", hicdc2hicoutput)
  if (hicdc2hicoutputdir == hicdc2hicoutput) {
    hicdc2hicoutputdir <- gsub("\\[^\\]+$", "", hicdc2hicoutput)
  }
  if (hicdc2hicoutputdir == hicdc2hicoutput) {
    hicdc2hicoutputdir <- gsub("\\\\[^\\\\]+$", "", hicdc2hicoutput)
  }
  if (!hicdc2hicoutputdir == hicdc2hicoutput & !dir.exists(hicdc2hicoutputdir)) {
    dir.create(hicdc2hicoutputdir, showWarnings = FALSE, 
               recursive = TRUE, mode = "0777")
  }
  jarpath<-"/athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar"
  ifelse(.Platform$OS.type == "windows" & Sys.getenv("R_ARCH") == 
           "/i386", min(memory, 2), memory)
  if (mode == "zvalue") {
    system2("java", args = c(paste0("-Xmx", as.character(memory), 
                                    "g"), "-jar", path.expand(jarpath), "pre", "-v", 
                             "-d", "-n", "-r", binsize, "-m", -2147400000, path.expand(tmpfile), 
                             path.expand(hicdc2hicoutput), gen_ver))
  }
  else {
    system2("java", args = c(paste0("-Xmx", as.character(memory), 
                                    "g"), "-jar", path.expand(jarpath), "pre", "-v", 
                             "-d", "-r", binsize, path.expand(tmpfile), path.expand(hicdc2hicoutput), 
                             gen_ver))
  }
  system2("rm", args = path.expand(tmpfile))
  return(hicdc2hicoutput)
}