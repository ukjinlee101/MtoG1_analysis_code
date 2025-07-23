library(stringr)

parseUCSCfiles = function(file) {
  #Load the files.txt file from ENCODE that annotates the
  #ChIP-seq experiments
  encodeTFBSannotation = fread(file, header=FALSE)
  #Parse the annotation into data.table format
  setnames(encodeTFBSannotation, "V1", "filename")
  ct = str_match(encodeTFBSannotation$V2, "cell=(.*?);")[, 2]
  tr = str_match(encodeTFBSannotation$V2, "treatment=(.*?);")[, 2]
  ab = str_match(encodeTFBSannotation$V2, "antibody=(.*?);")[, 2]
  encodeTFBSannotation
  encodeTFBSannotation[, V2:=NULL]
  encodeTFBSannotation[, cellType:=ct]
  encodeTFBSannotation[, treatment:=tr]
  encodeTFBSannotation[, antibody:=ab]
  encodeTFBSannotation[, description:=paste("ChIP", ct, ab)]
  encodeTFBSannotation[, filename:=sub(".gz", "", filename)]
  return(encodeTFBSannotation);
}
setwd("/data/groups/lab_bock/shared/resources/regions/LOLACore/mm10/encodeTFBSmm10/")
file = "files.txt"
anno = parseUCSCfiles(file)
write.table(anno, "index.txt",sep="\t",quote=FALSE,row.names=FALSE)
