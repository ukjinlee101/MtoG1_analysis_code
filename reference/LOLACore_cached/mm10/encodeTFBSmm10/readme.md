I downloaded all the files from UCSC, then used my R function parseUCSCfiles () in funcPublicData.R to parse the files.txt file into an index.txt file.[NS]
JK: pulled that function into new script parseUCSCfiles.R


```{r}
library(stringr)
utility("funcPublicData.R")
file = "/data/groups/lab_bock/shared/resources/regions/LOLACore/mm10/encodeTFBSmm10/files.txt"
anno = parseUCSCfiles(file)
write.tsv (anno, "index.txt")
```
