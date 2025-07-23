I downloaded all the files from UCSC, then used my R function parseUCSCfiles () in funcPublicData.R to parse the files.txt file into an index.txt file.

```{r}
library(stringr)
utility("funcPublicData.R")
file = "~/fhgfs/share/data/encodeTFBS/files.txt"
anno = parseUCSCfiles(file)
write.tsv (anno, "index.txt")
```
