These files were provided to me by Rebecca Hannah from CODEX:

wget http://hscl.cimr.cam.ac.uk/Files_for_transfer/Nathan/mm10_Archive.zip
wget http://hscl.cimr.cam.ac.uk/Files_for_transfer/Nathan/hg19_Archive.zip

After unarchiving them, I manually renamed the excel sheet she sent to get the index.txt file.

I had to manually remove the last line of this file:
GSM1122302_CBFB.bed

Because it was truncated:
hrM	0	400


I also did a search/replace for "N/A" to "NA" so that R will correctly recognize these as NAs, and provide generic descriptions if applicable.
