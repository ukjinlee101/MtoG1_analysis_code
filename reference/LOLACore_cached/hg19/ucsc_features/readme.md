Some notes on how I grabbed these tables from UCSC:

$dest = "../data/ucsc";
`mkdir $dest`;
`wget --directory-prefix=$dest http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz`;

* Grab a list of available UCSC tables:
```
echo "show tables;" | mysql --user=genome --host=genome-mysql.cse.ucsc.edu -D hg19 -A > ucsc.tables.txt
```

* Select the tables you want to grab; for example:

nestedRepeats
simpleRepeat
rmsk
microsat
genomicSuperDups
chainSelf
cpgIslandExt


* Then you can grab the columns for the tables (might have to check them, as they differ by table), with a command like this:

```
mysql --user=genome --host=genome-mysql.cse.ucsc.edu hg19 -A -sre 'SELECT chrom, start, end from nestedRepeats' > file.bed
```

* Now I created an index file with `ls *.bed > index.txt`

And manually named the description columns.

