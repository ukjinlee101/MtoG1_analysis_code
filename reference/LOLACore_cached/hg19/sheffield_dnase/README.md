To remove the first line of each file

for FILE in `ls *.bed`; do
	echo $FILE
	tail -n +2 "$FILE" > "$FILE.temp"
	rm $FILE
	mv $FILE.temp $FILE
done;

The first line was:

chr	start	stop	original_cluster	refined_cluster	original_distance	refined_distance

I can probably take all these other columns out completely to save space...
