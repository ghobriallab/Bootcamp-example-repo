#!/usr/bin/bash

gsutil ls -d gs://ukbb-pilot-april2024/cellranger/*/*/*TCR* | grep UKBBP >list.txt
readarray -t myArray < <(cat list.txt)
printf '%s\n' "${myArray[@]}"
for dir in "${myArray[@]}";
do
	var2="$(basename $dir)" #gets on the innermost directory name
	var3=${dir%/} #removes the last forward slash
	mkdir $var2
	mkdir $var2/outs

        gsutil -m cp -r  $var3/outs/*.csv $var2/outs/ ;
done