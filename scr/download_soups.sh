#! /usr/bin/bash

gsutil ls -d gs://ukbb-pilot-april2024/2_mismatch/souporcell/output/UKBBP* > soups_list.txt
readarray -t myArray < <(cat soups_list.txt)
printf '%s\n' "${myArray[@]}"
for dir in "${myArray[@]}";
do
        var2="$(basename $dir)" #gets on the innermost directory name
        var3=${dir%/} #removes the last forward slash
        # mkdir $var2

        gsutil -m cp -r  $var3 . ;
done
