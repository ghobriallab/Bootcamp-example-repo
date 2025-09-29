#! /bin/bash

input_dir=${1}

## Convert bed bim fam to vcf
mapfile -t chromosomes < <(ls   ${input_dir}| awk -F. '{print $1"."$2}' | uniq)

for chromosome in "${chromosomes[@]}";do
    out_chr=$(echo $chromosome | awk -F'_' '{print $2}')
    echo "Processing file "${chromosome}
    plink --bfile  ${input_dir}"/"${chromosome} --recode vcf --out  ${input_dir}"/h19_"${out_chr}
done


# move bed bim fam  and logs on separate directories
bed_dir=${input_dir}/beds/
log_dir=${input_dir}/logs/



if [ ! -d  $bed_dir ]; then
    mkdir $bed_dir
fi

mv ${input_dir}/*bed $bed_dir
mv ${input_dir}/*bim $bed_dir
mv ${input_dir}/*fam $bed_dir

if [ ! -d  $log_dir ]; then
    mkdir $log_dir
fi

mv ${input_dir}/*log $log_dir

## concat vcf files

bcftools concat ${input_dir}/*vcf -o ${input_dir}/h19_snps.vcf


## Convert from h19 to h38

picard LiftoverVcf -I ${input_dir}/h19_snps.vcf \
                    -O ${input_dir}/h38_snps.vcf \
                    -C ${chain} \
                    -R ${reference} \
                    -REJECT ${input_dir}/h38_snps_rejected.vcf 



# move h19 files into a separate directory
h19_dir=${input_dir}/h19/

if [ ! -d  $h19_dir ]; then
    mkdir $h19_dir
fi

mv ${input_dir}/h19_* $h19_dir
