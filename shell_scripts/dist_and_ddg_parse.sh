#!/bin/bash

# in
known_ddg_file=~/Downloads/known_mutations.mcsm.parsed.txt
unknown_ddg_file=~/Downloads/unknown_mutations.mcsm.parsed.txt

unknown_dist_file=~/Downloads/test.txt
known_dist_file=~/Downloads/known_mutations.txt

# out
ddg_file=~/Documents/comp_mut/results/ddg.txt
distances_file=~/Documents/comp_mut/results/distances.txt

rm ${ddg_file}
rm ${distances_file}

echo -e "change\tddg\tstatus" > ${ddg_file}
cat ${known_ddg_file} | sed "s/$/\tknown/" | tr " " "\t" >> ${ddg_file}
cat ${unknown_ddg_file}  | sed "s/$/\tunknown/" | tr " " "\t" >> ${ddg_file}

echo -e "dist\tstatus" > ${distances_file}
cat ${known_dist_file} | sed "s/ (alt loc A)//" | tail +4 | tr " " "\t" | cut -f10 | sed "s/$/\tknown/" >> ${distances_file}
cat ${unknown_dist_file} | sed "s/ (alt loc A)//" | tail +4 | tr " " "\t" | cut -f10 | sed "s/$/\tunknown/" >> ${distances_file}
