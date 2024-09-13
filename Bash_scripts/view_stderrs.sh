#! /bin/bash

filepath=$1

echo -n "Time taken: "
cat $filepath/j*.bed_checkRes* | grep 'Elapsed (wall clock) time (h:mm:ss or m:ss): '

echo -e "\nExit status: "
cat $filepath/j*.bed_checkRes* | grep 'Exit status' 

echo -e "\nPercent of CPU each job got: "
cat $filepath/j*.bed_checkRes* | grep 'Percent of CPU' 

echo -ne "\nRAM: "
cat $filepath/j*.bed_checkRes* | grep 'Maximum resident set size (kbytes): ' | awk '{sum+=$NF+0} END{print "Total RAM used is " sum " and average RAM used is " sum/NR}'
cat $filepath/j*.bed_checkRes* | grep 'Maximum resident set size (kbytes):'

echo -ne "\n Example of Command used: "
head -1 $filepath/j0.bed_checkRes.stderr

echo "----------End of report----------"
