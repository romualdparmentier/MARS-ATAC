#!/bin/bash
cd /home/romuaald/MARS-ATAC/data

# Set the mininum of reads in all bam files

unset min # Remove any pre-set variable named "min"

for file in $(ls -t | grep _cut.bam)
do
  name=$(basename -- $file)
  echo $name
  count=$(samtools view -c $file)

  if test -z "$min" #Test if the variable "min" is set
    then
    min=$count #If not (=TRUE) then assign min to the value of count
    echo first initialization, min = $min

  elif (("$count">"$min")) #If the read count number of the next file is superior to min,then min remains unchanged
    then
    echo count for $name = $count, superior to min = $min

  elif (("$count"<"$min")) #If the read count number is inferior to min, then min is assigned to the count value of the next file
    then
    min=$count
    echo new min is $min

  fi

done

## Calculate the ratio for downsampling thanks to the mininum

# First create an output folder named wih the date of day
dt=$(date '+%Y_%m_%d_%Hh_%Mmin/');
mkdir ~/Bureau/$dt

for file in $(ls -t | grep _cut.bam)
do
  name=$(basename -- $file | head -c 3) #Extract the sample sequencing code (3 first letter of thhe filename)
  prefix="/home/romuald/MARS-ATAC/exp/ATAC_bam_downsampling/"$dt #Path to store the outut in the previously created folder
  suffix="_cut_downsampled.bam"
  count=$(samtools view -c $file) #Count of the number of reads in the bam file
  echo "Downsampling value for $name = $min"
  downsampling_ratio=$(echo "$min/$count" | bc -l)
  samtools view -s 1$downsampling_ratio -b $file > "$prefix$name$suffix" # Downsampling (-s) of the file with the seed°1 (1$ratio)
done
