#!/bin/bash
cd /media/romuald/TOSHIBA_EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_bam_downsampling/2021_01_01

source ~/python_env/MACS2/bin/activate

for file in in $(ls -t | grep cut_downsampled)

do
  name=$(basename -- $file | head -c 3)
  dir=/media/romuald/TOSHIBA_EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_peak_calling/2021_01_01/
  macs2 callpeak -t $file -n $name\_cut --outdir $dir -f BAMPE -g hs -B --broad --broad-cutoff 0.1
done
