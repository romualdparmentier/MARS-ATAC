#!/bin/bash
cd /home/romuald/MARS-ATAC/exp/ATAC_bam_downsampling
cd $(ls -ltr |grep ^d |tail -1 | awk '{print $9}') # chande directory to the last folder created ($9 is the 9th field of the file metadata >> i.e complete file name)

source ~/python_env/MACS2/bin/activate

dt=$(date '+%Y_%m_%d_%Hh_%Mmin/')
mkdir /home/romuald/MARS-ATAC/exp/ATAC_peak_calling/$dt

for file in in $(ls -t | grep cut_downsampled)

do
  name=$(basename -- $file | head -c 3)
  dir=/home/romuald/MARS-ATAC/exp/ATAC_peak_calling/$dt
  macs2 callpeak -t $file -n $name\_cut --outdir $dir -f BAMPE -g hs -B --broad --broad-cutoff 0.1
done
