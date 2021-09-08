library(ezknitr)

# This script launches all the R markdown script (bin  ) and prodoced a report in each appropriate experiment (exp) folder

# 1 ATAC_peak_annotation

ezknit(
  file = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/bin/ATAC_peak_annotation.Rmd", 
  out_dir = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_peak_annotation/2021_01_01/",
  keep_md = F 
  )

# 2 ATAC_peak_readcount

ezknit(
  file = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/bin/ATAC_peak_readcount.Rmd", 
  out_dir = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_peak_readcount/2021_01_01/",
  keep_md = F 
  )



# 3 ATAC_peak_dynamics

ezknit(
  file = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/bin/ATAC_peak_dynamics.Rmd", 
  out_dir = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_peak_dynamics/2021_01_01/",
  keep_md = F 
  )

# 4 ATAC_differential_accessibility

ezknit(
  file = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/bin/ATAC_differential_accessibility.Rmd", 
  out_dir = "/media/romuald/TOSHIBA EXT/Romuald_NGS/MARS-ATAC/exp/ATAC_differential_accessibility/2021_01_01/",
  keep_md = F 
  )