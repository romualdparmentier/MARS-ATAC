#################################################################
###################### User inputs ##############################   # Part of the script to change for each analysis
#################################################################

# Working directory
the_dir = "~/Bureau/ATAC_2020/"
setwd(the_dir)

conditions_list = list("Xvivo"= list("00h" = c("LBL","LAS", "JLL")),
                       "MP" = list("05h" = c("QQY", "QR3", "QR8"),
                                   "24h" = c("QR0","QR5","QR9"),
                                   "48h" = c("QSJ","QR7","QFV")))
                       
# Define list of genes of interest
  list_gene = c("GATA2","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1","HHEX","TAL1","GATA1")
#list_gene = c("GATA2","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1") # Sans HEXX et GATA1 qui bug

# Liste à définir à la main : il faut faire des tests de graphiques pour la définir
max_read = c(25,20,15,10,30,15,20,20,40,40,10)
# max_read = c(25,20,15,10,30,15,20,20)# Sans HEXX et GATA1 qui bug

#################################################################
################### Libraries ###################################
#################################################################

library(Gviz)
library(dplyr)
library(GenomicRanges)

#################################################################
################## Functions and constants ######################
#################################################################

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#################################################################
############################# Script ############################
#################################################################

#*******
# Create saving directories
#*******

dir.create("results/Gene_Coverage/", showWarnings = FALSE )

#*******
# Loading and storing files - common for every loop turn
#*******

# Load promoters data base (from RAVI)
prom_gene_fantom_gr = loadRData(paste0(the_dir,"data/prom_gene_fantom_gr.rdata"))

# Adding FANTOM5 information for HHEX (lacking in original dataset, don't know why ?)

HHEX_prom = GRanges(seqnames = "chr10",
                    ranges = IRanges(start = c(94449649,94449675,94449703,94451574), end = c(94449664,94449694,94449718,94451587)), 
                    strand = "*")

mcols(HHEX_prom) = tibble(gene = "HHEX")
prom_gene_fantom_gr = c(prom_gene_fantom_gr, HHEX_prom)
prom_gene_fantom_gr = sort(prom_gene_fantom_gr)

# Extract peak control 00h directory and name
peak_control = dir(path = paste0(the_dir, "results/genomic_ranges/static_peaks"), pattern = "00h_Xvivo_gr.rda", full.names = TRUE)
peak_control_name = dir(path = paste0(the_dir, "results/genomic_ranges/static_peaks"), pattern = "00h_Xvivo_gr.rda")
peak_control_name = substr(peak_control_name, start = 20,(nchar(peak_control_name)-7))

# Extract a random bam control among the three files used for intersection
# random = sample(1:3, 1) # sélection aléatoire d'un BAM parmi les trois pour ce point de temps
# bam_chosen = conditions_list[["Xvivo"]][["00h"]][random]
# list_bam_files = list(paste0(the_dir,"results/downsampled_BAM/", bam_chosen, "_downsampled.bam"))

list_bam_files = paste0(the_dir,"results/downsampled_BAM/JLL_bedgraph.bedgraph")
names(list_bam_files) = substr(list_bam_files, nchar(list_bam_files)-20, nchar(list_bam_files)-18)

# Extract and organize gene information from genesymbol via Gviz
data(genesymbol, package = "biovizBase")
df_region = tibble()

for(i in 1:length(list_gene)) {
    region = genesymbol[list_gene[i]]
    region = range(region, ignore.strand = TRUE)
    region = region + 5000
    region = keepStandardChromosomes(region)
    df_region[i,"gene"] = list_gene[i]
    df_region[i,"chr"] = as.character(seqnames(region))
    df_region[i,"start"] = start(region)
    df_region[i,"end"] = end(region)
    df_region[i,"max_read"] = max_read[i]
    print(df_region[i,])
}

#*****
# Loop to draw a tracking plot for every condition
#*****

for (sublist in 1:length(conditions_list)) {

  if (names(conditions_list[sublist]) != "Xvivo") { # no drawing with only Xvivo point

        #########
        # Loading and storing files
        #########

        # Extract peak files directions             # double pattern => .+ pour &
        peak_files = dir(path = paste0(the_dir, "results/genomic_ranges/static_peaks"), pattern = names(conditions_list[sublist]), full.names = TRUE)
        peak_files = c(peak_control, peak_files)

        # Extract peak files names
        peak_files_name = dir(path = paste0(the_dir, "results/genomic_ranges/static_peaks"), pattern = names(conditions_list[sublist]))
        peak_files_name = substr(peak_files_name, start = 20,(nchar(peak_files_name)-7))
        peak_files_name = c(peak_control_name, peak_files_name)

        # Load peak files in a list and name it
        list_peaks_kept = lapply(X = peak_files, FUN = loadRData)
        names(list_peaks_kept) = peak_files_name

        # Extract random BAM files (among the 3 of the intersection) corresponding to peak files
         for (time_point in 1:length(conditions_list[[sublist]])) {
           random = sample(1:3, 1) # sélection aléatoire d'un BAM parmi les trois pour ce point de temps
           bam_chosen = conditions_list[[sublist]][[time_point]][random]
           bam_files = paste0(the_dir,"results/downsampled_BAM/", bam_chosen, "_downsampled.bam")
           list_bam_files = c(list_bam_files, bam_files)
         }
        
        names(list_bam_files)[-1] = substr(list_bam_files[-1], nchar(list_bam_files[-1])-18, nchar(list_bam_files[-1])-16)
        


       #########
       # Tracking plots creation for each gene
       #########

       for (i in 1:nrow(df_region)) { # one plot per gene of the list

         # Chromosom with gene illustration
         # ********************************
         axTrack <- GenomeAxisTrack() # scale
         idxTrack <- IdeogramTrack(genome="hg19", chromosome= df_region$chr[i]) # drawing


         # Isoform tracking for the studied gene and promoters
         # *****************    **********************************
         ensGenes <- UcscTrack(genome = "hg19", chromosome = df_region$chr[i],
                               background.title = "transparent",
                               alpha.title = 1, col.title = "black",
                               cex.title = 1.5, track = "ensGene",
                               from = df_region$start[i], to = df_region$end[i],
                               trackType = "GeneRegionTrack", rstarts = "exonStarts",
                               rends = "exonEnds", transcript = "name", strand = "strand",
                               fill = "lightgrey", name = df_region$gene[i])

         # recupère seulement ce qui est commun entre les promoteurs et les régions définies dans grange
         prom_peak = subsetByOverlaps(prom_gene_fantom_gr,
                                      GRanges(seqnames = df_region$chr[i],
                                              ranges= IRanges(start = df_region$start[i], end = df_region$end[i])))

         df_prom_retained = tibble("start" = start(prom_peak),
                                   "end" = end(prom_peak),
                                   "chr" = as.character(seqnames(prom_peak)))
         
         if(nrow(df_prom_retained) != 0){

         # faire apparaître les promoteurs sur le track des isoformes
         promoterTrack = HighlightTrack(trackList = ensGenes,
                               start = df_prom_retained$start,
                               end = df_prom_retained$end,
                               chromosome = df_region$chr[i],
                               col = "green", fill = "green")}


         # DataTracking : one per time_point
         # *********************************

         list_datatrack = list()

          for (time_point in 1:(length(conditions_list[[sublist]])+1)) { #nb de time point de la condition + Xvivo

            # ATTENTION VERIFIER QUE LE NOM correspond bien au range !!!
            dt = DataTrack(type = "histogram",
                         name = paste0(names(list_peaks_kept[time_point]), "_track"),
                         background.title = "transparent",
                         fill.histogram = "#0072B2",
                         col.histogram = "#0072B2",
                         alpha.title = 1,
                         col.title = "black",
                         col.axis = "black",
                         cex.title = 1,
                         range = unlist(list_bam_files[time_point]),
                         genome = "hg19",
                         ylim = c(0,df_region$max_read[i]),
                         window = -1,
                         chromosome = df_region$chr[i])

          list_datatrack = c(list_datatrack, dt)  #append ?

          }
        # list_datatrack[[]]$
        names(list_datatrack) = names(list_peaks_kept)

        #########
        # Assembling every tracking plot into one
        #########

        # Blank space to improve plot readibility
        blank_space = GenomeAxisTrack(col = "white",fill = "white", fontcolor = "white")

        # Gather plots in a list
        if (exists("promoterTrack") == TRUE){
          list_to_plot = list(idxTrack, axTrack, blank_space, promoterTrack)
        } else{
            list_to_plot = list(idxTrack, axTrack, blank_space, ensGenes)
            }

        # Create list of peaks detected in the gene region +/- 5kb

        for(time_point in 1:length(list_peaks_kept)) {

           # Determine if peaks are detected in this gene for this time_point
           peaks_in = subsetByOverlaps(list_peaks_kept[[time_point]],
                                       GRanges(seqnames = df_region$chr[i],
                                               ranges= IRanges(start = df_region$start[i], end = df_region$end[i])))

           # if not detected, no graph for this time_point for this gene
           if (length(peaks_in) == 0){ print("no_peaks_detected")}

           # if detected, graph with highligted peaks
           else {
             df_peaks_retained = tibble("start" = start(peaks_in),
                                        "end" = end(peaks_in),
                                        "chr" = as.character(seqnames(peaks_in)))

             ht = HighlightTrack(trackList = list_datatrack[time_point],
                                 start = df_peaks_retained$start,
                                 end = df_peaks_retained$end,
                                 chromosome = df_region$chr[i])

             list_to_plot = c(list_to_plot, blank_space, ht) # append ?
             }
          }
############# PROBLEME !!!!!!!!!!
          
           # Open a pdf file => enregistrement ??? demande d'ouverture ???
           pdf(file =  paste0(the_dir, "results/Gene_coverage/", df_region$gene[i], "_Coverage_plot.pdf"), width = 8, height = 10)

           # Create a plot
           plotTracks(list_to_plot, from = df_region$start[i], to = df_region$end[i], showTitle = TRUE)

           # Close the pdf file
           dev.off()

        }
        } # fin du if 
} # fin du for sur les conditions


  
