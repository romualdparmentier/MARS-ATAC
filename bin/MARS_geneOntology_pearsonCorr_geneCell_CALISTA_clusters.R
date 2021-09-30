
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyverse)

###########
# Fonctions 
###########

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#################################
### Gene Ontology ==> par cluster
#################################

# Paramètre pour charger les fichiers
clusters = 1:5
minCell = 5
COR_cutoff = 0.7

# GO Terms : ALL  

for(donor in 1:2){
  
  for (cluster_num in 1:length(clusters)){
    
    print(paste0("Donor_", donor,"_Cluster_", cluster_num))
    
    # Fichier contenant la liste des gènes à tester en GO
    df = read.table(paste0(
      "/home/romuald/Bureau/Analyses_MARS-ATAC/Diagrammes_corrélation/output/Cell", minCell,
      "_Cutoff0", COR_cutoff*100, # *100 parce que le nom du dossier est 070 ou 085
      "/list_gene_cluster_correlations/list_gene_cluster", cluster_num,
      "_Pearson_cutoff_", COR_cutoff,
      "_MinCell_", minCell,
      "_donor", donor,
      ".csv"))
    
    colnames(df) = "gene_name"
    
    # On stock le nom des gene (hgnc symbol) dans une variable
    gene_list = df$gene_name
    
    name_correspondace = bitr(df$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

    print(paste0("Donor_", donor,"_Cluster_", cluster_num, "_Entrez_ID_retrieval_Done"))

    colnames(name_correspondace) = c("gene_name","entrezgene_id")
    
    df_completed = dplyr::left_join(name_correspondace, df ,by = "gene_name")
    df_completed = drop_na(df_completed)
    
    print(paste0("Donor_",donor,"_Cluster_", cluster_num, "_Enrich_GO_onGoing..."))
    
    ego <- enrichGO(
      gene = df_completed$entrezgene_id,
      OrgDb = org.Hs.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE)
    
    print(paste0("Donor_",donor,"_Cluster_", cluster_num, "_Enrich_GO_done"))
    
    # Si on veut réduire la redondance des termes GO en les groupant par un score de similarité
    # ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    
    save(ego,
      file = paste0(
        "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
        "_Pearson", COR_cutoff,
        "/each_cluster/ego_all_cluster", cluster_num,
        "_donor", donor,
        ".rda")) # Sauvegarde le fichier R
    
    write.csv2(ego@result, 
      file = paste0(
        "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
        "_Pearson", COR_cutoff,
        "/each_cluster/Summary_ego_ALL_cluster", cluster_num,
        "_donor", donor,
        ".csv")) # Excel du GO terms
    
    png(file = paste0(
          "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
          "_Pearson", COR_cutoff,
          "/each_cluster/plot_ego_ALL_cluster", cluster_num,
          "_donor", donor,
          ".png"), 
        width = 1200, height = 800) # prépare un png
    
    print(dotplot(
      ego,
      showCategory = 50)) #Imprime le plot dans le png préparé
    
    dev.off()# s'assure que le png est bien fermé
    
        }
}


####################################################################################
# Compare GO terms avec la fonction CompareCluster implémentée dans Cluster Profiler
####################################################################################

# Paramètre pour charger les fichiers
clusters = 1:5
COR_cutoff = 0.70
minCell = 5


# Paramètre pour les boucles
list_entrez = list()
list_name = list()

for(donor in 1:2){
  
  print(paste("Analysis_Donor_",donor,"Ongoing..."))
  
  for (cluster_num in 1:length(clusters)){

    print(paste0("Donor_", donor,"_Cluster_", cluster_num))
    
    # Fichier contenant la liste des gènes à tester en GO
    df = read.table(paste0(
      "/home/romuald/Bureau/Analyses_MARS-ATAC/Diagrammes_corrélation/output/Cell", minCell,
      "_Cutoff0", COR_cutoff*100, # *100 parce que le nom du dossier est 070 ou 085
      "/list_gene_cluster_correlations/list_gene_cluster", cluster_num,
      "_Pearson_cutoff_", COR_cutoff,
      "_MinCell_", minCell,
      "_donor", donor,
      ".csv"))
    
    colnames(df) = "gene_name"
    
    # On stock le nom des gene (hgnc symbol) dans une variable
    gene_list = df$gene_name
    
    name_correspondace = bitr(df$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    colnames(name_correspondace) = c("gene_name","entrezgene_id")
    
    df_completed = df %>% left_join(name_correspondace, df ,by = "gene_name")
    df_completed = df_completed %>% drop_na()
    
    list_entrez[[cluster_num]] = df_completed$entrezgene_id
    list_name[[cluster_num]] = df_completed$gene_name
    names(list_entrez)[[cluster_num]] = paste("cluster_",cluster_num)
    
  }
  
  print("Comparison_GO_Terms_Ongoing...")
  
  # Compare GO terms entre les clusters
  comp_clust = compareCluster(geneClusters = list_entrez, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
  
  print("Comparison_GO_Terms_DONE")

  # Evite la redondance des GO terms
  comp_clust <- clusterProfiler::simplify(comp_clust, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  save(list_entrez,
      file = paste0(
        "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
        "_Pearson", COR_cutoff,
        "/compare_cluster/list_entrezID_all_cluster",
        "_donor", donor,
        ".rda")) # Sauvegarde la liste des entrez id pour chaque cluster
  
  save(comp_clust,
      file = paste0(
        "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
        "_Pearson", COR_cutoff,
        "/compare_cluster/Simplified_comp_cluster",
        "_donor", donor,
        ".rda")) # Sauvegarde le fichier R
    
    write.csv2(comp_clust@compareClusterResult, 
      file = paste0(
        "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
        "_Pearson", COR_cutoff,
        "/compare_cluster/Simplified_Summary_comp",
        "_donor", donor,
        ".csv")) # Excel du GO terms
    
    png(file = paste0(
          "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
          "_Pearson", COR_cutoff,
          "/compare_cluster/plot_simplified_comp_clust",
          "_donor", donor,
          ".png"), 
        width = 1200, height = 800) # prépare un png
    
    print(dotplot(
      comp_clust,
      showCategory = 50)) #Imprime le plot dans le png préparé
    
    dev.off()# s'assure que le png est bien fermé
}

########################################
### Compare cluster - All donor together
#######################################

# Paramètre pour charger les fichiers
clusters = 1:5
COR_cutoff = 0.70
minCell = 5

list_entrez = list()
list_name = list()
list_all = list()

for(donor in 1:2){
  
  print(paste("Analysis_Donor_",donor,"Ongoing..."))
  
  for (cluster_num in 1:length(clusters)){

    print(paste0("Donor_", donor,"_Cluster_", cluster_num))
    
    # Fichier contenant la liste des gènes à tester en GO
    df = read.table(paste0(
      "/home/romuald/Bureau/Analyses_MARS-ATAC/Diagrammes_corrélation/output/Cell", minCell,
      "_Cutoff0", COR_cutoff*100, # *100 parce que le nom du dossier est 070 ou 085
      "/list_gene_cluster_correlations/list_gene_cluster", cluster_num,
      "_Pearson_cutoff_", COR_cutoff,
      "_MinCell_", minCell,
      "_donor", donor,
      ".csv"))
    
    colnames(df) = "gene_name"
    
    # On stock le nom des gene (hgnc symbol) dans une variable
    gene_list = df$gene_name
    
    name_correspondace = bitr(df$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    colnames(name_correspondace) = c("gene_name","entrezgene_id")
    
    df_completed = df %>% left_join(name_correspondace, df ,by = "gene_name")
    df_completed = df_completed %>% drop_na()
    
    list_entrez[[cluster_num]] = df_completed$entrezgene_id
    list_name[[cluster_num]] = df_completed$gene_name
    names(list_entrez)[[cluster_num]] = paste0("#",cluster_num," d",donor)
    
  }
  
  list_all = c(list_all, list_entrez)
  
}

# Permet de lister tous les cluster un à côté des autres, le cluster 4 avant le 3
list_all = list_all[names(list_all)[c(1,6,2,7,4,9,3,8,5,10)]]
  
# Compare GO terms entre les clusters
comp_clust = compareCluster(geneClusters = list_all, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
  
# Evite la redondance des GO terms
comp_clust_simp <- clusterProfiler::simplify(comp_clust, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
save(list_entrez,
  file = paste0(
    "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
    "_Pearson", COR_cutoff,
    "/compare_cluster/Simplified_all_donors/list_entrezID_all_cluster_all_donor.rda")) # Sauvegarde la liste des entrez id pour chaque cluster
  
save(comp_clust_simp,
  file = paste0(
    "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
    "_Pearson", COR_cutoff,
    "/compare_cluster/Simplified_all_donors/Final_comp_cluster_all_donor.rda")) # Sauvegarde le fichier R
    
write.csv2(comp_clust_simp@compareClusterResult, 
  file = paste0(
  "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
  "_Pearson", COR_cutoff,
  "/compare_cluster/Simplified_all_donors/Final_Summary_comp_all_donor.csv")) # Excel du GO terms
    
png(file = paste0(
  "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/MinCell", minCell,
  "_Pearson", COR_cutoff,
  "/compare_cluster/Simplified_all_donors/Final_plot_comp_clust_all_donor.png"), 
  width = 1600, height = 1200) # prépare un png
    
print(dotplot(
  comp_clust_simp,
  showCategory = 50)) #Imprime le plot dans le png préparé
    
  dev.off()# s'assure que le png est bien fermé