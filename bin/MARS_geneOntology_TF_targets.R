#######################
## Gene Ontology ==> TF
#######################
 #Liste des gènes targetés par les TF (DE ou non DE) entre 05 et 24 ou 24 et 48
 TF_0524 = loadRData("/home/romuald/Bureau/Analyses_MARS-ATAC/TF_vs_Prom_vs_Expr_redo/intermediary_generated_files/list_target_genes_Change_NoChange_05_24.rda")
 TF_2448 = loadRData("/home/romuald/Bureau/Analyses_MARS-ATAC/TF_vs_Prom_vs_Expr_redo/intermediary_generated_files/list_target_genes_Change_NoChange_24_48.rda")

 list_TF_target = list(TF_0524, TF_2448)

 TF_change_nochange = c("target_TF_change", "target_TF_noChange")
 time_interval = c("05_vs_24","24_vs_48")

 for(time in 1:length(time_interval)){

   for (TF_statut in 1:length(TF_change_nochange)){

     print(paste0("Time_interval_",time_interval[time],"___",TF_change_nochange[TF_statut]))

     # On stock le nom des gene (hgnc symbol) dans une variable
     gene_list_name = list_TF_target[[time]][[TF_statut]]$Gene

     name_correspondace = bitr(gene_list_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
     colnames(name_correspondace) = c("gene_name","entrezgene_id")

     name_correspondace = name_correspondace %>% drop_na()

     # GO Terms : ALL

     print("Enrich_GOTerms_Ongoing...")

     ego <- enrichGO(gene = name_correspondace$entrezgene_id,
       OrgDb         = org.Hs.eg.db,
       ont           = "ALL",
       pAdjustMethod = "BH",
       pvalueCutoff  = 0.05,
       qvalueCutoff  = 0.2,
       readable      = TRUE)

     print("Enrich_GOTerms_DONE.")

     save(ego,
       file = paste0(
         "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/each_interval/ego_ALL_",TF_change_nochange[TF_statut],
         "_",time_interval[time],".rda")) # Sauvegarde le fichier R

     write.csv2(ego@result,
     file = paste0("/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/each_interval/summary_ego_ALL_",TF_change_nochange[TF_statut],
       "_",time_interval[time],".csv")) # Excel du GO terms

     png(file = paste0("/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/each_interval/ego_BP_plot_",TF_change_nochange[TF_statut],
       "_",time_interval[time],".png"),
       width = 1200, height = 800) # prépare un png

     print(dotplot(ego, showCategory = 50)) # Imprime le plot dans le png préparé

     dev.off()# s'assure que le png est bien fermé

   }

 }




 #############################################################
 ### Compare TF targets time interval
 #############################################################

 #Liste des gènes targetés par les TF (DE ou non DE) entre 05 et 24 ou 24 et 48
 TF_0524 = loadRData("/home/romuald/Bureau/Analyses_MARS-ATAC/TF_vs_Prom_vs_Expr_redo/intermediary_generated_files/list_target_genes_Change_NoChange_05_24.rda")
 TF_2448 = loadRData("/home/romuald/Bureau/Analyses_MARS-ATAC/TF_vs_Prom_vs_Expr_redo/intermediary_generated_files/list_target_genes_Change_NoChange_24_48.rda")

 list_TF_target = list(TF_0524, TF_2448)
 list_entrez = list()
 list_all = list()

 TF_change_nochange = c("target_TF_change", "target_TF_noChange")
 time_interval = c("05_vs_24","24_vs_48")

 for(time in 1:length(time_interval)){

   for (TF_statut in 1:length(TF_change_nochange)){

     print(paste0("Time_interval_",time_interval[time],"___",TF_change_nochange[TF_statut]))

     # On stock le nom des gene (hgnc symbol) dans une variable
     gene_list_name = list_TF_target[[time]][[TF_statut]]$Gene

     name_correspondace = bitr(gene_list_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
     colnames(name_correspondace) = c("gene_name","entrezgene_id")

     name_correspondace = name_correspondace %>% drop_na()

     list_entrez[[TF_statut]] = name_correspondace$entrezgene_id
     names(list_entrez)[TF_statut] = paste(time_interval[time],TF_change_nochange[TF_statut])

   }

   list_all = c(list_all, list_entrez)

 }

 comp_clust = compareCluster(geneClusters = list_all, fun = "enrichGO", OrgDb = "org.Hs.eg.db")

 # Evite la redondance des GOTerms
 comp_clust <- clusterProfiler::simplify(comp_clust, cutoff = 0.5, by = "p.adjust", select_fun = min)

   save(list_all,
       file = paste0(
         "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/compare_interval/list_entrezID_all_interval_statut.rda")) # Sauvegarde la liste des entrez id pour chaque cluster

   save(comp_clust,
       file = paste0(
         "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/compare_interval/Simplified_Comp_interval_statut.rda")) # Sauvegarde le fichier R

     write.csv2(comp_clust@compareClusterResult,
       file = paste0(
         "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/compare_interval/Simplified_Summary_comp_interval_statut.csv")) # Excel du GO terms

     png(file = paste0(
           "/home/romuald/Bureau/Analyses_MARS-ATAC/GO/output/TF/compare_interval/Simplified_plot_comp_interval_statut.png"),
         width = 1200, height = 800) # prépare un png

     print(dotplot(
       comp_clust,
       showCategory = 50)) #Imprime le plot dans le png préparé

     dev.off()# s'assure que le png est bien fermé
