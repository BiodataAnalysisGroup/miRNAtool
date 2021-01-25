library(multiMiR)
library(enrichR)
library(saveImageHigh)

functional_analysis <- function(sign.table.f, validated_or_predicted, kegg_enrich_criterion, go_criterion, output_dir){
  
  
  top_miRNAs<-sign.table.f$ID
  multimir_results <- NULL
  
  if (validated_or_predicted == 'validated'){
  #validated, predicted, disease.drug
    multimir_results<- get_multimir(org = 'hsa',
                                    mirna = top_miRNAs,
                                    table = 'validated',
                                    summary = TRUE)
    
    #multimir_results.validated@summary$target_symbol
  } else if (validated_or_predicted == 'predicted') {
    multimir_results <- get_multimir(org = 'hsa',
                                     mirna = top_miRNAs,
                                     table  = 'predicted',
                                     summary = TRUE)
    #length(multimir_results.predicted@summary[,1])
  } else {
    print('Wrong input value in validated_or_predicted variable!')
    return(NULL)
  }
  ##
  res<-as.data.frame(multimir_results@summary)
  res$barcode<- paste(res$mature_mirna_id, res$target_symbol)
  res.plot<-res[!duplicated(res$barcode), ]
  write.csv(res$target_symbol, paste('output/Tables/multimir_results_',validated_or_predicted,'_summary_target_symbol.csv', sep = ''))
  
  targetgenes_count<-res.plot %>% group_by(target_symbol) %>% summarise(count= n()) %>% filter(!is.na(target_symbol))
  
  gos <-targetgenes_count
  gos <- gos[order(-gos$count), ]  # sort
  gos <- gos[c(2:50),]
  gos$target_symbol <- factor(gos$target_symbol, levels=gos$target_symbol)
  head(gos)
  
  # Diverging Barcharts
  
  #png(file=paste(output_dir,'/diverging_barcharts_plate.png', sep = ''), width=900, height=600)
  myplot <- ggplot(gos, aes(x=target_symbol, y=as.numeric(count)) )+ 
    geom_bar(stat='identity', width=.2,position="dodge")  +
    coord_flip()
  save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/diverging_barcharts_plate.pdf', sep = ''), width = 6, height = 6)
  save_image({print(myplot)}, file.name = paste(output_dir,'/diverging_barcharts_plate.png', sep = ''), width = 6, height = 8)
  #print(myplot)
  #dev.off()
  
  genenames<- res.plot$target_symbol
  
  genels = c()
  strGenes <- as.character(genenames)  
  splitGenes = strsplit(strGenes, "[;]")
  splitGenes <-matrix(unlist(splitGenes)) 
  genesym <- unlist(splitGenes[!duplicated(splitGenes)])
  
  # KEGG enrichment and plot
  dbs <- list()
  dbs <- "KEGG_2019_Human"
  enriched <- enrichr(genesym, dbs)
  KEGG_enrich<- as.data.frame(enriched[["KEGG_2019_Human"]])
  KEGG_enrich.f<- subset(KEGG_enrich, Adjusted.P.value < kegg_enrich_criterion)
  write.csv(KEGG_enrich.f, paste('output/Tables/KEGG_enrich_f_',validated_or_predicted,'.csv', sep = ''))
  
  
  gos <- KEGG_enrich.f
  gos <- gos[order(-gos$P.value), ]  
  gos$Term <- factor(gos$Term, levels=gos$Term)
  
  #png(file=paste(output_dir,'/KEGG_enrichment',plate,'.png', sep = ''), width=1500, height=1500)
  myplot <- ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
    geom_bar(stat='identity', width=.2,position="dodge")  +
    coord_flip()
  save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/KEGG_enrichment',plate,'.pdf', sep = ''), width = 10, height = 14)
  save_image({print(myplot)}, file.name = paste(output_dir,'/KEGG_enrichment',plate,'.png', sep = ''), width = 10, height = 14)
  #print(myplot)
  #dev.off()
  
  #GO and plot
  dbs <- listEnrichrDbs()
  dbs <- "GO_Biological_Process_2018"
  enriched <- enrichr(genesym, dbs)
  GO_enrich<- as.data.frame(enriched[["GO_Biological_Process_2018"]])
  GO_enrich.f<- subset(GO_enrich, Adjusted.P.value < go_criterion)
  write.csv(KEGG_enrich.f, paste('output/Tables/GO_enrich_f_',validated_or_predicted,'.csv', sep = ''))
  
  gos <- GO_enrich.f 
  gos <- gos[order(-gos$Adjusted.P.value), ]  # sort
  gos$Term <- factor(gos$Term, levels=gos$Term)
  head(gos)
  
  # Diverging Barcharts
  #png(file=paste(output_dir,'/GO.png', sep = ''), width=1500, height=1500)
  myplot <- ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
    geom_bar(stat='identity', width=.1,position="dodge")  +
    coord_flip()
  save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/GO.pdf', sep = ''), width = 20, height = 20)
  save_image({print(myplot)}, file.name = paste(output_dir,'/GO.png', sep = ''), width = 20, height = 20)
  #print(myplot)
  #dev.off()
  
}