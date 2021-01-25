# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])
getwd()
analysis.path=getwd()  ;analysis.path
#analysis.path <- "D://exosomes//data" 

#setwd(analysis.path)
getwd() 

library(data.table)
library(tidyverse)
library(naniar)
library(prada)

list.files()

mirs<- fread("miRs_annotation_3plates.csv")
data<- fread("miRNome_data.csv")
meta<- fread("phenodata.csv")

head(meta)


data<- as.data.frame(data)
data$ID<- mirs$`miRNA ID`
data$plate<- mirs$Plate
head(data)

#i took only the first plate
#TO DO: QC and normalization analysis for each plate
data<- data[!duplicated(data$Well),]
data[,c(2:7)]<- apply(apply(data[,c(2:7)], 2, gsub, patt=",", replace="."), 2, as.numeric)

####QC

# undetermined plot
data.plot<-data[,c(2:7)]
myplot <- vis_miss(data.table(data.plot))
gg_miss_upset(data.table(data.plot))
# boxplot quality before normalization
# TO DO: exclude all the figures on the working directory (ask Nikos for his tool)
Mnemiopsis_count_data<- data.plot
pseudoCount = log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)

df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix 
ggplot(df, aes(x = Samples, y = count)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))+
  coord_flip()

# TO DO: exclude if there are many NAs

##normalization/ case
RTC<- filter(data, ID == "PPC" | ID == "miRTC")
RTC.c<- t(RTC[,c(5:7,8)]%>%  
            group_by(ID) %>%
            summarise_all("mean"))

RTC.n<- t(RTC[,c(2:4,8)]%>%  
            group_by(ID) %>%
            summarise_all("mean"))

RTC.c<- as.data.frame(RTC.c)
RTC.n<- as.data.frame(RTC.n)

names(RTC.c)<- RTC.c[1,]
names(RTC.n)<- RTC.n[1,]

RTC.c<- RTC.c[-1,]
RTC.n<- RTC.n[-1,]

RTC.c$criterion<- as.numeric(RTC.c$miRTC) - as.numeric(RTC.c$PPC)
RTC.n$criterion<- as.numeric(RTC.n$miRTC) - as.numeric(RTC.n$PPC)


RTC.all<- rbind(RTC.c, RTC.n)
rm(RTC.c);rm(RTC.n)
if (RTC.all$criterion < 5) {
  RTC.all$RTC = paste("Pass")
  
} else
{
  RTC.all$RTC = paste("Fail")
}

#TO DO exclude from the matrix the samples with "Fail"


#divide each target of interest by the mean Ct value of the Actin Normalizing gene

endogenous<- filter(data, ID == "SNORD61" | ID == "SNORD68" | 
                      ID == "SNORD72" | ID == "SNORD95" | ID=="RNU6-2")
exogenous<- filter(data, ID == "cel-miR-39")

gg_miss_upset(data.table(endogenous[,c(2:7)]))
gg_miss_upset(data.table(exogenous[,c(2:7)]))

endogenous<- colMeans(endogenous[,c(2:7)], na.rm = T)
exogenous<- colMeans(exogenous[,c(2:7)], na.rm = T)

normalized_data<- filter(data, ID != "SNORD61" & ID != "SNORD68" & 
                           ID != "SNORD72" & ID != "SNORD95" & 
                           ID!="RNU6-2" &  ID != "cel-miR-39" &
                           ID != "PPC" & ID != "miRTC")

#TO DO: the user could select the endogenous or exogenous normalazitaion
##normalize with endogenous
normalized_data[,c(2:7)]<-normalized_data[,c(2:7)]/endogenous

##normalized with exogenous

normalized_data[,c(2:7)]<-normalized_data[,c(2:7)]/exogenous

#log transform the data to develop normality in data

##boxplot with quality
#corrected data
Mnemiopsis_count_data<- normalized_data[,c(2:7)]
pseudoCount = log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)

df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix 
ggplot(df, aes(x = Samples, y = count)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))+
  coord_flip()

#TO DO: exclude a sample manually (because of a plot) and then continue to the the diff analysis
#TO DO: merge the three plates before the diff analysis

#DIff analysis

list<- as.character(meta$Sample_ID)
dfexp<-normalized_data[,list]

row.names(dfexp)<- normalized_data$ID

pheno <- factor(meta$Group)

phenoMat <- model.matrix(~pheno)
colnames(phenoMat) <- sub("^pheno","",colnames(phenoMat))
phenoMat;dim(phenoMat)

library(limma)


fit <- lmFit(object = dfexp,design = phenoMat)
gc()
set.seed(6)
fit <- eBayes(fit)

gc()
degCLLs <- topTable(fit,number =nrow(dfexp),adjust.method = "fdr",sort.by = "p")
head(degCLLs)


sign.table<- as.data.frame(degCLLs)
sign.table$ID<- row.names(sign.table)
sign.table.f<- sign.table[sign.table$adj.P.Val<= 0.01, ] #TO DO: user's selection

##Heatmap
data<-as.data.frame(normalized_data)
head(data)
#data$ID<- row.names(data)


data.all<-merge(data, sign.table.f, by.x = "ID", by.y = "ID")
colnames(data.all)
heat.data<- data.all[,c(3:8)]
row.names(heat.data)<- data.all$ID

z.mat <- t(scale(t(heat.data), center=TRUE, scale=TRUE))
head(z.mat)


library(ComplexHeatmap)
column_ha3 = HeatmapAnnotation(Groups = meta$Group,
                               na_col="white"
)

Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        top_annotation = column_ha3, row_names_gp = gpar(fontsize = 4), show_row_names = FALSE) #column_km = 2)


################
library("multiMiR")
top_miRNAs<-sign.table.f$ID

#validated, predicted, disease.drug
multimir_results.validated <- get_multimir(org     = 'hsa',
                                           mirna   = top_miRNAs,
                                           table   = 'validated',
                                           summary = TRUE)

multimir_results.validated@summary$target_symbol

multimir_results.predicted <- get_multimir(org     = 'hsa',
                                           mirna   = top_miRNAs,
                                           table   = 'predicted',
                                           summary = TRUE)
length(multimir_results@summary[,1])


##
res<-as.data.frame(multimir_results.validated@summary)
res$barcode<- paste(res$mature_mirna_id, res$target_symbol)
res.plot<-res[!duplicated(res$barcode), ]

targetgenes_count<-res.plot %>% group_by(target_symbol) %>% summarise(count= n()) %>% filter(!is.na(target_symbol))

gos <-targetgenes_count
gos <- gos[order(-gos$count), ]  # sort
gos <- gos[c(2:50),]
gos$target_symbol <- factor(gos$target_symbol, levels=gos$target_symbol)
head(gos)
# Diverging Barcharts
ggplot(gos, aes(x=target_symbol, y=as.numeric(count)) )+ 
  geom_bar(stat='identity', width=.2,position="dodge")  +
  coord_flip()


##
library(enrichR)
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
KEGG_enrich.f<- subset(KEGG_enrich, Adjusted.P.value< 0.05)

gos <- KEGG_enrich.f
gos <- gos[order(-gos$P.value), ]  
gos$Term <- factor(gos$Term, levels=gos$Term)
ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
  geom_bar(stat='identity', width=.2,position="dodge")  +
  coord_flip()


#GO and plot
dbs <- listEnrichrDbs()
dbs <- "GO_Biological_Process_2018"
enriched <- enrichr(genesym, dbs)
GO_enrich<- as.data.frame(enriched[["GO_Biological_Process_2018"]])
GO_enrich.f<- subset(GO_enrich, Adjusted.P.value< 0.05)

gos <- GO_enrich.f 
gos <- gos[order(-gos$Adjusted.P.value), ]  # sort
gos$Term <- factor(gos$Term, levels=gos$Term)
head(gos)
# Diverging Barcharts
ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
  geom_bar(stat='identity', width=.1,position="dodge")  +
  coord_flip()
