#Package readr is used for reading a tsv file.
#Package dplyr is used for doing some data manipulation functions.
library(readr)
library(dplyr)
#Place the disease’s name here. Note that the name must exact match the one shown in the original GWAS table.
wanted_disease <-  c("Crohn's disease")

GWAS<-read_tsv('gwas_catalog.tsv')
GWAS1<-filter(GWAS, GWAS$'REPLICATION SAMPLE DESCRIPTION'!="NA")
GWAS2<-filter(GWAS1,GWAS1$'P-VALUE'<=(10^(-7)),header=TRUE)
anth_traits<-read.csv('supplementary_table_1.csv')
GWAS3<- anti_join(GWAS2,anth_traits,by=c("DISEASE/TRAIT"="GWAS.Disease.Trait"))

Hugo_names <-read.csv('genes.csv',sep='\t', header=FALSE)
colnames(Hugo_names)='Approved.Symbol'
GWAS3_split <- strsplit(GWAS3$`REPORTED GENE(S)`,",")
Genes <- unique(unlist(GWAS3_split))
Genes_table <- data.frame(Genes = Genes)
match_Hugo <- semi_join(Genes_table,Hugo_names,by=c("Genes"="Approved.Symbol"))

Drug_Gene <- read.table('Drug_Gene_table.txt',sep='\t',header=TRUE)
match_DrugGene <- semi_join(Drug_Gene,match_Hugo,by=c("entrez_gene_symbol"="Genes"))

d <-  data.frame(wanted_disease)
wanted_genes_table <-  semi_join(GWAS3,d,by=c("DISEASE/TRAIT"="wanted_disease"))
all_genes <- unlist(strsplit(wanted_genes_table$`REPORTED GENE(S)`,","))
Disease_Genes <- data.frame(intersect(all_genes,match_Hugo$Genes))

combDrug = read_tsv('linkkCT.tsv')
CTD = subset(combDrug, CONDITION %in% c("Crohn Disease","Crohn's Disease","Crohns Disease"))
colnames(Disease_Genes) <- "names"
Potential_Drug <- merge(match_DrugGene,Disease_Genes, by.x="entrez_gene_symbol",by.y="names")

x = CTD$INTERVENTION_NAME
y = Potential_Drug$drug_name
Drug_Result = subset(y, !y%in%x)
