#Package readr is used for reading a tsv file.
#Package dplyr is used for doing some data manipulation functions.
library(readr)
library(dplyr)
#Place the disease’s name here. Note that the name must exact match the one shown in the original GWAS table.
wanted_disease <-  c("Crohn's disease")

#Started with 24455 rows of data from the original GWAS table, we first eliminated associations annotated as not replicated by filtering out all NA values in “replication sample description” column. Then we discarded the data with p-value greater than 10-7. After that, we also filtered out anthropometric traits, referred from the paper’s supplementary table 1, since they are not related to drugs discovery. The function used for this last step is called anti_join from dplyr package.

#After this step, we ended up with 3907 rows of data.
GWAS<-read_tsv('gwas_catalog.tsv')
GWAS1<-filter(GWAS, GWAS$'REPLICATION SAMPLE DESCRIPTION'!="NA")
GWAS2<-filter(GWAS1,GWAS1$'P-VALUE'<=(10^(-7)),header=TRUE)
anth_traits<-read.csv('supplementary_table_1.csv')
GWAS3<- anti_join(GWAS2,anth_traits,by=c("DISEASE/TRAIT"="GWAS.Disease.Trait"))

#The next step is to match gene names in our current table with the recognizable HUGO gene names. The function used in this step is called semi_join form dplyr package.

#After this step, we ended up with 2328 unique genes.
Hugo_names <-read.csv('genes.csv',sep='\t', header=FALSE)
colnames(Hugo_names)='Approved.Symbol'
GWAS3_split <- strsplit(GWAS3$`REPORTED GENE(S)`,",")
Genes <- unique(unlist(GWAS3_split))
Genes_table <- data.frame(Genes = Genes)
match_Hugo <- semi_join(Genes_table,Hugo_names,by=c("Genes"="Approved.Symbol"))

#Next, we obtained a drug-gene-interaction table from DGIdb. This table provides information of drugs that correspond to each gene. Among all the genes in this table, we are only interested in ones that are in our filtered GWAS table, so we filtered the rest out.

#After this step, we ended up with 671 unique genes that correspond to 5772 drugs.
Drug_Gene <- read.table('Drug_Gene_table.txt',sep='\t',header=TRUE)
match_DrugGene <- semi_join(Drug_Gene,match_Hugo,by=c("entrez_gene_symbol"="Genes"))

#This step is the one that we get the corresponding drugs for the specific disease we chose.

d <-  data.frame(wanted_disease)
wanted_genes_table <-  semi_join(GWAS3,d,by=c("DISEASE/TRAIT"="wanted_disease"))
all_genes <- unlist(strsplit(wanted_genes_table$`REPORTED GENE(S)`,","))
Disease_Genes <- data.frame(intersect(all_genes,match_Hugo$Genes))

#Finally, we compared the corresponding drugs we acquired from the previous step with drugs from ClinicalTrials database. Those that are not overlap between these two are our new potential drugs.

combDrug = read_tsv('linkkCT.tsv')
CTD = subset(combDrug, CONDITION %in% c("Crohn Disease","Crohn's Disease","Crohns Disease"))
colnames(Disease_Genes) <- "names"
Potential_Drug <- merge(match_DrugGene,Disease_Genes, by.x="entrez_gene_symbol",by.y="names")

x = CTD$INTERVENTION_NAME
y = Potential_Drug$drug_name
Drug_Result = subset(y, !y%in%x)
