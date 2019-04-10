# dm_loaddata.R
# Ken McGarry updated: 6/6/18
# load data and pre-process, removed some columns manually from excel files and some using R, also renamed variables

library(tidyverse)
memory.limit(1010241024*1024) # use more RAM memory (10 GBs)
cat("\nLoading files....")

#drug_targets <- load_drugtargets()  # loads in and preprocesses the drug.targets file.
#drug_targets <- drug_targets[!(duplicated(drug_targets[c("DrugName","Gene")]) | duplicated(drug_targets[c("DrugName","Gene")], fromLast = TRUE)), ]

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12K unique proteins with approx 60K interactions between them. http://hint.yulab.org/
# Downloaded on 6/6/18
ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2018.csv", header=TRUE,stringsAsFactors = FALSE)
# Convert gene names all to uppercase
ppi_hint <- mutate_all(ppi_hint, funs(toupper))
# get rid of duplicate A and B columns - NB for some reason duplicated() function cannot find them!!!
ppi_hint <- ppi_hint %>% 
  dplyr::filter(Gene_A != Gene_B)

# add bioplex ppi data to hint ppi data
# http://bioplex.hms.harvard.edu/
bioplex <- read.csv(file="C://R-files//proteins//BioPlex.csv", header=TRUE, sep=",")
ppi <- rbind(ppi_hint,bioplex) 
ppi <- ppi %>% 
  dplyr::filter(Gene_A != Gene_B)
ppi_hint <- ppi

# remove genes with multiple entries..
ppi_hint <- filter(ppi_hint, !grepl("\\|",Gene_A))
ppi_hint <- filter(ppi_hint, !grepl("\\|",Gene_B))
ppi_hint <- filter(ppi_hint, !grepl("-",Gene_A))
ppi_hint <- filter(ppi_hint, !grepl("-",Gene_B))

rm(bioplex,ppi)
# After, the paired data are sent to igraph - we end up with 15,554 proteins with 109,746 edges

disgene <- file.path('C://R-files//DiseaseModules//','gene_disease_associations.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
disgene <- data.frame(lapply(disgene, as.character), stringsAsFactors=FALSE) # needs to be strings NOT factors
disgene$diseaseId <- gsub("umls:","",disgene$diseaseId) # get rid of the bloody "umls:" from disgene for good!
names(disgene)[names(disgene)=="diseaseId"] <- "umls"
disgene <- filter(disgene,originalSource=="BEFREE")  # most reliable source of data
disgene <- disgene[!(duplicated(disgene[c("umls","geneSymbol")]) | duplicated(disgene[c("umls","geneSymbol")], fromLast = TRUE)), ]
disgene <- filter(disgene,score > 0.05)  # use only "high" confidence associations, majority are pitifully low
# Using settings of "BEFREE" and 0.1 we get 1,637 disease genes

indications  <- file.path('C://R-files//sider', 'indications.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
names(indications)[names(indications)=="umls_cui_from_label"] <- "umls"
indications <- indications[,c(1:2,5,7,9,10)]

meshtree <- file.path('C://R-files//disease//','meshtreefull.csv') %>% read.delim(na.strings='',sep=',',header=TRUE,comment.char="#")
meshtree <- data.frame(lapply(meshtree, as.character), stringsAsFactors=FALSE) # needs to be strings NOT factors
meshtree <- meshtree[,1:3]

meshtree <- meshtree %>% 
  filter(str_detect(MeSH, "C"))  # use stringr to seek first char "C" in each row ("C" is for diseases)

# create data frame containing top level disease categories (C1 - C26)
meshcats <- meshtree %>% 
  filter(str_length(MeSH) == 3)  # Top level categories have three chars

# need to convert between meshcodes and their ids, umls codes for diseases  
mappings <- file.path('C://R-files//disease//','map_umls_2_mesh.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#")
mappings <- data.frame(lapply(mappings, as.character), stringsAsFactors=FALSE) # needs to be strings NOT factors
mappings$umls <- gsub("umls:","",mappings$umls) # get rid of the bloody "umls:" from mappings for good!

# Convert id's ("D008105") into umls codes (umls:C0238065) using mappings structure
umls=""
for(i in 1:nrow(meshtree)){
  if(length(mappings[which(mappings[,3] == meshtree[i,2]),1])>0) # check for zero length returns when umls maps dont exist
    umls[i] <- mappings[which(mappings[,3] == meshtree[i,2]),1]
  else
    umls[i]="UNKNOWN" # 
}

maindata <- data.frame(umls,meshtree,stringsAsFactors=FALSE) # add the umls to mesh ids
maindata <- filter(maindata,umls != "UNKNOWN")  # remove entries without umls code
maindata$umls <- gsub("umls:","",maindata$umls) # get rid of the bloody "umls:" from maindata for good!

# Table of freq counts for each disease type
tempdata <- stringr::str_extract(maindata$MeSH,"^.{3}")  # str_extract doesnt work with pipes!
table(tempdata)

# loads in and preprocesses the drug.targets file.
drug_targets <- load_drugtargets()  
drug_targets <- drug_targets[!(duplicated(drug_targets[c("DrugName","Gene")]) | duplicated(drug_targets[c("DrugName","Gene")], fromLast = TRUE)), ]

# Remove small quantity proteins; Adhesion; Nuclear Other; Antibody; CD Molecules; Ribosomal; Cytokine; Surface Antigen; Membrane other
drug_targets <-  # Only keep protein target types with at least 50 occurences
  drug_targets %>%
  add_count(TargetClass,sort=TRUE) %>%
  filter(n > 50)

# read in protein classification based on pharos database (it has the protein families),
# assign them to ppi network for classification algorithims.
protein_class <- read.csv("c:\\R-files\\proteins\\pharos_v4.6.2.csv", header=TRUE,stringsAsFactors = FALSE,na.strings=c("", "NA"))
protein_class <- protein_class[,c(3,8)]
names(protein_class)[names(protein_class)=="DTO.Family"] <- "TargetClass"
names(protein_class)[names(protein_class)=="HGNC.Sym"] <- "Gene"
protein_class <- protein_class[!(duplicated(protein_class[c("TargetClass","Gene")]) | duplicated(protein_class[c("TargetClass","Gene")], fromLast = TRUE)), ]
protein_class <- na.omit(protein_class)

protein_class$TargetClass <- gsub('TF; Epigenetic', 'TF', protein_class$TargetClass)
protein_class$TargetClass <- gsub('oGPCR', 'GPCR', protein_class$TargetClass)




