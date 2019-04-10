# dm_functions.R
library(dplyr) 
library(ChemmineR)
library(ChemmineOB)
library(ape)
#library(sparcl)
library(cluster) # used for kmeans and silhoutte plot
library(xtable)
library(gplots) 
library(scatterplot3d) 
library(igraph)
library(ROCR)
library(VennDiagram)
library(ggplot2)
library(linkcomm)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library("AnnotationDbi")
library(GOplot)
library(scales)
library(rentrez)
library(stringr)
# huge amount of data gets loaded in here onwards!
library(ontologySimilarity)
library(ontologyIndex)
library(infotheo)
library(KEGGprofile)
library(KEGG.db)
library(Matrix)
data(go)
data(gene_GO_terms)
data(GO_IC)

# source("https://bioconductor.org/biocLite.R")
# biocLite("ReactomePA")

## --------------------- FUNCTION DEFINITIONS -----------------------

# Makes first letter of string uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# the cosine formula gives a similarity. It is 1 if the vectors are pointing in the same direction. 
# Distance measures need the value to be 0 when vectors are the same so 1 - similarity = distance. 
# Many uses need distance rather than similarity 
cosineDist1 <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}


# This vesrion uses the dot product matrix, to compute the cosine similarity matrix: 
# Input S is the matrix of dot product, dt is dataset. Divide the dot product by the norms of vectors. 
cosineDist2 <- function(S){
  doc_norm <- apply(as.matrix(dt),1,function(x) norm(as.matrix(x),"f")) 
  divide_one_norm <- S/doc_norm 
  cosine <- t(divide_one_norm)/doc_norm
  return (cosine)
}


# get_drug_names() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drug_names <- function(umls,rlist) {
  umls <- gsub("umls:","",umls) # get rid of the bloody "umls:" in the string
  cat("\nCode is..",umls)
  ilist <- filter(indications, umls_cui_from_meddra == umls) 
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  if(length(ilist) > 0){
    for (j in 1:length(ilist)){
      cat("\ndrug",j,"is", ilist[j])
    }
    return(ilist)
  }else{
    cat("\n","Sorry, no drugs found...check umls code is correct for your disease")
    return(NULL)}
}

# get_drugs_plus() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# The restricted list of drugs we cant use is passed to this function.
# It returns the drug names plus other information.
get_drugs_plus <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls) 
  mydrugs <- ilist # temp storage
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  
  if(length(ilist) == 0){  # stuck in a 2nd null test as zero length data creeping in and causing crashes.
    cat("\n","No drugs found...check umls code is correct for your disease or maybe no drugs for this disease")
    return(NULL)}
  
    options(warn=-1)
    ilist <- filter(mydrugs,drugbank_name == ilist) # causes warning message
    options(warn=0)
    ilist <- ilist[!duplicated(ilist$drugbank_name),]
    
  if(nrow(ilist) > 0){
    for (j in 1:nrow(ilist)){
      cat("\ndrug",j,"is", ilist[j,2])
    }
    ilist <- select(ilist,drugbank_id,drugbank_name,umls_cui_from_meddra,meddra_name)
    return(ilist)
  }else{
    cat("\n","No drugs found...check umls code is correct for your disease or maybe no drugs for this disease")
    #cat("\n i is...",j)
    return(NULL)}
}

# Convert the ID to umls code in order to access indications and obtain drug_id and drug_name, 
# digestive dataframe is passed as id; use ID to match with MeshID and get umls
id2umls <- function(id){
  x <- vector(mode="character",length=nrow(id))
  y <- mappings[1,] # instantiate a temporary vector, 
    
  for (i in 1:nrow(id)){
    x[i] <- id$ID[i] 
    tempy <- filter(mappings, meshId == x[i])  # if tempy has meshid then keep, else ignore and do save it
    if(nrow(tempy) > 0){
      y[i,] <- tempy[1,]
      y[i,]$umls <- gsub("umls:","",y[i,]$umls) # remove "umls:" from string
      #cat("\ni is ...",i)
    }
   
  }
  y <- na.omit(y)  # get rid of records with NA where no umls Id's exist for the Mesh id's
  return(y)
}


# Adds the MESH code into drugs_list dataframe, useful info when relating drugs to disease categories
# 
add_meshcode <- function(drugs){
  x <- vector(mode="character",length=nrow(drugs))
  y <- vector(mode="character",length=nrow(drugs))
  #drugplus <- 0
  
  for (i in 1:nrow(drugs)){
    tempx <- paste("umls:",drugs[i,3],sep="") # add the bloody "umls:" back in, and some magic numbers.
    tempy <- filter(mappings, umls == tempx) 
    tempz <- filter(digestive, ID == tempy$meshId)
    x[i] <- tempz$MeSH
    y[i] <- tempz$ID
  }
  
  drugsplus <- cbind(drugs,x)   # add the x or Mesh values
  drugsplus <- cbind(drugsplus,y)   # add the y or ID values
  colnames(drugsplus)[5] <- "MeSH"   # change from x and y to better names
  colnames(drugsplus)[6] <- "ID"
  
  return(drugsplus)
}


# R provides a tail and head command to view last six and first six elements, so why not the middle six?
middle <- function(mydata) {
  len <- nrow(mydata)
  startpoint <- round(len/2)
  endpoint <- startpoint+5
  mydata[startpoint:endpoint,]
  
}


# Supply get_disease_genes() with a "drug_list" which has "umls_cui_from_meddra" field to tie 
# with "diseaseId" in "disgene", For each disease see what genes are implicated and return dataframe
get_disease_genes <- function(mydrugs){
  implicated <-0 # instantiate, unfortunatly makes a zero entry, delete this at end of function.
  
  disgene$diseaseId <- gsub("umls:","",disgene$diseaseId) # get rid of the bloody "umls:" from disgene for good!
  mydrugs <- mydrugs[!duplicated(mydrugs[,c('umls_cui_from_meddra','meddra_name')]),]  # just keep unique diseases
  mydrugs <- select(mydrugs,umls_cui_from_meddra,meddra_name, MeSH,ID)  # keep name and ids etc
  
  for (i in 1:nrow(mydrugs)){
    tempx <- mydrugs[i,1]  # seek out the umls code
    tempy <- filter(disgene, diseaseId == tempx) 
    if(nrow(tempy) > 0){
      cat("\nFound gene(s)!")
      tempy <- select(tempy,diseaseId, geneName, diseaseName)
      implicated <- rbind(tempy,implicated)
    }
  }
  
  implicated <- implicated[-nrow(implicated),]            # last entry is zero so remove it
  return(implicated)
}


# Convert drugbank IDs (DB00035) to drugnames (Desmopressin)
ID2name <- function(DBid){
  thenames <- vector(mode="character",length=length(DBid))
  
  for (i in 1:length(DBid)){
    dname <- filter(indications, drugbank_id == DBid[i])
    dname <- dname[!duplicated(dname[,'drugbank_name']),]
    thenames[i] <- dname$drugbank_name
  }
  return(thenames)
}

# plot_chemsim() creates plots generated from the data computed by drugstructure_gi.R code .
plot_chemsim <- function(){
  plot.new()
  
  hc <- hclust(as.dist(1-simMA), method="complete") 
  #par(cex=0.1)
  heatmap.2((1-simMA), Rowv=as.dendrogram(hc), 
            Colv=as.dendrogram(hc), 
            #col=greenred(10),
            keysize = 2,
            key=TRUE,
            #col=bluered(256),
            col=colorpanel(40, "white","yellow", "darkblue"), 
            density.info="none", trace="none")
  
  # Creates the similarity score matrix and cluster them.  
  colnames(simMA)<-drugnames
  rownames(simMA)<-drugnames
  hc <- hclust(as.dist(1-simMA), method="complete") 
  plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=FALSE)
  
  cl <- kmeans(simMA,10,nstart=50) #cl <- kmeans(simMA,10,nstart=5)
  sk <- silhouette(cl$cl,dist(simMA))
  plot(sk)
  
  y <- cutree(hc,25) #10
  par(cex=0.8)
  ColorDendrogram(hc,y=y,labels=drugnames,branchlength = 0.7,cex = 0.7)  
  
}


# get_linked_diseases() that are not C06 but connected to C06 by shared genes. Modified to include all
# diseases associated with a disorder, 30/11/17 Need to inclued drugs associated with disorders.
get_linked_diseases <- function(dgenes){
  linked_diseases <- disgene[1,] # instantiate before use
  
  for (i in 1:length(dgenes)){
    #cat("\nlength of dgenes is ",length(dgenes))
    gene <- dgenes[i] # get genes individually and see what disorders they linked with
    glist <- filter(disgene, geneName == gene)  # This bit is OK
    if(nrow(glist) > 0){
      linked_diseases <- rbind(linked_diseases,glist) }
  }
  
  linked_diseases <- arrange(linked_diseases,diseaseName)  # sort alphabetically
  linked_diseases <- linked_diseases[,c(1,2,4,6)]            # keep only key variables
  # The complicated line below removes duplicate entries, there are quite a few and Im not sure how they got in.
  # but if diseasename and gene name in the same row occur then keep only one copy,
  linked_diseases  <- linked_diseases[!(duplicated(linked_diseases[c("diseaseName","geneName")]) | duplicated(linked_diseases[c("diseaseName","geneName")], fromLast = TRUE)), ]
  linked_diseases <- linked_diseases[-nrow(linked_diseases),]     # last entry is zero so remove it
  linked_diseases <- anti_join(linked_diseases, disease_umls, by="diseaseName") # If C06 disorders appear , remove them.
  
  return(linked_diseases)
}


# input  a gene or list of genes and get all diseases asscoiated with these genes 
# really best just using JUST one gene.
get_all_linked_diseases <- function(dgenes){
  all_diseases <- disgene[1,] # instantiate before use
  
  for (i in 1:length(dgenes)){
    gene <- dgenes[i]
    glist <- filter(disgene, geneName == gene)
    if(nrow(glist) > 0){
        all_diseases <- rbind(glist,all_diseases)
    }
  }
  
  all_diseases <- all_diseases[-nrow(all_diseases),]     # last entry is fake so remove it
  all_diseases <- all_diseases[!duplicated(all_diseases[,'diseaseName']),]   # get rid of the many duplicates
  all_diseases <- select(all_diseases,diseaseId,geneId,geneName,diseaseName)  # drop "score" variable
  all_diseases <- arrange(all_diseases,diseaseName)  # sort alphabetically
  return(all_diseases)
}

# print tables in LaTex format for inclusion into paper
print_tables <- function(){
  temp_table <- disgene[,c(1,4:6)]
  
  tli.table <- xtable(head(temp_table))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable((digestive[86:101,]))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(head(digestive))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(middle(drug_list))
  #digits(tli.table)[c(2,6)] <- 0
  print(tli.table,floating=FALSE)
 
  tli.table <- xtable(filter(drug_list, MeSH == "C06.405"))
  print(tli.table,floating=FALSE)
  
  tli.table <- xtable(filter(gene_list, diseaseId == "C0017178"))
  print(tli.table,floating=FALSE)
  
}

# goanalysis() will enrich a gene with GO terms
# depends on clusterprofiler library and several other things...
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis
go_analysis <- function(yourgenes,ontotype){
  cat("\n",yourgenes)
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ego <- enrichGO(gene          = eg[,2],
                  #universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = ontotype, # one of CC, BP or MF
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  return(ego)
}



# getDiseaseModules() pass it the linkcomm structure, the appropriate data will be
# converted into a dataframe - ready to be passed to GO enrichment functions before
# sending it to GObubble for analysis and display.
getDiseaseModules <- function(linkdata,batch){
  category <- c("MF","BP","CC")
  enrich <- c("GO:0017091", "AU-rich element binding","1/1", "23/16982", "0.001354375","RU12","FU",99)
  
  # remove modules with fewer than 20 genes - as per Menche 2015 paper
  linkdata$clusters <- Filter(function(x)length(x) > 20, linkdata$clusters)
  cat("\nFound ",length(linkdata$clusters), " usable modules.")

  if(is.null(batch)){
    cat("\nERROR: you need to enter ''all'' for all modules or a range in quotes e.g. ''77:89''")
    return(NULL)}
  
  if(batch =="all") {  # process all diseasemodules if batch ==all
    indexstart <- 1
    indexend <- length(linkdata$clusters)
  }else{
    tempy <- unlist(strsplit(batch,":"))  # its a range
    indexstart <- as.numeric(tempy[1])
    indexend <- as.numeric(tempy[2])
  }
  
  #for (i in 1:length(linkdata$clusters)){        # i=num of disease modules
  temp_i <- vector(mode="integer", length=indexend-indexstart);
  z <-1
  
  for (i in indexstart:indexend){        # i=num of disease modules
    items <- getNodesIn(linkdata, clusterids = i)
    temp_i[z] <- i
    z <- z +1
    for (k in 1:length(items)){              # k=num genes in disease module
      for (j in 1:length(category)){  # enrich from MF, BP and CC
        temp_enrich <- go_analysis(items[k],category[j])
        if(!is.null(temp_enrich) && nrow(temp_enrich)>0){
          temp_enrich <- temp_enrich[,c(1:5,8)]
          temp_enrich[,7] <- category[j]       # add the category e.g MF
          temp_enrich[,8] <- i              # add the disease module number (i.e. cluster number)
          enrich <- rbind(enrich,temp_enrich)}
      }
    }
  }
  
  # Fix the dataframe: add z-score, rename V7, add disease module number
  # To match what GOplot i.e. GObubble expects to find
  names(enrich)[names(enrich)=="Description"] <- "term"
  names(enrich)[names(enrich)=="pvalue"] <- "adj_pval"
  names(enrich)[names(enrich)=="V7"] <- "category"
  names(enrich)[names(enrich)=="V8"] <- "DiseaseModule"
  names(enrich)[names(enrich)=="geneID"] <- "genes"
  # convert p-value into "zscore", rename "pvalue" to "adj_pval" 
  namevector <- "zscore"
  enrich[ , namevector] <- qnorm(1 - as.numeric(enrich$adj_pval)/2)
  namevector <- "logFC"
  enrich[ , namevector] <- runif(nrow(enrich), -2, 3)#10#qnorm(1 - as.numeric(enrich$adj_pval)/2)
  namevector <- "count"
  enrich[ , namevector] <- 0  # Number of genes attached to this term.
  enrich$adj_pval <- as.numeric(enrich$adj_pval)
  enrich <- enrich[-1, ]     # 1st entry is rubbish so remove it
 
  # Set "count" for each term
  enrich <- setcount(enrich,temp_i)

  return(enrich)
}

# setcount() gets a count of the terms assigned to each disease module.
setcount <- function(dms,ind){
  #countn <- unique(dms$DiseaseModule) # How many disease modules are there?
  countn <- length(ind) # How many disease modules are there? based on index range?
  cat("\nThere are ",length(ind)," disease modules..numbered from",ind)
  for (j in 1:length(ind)){
    cat("\nJ is now",j)
    temp_dms <- filter(dms,DiseaseModule == (ind[j]))
    nterm <- unique(temp_dms$term) # How many unique terms do we have for this disease module?
    for (k in 1:length(nterm)){
      tcount <- nrow(filter(dms,term == nterm[k]))
      dms$count[dms$term == nterm[k] & dms$DiseaseModule == (ind[j])] <- tcount
      #cat("\nFor DM",ind[j]," tcount is ",tcount)
    }
  }
  
  return(dms)
} 

# This is a far quicker version of getDiseaseModules(). This version uses the ontologySimilarity packages
# by Daniel Green. Unfortunatley, Ive hand coded a lot hence its a big function.
createDiseaseModules <- function(linkdata){
  tempgenes <- names(gene_GO_terms)
  enrich <- data.frame(ID="GO:0000666", genes="RU12",DiseaseModule=666,adj_pval=0.001,zscore=6.001,logFC=0.2,
                       category="FU",term="Satanic like behaviour",stringsAsFactors=FALSE) #instantiate.
  # remove modules with fewer than 10 genes - as per Menche 2015 paper - 20 is too restrictive
  newclusters <- Filter(function(x)length(x) > 10, linkdata$clusters)
  cat("\nFound ",length(newclusters), " usable modules.")
  j<- 0  # set counter for clusters bigger than 10 
  
  for (i in 1:length(linkdata$clusters)){
    tempnodes <- getNodesIn(linkdata, clusterids = i,type="names")
    if(length(tempnodes) >= 10){
      j <- j+1
      #cat("\nj =...",j)
      tempnodes <- tempnodes[tempnodes %in% tempgenes]  # remove genes that do not exist in GO data 
      tempgo <- gene_GO_terms[tempnodes]
      cc <- go$id[go$name == "cellular_component"]
      bp <- go$id[go$name == "biological_process"]
      mf <- go$id[go$name == "molecular_function"] 
      temp_cc <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=cc, x))
      temp_bp <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=bp, x))
      temp_mf <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=mf, x))
      tmp_enrich_c <- data.frame(unlist(temp_cc),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...CC")
      tmp_enrich_c <- cbind(tmp_enrich_c,rownames(tmp_enrich_c),stringsAsFactors=FALSE) # gene names
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(j,length(unlist(temp_cc))))               # diseasemodule number
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(0.01,length(unlist(temp_cc))))            # adj_pval
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(3.2,length(unlist(temp_cc))))             # zscore
      tmp_enrich_c <- cbind(tmp_enrich_c,rep(0.2,length(unlist(temp_cc))))             # logFC
      tmp_enrich_c <- cbind(tmp_enrich_c,rep("CC",length(unlist(temp_cc))),stringsAsFactors=FALSE) # category
      tmp_enrich_c <- cbind(tmp_enrich_c,unname(go$name[unlist(temp_cc)]),stringsAsFactors=FALSE)
      
      tmp_enrich_b <- data.frame(unlist(temp_bp),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...BP")
      tmp_enrich_b <- cbind(tmp_enrich_b,rownames(tmp_enrich_b),stringsAsFactors=FALSE) # gene names
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(j,length(unlist(temp_bp))))               # diseasemodule number
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(0.01,length(unlist(temp_bp))))            # adj_pval
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(3.2,length(unlist(temp_bp))))             # zscore
      tmp_enrich_b <- cbind(tmp_enrich_b,rep(0.2,length(unlist(temp_bp))))             # logFC
      tmp_enrich_b <- cbind(tmp_enrich_b,rep("BP",length(unlist(temp_bp))),stringsAsFactors=FALSE) # category
      tmp_enrich_b <- cbind(tmp_enrich_b,unname(go$name[unlist(temp_bp)]),stringsAsFactors=FALSE)
      
      tmp_enrich_m <- data.frame(unlist(temp_mf),stringsAsFactors=FALSE)  # GO ID's
      #cat("\nCBIND...MF")
      tmp_enrich_m <- cbind(tmp_enrich_m,rownames(tmp_enrich_m),stringsAsFactors=FALSE) # gene names
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(j,length(unlist(temp_mf))))               # diseasemodule number
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(0.01,length(unlist(temp_mf))))            # adj_pval
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(3.2,length(unlist(temp_mf))))             # zscore
      tmp_enrich_m <- cbind(tmp_enrich_m,rep(0.2,length(unlist(temp_mf))))             # logFC
      tmp_enrich_m <- cbind(tmp_enrich_m,rep("MF",length(unlist(temp_mf))),stringsAsFactors=FALSE) # category
      tmp_enrich_m <- cbind(tmp_enrich_m,unname(go$name[unlist(temp_mf)]),stringsAsFactors=FALSE)
      
      colnames(tmp_enrich_b)[1] <- "ID"; 
      colnames(tmp_enrich_b)[2] <- "genes"; 
      colnames(tmp_enrich_b)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_b)[4] <- "adj_pval"; 
      colnames(tmp_enrich_b)[5] <- "zscore"; 
      colnames(tmp_enrich_b)[6] <- "logFC"; 
      colnames(tmp_enrich_b)[7] <- "category" 
      colnames(tmp_enrich_b)[8] <- "term" 
      rownames(tmp_enrich_b) <- c()
      
      colnames(tmp_enrich_c)[1] <- "ID"; 
      colnames(tmp_enrich_c)[2] <- "genes"; 
      colnames(tmp_enrich_c)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_c)[4] <- "adj_pval"; 
      colnames(tmp_enrich_c)[5] <- "zscore"; 
      colnames(tmp_enrich_c)[6] <- "logFC"; 
      colnames(tmp_enrich_c)[7] <- "category" 
      colnames(tmp_enrich_c)[8] <- "term" 
      rownames(tmp_enrich_c) <- c()
      
      colnames(tmp_enrich_m)[1] <- "ID"; 
      colnames(tmp_enrich_m)[2] <- "genes"; 
      colnames(tmp_enrich_m)[3] <- "DiseaseModule" 
      colnames(tmp_enrich_m)[4] <- "adj_pval"; 
      colnames(tmp_enrich_m)[5] <- "zscore"; 
      colnames(tmp_enrich_m)[6] <- "logFC"; 
      colnames(tmp_enrich_m)[7] <- "category" 
      colnames(tmp_enrich_m)[8] <- "term" 
      rownames(tmp_enrich_m) <- c()
      
      #cat("\nRBIND....CC , BF & MF to enrich")
      enrich <- rbind(enrich,tmp_enrich_c)
      enrich <- rbind(enrich,tmp_enrich_b)
    }
  }
  enrich <- enrich[-1, ]     # 1st entry is rubbish so remove it
  tmp_count <- table(enrich$ID)  # obtain a count of how many times each GO ID appears.
  idcount <- vector(mode="integer",length=nrow(enrich))
  
  idcount <- as.vector(tmp_count)[match(enrich$ID, names(tmp_count))] # almost too clever
  enrich <- cbind(enrich,idcount)
  colnames(enrich)[9] <- "count"   

  return(enrich)
}

# fix_C06() join C06 disease names with MeSH codes - its part manual process due to missing values
# C06 diseases joined with genes and drugs.
fix_C06 <- function(){
  C06Disease <- unique(gene_list$diseaseName)
  C06_ID <- vector(mode="character",length=length(C06Disease))
  
  for (i in 1:length(C06Disease)){
    tempy <- filter(digestive,Term == C06Disease[i])
    if(nrow(tempy) > 0){
      C06_ID[i] <-tempy$MeSH}
    if(nrow(tempy) == 0){
      C06_ID[i] <- "manual entry required"
    }
  }
  
  # Fix missing entries by hand
  C06_ID[9] <- "C06.301.623" # Liver Neoplasms
  C06_ID[2] <- "C06.552.630.400" # Liver Cirrhosis, Biliary 
  C06_ID[8] <- "C06.405.748.789"  # Malignant neoplasm of stomach Stomach
  C06_ID[22] <- "C06.405.205.731.249"  # Colitis, Ulcerative [C06.405.205.731.249]
  C06_ID[27] <- "C06.405.205" # 
  C06_ID[32] <- "C06.405.469.275.800.849"  # Stomach Ulcer [C06.405.469.275.800.849]
  C06_ID[34] <- "C06.405.469.637"  #Malabsorption Syndromes [C06.405.469.637]
  C06_ID[38] <- "C06.552" #liver diseases
  C06_ID[53] <- "C06.552.830.150" #     Porphyria, Acute Intermittent 
  C06_ID[17] <- "C06.405.117.119.500.484" # Gastroesophageal Reflux 
  C06_ID[10] <- "C06.552.697.160" #"Liver carcinoma"
  C06_ID[11] <- "C06.301.761.249.500" # Insulinoma 
  
  disease <- as.data.frame(cbind(C06Disease,C06_ID),stringsAsFactors=FALSE)
  return(disease)
}

# Disease ontology enrichment - based on the gene list it provides a list of
# diseases associated with these genes.
# https://bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/semanticAnalysis.html
DO_analysis <- function(mygenes){
  #data(geneList)
  mygenes <- bitr(mygenes,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- enrichDO(gene          = mygenes$ENTREZID,
              ont           = "DO",      pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",     #universe      = names(geneList),
              minGSSize     = 5,        maxGSSize     = 500,
              qvalueCutoff  = 0.05,     readable      = FALSE)
  head(x)
  upsetplot(x,n=20,main.bar.color="blue", sets.bar.color="gray",text.scale=1.5,order.by = c("freq", "degree"))
  enrichMap(x,n=50,fixed=FALSE)  # network view using TKPLOT()
  return(x)

}

# score_pathways(), scores the associated pathways from KEGG based on GeneRatio.
# must compute some sort of value to rank biological plausibility/importance
score_pathways <- function(dm){
  # How many disease mods do we have?
  
  kegpaths <- kegg_analysis(nonC06_nsc$geneName)
  score_path <- vector(mode="integer",length(nrow(dm))) #how may pathways do we have?
  for (i in 1:nrow(dm)){
    temp <- unlist(strsplit(dm[i]$GeneRatio,"/")) # split on the backslash symbol
    score_path[i] <- dm[i]$Count/as.numeric(temp[2])  # calculate the actual GeneRatio
  }
  return(score_path)
}

# KEGG over-representation test
kegg_analysis <- function(yourgenes){
  #cat("\nyourgenes are: ",yourgenes)
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db");
  eg<-eg[,2]
  kk <- enrichKEGG(gene= eg, organism= 'hsa', pvalueCutoff = 0.05)
  
  return(kk)
}

# Ranks the most salient disease modules based on GO annotations.
# 27/12/17
rank_alldm_go <- function(dm){
  dm <- dm[dm$ID %in% go$id,] # ensure missing GO terms are removed
  dm <- dm[dm$ID %in% attributes(GO_IC)$name,] # ensure missing IC terms are removed
  countdm <- length(unique(dm$DiseaseModule))
  listdm <- unique(dm$DiseaseModule)
  #cat("\nFound ",countdm,"Disease Modules.")
  for (i in 1:countdm){
    tempmod <- filter(dm,DiseaseModule == listdm[i])
    tempMF <- filter(tempmod,category =="MF")
    #cat("\ndismod has ",nrow(tempMF)," MF")
    cat("\ndismod",listdm[i], "has",length(unique(tempmod$genes))," genes and ",(table(tempmod$category))," GO annotations")
  }
  terms_by_disease_module <- split(dm$ID,dm$DiseaseModule)  # do split by disease module
  terms_by_disease_module <- unname(terms_by_disease_module)   # Remove names for the moment
  sim_matrix <- get_sim_grid(ontology=go,information_content=GO_IC,term_sets=terms_by_disease_module)
  # see how the disease modules cluster
  dist_mat <- max(sim_matrix) - sim_matrix  # need a distance matrix, not a similarity matrix
  #clusterdetails <- hclust(as.dist(dist_mat),"ave")
  #plot(hclust(as.dist(dist_mat)))
  return(sim_matrix)
}

# rank_alldm_pathways(), scores the associated pathways from KEGG based on GeneRatio.
# 28/12/17
# Barrett Esophagus; Adenomatous Polyposis Coli; Gastrointestinal Stromal Tumors; 
# Cholecystitis ; Cholelithiasis; Primary biliary cirrhosis; Cholestasis
rank_alldm_pathways <- function(dm){
  nmods <- length(unique(dm$DiseaseModule)) 
  listdm <- unique(dm$DiseaseModule)
  score_path <- data.frame(score=rep(0,nmods),module=rep("disease",nmods),stringsAsFactors=FALSE); 
  cat("\nWe have ",nmods,"disease modules")# How many disease mods do we have?
  for (j in 1:nmods){
    cat("\nprocessing module ",j," which is...",listdm[j])
    temp <- unlist(strsplit(listdm[j],"_")) # split on the 'underscore' symbol
    tmpmod <- filter(disgene1, grepl(temp[1],diseaseName))
    thegenes <- unique(tmpmod$geneName)
      if(length(thegenes) > 1){
       subsetgenes <- sample(thegenes,0.7*length(thegenes))}else{
       subsetgenes <- thegenes
      }
   
    kegpaths <- kegg_analysis(subsetgenes)
    score_path[j,2] <- listdm[j]
    if(!is.null(kegpaths)){
      tmpscore <- rep(0,nrow(kegpaths))
      for (i in 1:nrow(kegpaths)){
        temp <- unlist(strsplit(kegpaths[i]$GeneRatio,"/")) # split on the backslash symbol
        tmpscore[i] <- kegpaths[i]$Count/as.numeric(temp[2])  # calculate the actual GeneRatio
        #cat("\ntmpscore= ",tmpscore)
        #tmpscore[i] <- nrow(kegpaths)
      } 
      score_path[j,1] <- sum(tmpscore)
      score_path[j,2] <- listdm[j]
    }
   # score_path <- rbind(tmpscore,score_path)
  }
  score_path[is.na(score_path)] <- 0.00   # remove NA where no keggpath could be found
  return(score_path)
}


# 28/12/17
# module_overlap() provides percentage similarity values where the key modules overlap  with other 
# modules from different diseases. Calls up hyper_matrix(). We need lists() of genes based around 
# disease modules. Choose only a select few modules, otherwise table will be unreadable.
module_overlap <- function(){
  
  # Compare non-C06 diseases between themselves for shared genes
  nonC06.list <- list(Alz=nonC06_alz$geneName,
                    Asth=nonC06_asth$geneName,
                    Diab=nonC06_dia$geneName,
                    Hypo=nonC06_hyp$geneName,
                    Scc=nonC06_nsc$geneName,
                    RH=nonC06_ra$geneName,
                    Sch=nonC06_sch$geneName,
                    Park=nonC06_park$geneName,
                    Aut=nonC06_aut$geneName,
                    Obes=nonC06_obs$geneName)
  
  # Build up gene lists for separate C06 diseases
  tmpmod <- filter(disgene1,diseaseName == "Primary biliary cirrhosis")
  pbc_genes <- unique(tmpmod$geneName)
  tmpmod <- filter(disgene1,diseaseName == "Gastrointestinal Diseases")
  gd_genes <- unique(tmpmod$geneName)  
  tmpmod <- filter(disgene1,diseaseName == "Cholecystitis")
  cho_genes <- unique(tmpmod$geneName)  
  tmpmod <- filter(disgene1,diseaseName == "Gastrointestinal Stromal Tumors")
  gst_genes <- unique(tmpmod$geneName)
  tmpmod <- filter(disgene1,diseaseName == "Barrett Esophagus")
  bar_genes <- unique(tmpmod$geneName)  
  tmpmod <- filter(disgene1,diseaseName == "Adenomatous Polyposis Coli")
  adc_genes <- unique(tmpmod$geneName)
  tmpmod <- filter(disgene1,diseaseName == "Inflammatory Bowel Diseases")
  ibs_genes <- unique(tmpmod$geneName)  
  tmpmod <- filter(disgene1,diseaseName == "Celiac Disease")
  cel_genes <- unique(tmpmod$geneName)
  tmpmod <- filter(disgene1,diseaseName == "Crohn Disease")
  cro_genes <- unique(tmpmod$geneName)  
  tmpmod <- filter(disgene1,diseaseName == "Liver carcinoma")
  liv_genes <- unique(tmpmod$geneName)
  
  # Compare shared genes between main C06 diseases
  C06.list <- list(Cirr=pbc_genes,    # "Primary biliary cirrhosis"
                  Gast=gd_genes,         # "Gastrointestinal Diseases"
                  Chol=cho_genes,         # "Cholecystitis" 
                  Gtum=gst_genes,         # "Gastrointestinal Stromal Tumors"
                  Barr=bar_genes,       # "Barrett Esophagus"
                  Adeno=adc_genes,   # "Adenomatous Polyposis Coli"
                  IBS=ibs_genes,   # "Inflammatory Bowel Diseases"
                  Celi=cel_genes,  # "Celiac Disease"
                  Croh=cro_genes,  #"Crohn Disease"   
                  Liver=liv_genes)  # Liver carcinoma
  
  # gene.list <- list() # my version of list
  universeOfGenes <- 1000; # how many unique genes in total set?
  mat_table <- hyper_matrix(nonC06.list, universeOfGenes)
  mat_table <- hyper_matrix(C06.list, universeOfGenes)
  
  # Now compare the key C06 and nonC06 diseases for shared genes
  Compare.list <- list(Barr=bar_genes,       # "Barrett Esophagus"
                       Adeno=adc_genes,   # "Adenomatous Polyposis Coli"
                       IBS=ibs_genes,   # "Inflammatory Bowel Diseases"
                       Celi=cel_genes,  # "Celiac Disease"
                       Croh=cro_genes, 
                       RH=nonC06_ra$geneName,
                       Sch=nonC06_sch$geneName,
                       Park=nonC06_park$geneName,
                       Aut=nonC06_aut$geneName,
                       Obes=nonC06_obs$geneName)
  
  mat_table <- hyper_matrix(Compare.list, universeOfGenes)
  
  return(mat_table)
}


# score_go() give a score to each disease module based on mutual information from similarity matrix
# derived from the GO annotation.
score_go <- function(dm,disease){
  dm <- dm[dm$ID %in% go$id,] # ensure missing GO terms are removed
  dm <- dm[dm$ID %in% attributes(GO_IC)$name,] # ensure missing IC terms are removed
  terms_by_disease_module <- split(dm$ID,dm$DiseaseModule)  # do split by disease module
  terms_by_disease_module <- unname(terms_by_disease_module)   # Remove names for the moment
  sim_matrix <- get_sim_grid(ontology=go,information_content=GO_IC,term_sets=terms_by_disease_module)
# Calculate mutual information from the similarity matrix, provides a score of sorts for each disease module
  nbins <- sqrt(NROW(sim_matrix))
  dat <- infotheo::discretize(sim_matrix,"equalwidth", nbins) # use full package extension
  IXY <- infotheo::mutinformation(dat,method= "emp")
  IXY2 <-infotheo::mutinformation(dat[,1],dat[,2])
  H <- infotheo::entropy(infotheo::discretize(sim_matrix[1,]),method="shrink")

  # see how the disease modules cluster
  dist_mat <- max(sim_matrix) - sim_matrix
  plot(hclust(as.dist(dist_mat)))
  
  for (i in 1:nrow(sim_matrix)){
    cat("\nDisease is", disease, "Module[",i,"] biological value = ",IXY[i])  }
  
  return(IXY)
}


# Group the most salient disease modules and cluster them, this is for the table
# in Latex file. need to get either disease name or MESH id for dendro plot.
score_alldm_go <- function(dm){
  dm <- dm[dm$ID %in% go$id,] # ensure missing GO terms are removed
  dm <- dm[dm$ID %in% attributes(GO_IC)$name,] # ensure missing IC terms are removed
  countdm <- length(unique(dm$DiseaseModule))
  cat("\nFound ",countdm,"Disease Modules.")
  terms_by_disease_module <- split(dm$ID,dm$DiseaseModule)  # do split by disease module
  #terms_by_disease_module <- unname(terms_by_disease_module)   # Remove names for the moment
  sim_matrix <- get_sim_grid(ontology=go,information_content=GO_IC,term_sets=terms_by_disease_module)
  # see how the disease modules cluster
  dist_mat <- max(sim_matrix) - sim_matrix  # need a distance matrix, not a similarity matrix
  #clusterdetails <- hclust(as.dist(dist_mat),"ave")
  #plot(hclust(as.dist(dist_mat)))

  #hc <- hclust(as.dist(1-dist_mat), method="complete") 
  #plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=FALSE)
  #KM <- kmeans(as.dist(1-dist_mat), 15, nstart = 20)
  #optimum_clusters((dist_mat))
  
  return(dist_mat)
}




# Ranks the NEW disease modules, this is for the table in Latex file. 
# 
rank_group_pathways <- function(dm){
  dm <- dm[dm$ID %in% go$id,] # ensure missing GO terms are removed
  dm <- dm[dm$ID %in% attributes(GO_IC)$name,] # ensure missing IC terms are removed
  countdm <- length(unique(dm$newgroup))
  cat("\nFound ",countdm,"New Disease Modules to Rank.")
  terms_by_disease_module <- split(dm$ID,dm$newgroup)  # do split by disease module
  #terms_by_disease_module <- unname(terms_by_disease_module)   # Remove names for the moment
  sim_matrix <- get_sim_grid(ontology=go,information_content=GO_IC,term_sets=terms_by_disease_module)

  # Calculate mutual information from the similarity matrix, provides a score of sorts for each disease module
  nbins <- sqrt(NROW(sim_matrix))
  dat <- infotheo::discretize(sim_matrix,"equalwidth", nbins) # use full package extension
  IXY <- infotheo::mutinformation(dat,method= "emp")
  IXY2 <-infotheo::mutinformation(dat[,1],dat[,2])
  H <- infotheo::entropy(infotheo::discretize(sim_matrix[1,]),method="shrink")
  
  # ADJUST RANK WITH BIOLOGICAL INPUT FROM KEGG
  collate_scores <- jaccard(IXY)  # ensure Matrix library is loaded.
  collate_scores <- diag(collate_scores)
  
  for (i in 1:nrow(sim_matrix)){
    cat("\nDiseaseModule is",  "Module[",i,"] biological value = ",IXY[i])  }
  rankedmods <- IXY
  
  return(rankedmods)
}
  
# getdrugs() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drugs <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls)
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  if(length(ilist) > 0){
    for (j in 1:length(ilist)){
      cat("\ndrug",j,"is", ilist[j])
    }
    return(ilist)
  }else{
    cat("\n","Sorry, no drugs found...check umls code is correct for your disease")
    return(NULL)}
}


# print_dm_table() will generate the latex stuff based on annoations and ranking 
# methods to create the disease module table for the paper, containing:
#   C06/DX0 numbers
#   GO enrichment counts
#   KEGG enrichment counts
#   Biological plausibility ranking
#   Current drugs / DX0 drug reposition candidates
#   components/complexity

print_dm_table <- function(yourtable){

  tablehead <- xtable(head(yourtable))
  tablemiddle <- xtable(middle(yourtable))
  tabletail <- xtable(tail(yourtable))
  dm.table <- rbind(tablehead,tablemiddle,tabletail)
  #digits(tli.table)[c(2,6)] <- 0
  print(dm.table,floating=FALSE)
  
}


# make_C06_mods() is used by join_dm() to build complete list of disease modules (C06 and non-C06)
# The end results is then clustered to see similar/dissimilar modules and more importanatly the 
# degree of module overlap. dplyr::select(dmname_enrich,category,ID,term,genes,DiseaseModule)
# Join C06 structure with appropriate genes and re-annotate then merge.....depth of Cholestasis is C06.130.120
make_C06_mods <- function(){
  start.time <- Sys.time() # start the clock
  C06mods <- data.frame(category="FAKE2",ID="GO:0000666", term="happy behaviour",genes="gene12",DiseaseModule=66,
                        stringsAsFactors=FALSE) #instantiate.
  # break the 55 diseases into the seven broad C06.XXX groups-else we do not have enough genes to form viable communities!
  n <- nrow(C06)
  C06code <- C06$C06_ID
  for (j in 1:n){
    C06code[j] <- substr(C06$C06_ID[j],1,7)}
  C06 <- cbind(C06,C06code)
  
  tempgene <- disgene[,c(1,2,4,6)]  # use a temp variable
  C06code <- vector(mode="character", length=nrow(tempgene))
  tempgene <- cbind(tempgene,C06code) 
  tempgene <- data.frame(lapply(tempgene, as.character), stringsAsFactors=FALSE)
  C06 <- data.frame(lapply(C06,as.character),stringsAsFactors=FALSE)
  for (i in 1:nrow(C06)){
    tempgene[tempgene$diseaseName == C06$C06Disease[i], "C06code"] <- C06$C06code[i]
  }

  C06[] <- lapply(C06, as.character) # keep as strings NOT factors
  Disease <- unique(C06$C06code)
  for (i in 1:length(Disease)){     # for every C06 disease type see what modules they form.
    dg <- dplyr::filter(tempgene,C06code==Disease[i])
    cat("\ni=",i,"  ",Disease[i]," with ", length(unique(dg$geneName)), " genes.")
    usethese <- unique(dg$geneName)
    if(length(usethese) > 20){   # use no more than 100 implicated genes
      usethese <- usethese[sample(1:length(usethese), 20,replace=FALSE)] }
    cat("\nUsing ",length(usethese)," genes.")
    if(i == 7){usethese <- c("C5","GSK3B","IFNG","IL18","HMGB1","NLRP1","NLRP3")} # provide Peritonitis with genes not in MY database
    usethese[] <- lapply(usethese,as.character)
    usethese <- str_to_upper(usethese); usethese <- unique(usethese)
    tempinteractions <- use_rentrez(usethese)
    tempinteractions[] <- lapply(tempinteractions, as.character) # keep as strings NOT factors
    tempinteractions[,1] <- str_to_upper(tempinteractions[,1])
 
    
    lcom <- getLinkCommunities(tempinteractions, hcmethod = "single")  # consider cutting density partition manually
    lmods <- createDiseaseModules(lcom)  
    if(length(lcom$clusters) < 10){  # if fewer than 10 modules try again - cut 
      lcom <- newLinkCommsAt(lcom, cutat = 0.5) # cut it at 0.5 or 0.6
      lmods <- createDiseaseModules(lcom) 
    }
    lmods <- dplyr::select(lmods,category,ID,term,genes,DiseaseModule)
    if(nrow(lmods) > 4000){
      lmods <- sample_n(lmods,4000)}
    lmods$DiseaseModule <- paste(C06$C06code[i],lmods$DiseaseModule,sep="_")
    C06mods <- rbind(C06mods,lmods)
  }
  C06mods <- C06mods[-1, ]     # 1st entry is rubbish so remove it
  
  end.time <- Sys.time() # stop clock and figure out time taken
  time.taken <- end.time - start.time
  cat("\n Execution took - ",time.taken)
  return(C06mods)
}


# join_dm() will concatenate the C06 diseases with the non C06 diseases
# afterwards, clustering can be used on their GO to determine similarities/differences.
# They are: alzmods; asthmods; autmods; diamods; hypmods; nscmods; obsmods; parkmods; ramods; schmods; 
# At the moment the above are Global data structures
join_dm <- function(){
  # WARNING DO ONLY ONCE!!!!
  alzmods_enrich$DiseaseModule <- paste("Alzheimers",alzmods_enrich$DiseaseModule,sep="_")
  autmods_enrich$DiseaseModule <- paste("Autism",autmods_enrich$DiseaseModule,sep="_")
  asthmods_enrich$DiseaseModule <- paste("Asthma",asthmods_enrich$DiseaseModule,sep="_")
  diamods_enrich$DiseaseModule <- paste("Diabetes",diamods_enrich$DiseaseModule,sep="_")
  hypmods_enrich$DiseaseModule <- paste("Hypertension",hypmods_enrich$DiseaseModule,sep="_")
  nscmods_enrich$DiseaseModule <- paste("non-small-cell-carcinoma",nscmods_enrich$DiseaseModule,sep="_")
  obsmods_enrich$DiseaseModule <- paste("Obesity",obsmods_enrich$DiseaseModule,sep="_")
  parkmods_enrich$DiseaseModule <- paste("Parkinson",parkmods_enrich$DiseaseModule,sep="_")
  ramods_enrich$DiseaseModule <- paste("RheumatoidArthritis",ramods_enrich$DiseaseModule,sep="_")
  schmods_enrich$DiseaseModule <- paste("Schizophrenia",schmods_enrich$DiseaseModule,sep="_")
  
  allmods <- rbind(dplyr::select(diamods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(autmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(asthmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(hypmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(nscmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(obsmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(parkmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(ramods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(schmods_enrich,category,ID,term,genes,DiseaseModule),
                   dplyr::select(alzmods_enrich,category,ID,term,genes,DiseaseModule))
  
  C06mods <- make_C06_mods()
  C06mods$genes
  allmods <- rbind(allmods,C06mods)
  return(allmods)
}

#save(C06mods,allmods, file = "C06-20thDec-2017.RData")


# merge_dm() will merge modules based on GO biological similarity. Where "dm" is
# a dataframe of modules (small_data). Then create an n x m matrix of dismods and genes, 
# then cluster it.
merge_dm <- function(dm){
  umods <- unique(dm$DiseaseModule)
  nmods <- length(umods)
  elist <- data.frame(genes="RU12",disease="Pox_1",stringsAsFactors = FALSE)
  
  for (i in 1:nmods){  # create an edgelist
    cat("\nCreating DM ",i, " of ",nmods)
    tempstuff <- filter(dm,DiseaseModule == umods[i])  # work way thru all disease modules
    tempstuff <- tempstuff %>% distinct(genes,DiseaseModule, .keep_all = TRUE)
    elist_temp <- cbind(genes=tempstuff$genes,disease=tempstuff$DiseaseModule)
    elist <- rbind(elist,elist_temp)
  }
  
  elist <- elist[-1,]   # remove silly 1st entry
  mat <- matrix(0, length(unique(elist[,1])), length(unique(elist[,2])))
  rownames(mat) <- unique(elist[,1])
  colnames(mat) <- unique(elist[,2])

  len_el <- nrow(elist)
  for(i in 1:len_el){ 
    matcols <- as.character(elist[i,]$genes)
    matrows <- as.character(elist[i,]$disease)
    cat("\n - ",i," ", matcols," ", matrows)
    mat[matcols,matrows] <- 1 
    
  }

  rownames(mat) <- unique(elist[,1])
  colnames(mat) <- unique(elist[,2])
  
  #d = dist(mat, method = "binary")
  #hc = hclust(d, method="ward")
  #plot(hc)
  #cluster.means = aggregate(mat,by=list(cutree(hc, k = 6)), mean)
  
  testcluster <- biclust(x = mat, method=BCBimax(),number=50,minr=5,minc=3) # needs biclust library
  testcluster
  #testcluster <- biclust(x = disma, method=BCBimax()) # needs biclust library
  
  #plotclust(testcluster,mat,noC=12)
  #drawHeatmap(x = mat, bicResult = testcluster, number = 2)
  #drawHeatmap2(x = mat, bicResult = testcluster, number = 1) 
  #drawHeatmap2(x = mat, bicResult = testcluster,number=10)
  
  testcluster <- biclust(mat,method=BCSpectral(),numberOfEigenvalues=1,withinVar=100)
  testcluster
  drawHeatmap2(x = mat, bicResult = testcluster,number=4)
  
  writeBiclusterResults("BI-results-21.txt", testcluster,"Disease modules", dimnames(mat)[1][[1]],
                        dimnames(mat)[2][[1]])
  
  #xmotif<-biclust(x=mat, method=BCXmotifs(), number=50, alpha=0.05,nd=20, ns=20, sd=5)
  #drawHeatmap2(x = mat, bicResult = xmotif,number=1)
  
  return(mat)
  
}


# numbers_only() does a string contain numbers only? I used this to debug the NCBI gene server, as 
# it appears to be sending me ENTREZID id's instead of SYMBOL names.
numbers_only <- function(x) !grepl("\\D", x)

# obtained from stackoverflow question.
jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}

# NOT   working as fucntion - but hand pasting code into console works
# How much RAM is taken up by current R session? Obtained from StackOverFlow
# https://stackoverflow.com/questions/1395270/determining-memory-usage-of-objects
ram_used <- function(){
   #sort(sapply(ls(), function(x) format(object.size(get(x)), unit = 'auto')))
   object.size(x=lapply(ls(), get))
   print(object.size(x=lapply(ls(), get)), units="Mb")
   
   Mb <- ls() %>% sapply(. %>% get %>% object.size %>% '/'(10^6))
   cbind(Mb, "Mb") %>% as.data.frame
   sum(Mb)
   }


# Another stackoverflow soloution by Yan
# https://stats.stackexchange.com/questions/95523/significance-of-overlap-between-multiple-lists
# gene.list <- list(listA=paste0("gene",c(1,2,3,4,5,6,7,8,9)),
# where the upper triangle on the right is the lengths of overlap of each pair, and the bottom triangle 
# on the left is the significance of the overlap by hypergeometric test. Here in this toy example, the 
# overlap between listA and listB is significant in a world a 14 genes if you choose 0.05 as your p-value 
# cutoff. Any other pair is not significantly overlapping.


hyper_matrix <- function(gene.list, background){
  # generate every combinations of two gene lists
  combination <- expand.grid(names(gene.list),names(gene.list))
  combination$values <- rep(NA, times=nrow(combination))
  
  # convert long table into wide
  combination <- reshape(combination, idvar="Var1", timevar="Var2", direction="wide")
  rownames(combination) <- combination$Var1
  combination <- combination[,-1]
  colnames(combination) <- gsub("values.", "", colnames(combination))
  
  # calculate the length of overlap of each pair
  for(i in colnames(combination)){
    for(j in rownames(combination)){
      combination[j,i]<-length(intersect(gene.list[[j]],gene.list[[i]]))
    }
  }
  
  # calculate the significance of the overlap of each pair
  for(m in 1:length(gene.list)){
    for(n in 1:length(gene.list)){
      if(n>m){
        combination[n,m] <- phyper(combination[m,n]-1, length(gene.list[[m]]), background-length(gene.list[[m]]), length(gene.list[[n]]), lower.tail=F)
        # note that the phyper function (lower.tail=F) give the probability of P[X>x], so the the overlap length should subtract 1 to get a P[X>=x].
      }
    }
  }
  # round to 2 digit.
  return(round(combination,2))
}

optimum_clusters <- function(df){
  # Elbow method
  fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method")
  
  # Silhouette method
  fviz_nbclust(df, kmeans, method = "silhouette") +
    labs(subtitle = "Silhouette method")

  # Gap statistic
  # nboot = 50 to keep the function speedy. recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  set.seed(123)
  fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50) +
   labs(subtitle = "Gap statistic method")

}

