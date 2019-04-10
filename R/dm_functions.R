# dm_functions.R
library(dplyr) 
#library(ChemmineR)
#library(ChemmineOB)
#library(ape)
#library(sparcl)
#library(cluster) # used for kmeans and silhoutte plot
#library(xtable)
library(gplots) 
#library(scatterplot3d) 
library(igraph)
library(diffusr)
#library(ROCR)
#library(VennDiagram)
library(ggplot2)
#library(linkcomm)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#library("AnnotationDbi")
#library(GOplot)
#library(scales)
library(NCBI2R)
library(rentrez)
#library(stringr)
# huge amount of data gets loaded in here onwards!
#library(ontologySimilarity)
#library(ontologyIndex)
#library(infotheo)
#library(KEGGprofile)
#library(KEGG.db)
#library(Matrix)
#data(go)
#data(gene_GO_terms)
#data(GO_IC)

# source("https://bioconductor.org/biocLite.R")
# biocLite("ReactomePA")

## --------------------- FUNCTION DEFINITIONS -----------------------

# Makes first letter of string uppercase
uppercase <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# R already provides a tail and head command to view last six and first six elements, so 
# why not the middle six?, also supply x if you want more/less than the usual 5 records.
middle <- function(mydata,x=NULL) {
  len <- nrow(mydata)
  startpoint <- round(len/2)
  endpoint <- startpoint+5
  if(!is.null(x)){
    endpoint <- startpoint+(x-1)
  }
  mydata[startpoint:endpoint,]
  
}


build_network <- function(ppi_data){
  # assign 1=disease; 0=non-disease to each protein
  un_disease <- (unique(disgene$geneSymbol))        # XXX unique disease genes
  length(un_disease)
  un_ppi <- (unique(c(ppi_data$Gene_A,ppi_data$Gene_B)))      # XXX unique general proteins in ppi
  length(un_ppi)
  
  joint_ppi <- un_disease[un_disease %in% un_ppi]  # 
  not_ppi <- un_disease[!un_disease %in% un_ppi]  # 
  
  # dataframe containing disease proteins and non-disease proteins. Annotate with:
  # 1. is it disease implicated; 2. is it a hub? 3. is it a drug target?
  # create ppi network (igraph object) and annotate with target or not target, disease or not disease
  ppi_net <- graph.data.frame(ppi_data,directed = FALSE)
  #ppi_net <- graph_from_data_frame(ppi_data,directed = FALSE,vertices = NULL)
  
  #ppi_net <- as.undirected(ppi_net); 
  ppi_net <- igraph::simplify(ppi_net)  # remove duplicates and self-loops
  ppi_net <- delete_isolates(ppi_net)
  delete.vertices(igraph::simplify(ppi_net), degree(ppi_net)==0)  # remove proteins with no partners
  
  V(ppi_net)[1:vcount(ppi_net)]$disease <- 0   # Intialise all to zeros
  V(ppi_net)[1:vcount(ppi_net)]$hub <- 0   # Intialise all to zeros
  V(ppi_net)[1:vcount(ppi_net)]$ptype <- "unknown"   # Intialise protein "type" to unknown
  V(ppi_net)[1:vcount(ppi_net)]$dislist <-list(c("unknown"))   # Intialise list of diseases associated with this protein
  V(ppi_net)[1:vcount(ppi_net)]$drugs   <-list(c("unknown"))   # Intialise list of drugs associated with this protein
  
  
  # get main component only - ignore lessor weakly connected groups
  V(ppi_net)$comp <- components(ppi_net)$membership
  ppi_net <- induced_subgraph(ppi_net,V(ppi_net)$comp==1)
  
  # remove from joint_ppi the lost nodes 
  survivors <- V(ppi_net)$name
  joint_ppi <- un_disease[un_disease %in% survivors] 
  
  ppi_net <- set_vertex_attr(ppi_net,"disease",joint_ppi,1) # Now assign "1" if protein is a disease protein (very neat!)
  
  V(ppi_net)$label <- V(ppi_net)$name  # names are not plotted in igraph, but labels are
  
  return(ppi_net)
}

# remove anything not connected to the giant central component (GCC)
delete_isolates <- function(gt) {
  isol <- V(gt)[degree(gt)==0]
  gt <- delete.vertices(gt, isol)
  return(gt)
  
}


# Calculate some statistics about the disease gene network
# returns a list: net and nodes
get_gstatistics <- function(gt) {
  net <- data.frame( 
    modu=igraph::modularity(gt, membership(cluster_walktrap(gt))),
    avepath=igraph::average.path.length(gt),
    nedges=igraph::ecount(gt),
    nverts=igraph::vcount(gt),
    transit=igraph::transitivity(gt),
    diam=igraph::diameter(gt,weights=NA),
    connect=igraph::is.connected(gt))
  
  nodes <- data.frame(   
    closeness=igraph::estimate_closeness(gt,mode="all",cutoff=3),
    degree=(igraph::degree(gt)),
    betweenness=igraph::estimate_betweenness(gt,directed=FALSE,cutoff=3),
    hubness=igraph::hub_score(gt)$vector,
    central=vector(mode="integer", length=net$nverts),
    comm=vector(mode="integer", length=net$nverts))
  
  tmp <- igraph::cluster_walktrap(gt)
  nodes$comm <- as.vector(membership(tmp))
  alpha <- igraph::alpha_centrality(ppi_net,alpha=0.1)  
  nodes$central <- as.vector(alpha)
  
  cat("\nOverall network statistics:")
  cat("\n   Modularity ",net$modu)
  cat("\n   Average path ",net$avepath)
  cat("\n   N edges ",net$nedges)
  cat("\n   N vertices ",net$nverts)
  cat("\n   Transitivity ",net$transit)
  cat("\n   Diameter ",net$diam)
  cat("\n   Is connected? ",net$connect)
  gstats = list(net=net, nodes=nodes)
  return(gstats)
}


# Drug targets of the various drugs - usefully contains protein type (e.g. GPCR) as well.
# from http://drugcentral.org/download
# DrugCentral is a comprehensive drug information resource for FDA drugs and drugs approved outside USA. The 
# resources can be searched using: drug, target, disease, pharmacologic action, terms. 
load_drugtargets <- function(){
  drug_targets <- read.csv(file="C://R-files//disease//drug.target.interaction.tsv", header=TRUE, sep="\t",stringsAsFactors = FALSE)
  names(drug_targets)[names(drug_targets)=="DRUG_NAME"] <- "DrugName"
  names(drug_targets)[names(drug_targets)=="TARGET_CLASS"] <- "TargetClass"
  names(drug_targets)[names(drug_targets)=="GENE"] <- "Gene"
  drug_targets <- drug_targets[,c(1,4,6)]  # Keep only need three variables
  drug_targets$DrugName <- uppercase(drug_targets$DrugName)   # convert first letter to uppercase to match existing data
  drug_targets <- na.omit(drug_targets)# remove NA's
  
  # now unlist special entries, I edited the original file and replaced "|" with "/"
  drug_targets<-
    drug_targets %>% 
    mutate(Gene=strsplit(as.character(Gene), "/")) %>%   # symbols=Gene
    unnest(Gene)
  drug_targets$Gene <- toupper(drug_targets$Gene)  # all to uppercase
  
  # shorten some names, for ease printing etc
  drug_targets$TargetClass <- gsub('Nuclear hormone receptor', 'NR', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Transcription factor', 'TF', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Membrane receptor', 'Membrane', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Ion channel', 'IC', drug_targets$TargetClass)
  return(drug_targets)
}



# Every ppi_net protein needs annotating with implicated diseases (if any)
# Each disease implicated protein may have more than one disease associated with it.
# Majority of proteins will have none or at least "unknown"
annotate_ppi_diseases <- function(ppi_net){
  cat("\nAnnotating ppi with known diseases.....")
  glist <- V(ppi_net)$name

  for (i in 1:length(glist)){
    templist <- filter(disgene,geneSymbol == glist[i])
    if(nrow(templist)!=0){
      #cat("\nFound a disease for ",glist[i])
      V(ppi_net)$dislist[i] <- list(templist$diseaseName)
    }
  }
  return(ppi_net)
}


# Every ppi_net protein needs annotating with the drugs targetted at it (if any)
# Each target protein may have more than one drug aimed at it.
# Majority of proteins will have none 
annotate_ppi_drugs <- function(ppi_net){
  cat("\nAnnotating ppi with known drugs.....")
  glist <- V(ppi_net)$name
 
  for (i in 1:length(glist)){
    templist <- filter(drug_targets,Gene == glist[i])
    if(nrow(templist)!=0){
      #cat("\nFound a disease for ",glist[i])
      V(ppi_net)$drugs[i] <- list(templist$DrugName)
    }
  }
  return(ppi_net)
}

# For each protein in the list findout how many publications it has been mentioned in
# break down into blocks of BLOCKSIZE. Otherwise server error occurs.
get_pubs <- function(protein_list){
  total.time <- Sys.time()
  cat("\nObtaining a publication count for each protein in your list - may take some time.....")
  BLOCKSIZE <- 250
  npubs <- NULL
  len <- length(protein_list)
  if(len > BLOCKSIZE){
    nblocks <- len %/% BLOCKSIZE  # how many blocks?
    remainder <- len %% BLOCKSIZE  # remainder, to be mopped up with one final block
    bstart <- 1; bfinish <- BLOCKSIZE
    if(remainder > 0){nblocks <- nblocks+1}  # some left over, so another final block needed
    for(i in 1:nblocks){
      start.time <- Sys.time()
      if(i==1){bstart<-1;bfinish <- BLOCKSIZE} # 1st time around 
      if(i>1){bstart <- bstart+BLOCKSIZE; bfinish <- bfinish + BLOCKSIZE}
      if(i== nblocks){bfinish <- len}
      #cat("\nbstart=",bstart," bfinish=",bfinish)  #debug
      ptemp <- count_articles(protein_list[bstart:bfinish],BLOCKSIZE)
      npubs <- rbind(npubs,ptemp)
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      cat(" took ",time.taken)
      if(bstart==bfinish){end.time <- Sys.time(); time.taken <- end.time - start.time;return(npubs)}
    }
  }
  
  if(len < (BLOCKSIZE-1)){npubs <- count_articles(protein_list,BLOCKSIZE)} #if num proteins < blocksiez then only one call needed.
  end.time <- Sys.time()
  time.taken <- end.time - total.time
  cat(": total time took ",time.taken)
  return(npubs)
}

# See how many research articles are written about our proteins. Uses rentrez package.
count_articles <- function (protein_block,bsize){
  cat("\n     Obtaining block of up to ",bsize," proteins.....")
  for (i in 1:length(protein_block)){
    pname <- paste(protein_block[i],'[GENE]) AND (Homo sapiens[ORGN])',sep="")
    ids <- entrez_search(db="pubmed", term=pname,retmax=10000)
    
    atemp <- cbind(protein_block[i],length(ids$ids))
    
    if(i!=1){
      articles <- rbind(articles,atemp)} 
    else{
      articles <- atemp}
  }
  return(articles)
}






