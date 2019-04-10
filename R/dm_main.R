# dm_main.R
# Search for disease modules from igraph object based on known diseased implicated proteins
# A disease module may exist across functional and topological proteins structures.
# Ken McGarry updated: 23/6/18

setwd("C:/R-files/DiseaseModules")    # point to where my code lives
source("dm_functions.R")
source("dm_loaddata.R")
load("C:/R-files/DiseaseModules/dm_19thJune.RData")  # reload precomputed data to save time

# ---------------- Stages:  --------------------
# 1. Create a list of the 26 disease categories using mesh datafile
# 2. Build ppi network using pairwise protein connections from the several databases
# 3. Get genes associated with each disease from "disgene" dataframe.
# 4. Get the drugs associated with each disease.
# 5. Compute statistics relating to diseases and drugs
# 6. Determine topological modules, functional modules and disease modules (hard bit!)

### Stage 1: dm_load_data.R performs stage 1, simply list the mesh categories
meshcats

### Stage 2: Build ppi_net complex network and obtain statistics, we end up with 15,554 proteins 
# with 109,746 edges. build_network.R annotates each gene with its: disease gene status; if targeted by drug;  
ppi_net <- build_network(ppi_hint)
#net_stats <- get_gstatistics(ppi_net)  # statistics takes 10-15 minutes to compute on my laptop.


### Stage 3: every ppi_net protein need annotating with drugs targetted at it (if any)
ppi_net <- annotate_ppi_drugs(ppi_net)


### Stage 4: every ppi_net protein need annotating with implicated diseases (if any)
ppi_net <- annotate_ppi_diseases(ppi_net)


### Stage 5: Compute statistics relating to diseases and drugs
# promiscuity of drugs to proteins and vice-versa, figure out level of bias: publication count to protein interacts
plist <- V(ppi_net)$name

# manual access - Takes about 10 minutes for every 1,000 protein names- Careful :fuckup from 10000 to 10500
##count_temp <- get_pubs(plist[15001:15554])   # Ive had to manually change numbers in square brackets, keeping eye on errors
##count_pubs <- rbind(count_pubs,count_temp)  # concatenate the returned pubs

# sort out data frame for degree,
tempdata <- net_stats[[2]]
rnames <- row.names(tempdata)
tempdata <- cbind(rnames,tempdata$degree) # Use degree as a count of interactions
tempdata <- as.data.frame(tempdata,stringsAsFactors=FALSE)
tempdata$V2 <- as.numeric(tempdata$V2)
colnames(tempdata)<- c("protein","degree")

# sort out dataframe for publications,
count_pubs <- as.data.frame(count_pubs,stringsAsFactors=FALSE)
count_pubs$V2 <- as.numeric(count_pubs$V2)
colnames(count_pubs)<- c("protein","publications")

# join count_pubs with tempdata with "protein" as key. Data ready for plotting
plot_data1 <- count_pubs %>%
  left_join(tempdata, by="protein")

plot_data1 <- filter(plot_data1,degree < 250)  # get rid of five outliers with degree more than 300 (prettier plot)
plot_data1 <- filter(plot_data1,publications != 0)  # get rid of proteins with zero publications (more sensible plot)

# do a nice ggplot scatterplot of publications Versus degree
p <- ggplot(plot_data1, aes(publications,degree)) +
  geom_point(alpha = 1/6,colour = "blue", size = 2) +
  scale_x_continuous(breaks=seq(0, 10000, 1000))  +
  scale_y_continuous(breaks=seq(0, 250, 50))
  
p <- p + theme(axis.title.y = element_text(size = rel(1.8), angle = 90))
p <- p + theme(axis.title.x = element_text(size = rel(1.8), angle = 00))
p + theme(axis.text.x = element_text(face="bold", size=14, angle=0),
          axis.text.y = element_text(face="bold", size=14, angle=0))
p

# little point in plotting No drugs versus No proteins
#count_drugs <- get_drugs(plist,drug_targets)




#### bits of code examples of accessing attributes from igraph objects ####
#vertex_attr_names(ppi_net)
#get.vertex.attribute(ppi_net,"dislist")
#table(get.vertex.attribute(ppi_net,"disease"))
# get node index by name
# V(ppi_net)$name  # gets all the gene names
# V(ppi_net)$name[1]   # get a specific protein name by index 
#nindex <- which(V(ppi_net)$name == "TP53")   # search for a particular protein index by name
#get.vertex.attribute(ppi_net, "dislist", index=nindex)  # obtain an attribute by the index
# Each disease implicated protein may have more than one disease associated with it.
# Majority of proteins will have none or at least "unknown"
#V(ppi_net)$dislist[2602] <- list(c("headache","upset tummy"))   # set attribute individually, just an example
# Each target protein may have more than one drug aimed at it.
# Majority of proteins will have none or at least "none"
#V(ppi_net)$drugs[2602] <- list(c("cocaine","marijuana"))   # set attribute individually, just an example
#query <- get.vertex.attribute(ppi_net, index = nindex)  # if attribute name is missing like here then all attributes are returned


