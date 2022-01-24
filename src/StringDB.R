#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")
install.packages("igraph")

#Load
library(STRINGdb)
library(igraph) #to manipulate and plot data


string_db <- STRINGdb$new( version="11", species=9606, 
                           score_threshold=200, input_directory="") # 9606 for Human


full.graph <- string_db$get_graph()

# see how many proteins do you have    
vcount(full.graph)

# make a data frame 
df_edges <- as.data.frame(get.edgelist(full.graph))


# find top 200 proteins with the highest degree
top.degree.verticies <- names(tail(sort(degree(full.graph)), 200))


# extract the relevant subgraph
top.subgraph <- induced_subgraph(full.graph, top.degree.verticies)

plot(top.subgraph)

# count the number of proteins in it
vcount(top.subgraph)
