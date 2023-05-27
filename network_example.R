# Network example

setwd("/home/samuel/Documentos/U/tesis/")
rm(list=ls())

library(igraph)
library(paletteer) 
my_palette <- colorRampPalette(c("lightblue","lightgreen"))
cols_grid <- my_palette(100)
library(fields)

minmax_scaler <- function(x){
  return((x - min(x))/(max(x)-min(x)))
}


# help(package = "igraph")

# data from: https://github.com/bansallab/asnr/blob/master/Networks/Mammalia/raccoon_proximity_weighted/weighted_raccoon_matrix_3.graphml
igraph::read_graph("datasets/asnr/weighted_raccoon_matrix_3.graphml",
                   format = "graphml") -> G

# Articulo original: 
# https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2656.12422


# Plot options
igraph_options(vertex.size = 10, vertex.color = "lightblue",
               vertex.label.color = "black")


{
  pdf("raccoons_network.pdf",
      height = 6, width = 6)
  set.seed(1)
  par(mar = c(0,0,0,0),mfrow = c(1,1))
  plot(G,
       layout = layout_with_dh,
       vertex.label.color = "red")
  dev.off()
}


# Subgrafo
G_sub <- subgraph(G,c(19,8,2,18,1,21))
{
  pdf("raccoons_subnetwork.pdf",
      height = 6, width = 6)
  set.seed(1)
  par(mar = c(0,0,0,0),mfrow = c(1,1))
  plot(G_sub,layout = layout_with_dh,
       vertex.label = c(19,8,2,18,1,21),
       vertex.label.color = "red"
  )
  dev.off()
}

# Adj Matrix
A <- as_adj(G,attr = "weight")
A

xtable::xtable(as.matrix(A),digits = 0)

# Orden
vcount(G)

# Tamaño
ecount(G)

# Densidad
edge_density(G)

# k-estrella
## Gsub es 5-estrella con centro en 21

# Transitividad
transitivity(G)

# Distancia media
mean(distances(G))

# Asortatividad
assortativity_degree(G)

# Grado
degs <- degree(G)

# Vecinos
## Los vecinos del 11 son 3, 6 y 10

# Centralidad propia
eig_cent <- eigen_centrality(G)$vector

# Centralidad de intermediación
bet_cent <- betweenness(G)

{
  pdf("raccoons_centrality.pdf",
      height = 5,width = 8)
  set.seed(1)
  par(mar = c(4,0,1,0),mfrow = c(1,2))
  plot(G,layout = layout_with_dh,
       vertex.label.cex = 0.7,
       vertex.color = cols_grid[1+ceiling(99*minmax_scaler(eig_cent))],
       main = "Centralidad propia"
  )
  image.plot(legend.only = T,zlim=range(eig_cent),
             col = cols_grid,add = TRUE,
             legend.shrink = 0.8,
             horizontal = T)
  set.seed(1)
  plot(G,layout = layout_with_dh,
       vertex.label.cex = 0.7,
       vertex.color = cols_grid[1+ceiling(99*minmax_scaler(bet_cent))],
       main = "Centralidad de intermediación"
  )
  image.plot(legend.only = T,zlim=range(bet_cent),
             col = cols_grid,add = TRUE,
             legend.shrink = 0.8,
             horizontal = T)
  
  dev.off()
}

node_stats <- cbind(degs,eig_cent,bet_cent)
rownames(node_stats) <- 1:vcount(G)
colnames(node_stats) <- c("Grado","C. Propia","C. Intermediación")
xtable::xtable(node_stats,digits = 3)

# CFG
cols <- c("lightblue","lightgreen","lightpink")
cl <- cluster_fast_greedy(G)
modularity(cl)
clusters_CFG <- cut_at(cl,3)
cl <- cluster_walktrap(G)
modularity(cl)
clusters_WT <- cut_at(cl,3)
{
  pdf("raccoons_cluster.pdf",
      height = 5,width = 8)
  set.seed(1)
  par(mar = c(0,0,2,0),mfrow = c(1,2))
  plot(G,layout = layout_with_dh,
       vertex.label.cex = 0.7,
       vertex.color = cols[clusters_CFG],
       main = "Agrupamiento rápido y codicioso"
  )
  set.seed(1)
  plot(G,layout = layout_with_dh,
       vertex.label.cex = 0.7,
       vertex.color = cols[clusters_WT],
       main = "Agrupamiento Walktrap"
  )
  dev.off()
}

