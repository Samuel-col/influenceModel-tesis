
# Header ----
library(dplyr)
library(ggplot2)
library(igraph)

rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")
dir("datasets/hoaxy")

# Read Data ----
df_raw <- read.csv("datasets/hoaxy/reforma_tributaria_04_11_2022_5pm_mixed_spanish.csv")

df_raw %>% tibble %>% filter(tweet_type == "retweet") %>% 
  select(from_user_id,to_user_id,
         from_user_screen_name,to_user_screen_name,
         from_user_botscore,to_user_botscore) -> df

G <- graph_from_data_frame(select(df,from_user_id,to_user_id),directed = T)
G <- simplify(G)

save(G,df,file = "out_refTrib.RData")

class(G)
ecount(G)
vcount(G)

## Plots ----
igraph_options(vertex.label = "",edge.arrow.size = 0.5)
{
x11()
par(mar = c(0,0,0,0))
plot(G,vertex.size = log(1+degree(G)),layout = layout_with_fr)
}

{
x11()
dd <- degree_distribution(G)
plot(0:(length(dd)-1),dd, log = "xy",type = "b", col = "blueviolet", lwd = 2)
grid()
}

# Giant component ----
comp <- components(G)
comp$csize[1]/vcount(G)
G_gc <- subgraph(G,which(comp$membership==1))
df_gc <- df_raw[comp$membership==1,]

## Plot ----

N <- vcount(G_gc)
n_significativos <- 9
  
ifelse(rank(degree(G_gc,mode = "out"),
                   ties.method = "last")>N-n_significativos,
       N - rank(degree(G_gc,mode = "out"),
                ties.method = "last") + 1,NA) -> plt_vertex.label

{
pdf("data_influence_network_base_plot.pdf",height = 6,width = 6)
par(mar = c(0,0,0,0))
set.seed(2023)
plot(G_gc,vertex.size = log(1+degree(G_gc)),
     layout = layout_with_fr,
     vertex.label = plt_vertex.label,
     vertex.label.color = "red",
     vertex.color = "lightblue",
     vertex.label.cex = 0.4,
     edge.width = 0.5,
     edge.color = "#a0a0a060",
     edge.arrow.size = 0.2)
dev.off()
}

# tabla de más retuiteados
users <- df_gc %>% select(id = from_user_id,
                          name = from_user_screen_name) %>% 
  unique()
users <- rbind(users,
               df_gc %>% 
                 select(id = to_user_id,
                        name = to_user_screen_name) %>% 
                 unique()
               ) %>% unique()

rownames(users) <- NULL

d_out <- degree(G_gc,mode = "out")
tibble(id = as.numeric(names(d_out)),d_out = d_out) %>% 
  left_join(users,by = join_by(id == id)) %>% 
  arrange(desc(d_out)) -> df_by_out_degree

# df_by_out_degree %>% select(-id) %>% 
#   mutate(d_out = as.integer(d_out)) %>% 
#   head(10) -> tmp_output

xtable::xtable(tmp_output)

A <- as_adjacency_matrix(G_gc)
A <- as.matrix(A)
N <- vcount(G_gc)
K <- 2


# ¿Quién conecta los dos clusters?
which(degree(G_gc) == 12)

users %>% filter(
  id == 25490757 
  ) %>% select(id) -> my_id

plt_vertex.label <- ifelse(as.numeric(names(degree(G_gc))) == my_id[[1]][1],
                           "*",NA)

par(mar = c(0,0,0,0))
set.seed(2023)
plot(G_gc,vertex.size = log(1+degree(G_gc)),
     layout = layout_with_fr,
     vertex.label = plt_vertex.label,
     vertex.label.color = "red",
     vertex.color = "lightblue",
     vertex.label.cex = 4,
     edge.arrow.size = 0.3)

users %>% filter(
  id == 25490757 
) %>% select(name)


