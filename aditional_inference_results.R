library(paletteer) 
my.pal <- paletteer_d("jcolors::pal2")
library(igraph)
library(dplyr)
library(fields)

aw <- 1; bw <- 1; as <- 1; bs <- 1 # caso 1
path_code <- paste0("out_aw_",aw,"_bw_",bw,"_as_",as,"_bs_",bs,"__")

setwd("/home/samuel/Documentos/U/tesis/")

df_raw <- read.csv("datasets/hoaxy/reforma_tributaria_04_11_2022_5pm_mixed_spanish.csv")

df_raw %>% tibble %>% filter(tweet_type == "retweet") %>% 
  select(from_user_id,to_user_id,
         from_user_screen_name,to_user_screen_name,
         from_user_botscore,to_user_botscore) -> df

G <- graph_from_data_frame(select(df,from_user_id,to_user_id),directed = T)
G <- simplify(G)

comp <- components(G)
comp$csize[1]/vcount(G)
G_gc <- subgraph(G,which(comp$membership==1))
df_gc <- df_raw[comp$membership==1,]

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

# load("out_refTrib.RData")
load(paste0(path_code,"model_posterior_processed_samples.RData"))

d_out <- degree(G_gc,mode = "out")
tibble(id = as.numeric(names(d_out)),d_out = d_out,
       d_in = degree(G_gc,mode = "in"),
       I_hat = I_hat,
       O_hat = O_hat) %>% 
  left_join(users,by = join_by(id == id)) -> users_covs 

# Tabla - O ----
users_covs %>% arrange(desc(O_hat)) %>% 
  select(O_hat,d_out,name) %>% head(13) %>% 
  xtable::xtable()

# Color - O ----

my_palette <- colorRampPalette(my.pal[5:4])
cols_grid <- my_palette(100)
cols <- cols_grid[round((O_hat-min(O_hat))/(max(O_hat)-min(O_hat))*100)]

{
pdf("data_influence_network_O_color.pdf",height = 6,width = 6)
par(mar = c(4,0,0,0))
set.seed(2023)
plot(G_gc,vertex.size = log(1+degree(G_gc)),
     layout = layout_with_fr,
     vertex.label = NA,
     vertex.color = cols,
     vertex.label.cex = 0.4,
     edge.width = 0.5,
     edge.color = "#a0a0a060",
     edge.arrow.size = 0.2)
image.plot(legend.only = T,zlim=range(O_hat),
           col = cols_grid,add = TRUE,
           legend.shrink = 0.8,
           horizontal = T,
           legend.width = 0.5)
dev.off()
}

# Tabla - I ----

user_covs %>% arrange(desc(I_hat)) %>% 
  select(I_hat,d_in,name) %>% head(10) %>% 
  xtable::xtable()

# Color - I ----


my_palette <- colorRampPalette(my.pal[5:4])
cols_grid <- my_palette(100)
cols <- cols_grid[round((I_hat-min(I_hat))/(max(I_hat)-min(I_hat))*100)]
shapes <- c("square","circle")[1 + (O_hat < (-15))]

{
pdf("data_influence_network_I_color.pdf",height = 6,width = 6)
par(mar = c(0,0,0,0))
set.seed(2023)
plot(G_gc,vertex.size = log(1+degree(G_gc)),
     layout = layout_with_fr,
     vertex.label = NA,
     vertex.color = cols,
     vertex.shape = shapes,
     vertex.label.cex = 0.4,
     edge.arrow.size = 0.3)
dev.off()
}


