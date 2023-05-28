# JAGS implementation

# Header ----
rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")

library(R2jags)
library(igraph)
library(paletteer) 
my.pal <- paletteer_d("jcolors::pal2")

# Pre-model ----

load("out_refTrib.RData")
comp <- components(G)
comp$csize[1]/vcount(G)
G <- subgraph(G,which(comp$membership==1))

model <- function(){
  # likelihood
  for(i in 1:N){
    for(j in 1:N){
      y[i,j]        ~ dbern(ilogit(linPred[i,j]))
      linPred[i,j] <- O[i] + inprod(u[i,],u[j,])/sqrt(inprod(u[i,],u[i,]))
    }
  }
  
  # prior
  O[i] ~ dnorm(0,iw)
  for(i in 1:N){
    u[i,1:p] ~ dmnorm(u0,is*I)
  }
  iw ~ dgamma(aw,bw)
  is ~ dgamma(as,bs)
}

# prior parameters
p <- 2
N <- vcount(G)
u0 <- c(rep(0,p))
I <- diag(p)
aw <- 1; bw <- 1; as <- 1; bs <- 1 # caso 1
# aw <- 2; bw <- 1; as <- 2; bs <- 1 # caso 2
# aw <- 3; bw <- 2; as <- 3; bs <- 2 # caso 3
path_code <- paste0("out_aw_",aw,"_bw_",bw,"_as_",as,"_bs_",bs,"__")

# input
as.matrix(as_adjacency_matrix(G)) -> y
y <- 1*(y>0)
rownames(y) <- colnames(y) <- NULL
model_data <- list(
  y = y,
  p = p,
  N = N,
  u0 = u0,
  I = I,
  aw = aw,
  bw = bw,
  as = as,
  bs = bs
)

# parameters
model_parameters <- c("O","u","iw","is")

# initial values
set.seed(2023)
initial_values <- list(
  list(
    "u" = matrix(rnorm(N*p),N,p),
    "O" = c(rep(0,N))
  )
)

# mcmc settings
nburn  <- 0 # 5000
nthin  <- 1 # 10
niter  <- 55000
nchain <- length(initial_values)

# MCMC ----
set.seed(2023)
# t <- Sys.time()
# fit <- jags(data = model_data, inits = initial_values, 
#             parameters.to.save = model_parameters, model.file = model, 
#             n.chains = nchain, n.iter = niter, n.thin = nthin, n.burnin = nburn)
# Sys.time() - t

# save(fit,file = paste0(path_code,"refTribMod.RData"))
load(paste0(path_code,"refTribMod.RData"))
(fit$BUGSoutput$DIC -> DIC)

## Thin and burn ----

### convergence  functions ----

library(coda)

ESS <- function(x){
  effectiveSize(x)
}

ME <- function(x){
  sd(x)/sqrt(ESS(x))
}

CV <- function(x){ # correr
  sd(x)/abs(mean(x))
}

check_convergence <- function(x,lag.max = 4){
  if(is.null(dim(x))){
    
    res <- c(ESS(x),ME(x),CV(x),acf(x,lag.max = lag.max,plot = F)$acf[-1])
    names(res) <- c("TEM","EM","CV",paste0("AC",1:lag.max))
    return(res)
  
  }else{
    
    res <- apply(x,MARGIN = 2,FUN = check_convergence)
    return(t(res))
    
  }
}

plot_convergence <- function(x){
  if(is.null(dim(x))){
    stop("x must contain samples of several parameters.")
  }else{
    conv <- check_convergence(x)
    par(mfrow=c(2,3),bty = "l")
    for(i in 1:6){
      hist(conv[,i],main = colnames(conv)[i],
           freq = F, col = "#8a90fd",
           border = "#8a90fd",xlab = "",
           ylab = "")
      # grid()
      # plot(density(conv[,i]),main = colnames(conv)[i],
      #      col = "#8a9090",xlab = "")
      # grid()
      # polygon(density(conv[,i]),col = "#8a90fd")
    }
    
  }
}

### Burn ----

B <- dim(fit$BUGSoutput$sims.array)[1]
dev_samples <- fit$BUGSoutput$sims.array[,1,(N+1)]
# plot(dev_samples,type = "l",col = "#bbbbbb")

check_convergence(dev_samples,lag.max = 10)
{
  pdf(paste0(path_code,"brute_sample_loglike_chain.pdf"),
      height = 5,width = 8)
  par(bty = "l")
  plot(-0.5*dev_samples,type = "l",col = "#bbbbbb",
       ylab = "Log-verosimilitud",
       xlab = "Muestra")
  dev.off()
}
# Quemamos las primeras 5k

fit$BUGSoutput$sims.array <- fit$BUGSoutput$sims.array[-(1:5000),,]

### Thin ----

B <- dim(fit$BUGSoutput$sims.array)[1]
dev_samples <- fit$BUGSoutput$sims.array[,(N+1)]
check_convergence(dev_samples,lag.max = 10)
# plot(dev_samples[100:200],type = "l",col = "#bbbbbb")

# Adelgazamos cada 10
new_indexes <- seq(1,50000,by = 10)
fit$BUGSoutput$sims.array <- fit$BUGSoutput$sims.array[new_indexes,]
dev_samples <- fit$BUGSoutput$sims.array[,(N+1)]
check_convergence(dev_samples,lag.max = 10)
# plot(dev_samples[100:200],type = "l",col = "#bbbbbb")
{
  pdf(paste0(path_code,"final_sample_loglike_chain.pdf"),
      height = 5,width = 8)
  par(bty = "l")
  plot(-0.5*dev_samples,type = "l",col = "#bbbbbb",
       ylab = "Log-verosimilitud",
       xlab = "Muestra")
  dev.off()
}

# post-model ----

# extract samples
B <- dim(fit$BUGSoutput$sims.array)[1]

## deviance ----
dev_samples <- fit$BUGSoutput$sims.array[,(N+1)]


write.csv(check_convergence(dev_samples),paste0(path_code,"deviance_convergence.txt"))


## O ----
fit$BUGSoutput$sims.array[1:B,1:N] -> O_samples

pdf(paste0(path_code,"O_convergence.pdf"),height = 6,width = 9)
plot_convergence(O_samples)
dev.off()



### Posterior mean ----
O_hat <- colMeans(O_samples)
names(O_hat) <- NULL


## Identifiability ----

### MLE ----
# Finding MLE between sampled Us
which.min(dev_samples) -> MLE_id
U_MLE <- matrix(fit$BUGSoutput$sims.array[MLE_id,(N+4):(3*N + 3)],ncol = 2)

### Procrustean transformation ----

lapply(1:B,function(i){
  u <- matrix(fit$BUGSoutput$sims.array[i,(N+4):(3*N+3)],ncol = 2)
  return(u)
}) -> u_samples

#### Hoffs formula ----

invsqrtm <- function(A,p){
  A.eigen <- eigen(A)
  A.v <- A.eigen$vectors[,1:p]
  A.l <- A.eigen$values[1:p]
  return(A.v%*%diag(1/sqrt(A.l))%*%t(A.v))
}

procrus <- function(Z,Z0,p){
  A <- t(Z)%*%Z0%*%t(Z0)%*%Z
  return(Z%*%invsqrtm(A,p)%*%t(Z)%*%Z0)
}

# u_transformed <- array(NA,dim = c(N,p,B))
# 
# {
#   t  <- Sys.time()
#   tt <- Sys.time()
#   for(i in 1:B){
#     u_transformed[,,i] <- procrus(u_samples[[i]],U_MLE,p)
#     if(Sys.time()-tt>10){
#       cat(round(i/B*100),"% completado.\n",sep = "")
#       tt <- Sys.time()
#     }
#   }
#   cat("Procrustean transformation finished in ",Sys.time()-t,"\n")
# }
# save(u_transformed,file = paste0(path_code,"model_procrustean_transformed_samples.RData"))


### posterior mean ----
load(file = paste0(path_code,"model_procrustean_transformed_samples.RData"))
apply(u_transformed,MARGIN = c(1,2),FUN = mean) -> u_mean

pdf(paste0(path_code,"U_convergence.pdf"),height = 6,width = 9)
plot_convergence(apply(u_transformed,MARGIN = 3,FUN = c))
dev.off()

## I ----
sapply(1:B,function(i){
  u <- matrix(fit$BUGSoutput$sims.array[i,(N+4):(3*N + 3)],ncol = 2)
  I <- sqrt(diag(u%*%t(u)))
  return(I)
}) -> I_samples
I_samples <- t(I_samples)

pdf(paste0(path_code,"I_convergence.pdf"),height = 6,width = 9)
plot_convergence(I_samples)
dev.off()

I_hat <- colMeans(I_samples)


### Correlaciones ----

dim(O_samples)
mapply(cor,
       x = as.list(as.data.frame(t(O_samples))),
       y = as.list(as.data.frame(t(I_samples)))
) -> cor_O_I

{
  pdf(paste0(path_code,"O_I_corr_density.pdf"),
      height = 4,width = 5)
  par(bty = "l", mar = c(5,4,2,2))
  plot(density(cor_O_I),main = NA,
       xlab = expression(r[O][I]),
       ylab = "Densidad")
  grid()
  polygon(density(cor_O_I),
          col = my.pal[1],
          border = my.pal[1])
  dev.off()
}

quantile(tmp, probs = c(0.025,0.975))

corr_data <- cbind(O_hat,I_hat,degree(G,mode = "out"),
                   degree(G,mode = "in"),
                   # eigen_centrality(G,directed = F)$vector)
                   betweenness(G))
colnames(corr_data) <- c("O","I","d_out","d_in","betw_c")
write.csv(cor(corr_data),paste0(path_code,"correlations.txt"))


## tau ----
t <- Sys.time()
tau_sum <- matrix(0,N,N)
for(i in 1:B){
  u <- u_transformed[,,i]
  I <- I_samples[i,]
  tau_sum <- tau_sum + (u%*%t(u))/(I%*%t(I))
  if(Sys.time()-t>5){
    t <- Sys.time()
    cat(round(i/B*100,2),"% completed.\n",sep = "")
  }
}
tau_mean <- tau_sum/B




## Circular space ----

rowsNorms <- function(X){
  return(sqrt(diag(X%*%t(X))))
}
normalizeRows <- function(X){
  X_norms     <- rowsNorms(X)
  X_norms_mat <- X_norms%*%t(rep(1,ncol(X)))
  return(X/X_norms_mat)
}


# This posterior mean does not falls in the unit circle:
the_sum <- matrix(0,N,p)

t <- Sys.time()
for(i in 1:B){
  ui <- u_transformed[,,i]
  the_sum <- the_sum + normalizeRows(ui)
  if(Sys.time() - t > 2){
    t <- Sys.time()
    cat(round(i/B*100),"% completado.\n",sep="")
  }
}

the_pos1 <- the_sum/B
the_pos2 <- normalizeRows(the_pos1)

## Otros Parametros ----

sigma2_samples <- 1/fit$BUGSoutput$sims.array[,N+2]
omega2_samples <- 1/fit$BUGSoutput$sims.array[,N+3]

sigma2_mean <- mean(sigma2_samples)
omega2_mean <- mean(omega2_samples)


## GOF ----

graph_summary <- function(G){
  res <- c()
  res[1] <- vcount(G)
  res[2] <- edge_density(G)
  res[3] <- transitivity(G,type = "global")
  res[4] <- assortativity_degree(G,directed = F)
  res[5] <- mean_distance(G,directed = F)
  res[6] <- mean(degree(G,mode = "all"))
  res[7] <- sd(degree(G,mode = "all"))
  # Clustering modularity
  Gu     <- as.undirected(G,mode = "collapse")
  cl     <- cluster_fast_greedy(Gu)
  res[8] <- modularity(cl)
  
  
  names(res) <- c("Size","Density","Transitivity",
                  "Assortativity","MeanDistance",
                  "MeanDegree","SdtDegree",
                  "ClModularity")
  return(res)
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

posterior_stats <- NULL
options("width")$width - 5 -> bar_width
t <- Sys.time()
for(i in 1:B){
  # Linear predictor
  Oi              <- O_samples[i,]%*%t(rep(1,N))
  utu             <- u_transformed[,,i]%*%t(u_transformed[,,i])
  u_norm          <- (sqrt(diag(u_transformed[,,i]%*%t(u_transformed[,,i])))%*%t(rep(1,N)))
  linPred         <- Oi +  utu / u_norm 
  probs           <- ilogit(linPred)
  # Sample network
  y_new           <- matrix(rbinom(N*N,1,probs),ncol = N) 
  diag(y_new)     <- 0 
  G_new           <- graph_from_adjacency_matrix(y_new,mode = "directed")
  # Statistics
  st_new          <- graph_summary(G_new)
  posterior_stats <- rbind(posterior_stats,st_new)
  # Report
  if(Sys.time()-t>5){
    t <- Sys.time()
    cat("[",strrep("=",floor(i/B*bar_width)),
        strrep(" ",bar_width - floor(i/B*bar_width)),
        "]\r",sep = "")
  }
}
cat("\n")

rownames(posterior_stats) <- NULL
st0 <- graph_summary(G)


# Save ----
# save(u_mean,tau_mean,I_hat,O_hat,dev_samples,U_MLE,B,N,G,posterior_stats,st0,
#      O_samples,I_samples,the_pos1,the_pos2,p,sigma2_mean,
#      omega2_mean,sigma2_samples,omega2_samples,
#      file = paste0(path_code,"model_posterior_processed_samples.RData"))





# Plots ----
library(paletteer) 
my.pal <- paletteer_d("jcolors::pal2")
aw <- 1; bw <- 1; as <- 1; bs <- 1 # caso 1
# aw <- 2; bw <- 1; as <- 2; bs <- 1 # caso 2
# aw <- 3; bw <- 2; as <- 3; bs <- 2 # caso 3
path_code <- paste0("out_aw_",aw,"_bw_",bw,"_as_",as,"_bs_",bs,"__")

# rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")
library(igraph)
load("out_refTrib.RData")
load(paste0(path_code,"model_posterior_processed_samples.RData"))

## deviance ----
pdf(paste0(path_code,"loglike_trace.pdf"),height = 7, width = 9)
plot(-0.5*dev_samples,type = "l",ylab = "Log-Likelihood",
     xlab = "Iteration")
dev.off()

## O ----
plot(O_samples[,1],type = "l")
plot(density(O_samples[,1]),xlab = expression(O[1]),ylab = "Density")

# Colors
my_palette <- colorRampPalette(c("blue","orange"))
my_palette <- colorRampPalette(my.pal[5:4])
cols_grid <- my_palette(100)
cols <- cols_grid[round((O_hat-min(O_hat))/(max(O_hat)-min(O_hat))*100)]
O_cols <- cols


n_extremes <- 3
extremes <- (rank(O_hat) <= n_extremes) | (rank(O_hat) > N - n_extremes)

pdf(paste0(path_code,"compare_O_posterior_trace.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,4))
graphics::matplot(O_samples[,extremes],lty = 1,
        type = "l",col = scales::alpha(cols[extremes],0.5),
        ylab = expression(O[i]),
        xlab = "Iteración")
image.plot(legend.only = T,zlim=range(O_hat),
            col = cols_grid,add = TRUE,
            legend.shrink = 0.8,
            horizontal = F)
dev.off()

pdf(paste0(path_code,"compare_O_posterior_density.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,4))
plot(NULL,NULL,xlim = c(-100,0),ylim = c(0,3.5),
     xlab = expression(O[i]),ylab = "Densidad")
for(ex in which(extremes)){
  lines(density(O_samples[,ex]),lwd = 1, col = cols[ex])
}
image.plot(legend.only = T,zlim=range(O_hat),
           col = cols_grid,add = TRUE,
           legend.shrink = 0.8,
           horizontal = F)
dev.off()




### distrbutions of posterior means ----

par(bty = "l")
plot(density(O_hat),xlab = ,
     ylab = "Densidad")

{
  pdf(paste0(path_code,"O_posterior_mean.pdf"),height = 3,width = 4)
  par(bty = "l",mar = c(4,4,2,2))
  plot(density(O_hat),xlab = expression(paste(widehat(E)(O[i],"|",bold(Y)))),
       ylab = "Densidad",main = NA)
  grid()
  polygon(density(O_hat),border = my.pal[1],
          col = my.pal[1])
  dev.off()
}



### FDA ----

dens.list <- NULL
for(i in 1:ncol(O_samples)){
  dens.list[[i]] <- density(O_samples[,i],
                            from = min(O_samples),
                            to = max(O_samples),
                            n = 1024)
}

fd.grid   <- dens.list[[1]]$x
fd.values <- t(sapply(dens.list,function(l)l$y))

library(fda)

deriv_pen <- int2Lfd(0)
base      <- create.bspline.basis(range(fd.grid),nbasis = 50,norder = 5)
dat_fdpar <- fdPar(base,deriv_pen,lambda = 1)
O_fd      <- smooth.basis(argvals=fd.grid,y=t(fd.values),fdParobj=dat_fdpar)
O_fd      <- O_fd$fd

plot(O_fd)
methods(class = class(O_fd))

O_pca     <- pca.fd(O_fd,nharm = 10)
cumsum(O_pca$varprop)

O_fd_z    <- O_pca$scores[,1:3]
O_fd_lamb <- O_pca$values[1:3]
O_fd_w    <- O_fd_z/(rep(1,N)%*%t(sqrt(O_fd_lamb)))

plot(O_fd_z,pch = 20)
which.max(O_fd_z[,2])

O_fd_dist <- matrix(0,N,N)
for(i in 1:N){
  for(j in 1:(i-1)){
    w_w <- O_fd_w[i,] - O_fd_w[j,]
    O_fd_dist[i,j] <- O_fd_dist[j,i] <- sqrt(t(w_w)%*%w_w)
  }
}

library(cluster)
pam(O_fd_dist,2) -> PAM_fd
PAM_fd$clustering -> cluster_fd_O

library(dplyr)
tibble(functional_cluster = cluster_fd_O,
       I_hat      = I_hat,
       O_hat      = O_hat,
       in_degree  = degree(G,mode = "in"),
       out_degree = degree(G,mode = "out")) -> vertex_by_cluster

vertex_by_cluster %>% 
  group_by(functional_cluster) %>% 
  summarize(I_hat      = mean(I_hat),
            O_hat      = mean(O_hat),
            in_degree  = mean(in_degree),
            out_degree = mean(out_degree))

wilcox.test(I_hat ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(O_hat ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(in_degree ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(out_degree ~ functional_cluster, data = vertex_by_cluster)


hist(O_hat)
cor(O_hat>-20,cluster_fd_O,method = "spearman")

pdf(paste0(path_code,"network_by_influence_cluster.pdf"),height = 8,width = 8)
par(mfrow = c(1,1),mar = c(0,0,2,0))
set.seed(1)
plot(G,vertex.size = 2, vertex.label = "",
     vertex.color = cluster_fd_O,
     main = "Network by Influence Cluster")
dev.off()

## Indentifiability ----

### MLE ----
plot(U_MLE,pch = 20,main = "Social space pseudo-MLE")

### U ----
apply(u_mean,2,range)
{
  pdf(paste0(path_code,"latent_space.pdf"),height = 6, width = 7)
  plot(u_mean,pch = 20,
       xlab = "1ra dimensión", ylab = "2da dimensión",
       main = expression(paste("Media posterior de ",
                               bold(u)[i],"  ",i== 1,ldots,N)),
       xlim = c(-3,3),ylim  = c(-3,3))
  # points(U_MLE,pch = 20, col = "#99999960")
  abline(h = 0, v = 0, col = "#dddddd")
  dev.off()
}

#### Posterior social space by out score ----
pdf(paste0(path_code,"latent_space_by_O.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,5))
plot(u_mean,pch = 20,col = cols,
     main = expression(paste(widehat(E)(bold(u)[i],"|", bold(Y)),"  ",i==1,",",ldots,",",N)),
     #sub = expression(paste(widehat(E)(O[i],"|",bold(Y)),": Bajo = Verde, Alto = Naranja")),
     xlab = "1ra dimensión",ylab = "2da dimensión")
image.plot(legend.only = T,zlim=range(O_hat),
           col = cols_grid,add = TRUE,
           legend.shrink = 0.8,
           horizontal = F)
grid()
abline(h = 0, v = 0, col = "#dddddd")
dev.off()


## I ----
plot(I_samples[,1],type = "l",main = expression(paste(I[1]," trace")))
plot(density(I_samples[,1]))
matplot(I_samples[,1:4],type = "l")

index_cols <- ceiling((I_hat-min(I_hat))/(max(I_hat)-min(I_hat))*100)
cols <- cols_grid[ifelse(index_cols==0,1,index_cols)]

n_extremes <- 3
extremes <- (rank(I_hat) <= n_extremes) | (rank(I_hat) > N - n_extremes)

pdf(paste0(path_code,"compare_I_posterior_trace.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,4))
graphics::matplot(I_samples[,extremes],lty = 1,
        type = "l",col = scales::alpha(cols[extremes],0.2),
        ylab = expression(I[i]),
        xlab = "Iteracion")
image.plot(legend.only = T,zlim=range(I_hat),
            col = cols_grid,add = TRUE,
            legend.shrink = 0.8,
            horizontal = F)
dev.off()

pdf(paste0(path_code,"compare_I_posterior_density.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,4))
plot(NULL,NULL,xlim = c(0,6),ylim = c(0,1),
     xlab = expression(I[i]),ylab = "Densidad")
for(ex in which(extremes)){
  lines(density(I_samples[,ex]),lwd = 1, col = cols[ex])
}
image.plot(legend.only = T,zlim=range(I_hat),
           col = cols_grid,add = TRUE,
           legend.shrink = 0.8,
           horizontal = F)
dev.off()

### distributions of posterior means ----
plot(density(I_hat))
{
  pdf(paste0(path_code,"I_posterior_mean.pdf"),height = 3,width = 4)
  par(bty = "l",mar = c(4,4,2,2))
  plot(density(I_hat),xlab = expression(paste(widehat(E)(O[i],"|",bold(Y)))),
       ylab = "Densidad",main = NA)
  grid()
  polygon(density(I_hat),border = my.pal[2],
          col = my.pal[2])
  dev.off()
}

### FDA ----
dens.list <- NULL
for(i in 1:ncol(I_samples)){
  dens.list[[i]] <- density(I_samples[,i],
                            from = min(I_samples),
                            to = max(I_samples),
                            n = 1024)
}

fd.grid   <- dens.list[[1]]$x
fd.values <- t(sapply(dens.list,function(l)l$y))

# library(fda)

deriv_pen <- int2Lfd(0)
base      <- create.bspline.basis(range(fd.grid),nbasis = 50,norder = 5)
dat_fdpar <- fdPar(base,deriv_pen,lambda = 1)
I_fd      <- smooth.basis(argvals=fd.grid,y=t(fd.values),fdParobj=dat_fdpar)
I_fd      <- I_fd$fd

plot(I_fd)
methods(class = class(I_fd))

pdf(paste0(path_code,"I_fd_boxplot.pdf"),height = 5,width=7)
boxplot(I_fd,xlab = "I",ylab = "Density",
        main = "Functional boxplot of I posterior density",
        )#ylim = c(0,1))
abline(h = 0, col = "#dddddd")
dev.off()


I_pca     <- pca.fd(I_fd,nharm = 10)
cumsum(I_pca$varprop)

I_fd_z    <- I_pca$scores[,1:2]
I_fd_lamb <- I_pca$values[1:2]
I_fd_w    <- I_fd_z/(rep(1,N)%*%t(sqrt(I_fd_lamb)))

plot(I_fd_z,pch = 20)
which.max(I_fd_z[,2])

I_fd_dist <- matrix(0,N,N)
for(i in 1:N){
  for(j in 1:(i-1)){
    w_w <- I_fd_w[i,] - I_fd_w[j,]
    I_fd_dist[i,j] <- I_fd_dist[j,i] <- sqrt(t(w_w)%*%w_w)
  }
}

library(cluster)
pam(I_fd_dist,2) -> PAM_fd
PAM_fd$clustering -> cluster_fd_I

library(dplyr)
tibble(functional_cluster = cluster_fd_I,
       I_hat      = I_hat,
       O_hat      = O_hat,
       in_degree  = degree(G,mode = "in"),
       out_degree = degree(G,mode = "out")) -> vertex_by_cluster

vertex_by_cluster %>% 
  group_by(functional_cluster) %>% 
  summarize(I_hat      = mean(I_hat),
            O_hat      = mean(O_hat),
            in_degree  = mean(in_degree),
            out_degree = mean(out_degree))

wilcox.test(I_hat ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(O_hat ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(in_degree ~ functional_cluster, data = vertex_by_cluster)
wilcox.test(out_degree ~ functional_cluster, data = vertex_by_cluster)


table(cluster_fd_I,cluster_fd_O)

hist(I_hat)
cor(I_hat,cluster_fd_I)

pdf(paste0(path_code,"network_by_influenciability_cluster.pdf"),height = 8,width = 8)
par(mfrow = c(1,1),mar = c(0,0,2,0))
set.seed(1)
plot(G,vertex.size = 2, vertex.label = "",
     vertex.color = cluster_fd_I,
     main = "Network by Influenciability Cluster")
dev.off()

## tau ----

library(fields)
cols <- designer.colors(64,c("blue","white","red"))
cols <- designer.colors(64,c(my.pal[5],"white",my.pal[4]))
pdf(paste0(path_code,"tau.pdf"),width = 8,height = 7)
image.plot(1:N,1:N,tau_mean,col = cols,zlim = c(-1,1),
           xlab = expression(i),ylab = expression(j))
dev.off()

### clustering ----

library(cluster)
k <- 2

#### k-Medoids ----
pam(0.5-0.5*tau_mean,k = k) -> PAM
PAM$clustering -> id_groups_PAM

#### k-Means ----
kmeans(u_mean,k) -> KM
KM$cluster -> id_groups_KM

#### Walktrap ----
# https://igraph.org/r/html/latest/cluster_walktrap.html
# plot(cluster_walktrap(G),G,vertex.size = 2,vertex.label = "")
cluster_walktrap(G) -> WT
is_hierarchical(WT)
cut_at(WT,no = 4) -> id_groups_WT

#### Edge betweenness ----
# https://igraph.org/r/html/latest/cluster_edge_betweenness.html
# plot(cluster_edge_betweenness(G),G,vertex.size = 2,vertex.label = "")
# cluster_walktrap(G) -> EB
# is_hierarchical(EB)
# cut_at(EB,no = 4) -> id_groups_EB

#### Fast greedy ----
cluster_fast_greedy(as.undirected(G,mode = "collapse")) -> FG
cut_at(FG,no = 2)-> id_groups_FG

table(id_groups_PAM,id_groups_KM)
# table(id_groups_EB,id_groups_WT)
table(id_groups_WT,id_groups_FG)
table(id_groups_KM,id_groups_FG)

{
  pdf(paste0(path_code,"network_clusters.pdf"),height = 8,width = 8)
  par(mfrow = c(2,2),mar = c(0,0,2,0))
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_PAM,
       main = expression(paste("k-medoids en ",tau[ij],"  ",list(i,j==1,ldots,N))),edge.arrow.size = 0.3)
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_KM,
       main = expression(paste("k-means en ",bold(u)[i],"  ",list(i==1,ldots,N))),edge.arrow.size = 0.3)
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_WT,
       main = expression(paste("Walktrap en ",G)),edge.arrow.size = 0.3)
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_FG,
       main = expression(paste("Fast Greedy en ",G)),edge.arrow.size = 0.3)
  dev.off()
}


id_groups <- id_groups_KM
id_groups <- id_groups_FG


pdf(paste0(path_code,"tau_by_clusters.pdf"),width = 8,height = 7)
resort_index <- unlist(lapply(1:k,function(x)which(id_groups==x)))
image.plot(1:N,1:N,tau_mean[resort_index,resort_index],
           col = cols,xlab = expression(i),
           ylab = expression(j),zlim = c(-1,1),
           xaxt = "n",yaxt = "n")
mtext(resort_index[(1:10)*60],side = 1,at = (1:10)*60)
mtext(resort_index[(1:10)*60],side = 2,at = (1:10)*60)
dev.off()


### I by cluster ----

wilcox.test(I_hat[id_groups==1],
            I_hat[id_groups==2])

tapply(I_hat,factor(id_groups),mean)

boxplot(I ~ cl,data = data.frame(
  I = I_hat, cl = factor(id_groups)
))

### Posterior social space by group ----
pdf(paste0(path_code,"latent_space_by_cluster.pdf"),height = 6,width = 7)
plot(u_mean,pch = 20,col = my.pal[3+id_groups],
     main = expression(paste(widehat(E)(bold(u)[i],"|", bold(Y)),"  ",i==1,",",ldots,",",N)),
     xlab = "1ra dimensión",ylab = "2da dimensión")
grid()
abline(h = 0, v = 0, col = "#dddddd")
dev.off()

### Posterior social space by group and O ----
pdf(paste0(path_code,"latent_space_by_O_and_cluster.pdf"),height = 6,width = 7)
par(bty = "l",mar = c(5,4,4,5))
plot(u_mean,pch = c(15,17)[id_groups],
     col = scales::alpha(O_cols,0.3),
     main = expression(paste(widehat(E)(bold(u)[i],"|", bold(Y)),"  ",i==1,",",ldots,",",N)),
     xlab = "1ra dimensión",ylab = "2da dimensión")
image.plot(legend.only = T,zlim=range(O_hat),
           col = cols_grid,add = TRUE,
           legend.shrink = 0.8,
           horizontal = F)
grid()
abline(h = 0, v = 0, col = "#dddddd")
dev.off()



## circular space ----

{
  pdf(paste0(path_code,"political_space_1.pdf"),width = 5,height = 5)
  xgrid <- seq(-1,1,length = 300)
  plot(xgrid,sqrt(1-xgrid^2),type = "l",col = "gray",
       xlim = c(-1,1),ylim = c(-1,1),
       xlab = "1ra dimensión",ylab = "2da dimensión")
  abline(h = 0, v = 0, col = "#dddddd")
  lines(xgrid,-sqrt(1-xgrid^2),col = "gray")
  points(the_pos1,pch = 20, col = my.pal[3+id_groups])
  dev.off()
}

# plot(normalizeRows(the_pos1))
the_pos2 <- normalizeRows(the_pos1)
the_pos2

infl_id <- which(635 - rank(O_hat) <= 32)
infl <- the_pos2[infl_id,]

infl_diff_id <- which( ((infl[,2]>0)*(infl[,1] < 0.1))==1)

infl_id[infl_diff_id]


{
  pdf(paste0(path_code,"political_space_2.pdf"),width = 5,height = 5)
  xgrid <- seq(-1,1,length = 300)
  plot(xgrid,sqrt(1-xgrid^2),type = "l",col = "gray",
       xlim = c(-1,1),ylim = c(-1,1),
       xlab = "1ra dimensión",ylab = "2da dimensión")
  abline(h = 0, v = 0, col = "#dddddd")
  lines(xgrid,-sqrt(1-xgrid^2),col = "gray")
  points(the_pos2,pch = 20,
         col = scales::alpha(my.pal[3+id_groups],0.55))
  dev.off()
}

plot_influence_paths <- function(G,the_pos,O_hat,id_groups = 1:length(O_hat),n = 4,
                                 letters = F){
  letras <- c("a","b","c","d","e","f",
              "g","h","i","j","k","l",
              "m","n","o","p","q","r",
              "s","t","u","v","w","x",
              "y","z")
  # Indetify relationships
  top <- length(O_hat) + 1 - rank(O_hat)
  influencers <- c()
  for(i in 1:n){
    influencers <- c(influencers,which(top == i))
  }
  influenced <- NULL
  for(i in 1:n){
    influenced[[i]] <- as.numeric(neighbors(G,influencers[i],mode = "out"))
  }
  
  # Prepare plot
  xgrid_outer <- seq(-1,1,length = 300)
  xgrid_inner <- seq(-0.5,0.5,length = 300)
  plot(xgrid_outer,sqrt(1-xgrid_outer^2),type = "l",col = "gray",
       xlim = c(-1,1),ylim = c(-1,1),
       xlab = "1ra dimensión",ylab = "2da dimensión",
       main = "Posición de los influenciadores")
  # abline(h = 0, v = 0, col = "#dddddd")
  lines(xgrid_outer,-sqrt(1-xgrid_outer^2),col = "gray")
  lines(xgrid_inner,sqrt(0.25-xgrid_inner^2),col = "gray")
  lines(xgrid_inner,-sqrt(0.25-xgrid_inner^2),col = "gray")
  
  # plot relations
  for(i in 1:n){
    for(j in 1:length(influenced[[i]])){
      segments(x0 = the_pos[influencers[i],1]/2,y0 = the_pos[influencers[i],2]/2,
               x1 = the_pos[influenced[[i]][j],1], y1 = the_pos[influenced[[i]][j],2],
               col = "#bbbbbb70",lty = 2,lwd = 0.5)
    }
  }
  
  # plot individuals
  if(n<=length(letras)){
    influencers_pch <- letras[1:n]
  }else{
    influencers_pch <- rep(4,n)
  }
  points(the_pos[unlist(influenced),],pch = 6, col = my.pal[3+id_groups[unlist(influenced)]])
  points(the_pos[influencers,]/2,pch = influencers_pch,
         col = my.pal[3+id_groups[influencers]])
  
  
  cat("Top ",n," influentials reach ",round(length(unique(unlist(influenced)))/vcount(G)*100,2),"% of the network.\n",sep = "")
  
}

# Influentials
pdf(paste0(path_code,"influence_paths_influentials.pdf"),height = 8,width = 9)
n_influentials <- min(table(kmeans(O_hat,centers = 2)$cluster))
plot_influence_paths(G,the_pos2,O_hat,id_groups,n_influentials)
dev.off()
# Top must influentials
pdf(paste0(path_code,"influence_paths_top_5perc.pdf"),height = 8,width = 9)
plot_influence_paths(G,the_pos2,O_hat,id_groups,round(0.05*N))
dev.off()

pdf(paste0(path_code,"influence_paths_top_2perc.pdf"),height = 8,width = 9)
plot_influence_paths(G,the_pos2,O_hat,id_groups,round(0.02*N),letters = T)
dev.off()


## GOF ----

plots_names <- c(NA,
                 "Densidad",
                 "Transitividad",
                 "Asortatividad",
                 "Distancia Media",
                 NA,"Desviación Estándar del Grado",
                 "Modularidad de Agrupamiento")
pdf(paste0(path_code,"GOF.pdf"),height = 8,width = 5)
par(mfrow = c(3,2),mar = c(4,4,4,2))
for(k in c(2,3,4,5,7,8)){
  R <- range(c(quantile(posterior_stats[,k],probs = c(0.001,0.999)),st0[k]))
  plot(density(posterior_stats[,k]),
       main = plots_names[k],
       xlim = R, sub = "", lwd = 2,xlab = "",
       ylab = "Densidad")
  grid()
  polygon(density(posterior_stats[,k]),col = "#dddddda0",
          border = "#00000000")
  abline(v = st0[k], lwd = 2, col = "red")
}
dev.off()
