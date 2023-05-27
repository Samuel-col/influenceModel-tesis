# JAGS implementation (null model)

# Header ----
rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")

library(R2jags)
library(igraph)

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
      linPred[i,j] <- O + inprod(u[i,],u[j,])/sqrt(inprod(u[i,],u[i,]))
    }
  }
  
  # prior
  O ~ dnorm(0,iw)
  for(i in 1:N){
    u[i,1:p] ~ dmnorm(u0,is*I)
  }
  iw ~ dgamma(aw,bw)
  is ~ dgamma(as,bs)
}

# prior parameters
p  <- 2
N  <- vcount(G)
u0 <- c(rep(0,p))
I  <- diag(p)
aw <- 1
bw <- 1
as <- 1
bs <- 1
path_code <- paste0("out_NULL_aw_",aw,"_bw_",bw,"_as_",as,"_bs_",bs,"__")

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
    "O" = 0
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

fit$BUGSoutput$DIC -> DIC
DIC

## Thin and burn ----

### convergence  functions ----

library(coda)

ESS <- function(x){
  effectiveSize(x)
}

ME <- function(x){
  sd(x)/sqrt(ESS(x))
}

CV <- function(x){
  sd(x)/abs(mean(x))
}

check_convergence <- function(x,lag.max = 4){
  if(is.null(dim(x))){
    
    res <- c(ESS(x),ME(x),CV(x),acf(x,lag.max = lag.max,plot = F)$acf[-1])
    names(res) <- c("ESS","ME","CV",paste0("acf.lag",1:lag.max))
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
    par(mfrow=c(2,3))
    for(i in 1:6){
      hist(conv[,i],main = colnames(conv)[i],
           freq = F, col = "#8a90fd",
           border = "#8a90fd",xlab = "")
    }
    
  }
}

### Burn ----

B <- dim(fit$BUGSoutput$sims.array)[1]
dev_samples <- fit$BUGSoutput$sims.array[,1,2]
# plot(dev_samples,type = "l",col = "#bbbbbb")

# Quemamos las primeras 5k

fit$BUGSoutput$sims.array <- fit$BUGSoutput$sims.array[-(1:5000),1,]

### Thin ----

B <- dim(fit$BUGSoutput$sims.array)[1]
dev_samples <- fit$BUGSoutput$sims.array[,2]
check_convergence(dev_samples,lag.max = 10)

# Adelgazamos cada 5
new_indexes <- seq(1,50000,by = 10)
fit$BUGSoutput$sims.array <- fit$BUGSoutput$sims.array[new_indexes,]
dev_samples <- fit$BUGSoutput$sims.array[,2]
check_convergence(dev_samples,lag.max = 10)

# post-model ----

# extract samples
B <- dim(fit$BUGSoutput$sims.array)[1]

## deviance ----
dev_samples <- fit$BUGSoutput$sims.array[,2]


write.csv(check_convergence(dev_samples),paste0(path_code,"deviance_convergence.txt"))


## O ----
fit$BUGSoutput$sims.array[,1] -> O_samples

check_convergence(O_samples)


### Posterior mean ----
O_hat <- mean(O_samples)


## Identifiability ----

### MLE ----
# Finding MLE between sampled Us
which.min(dev_samples) -> MLE_id
U_MLE <- matrix(fit$BUGSoutput$sims.array[MLE_id,5:(p*N + 4)],ncol = 2)

### Procrustean transformation ----

lapply(1:B,function(i){
  u <- matrix(fit$BUGSoutput$sims.array[i,5:(2*N+4)],ncol = 2)
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
  u <- matrix(fit$BUGSoutput$sims.array[i,5:(p*N + 4)],ncol = 2)
  I <- sqrt(diag(u%*%t(u)))
  return(I)
}) -> I_samples
I_samples <- t(I_samples)

pdf(paste0(path_code,"I_convergence.pdf"),height = 6,width = 9)
plot_convergence(I_samples)
dev.off()

I_hat <- colMeans(I_samples)


corr_data <- cbind(I_hat,degree(G,mode = "out"),
                   degree(G,mode = "in"),
                   # eigen_centrality(G,directed = F)$vector)
                   betweenness(G))
colnames(corr_data) <- c("I","d_out","d_in","betw_c")
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
  Oi              <- O_samples[i]*matrix(1,N,N)
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
save(u_mean,tau_mean,I_hat,O_hat,dev_samples,U_MLE,B,N,G,posterior_stats,st0,
     O_samples,I_samples,the_pos1,the_pos2,p,
     file = paste0(path_code,"model_posterior_processed_samples.RData"))



# Plots ----

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

## Indentifiability ----

### MLE ----
plot(U_MLE,pch = 20,main = "Social space pseudo-MLE")

### U ----
plot(u_mean,pch = 20,
     xlab = "dim 1", ylab = "2 dim",
     main = "Social sapace bayesian estimator")
# points(U_MLE,pch = 20, col = "#99999960")
abline(h = 0, v = 0, col = "#dddddd")


## I ----
plot(I_samples[,1],type = "l",main = expression(paste(I[1]," trace")))
plot(density(I_samples[,1]))
matplot(I_samples[,1:4],type = "l")

n_extremes <- 3
extremes <- (rank(I_hat) <= n_extremes) | (rank(I_hat) > N - n_extremes)

pdf(paste0(path_code,"compare_I_posterior_trace.pdf"),height = 6,width = 7)
matplot(I_samples[,extremes],lty = 1,
        type = "l",col = scales::alpha(cols[extremes],0.5),
        ylab = expression(I[i]),
        xlab = "Iteration")
dev.off()

pdf(paste0(path_code,"compare_I_posterior_density.pdf"),height = 6,width = 7)
plot(NULL,NULL,xlim = c(0,6),ylim = c(0,1),
     xlab = expression(I[i]),ylab = "Density")
for(ex in which(extremes)){
  lines(density(I_samples[,ex]),lwd = 1, col = cols[ex])
}
dev.off()

### distributions of posterior means ----
plot(density(I_hat))
hist(I_hat,breaks = 15,freq = F)

## tau ----

library(fields)
cols <- designer.colors(64,c("blue","white","red"))
pdf(paste0(path_code,"tau.pdf"),width = 8,height = 7)
image.plot(1:N,1:N,tau_mean,col = cols,zlim = c(-1,1))
dev.off()

### clustering ----

library(cluster)
k <- 2

#### k-Medoids ----
pam(tau_mean,k = k) -> PAM
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
cluster_walktrap(G) -> EB
is_hierarchical(EB)
cut_at(EB,no = 4) -> id_groups_EB

#### Fast greedy ----
cluster_fast_greedy(as.undirected(G,mode = "collapse")) -> FG
cut_at(FG,no = 2)-> id_groups_FG

table(id_groups_PAM,id_groups_KM)
table(id_groups_EB,id_groups_WT)
table(id_groups_WT,id_groups_FG)
table(id_groups_KM,id_groups_FG)

{
  pdf(paste0(path_code,"network_clusters.pdf"),height = 8,width = 8)
  par(mfrow = c(2,2),mar = c(0,0,2,0))
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_PAM,
       main = "k-medoids on tau")
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_KM,
       main = "k-means on U")
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_WT,
       main = "Walktrap")
  set.seed(1)
  plot(G,vertex.size = 2, vertex.label = "",vertex.color = id_groups_FG,
       main = "Fast Greedy")
  dev.off()
}


id_groups <- id_groups_KM


pdf(paste0(path_code,"tau_by_clusters.pdf"),width = 8,height = 7)
resort_index <- unlist(lapply(1:k,function(x)which(id_groups==x)))
image.plot(1:N,1:N,tau_mean[resort_index,resort_index],
           col = cols,xaxt = "n",yaxt = "n",zlim = c(-1,1))
dev.off()


### Posterior social space by group ----
plot(u_mean,pch = 20,col = id_groups,
     main = "Social sapace bayesian estimator by groups")
abline(h = 0, v = 0, col = "#dddddd")

plot(u_mean,pch = 20,col = id_groups_WT,
     main = "Social sapace bayesian estimator by groups")
abline(h = 0, v = 0, col = "#dddddd")




## circular space ----

{
  pdf(paste0(path_code,"social_space.pdf"),width = 5,height = 5)
  xgrid <- seq(-1,1,length = 300)
  plot(xgrid,sqrt(1-xgrid^2),type = "l",col = "gray",
       xlim = c(-1,1),ylim = c(-1,1),
       xlab = "1 dim",ylab = "2 dim")
  abline(h = 0, v = 0, col = "#dddddd")
  lines(xgrid,-sqrt(1-xgrid^2),col = "gray")
  points(the_pos,pch = 20, col = id_groups)
  dev.off()
}


## GOF ----

plots_names <- c(NA,
                 "Densidad",
                 "Transitividad",
                 "Asortativdad",
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
  cat(quantile(posterior_stats[,k],probs = c(0.025,0.975)),"  : ")
  cat(st0[k],"\n")
}
dev.off()
