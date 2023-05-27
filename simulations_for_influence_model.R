# Simulations for the influence model
rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")

# Packages ----

## Plots ----
library(paletteer) 
my.pal <- paletteer_d("jcolors::pal2")
## C++ ----
library(Rcpp)
library(RcppArmadillo)
## Networks ----
library(igraph)
## Modeling ----
library(R2jags)


# Functions ----

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

ESS <- function(x){
  effectiveSize(x)
}

ME <- function(x){
  sd(x)/sqrt(ESS(x))
}

CV <- function(x){ # correr
  sd(x)/abs(mean(x))
  # ME(x)/abs(mean(x))
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

giant_comp <- function(G){
  comp <- components(G)
  comp$csize[1]/vcount(G)
  G <- subgraph(G,which(comp$membership==1))
  return(G)
}


model <- function(){
  # likelihood
  for(i in 1:N){
    for(j in 1:N){
      y[i,j]        ~ dbern(ilogit(linPred[i,j]))
      linPred[i,j] <- O[i] + inprod(u[i,],u[j,])/sqrt(inprod(u[i,],u[i,]))
    }
  }
  
  # prior
  for(i in 1:N){
    O[i]  ~ dnorm(0,iw)
    u[i,1:p] ~ dmnorm(u0,is*I)
  }
  iw ~ dgamma(aw,bw)
  is ~ dgamma(as,bs)
}

graph_summary <- function(G){
  res <- c()
  res[1] <- round(vcount(G))
  res[2] <- edge_density(G)
  res[3] <- transitivity(G,type = "global")
  res[4] <- assortativity_degree(G,directed = F)
  res[5] <- mean_distance(G,directed = F)
  res[6] <- mean(degree(G,mode = "all"))
  res[7] <- sd(degree(G,mode = "all"))
  # Clustering Stats
  Gu     <- as.undirected(G,mode = "collapse")
  cl     <- cluster_fast_greedy(Gu)
  res[8] <- modularity(cl)
  
  
  names(res) <- c("Order","Density","Transitivity",
                  "Assortativity","MeanDistance",
                  "MeanDegree","SdtDegree",
                  "ClModularity")
  return(res)
}

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

multivariate_CI <- function(U,alpha){
  Ubar   <- colMeans(U)
  USigma <- cov(U)
  dists  <- mahalanobis(U,Ubar,USigma)
  quant  <- quantile(dists,probs = 1-alpha) 
  return(quant)
}



# Data ----

aw <- 1; bw <- 1; as <- 1; bs <- 1 # caso 1
# aw <- 2; bw <- 1; as <- 2; bs <- 1 # caso 2
# aw <- 3; bw <- 2; as <- 3; bs <- 2 # caso 3
path_code <- paste0("out_aw_",aw,"_bw_",bw,"_as_",as,"_bs_",bs,"__")

load(paste0(path_code,"model_posterior_processed_samples.RData"))

# Simulate network ----

set.seed(2023)
# Linear predictor
Oi              <- O_hat%*%t(rep(1,N))
utu             <- u_mean%*%t(u_mean)
u_norm          <- (sqrt(diag(u_mean%*%t(u_mean)))%*%t(rep(1,N)))
linPred         <- Oi +  utu / u_norm 
probs           <- ilogit(linPred)
# Sample network
y               <- matrix(rbinom(N*N,1,probs),ncol = N) 
diag(y)         <- 0 
G               <- graph_from_adjacency_matrix(y,mode = "directed")
vcount(G);ecount(G)
# Giant component
# G               <- giant_comp(G)
# vcount(G);ecount(G)

graph_summary(G)


## Plot network ----

{
pdf("simulated_network_base_plot.pdf",height = 6,width = 6)
par(mar = c(0,0,0,0))
set.seed(2023)
plot(G,vertex.size = log(1+degree(G)),
     layout = layout_with_fr,
     vertex.label = NA,
     vertex.color = "lightblue",
     edge.width = 0.5,
     edge.color = "#a0a0a080",
     edge.arrow.size = 0.2)
dev.off()
}


# Model ----

## Prior parameters ----
p <- 2
N <- vcount(G)
u0 <- c(rep(0,p))
I <- diag(p)
aw <- 1; bw <- 1; as <- 1; bs <- 1 # caso 1

## Input ----
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

## Parameters ----
model_parameters <- c("O","u","iw","is")

## Initial values ----
set.seed(2023)
initial_values <- list(
  list(
    "u" = matrix(rnorm(N*p),N,p),
    "O" = c(rep(0,N))
  )
)

## MCMC settings ----
nburn  <- 0 # 5000
nthin  <- 1 # 10
niter  <- 55000
nchain <- length(initial_values)

## MCMC ----
set.seed(2023)
# t <- Sys.time()
# fit <- jags(data = model_data, inits = initial_values, 
#             parameters.to.save = model_parameters, model.file = model, 
#             n.chains = nchain, n.iter = niter, n.thin = nthin, n.burnin = nburn)
# Sys.time() - t
# Time difference of 20.67402 hours

## Save/Load ----
# save(fit,file = paste0(path_code,"simulationModel.RData"))
load(paste0(path_code,"simulationModel.RData"))
(fit$BUGSoutput$DIC -> DIC)

# Process samples ----

SAMPS <- fit$BUGSoutput$sims.array[-(1:5000),1,]
new_indexes <- seq(1,50000,by = 20)
SAMPS <- SAMPS[new_indexes,]
B <- dim(SAMPS)[1]

## Log-Likelihood ----

dev_samples <- SAMPS[,(N+1)]

check_convergence(dev_samples)

{
  pdf(paste0(path_code,"simulated_network_loglike_chain.pdf"),
      height = 5,width = 8)
  par(bty = "l",mar = c(5,4,2,2))
  plot(-0.5*dev_samples,type = "l",col = "#bbbbbb",
       ylab = "Log-verosimilitud",
       xlab = "Muestra")
  dev.off()
}

## O ----
SAMPS[,1:N] -> O_samples

## U ----
lapply(1:B,function(i){
  u <- matrix(SAMPS[i,(N+4):((p+1)*N+3)],ncol = 2)
  return(u)
}) -> u_samples

u_refference <- u_mean
u_transformed <- array(NA,dim = c(N,p,B))
{
  t  <- Sys.time()
  tt <- Sys.time()
  for(i in 1:B){
    u_transformed[,,i] <- procrus(u_samples[[i]],u_refference,p)
    if(Sys.time()-tt>10){
      cat(round(i/B*100),"% completado.\n",sep = "")
      tt <- Sys.time()
    }
  }
  cat("Procrustean transformation finished in ",Sys.time()-t,"\n")
}

## Variances ----
sigma2_samples <- 1/SAMPS[,N+2]
omega2_samples <- 1/SAMPS[,N+3]

# Credibility intervals ----

alpha <- 0.05

## O ----
O_CI <- apply(O_samples,MARGIN = 2,
              FUN = function(samps){
                quantile(samps,
                         probs = c(alpha/2,
                                   1-alpha/2))
              })
O_CI <- t(O_CI)

## U ----
U_CI <- apply(u_transformed,MARGIN = 1,
              FUN = function(Ut){
                multivariate_CI(t(Ut),alpha)
              })

## Variances ----
sigma2_CI <- quantile(sigma2_samples,
                      probs = c(alpha/2,1-alpha/2))

omega2_CI <- quantile(omega2_samples,
                      probs = c(alpha/2,1-alpha/2))

# Save/Load ----
# save(O_CI,U_CI,sigma2_CI,omega2_mean,
#      file = paste0(path_code,"simulationCI.RData"))
load(paste0(path_code,"simulationCI.RData"))

# Contenencia ----

## O ----

sapply(1:N,FUN = function(i){
  (O_CI[i,1] < O_hat[i]) && (O_hat[i] < O_CI[i,2])
}) -> O_in
mean(O_in)
sum(!O_in)

## U ----

sapply(1:N,FUN = function(i){
  Ui <- t(u_transformed[i,,])
  mahalanobis(u_mean[i,],colMeans(Ui),cov(Ui))
}) -> dists

U_in <- dists < U_CI
mean(U_in)
sum(!U_in)

# Gráficos ----

## O ----

tmp_N <- N
{
  pdf("simulated_network_O_intervals.pdf",
      height = 6,width = 3)
  par(bty = "l",mar = c(5,4,2,2))
  plot(NULL,NULL,ylim = c(1,tmp_N),
       xlim = c(min(O_CI[,1]),
                max(O_CI[,2])),
       xlab = expression(O[i]),
       ylab = expression(i))
  segments(x0 = O_CI[1:tmp_N,1],
           y0 = 1:tmp_N,
           x1 = O_CI[1:tmp_N,2],
           y1 = 1:tmp_N,
           col = "#b0b0b080")
  points(x = O_hat[1:tmp_N],y = 1:tmp_N,
         col = "red",pch = '.')
  dev.off()
}

tmp_N <- 100
{
  pdf("simulated_network_O_intervals_zoom.pdf",
      height = 6,width = 3)
  par(bty = "l",mar = c(5,4,2,2))
  plot(NULL,NULL,ylim = c(1,tmp_N),
       xlim = c(min(O_CI[,1]),
                max(O_CI[,2])),
       xlab = expression(O[i]),
       ylab = expression(i))
  segments(x0 = O_CI[1:tmp_N,1],
           y0 = 1:tmp_N,
           x1 = O_CI[1:tmp_N,2],
           y1 = 1:tmp_N,
           col = "#b0b0b080")
  points(x = O_hat[1:tmp_N],y = 1:tmp_N,
         col = "red",pch = '.')
  dev.off()
}

## U ----

tmp_N <- N
{
  pdf("simulated_network_U_intervals.pdf",
      height = 6,width = 3)
  par(bty = "l",mar = c(5,4,2,2))
  plot(NULL,NULL,ylim = c(1,tmp_N),
       xlim = c(0,max(U_CI)),
       xlab = expression(m[i]),
       ylab = expression(i))
  segments(x0 = rep(0,tmp_N),
           y0 = 1:tmp_N,
           x1 = U_CI[1:tmp_N],
           y1 = 1:tmp_N,
           col = "#b0b0b080")
  points(x = dists[1:tmp_N],y = 1:tmp_N,
         col = "red",pch = '.')
  dev.off()
}
