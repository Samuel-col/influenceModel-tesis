# Simulate from influence model
rm(list = ls())
setwd("/home/samuel/Documentos/U/tesis/")

# Packages ----

## Plots ----
library(fields)
library(paletteer) 
my.pal <- paletteer_d("jcolors::pal2")
## C++ ----
library(Rcpp)
library(RcppArmadillo)
sourceCpp("cascade.cpp")
sourceCpp("optimize.cpp")

## Networks ----
library(igraph)


# Functions ----

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

giant_comp <- function(G){
  comp <- components(G)
  comp$csize[1]/vcount(G)
  G <- subgraph(G,which(comp$membership==1))
  return(G)
}


sample_covariables <- function(N,cl1_ids = 1:round(N/2),
                               seed = 2023,
                               mode = "random",
                               type = c("realistic",
                                        "realistic"),
                               shift = 0){
  
  cl2_ids <- setdiff(1:N,cl1_ids)
  set.seed(seed)
  if(type[1] == "realistic"){
    I <- rgamma(N,shape = 14/5,rate = 7/5)
  }else if(type[1] == "constant"){
    I <- rep((14/5)/(7/5),N)
  }else{
    stop()
  }
  if(type[2] == "realistic"){
    O <- rgamma(N,shape = 3,rate = 3/4) + shift
  }else if(type[2] == "constant"){
    O <- rep((3)/(3/4) + shift,N)
  }else{
    stop()
  }
  X <- rep(0,N)
  if(mode == "random"){
    set.seed(seed)
    X[sample(1:N,2)] <- 2:3
  }else if(mode == "random clust"){
    set.seed(seed)
    X[sample(cl1_ids,1)] <- 2
    X[sample(cl2_ids,1)] <- 3
  }else if(mode == "best"){
    id.1 <- which.max(O)
    tmp_m <- max(O[-id.1])
    id.2 <- which(O==tmp_m)
    X[id.1] <- 3
    X[id.2] <- 2
    cl.g.ids <- NULL
  }else if(mode == "best clust"){
    tmp_m <- max(O[cl1_ids])
    id.1 <- which(O == tmp_m)
    tmp_m <- max(O[cl2_ids])
    id.2 <- which(O == tmp_m)
    X[id.1] <- 3
    X[id.2] <- 2
  }else{
    stop()
  }
  return(list(
    X = X, O = O, I = I
  ))
}

sample_tau <- function(N,cl1_ids = 1:round(N/2),
                         cl = 1/2,seed = 2023){
  set.seed(seed)
  n1 <- length(cl1_ids)
  n2 <- N - n1
  tau <- matrix(0,N,N)
  tau[cl1_ids,cl1_ids]   <- 1 - rbeta(n1*n1,1-cl,1)
  tau[-cl1_ids,-cl1_ids] <- 1 - rbeta(n2*n2,1-cl,1)
  tau[cl1_ids,-cl1_ids]  <- rbeta(n1*n2,1-cl,1)
  tau[-cl1_ids,cl1_ids]  <- rbeta(n1*n2,1-cl,1)
  return(2*tau - 1)
}

sample_network <- function(N,cl1_ids = 1:round(N/2),
                           seed = 2023,mode = "random",
                           cl = 1/2,
                           type = c("realistic",
                                    "realistic"),
                           shift = 0){
  # Sample parameters
  tau  <- sample_tau(N = N,cl1_ids = cl1_ids,
                     cl = cl,seed = seed)
  covs <- sample_covariables(N = N, cl1_ids = cl1_ids,
                             seed = seed, mode = mode,
                             type = type,shift = shift)
  O <- covs$O
  I <- covs$I
  # Linear predictor
  O_mat    <- rep(1,N)%*%t(O)
  I_mat    <- I%*%t(rep(1,N))
  lin_pred <- O_mat + tau*I_mat
  # Sample network
  probs       <- ilogit(lin_pred)
  y_new       <- matrix(rbinom(N*N,1,probs),ncol = N)
  diag(y_new) <- 0 
  G_new       <- graph_from_adjacency_matrix(y_new,mode = "directed")
  
  obj <- list(G = G_new,X = covs$X,y = y_new,
              O = O, I = I,tau = tau)
  class(obj) <- c("influenceNetwork","list")
  return(obj)
}

density.influenceNetwork <- function(iN){
  edge_density(iN$G)
}

plot.influenceNetwork <- function(sn,
                                 vertex.size = 3,
                                 vertex.label = NA,
                                 edge.arrow.size = 0.3,
                                 vertex.color = my.pal[sn$X + 2],
                                 ...){
  plot(sn$G,vertex.size = vertex.size,
       vertex.label = vertex.label,
       edge.arrow.size = edge.arrow.size,
       vertex.color = vertex.color,...)
}


shift_to_meanDeg <- function(shift,N,cl1_ids = 1:round(N/2),
                         mode,cl,type,k = 5){
  sum_Den <- 0
  mean_Den <- median(sapply(1:k,FUN = function(j){
    density(sample_network(N,seed = 2023 + j,shift = shift,
                           type = type,cl = cl, mode = mode,
                           cl1_ids = cl1_ids))
  }))
  return(mean_Den*(N-1))
}


run_simulation <- function(N,cl1_ids = 1:round(N/2),
                           mode = "random",cl = 1/2,
                           type = c("realistic",
                                    "realistic")){
  cat("\n")
  cat("Clustering degree:",cl,"\n")
  cat("Initiators selection method:",mode,"\n")
  cat("I and O sampling method:",type[1],"and",type[2],"\n")
  
  cpp_optimize(function(shift){
    shift_to_meanDeg(shift,N = N, mode = mode,
                     cl = cl, type = type,
                     k = 7)
  },10)$shift -> shift
  
  # Simulate networks
  SN <- list()
  for(i in 1:4){
  SN[[i]] <- sample_network(N = N, cl1_ids = cl1_ids,
                            seed = i, mode = mode,
                            cl = cl,type = type,
                            shift = shift)
  }
  # Run cascade
  times <- c()
  cascades <- list()
  for(i in 1:4){
    t <- Sys.time()
    y <- as.matrix(get.adjacency(SN[[i]]$G))
    cascade_simulator(y,SN[[i]]$X,
                SN[[i]]$O,SN[[i]]$I,
                SN[[i]]$tau,
                seed = i)
    dif_time <- Sys.time() - t
    units(dif_time) <- "secs"
    times <- c(times,as.numeric(dif_time)) 
    cascades[[i]] <- read.csv(".tmp_results.txt")
  }
  sim_data <- list(networks = SN, cascades = cascades,
                   times = times, N = N,mode = mode,
                   cl = cl,type = type,shift = shift)
  
  class(sim_data) <- c("influenceDiffusion","list")
  return(sim_data)
}

cascade_reach <- function(cascade){
  IURS <- cascade[,5:8]  
  N <- sum(IURS[1,])
  SvR <- sum(IURS[nrow(IURS),3:4])
  return(SvR/N)
}

print.influenceDiffusion <- function(iD){
  cat("\n")
  cat("Graph order = ",iD$N,"\n",sep = "")
  cat("Clustering degree:",iD$cl,"\n")
  cat("Initiators selection method:",iD$mode,"\n")
  cat("I and O sampling method:",iD$type[1],"and",iD$type[2],"\n")
  cat("Run\tDensity\tReach\tTime\n",sep = "")
  for(i in 1:4){
    cat(i,"\t",
        round(edge_density(iD$networks[[i]]$G),4),"\t",
        round(cascade_reach(iD$cascades[[i]]),5),"\t",
        iD$times[i],"\n",sep = "")
  }
}

plot_cascade <- function(cascade,col = my.pal[2:5],...){
  IURS <- cascade[,5:8]  
  N <- sum(IURS[1,])
  matplot(IURS,type = "l",col= col,lty = 1,lwd = 2,log = "y",
          xlab = "Jumps",ylab = "Individuals",#main = "States flow",
          ylim = c(1,N*2),...)
  grid()
  legend("top",legend = c("Ignores","Undecided","Rejects","Supports"),
         lty = 1, col = col,lwd = 3,bg = "#ffffffb0",
         horiz = TRUE)
}

plot_jumps <- function(cascade,col = my.pal,...){
  plot(cascade[,4],log = "y",type = "l",pch = 20,
       main = "Jump Times",col = "#bbbbbb",
       xlab = "Jumps",ylab = "Jump length (t)",
       ...)
  points(cascade[,4],col = col[cascade[,3]+4],
         pch = 20)
  legend("bottom",legend = c("Information jump","Influence jump"),
         pch = 20,bg = "#ffffffb0",
         col = col[4:5],horiz = TRUE)
  grid()
}

plot.influenceDiffusion <- function(iD,pal = my.pal){
  
  cl   <- iD$cl
  mode <- iD$mode
  type <- iD$type
  N    <- iD$N
  file.name <- paste0("sim_","N_",N,"_cl_",round(cl,2),"_mode_",mode,"_I_",type[1],"_O_",type[2],"_")
  # States flow
  pdf(paste0(file.name,"cascade.pdf"),width = 14,height = 8)
  par(mfrow=c(2,2),bty = "l")   
  for(i in 1:4){
    plot_cascade(iD$cascades[[i]])
  }
  dev.off()
  
  # Jumps length
  upper_l <- max(sapply(iD$cascades,FUN = function(df){max(df[,4])}))
  lower_l <- min(sapply(iD$cascades,FUN = function(df){min(df[,4])}))
  pdf(paste0(file.name,"jumps.pdf"),width = 14,height = 8)
  par(mfrow=c(2,2),bty = "l")   
  for(i in 1:4){
    plot_jumps(iD$cascades[[i]],ylim = c(lower_l/2,upper_l))
  }
  dev.off()

  # Border size
  upper_l <- max(sapply(iD$cascades,FUN = function(df){max(df[,9])}))
  lower_l <- min(sapply(iD$cascades,FUN = function(df){min(df[,9])}))
  max_n_jumps <- max(sapply(iD$cascades,FUN = nrow))
  pdf(paste0(file.name,"border.pdf"),width = 7,height = 5)
  plot(iD$cascades[[1]][,9],type = "l",
       main = "Border size",
       xlab = "Jumps",ylab = "",ylim = c(lower_l,upper_l*1.3),
       lty = 1,xlim = c(1,max_n_jumps),col = pal[2])
  grid()
  for(i in 2:4){
    lines(iD$cascades[[i]][,9],lty = 1,col = pal[1+i])
  }
  legend("top",legend = 1:4,lty = 1, col = pal[2:5],
         title = "Run",bg = "white",horiz = TRUE)
  dev.off()
  
  # Assortativity
  upper_l <- max(sapply(iD$cascades,FUN = function(df){max(df[,10],na.rm = T)}))
  lower_l <- min(sapply(iD$cascades,FUN = function(df){min(df[,10],na.rm = T)}))
  max_n_jumps <- max(sapply(iD$cascades,FUN = nrow))
  pdf(paste0(file.name,"assortativity.pdf"),width = 7,height = 5)
  plot(iD$cascades[[1]][-nrow(iD$cascades[[1]]),10],type = "l",
       main = "Assortativity",
       xlab = "Jumps",ylab = "",ylim = c(lower_l,upper_l*1.3),
       lty = 1,xlim = c(1,max_n_jumps),col = pal[2])
  grid()
  for(i in 2:4){
    lines(iD$cascades[[i]][-nrow(iD$cascades[[i]]),10],lty = 1,col = pal[1+i])
  }
  legend("top",legend = 1:4,lty = 1, col = pal[2:5],
         title = "Run",bg = "white",horiz = TRUE)
  dev.off()

}

plot_network_cascade <- function(iD,run,steps = 4L,
                                 pal = my.pal[2:5],
                                 vertex.label = "",
                                 edge.arrow.size = 0.3,
                                 vertex.size = 4,
                                 ...){

  G <- iD$networks[[run]]$G
  X_init <- iD$networks[[run]]$X
  
  comp   <- components(G)
  gc_ids <- comp$membership == 1
  G      <- subgraph(G,gc_ids)
  
  if(is.integer(steps) && length(steps) == 1){
    n_steps <- steps
    n_jumps <- nrow(iD$cascades[[run]])
    delta <- log2(n_jumps)/(n_steps-1)
    steps <- floor(2^( (0:(n_steps-1))*delta ))
  }else{
    n_steps <- as.integer(length(steps))
  }
  if(n_steps==1){
    mfrow <- c(1,1)
  }else if(n_steps == 2){
    mfrow <- c(1,2)
  }else if(n_steps %in% 3:4){
    mfrow <- c(2,2)
  }else if(n_steps %in% 5:6){
    mfrow <- c(3,2)
  }else if(n_steps %in% 7:9){
    mfrow <- c(3,3)
  }else if(n_steps %in% 10:12){
    mfrow <- c(4,3)
  }else{
    stop("This function can only plot until 12 steps.")
  }
  
  pdf(paste0("sim_N_",iD$N,"_cl_",iD$cl,"_mode_",iD$mode,"_I_",iD$type[1],"_O_",iD$type[2],"_run_",run,"_network_cascade.pdf"),height = 7,width = 10)
  par(mfrow=mfrow,mar = c(0,0,1,0))
  for(st in steps){
    cat(st," ")
    X_t <- update_X_until(X_init,iD$cascades[[run]],st)
    set.seed(1)
    plot(G,vertex.color = pal[X_t[gc_ids]+1],main = st,
         vertex.label = vertex.label,
         vertex.size = vertex.size,
         edge.arrow.size = edge.arrow.size,...)
  }
  cat("\n")
  dev.off()

}

njumps <- function(iD){
  sapply(iD$cascades,FUN = nrow)
}


ttotal <- function(iD){
  sapply(iD$cascades,FUN = function(tb){
    sum(tb[,4])
  })
}


IURSprop <- function(iD,state = "I"){
  
  if(state == "I"){
    st <- 0
  }else if(state == "U"){
    st <- 1
  }else if(state == "R"){
    st <- 2
  }else if(state == "S"){
    st <- 3
  }else{
    stop("El estado ingresado no existe.")
  }
  
  sapply(iD$cascades,FUN = function(tb){
    njp <- nrow(tb)
    tot <- tb[njp,st+5]
    return(tot/iD$N)
  })
}

reach <- function(iD){
  Rprop <- IURSprop(iD,state = "R")
  Sprop <- IURSprop(iD,state = "S")
  maxReach <- ifelse(Rprop>Sprop,Rprop,Sprop)
  return(maxReach)
}



# Test ----

## tau ----

pdf("sims_kappa_taus_image.pdf",
    height = 6,width = 7)
N <- 200
cols <- designer.colors(64,c(my.pal[5],"white",my.pal[4]))
par(mfrow = c(2,2),mar = c(4,4,4,4),
    oma = c(1,1,0.5,1),cex = 0.7)
for(cl in c(0,0.5,0.7,0.9)){
  tmp_tau <- sample_tau(N,cl = cl)
  diag(tmp_tau) <- 1
  image.plot(1:N,1:N,tmp_tau,
             main = paste("k = ",cl),
             col = cols, zlim = c(-1,1),
             xlab = expression(i),
             ylab = expression(j))
}
dev.off()


cls   <- seq(0,1,by = 0.1)
means <- c()
cvs   <- c()
k <- 10
N <- 400
ids <- rep(1:2,c(round(N/2),N - round(N/2)))
for(cl in cls){
  cat(cl,"\n")
  mods <- c()
  for(i in 1:k){
    iN <- sample_network(N = N,seed = i,cl = cl)
    mods[i] <- modularity(iN$G,ids)
  }
  means <- c(means,mean(mods))
  cvs <- c(cvs,sd(mods)/abs(mean(mods)))
}

data.frame(cl = cls,AVGmod = means,CVmod = cvs) -> cl_to_mod
xtable::xtable(cl_to_mod,rownames = F)
write.csv(cl_to_mod,file = "cl_to_modularity.csv",
          row.names = F)

## shift ----

## Sampler ----
for(cl in seq(0,1,by = 0.1)){
  sn_tmp <- sample_network(N = 400,cl = cl,
                           type = c("realistic",
                                    "constant"))
  cat(cl,":",edge_density(sn_tmp$G),"\n")
}

{
  par(mfrow = c(2,2),mar = c(0,0,1,0))
  # ER - rand
  sn_tmp <- sample_network(N = 400,cl = 0.2,
                           mode = "random")
  plot(sn_tmp,main = "ER - rand")
  # ER - best
  sn_tmp <- sample_network(N = 400,cl = 0.2,
                           mode = "best")
  plot(sn_tmp,main = "ER - best")
  # CL - rand
  sn_tmp <- sample_network(N = 400,cl = 0.9,
                           mode = "random clust")
  plot(sn_tmp,main = "CL - rand")
  # CL - best
  sn_tmp <- sample_network(N = 400,cl = 0.9,
                           mode = "best clust")
  plot(sn_tmp,main = "CL - best")
}

## Cascades ----
tmp_sim <- run_simulation(400,cl = 0.2,mode = "random",
                          type = c("constant","realistic"))
plot(tmp_sim$networks[[4]])
tmp_sim

plot_cascade(tmp_sim$cascades[[1]])
plot_jumps(tmp_sim$cascades[[1]])

plot(tmp_sim)
njumps(tmp_sim)
plot_network_cascade(tmp_sim,run = 3,steps = c(1,100,1000,1500,2360))


sample_network(400,cl = 0.2,mode = "random",
               type = c("constant","constant")) -> iN

# Simulations ----

types <- c("realistic","constant")
cls <- c(0,0.9)
#N <- 100 #144 seg
#N <- 200 #847 seg
#N <- 300 #1365 seg
N <- 1000
modes <- c("random","best","random clust","best clust")


cpp_optimize(function(shift){
  shift_to_meanDeg(shift,N = N, mode = "random",
                   cl = 0, type = c("constant","constant"),
                   k = 7)
},10) -> tmp


{
  t <- Sys.time()
  n_sim <- 0
  for(I_ty in types){
    for(O_ty in types){
      for(cl in cls){
        if(cl < 0.5){
          for(mod in modes[1:2]){
            if(O_ty == "constant" && mod == "best") next
            tmp_sim <- run_simulation(N,cl = cl,
                                      mode = mod,
                                      type = c(I_ty,O_ty))
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            save(tmp_sim,file = file.name)
            #
            n_sim <- n_sim + 1
            cat(round(n_sim*100/12,2),"% completado.\n")
          }
        }else{
          for(mod in modes[3:4]){
            if(O_ty == "constant" && mod == "best clust") next
            tmp_sim <- run_simulation(N,cl = cl,
                                      mode = mod,
                                      type = c(I_ty,O_ty))
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            save(tmp_sim,file = file.name)
            #
            n_sim <- n_sim + 1
            cat(round(n_sim*100/12,2),"% completado.\n")
          }
        }
      }
    }
  }
  Sys.time() - t
}

load("chain_N_1000_cl_0_mode_random_I_constant_O_realistic.RData")

plot(tmp_sim)

# Plots ----

{
  t <- Sys.time()
  n_sim <- 0
  for(I_ty in types){
    for(O_ty in types){
      for(cl in cls){
        if(cl < 0.5){
          for(mod in modes[1:2]){
            if(O_ty == "constant" && mod == "best") next
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            load(file = file.name)
            plot(tmp_sim)
            for(run in 1:4){
              plot_network_cascade(tmp_sim,run,steps = 9L,
                                   edge.width = 0.5,
                                   edge.color = "#a0a0a080",
                                   edge.arrow.size = 0.2,
                                   layout = layout_with_kk)
            }
            #
            n_sim <- n_sim + 1
            cat(round(n_sim*100/12,2),"% completado.\n")
          }
        }else{
          for(mod in modes[3:4]){
            if(O_ty == "constant" && mod == "best clust") next
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            load(file = file.name)
            plot(tmp_sim)
            for(run in 1:4){
              plot_network_cascade(tmp_sim,run,steps = 9L,
                                   edge.width = 0.5,
                                   edge.color = "#a0a0a080",
                                   edge.arrow.size = 0.2,
                                   layout = layout_with_mds)
            }
            #
            n_sim <- n_sim + 1
            cat(round(n_sim*100/12,2),"% completado.\n")
          }
        }
      }
    }
  }
  Sys.time() - t
}

## Specific plot ----
N <- 1000;cl <- 0.9;mod <- "random clust";I <- "constant";O <- "constant"
file.name <- paste0("chain_N_",N,"_cl_",
                    round(cl,2),"_mode_",
                    mod,"_I_",I,
                    "_O_",O,".RData")
load(file = file.name)
plot_network_cascade(tmp_sim,3,steps = 9L,
                     edge.width = 0.5,
                     vertex.size = 3,
                     edge.color = "#a0a0a010",
                     edge.arrow.size = 0.2,
                     layout = layout_with_mds)


# Tabla Resultados ----

types_L <- list(realistic = "No constante",
                constant = "Constante")
mod_L <- list("random" = "MAS",
              "best" = "Mejor",
              "random clust" = "MAS estratificado",
              "best clust" = "Mejor por grupo")

res_tab <- NULL
full_res_tab <- NULL
{
  t <- Sys.time()
  n_sim <- 0
  for(O_ty in types){
    for(I_ty in types){
      for(cl in cls){
        if(cl < 0.5){
          for(mod in modes[1:2]){
            if(O_ty == "constant" && mod == "best") next
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            load(file = file.name)
            res <- c(
              types_L[O_ty],
              types_L[I_ty],
              cl,
              mod_L[mod],
              max(njumps(tmp_sim)),
              round(mean(ttotal(tmp_sim)),2),
              round(mean(reach(tmp_sim)),2)
            )
            res_tab <- rbind(res_tab,res)
            for(ii in 1:4){
              full_res <- c(
                types_L[O_ty],
                types_L[I_ty],
                cl,
                mod_L[mod],
                njumps(tmp_sim)[ii],
                round(ttotal(tmp_sim)[ii],2),
                round(reach(tmp_sim)[ii],2)
              )
              full_res_tab <- rbind(full_res_tab,full_res)
            }
            #
          }
        }else{
          for(mod in modes[3:4]){
            if(O_ty == "constant" && mod == "best clust") next
            file.name <- paste0("chain_N_",N,"_cl_",
                                round(cl,2),"_mode_",
                                mod,"_I_",I_ty,
                                "_O_",O_ty,".RData")
            load(file = file.name)
            res <- c(
              types_L[O_ty],
              types_L[I_ty],
              cl,
              mod_L[mod],
              max(njumps(tmp_sim)),
              round(mean(ttotal(tmp_sim)),2),
              round(mean(reach(tmp_sim)),2)
            )
            res_tab <- rbind(res_tab,res)
            for(ii in 1:4){
              full_res <- c(
                types_L[O_ty],
                types_L[I_ty],
                cl,
                mod_L[mod],
                njumps(tmp_sim)[ii],
                round(ttotal(tmp_sim)[ii],2),
                round(reach(tmp_sim)[ii],2)
              )
              full_res_tab <- rbind(full_res_tab,full_res)
            }
            #
          }
        }
      }
    }
  }
  Sys.time() - t
}

colnames(res_tab) <- c("O","I","k","Muestreo",
                       "N° Saltos","Tiempo Total",
                       "Alcance")
rownames(res_tab) <- NULL

colnames(full_res_tab) <- c("O","I","k","Muestreo",
                            "N° Saltos","Tiempo Total",
                            "Alcance")
rownames(full_res_tab) <- NULL

View(res_tab)
View(full_res_tab)

xtable::xtable(res_tab)

res_df <- data.frame(
  O = Reduce('c',res_tab[,"O"]),
  I = Reduce('c',res_tab[,"I"]),
  kap = factor(Reduce('c',res_tab[,"k"])),
  Muestreo = Reduce('c',res_tab[,"Muestreo"]),
  totTime = Reduce('c',res_tab[,"Tiempo Total"]),
  Alcance = Reduce('c',res_tab[,"Alcance"])
)

full_res_df <- data.frame(
  O = factor(Reduce('c',full_res_tab[,"O"])),
  I = factor(Reduce('c',full_res_tab[,"I"])),
  kap = factor(Reduce('c',full_res_tab[,"k"])),
  Muestreo = factor(Reduce('c',full_res_tab[,"Muestreo"])),
  totTime = Reduce('c',full_res_tab[,"Tiempo Total"]),
  Alcance = Reduce('c',full_res_tab[,"Alcance"])
)

# ANOVA ----

## Modelos iniciales ----
model_reach <- aov(I(log(Alcance)) ~ O + I + kap + kap:Muestreo, 
                   data = full_res_df)
anova(model_reach)
{
  par(mfrow=c(2,2))
  plot(model_reach)
}
{
  par(mfrow = c(1,1))
  car::qqPlot(residuals(model_reach))
}

model_time <- aov(I(log(totTime)) ~ O + I + kap + kap:Muestreo, 
                  data = full_res_df)
anova(model_time)
{
  par(mfrow=c(2,2))
  plot(model_time)
}
{
  par(mfrow = c(1,1))
  car::qqPlot(residuals(model_time))
}


## Quitar outliers ----
model_reach <- aov(I(log(Alcance)) ~ O + I + kap + kap:Muestreo, 
                   data = full_res_df,
                   subset = -c(30))
anova(model_reach)
{
  par(mfrow=c(2,2))
  plot(model_reach)
}
{
  par(mfrow = c(1,1))
  car::qqPlot(residuals(model_reach),
              envelope=.99)
}
shapiro.test(residuals(model_reach))
tseries::jarque.bera.test(residuals(model_reach))

xtable::xtable(anova(model_reach))


model_time <- aov(I(log(totTime)) ~ O + I + kap + kap:Muestreo, 
                  data = full_res_df,
                  subset = -c(30))
anova(model_time)
{
  par(mfrow=c(2,2))
  plot(model_time)
}
{
  par(mfrow = c(1,1))
  car::qqPlot(residuals(model_time),
              envelope=.99)
}
shapiro.test(residuals(model_time))
tseries::jarque.bera.test(residuals(model_time))

xtable::xtable(anova(model_time))

{
  pdf("sims_qqplot_anovas.pdf",
      height = 5, width = 10)
  par(mfrow=c(1,2),
      mar = c(5,4,2,1))
  car::qqPlot(residuals(model_time),
              xlab = "Cuantiles teóricos",
              ylab = "Cuantiles empíricos")
  car::qqPlot(residuals(model_reach),
              xlab = "Cuantiles teóricos",
              ylab = "Cuantiles empíricos")
  dev.off()
}

## Homocedasticidad ----


library(dplyr)

full_res_df %>% tibble %>% 
  mutate(Tratamiento = factor(paste(O,I,kap,Muestreo,sep = "-"))) %>% 
  select(Tratamiento,totTime,Alcance) -> full_res_df_trat

car::leveneTest(I(log(totTime)) ~ Tratamiento,
                data = full_res_df_trat,
                subset = -c(30))

car::leveneTest(I(log(Alcance)) ~ Tratamiento,
                data = full_res_df_trat,
                subset = -c(30))

