# Gráficos de la función de distribución Gamma
setwd("/home/samuel/Documentos/U/tesis/")

## Oi


pdf("plot_gamma_Oi_Ii.pdf",height = 5, width = 8)
par(bty = "l",mfrow=c(1,2),mar = c(5,4,2,2))

shape <- 3
rate <- 3/4
shift <- 30
my.grid <- seq(-shift,-10,length = 300)
n.grid <- length(my.grid)
my.d <- dgamma(my.grid+shift,shape,rate)

curve(dgamma(x+shift,shape,rate),-shift,-10,
      #main = expression(paste("Función de densidad para ",O[i])),
      xlab = expression(O[i]),ylab = "Densidad",
      ylim = c(0,0.4))
grid()
polygon(c(my.grid,my.grid[n.grid:1]),c(my.d,rep(0,n.grid)),
        col = "#4090aab0")

#--

shape <- 21/5
rate <- 7/5
shift <- 0
my.grid <- seq(-shift,8,length = 300)
n.grid <- length(my.grid)
my.d <- dgamma(my.grid+shift,shape,rate)

curve(dgamma(x+shift,shape,rate),-shift,8,
      #main = expression(paste("Función de densidad para ",I[i])),
      xlab = expression(I[i]),ylab = "Densidad",
      ylim = c(0,0.4))
grid()
polygon(c(my.grid,my.grid[n.grid:1]),c(my.d,rep(0,n.grid)),
        col = "#b02070b0")


dev.off()


