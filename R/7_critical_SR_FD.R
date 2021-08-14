## funcao para calcular o numero de comunidades acima e abaixo de um limiar critico de riqueza

load("all_results_tb.RData")

#   100 traps

lapply (seq(1,length(pools)), function (i) {
## considerando somente deltaSR, por enquanto
 habitatSR_100 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_100 >=4)] [which (riqueza_100_4 >=5)])
 levels (habitatSR_100) [which (levels (habitatSR_100) == "forest.0.0.0")] <- "Natural"
 levels (habitatSR_100) [which (levels (habitatSR_100) == "grassland.0.0.0")] <- "Natural"
 levels (habitatSR_100) [which (levels (habitatSR_100) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaSR <- seq (range(deltaSR_100)[1], range(deltaSR_100)[2], 0.5)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit <- lapply (crit_deltaSR, function (limiar) {
 
   (table (habitatSR_100 [which (deltaSR_100 <  limiar)]))/length(deltaSR_100)
   
   }
)

# desmanchar a lista = 1 valor para cada limiar
delta_SR_crit <- do.call (rbind , delta_SR_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_100 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_100 >=4)] [which (riqueza_100_4 >=5)])
#levels (habitatSR_100) [which (levels (habitatSR_100) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_100) [which (levels (habitatSR_100) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_100) [which (levels (habitatSR_100) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit_cada_habitat <- lapply (crit_deltaSR, function (limiar) {
 
  (table (cada_habitatSR_100 [which (deltaSR_100 <  limiar)]) )/length(deltaSR_100)
 
   }
)

# desmanchar a lista
delta_SR_crit_cada_habitat <- do.call (rbind, delta_SR_crit_cada_habitat)
# remover ultima coluna
delta_SR_crit_cada_habitat <- delta_SR_crit_cada_habitat [,-which(colnames(delta_SR_crit_cada_habitat) == "edge.0.0.1")]
delta_SR_crit_cada_habitat <- as.data.frame(delta_SR_crit_cada_habitat)


png(file=paste ("criticalVal_",pools[i],".png"), width=22, height=12, units="cm", res=600, family="serif")
par(mfrow=c(1,2))
# add plot
plot(NA, xlim=c(range(crit_deltaSR)[1], range(crit_deltaSR) [2]),ylim = c (0,1),
    xlab = expression (paste ("Critical values of ",Delta,"SR",sep="")),
    ylab="Proportion of communities",
    main=expression (paste (Delta,"SR",sep="")))

# natural
points (crit_deltaSR, delta_SR_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaSR, delta_SR_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaSR, delta_SR_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaSR, delta_SR_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaSR, delta_SR_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaSR, delta_SR_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaSR, delta_SR_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_100),lwd=2,lty=2)

### considerando  SES delta FD

## considerando somente deltaSR, por enquanto

habitatSR_100 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_100 >=4)] [which (riqueza_100_4 >=5)])
levels (habitatSR_100) [which (levels (habitatSR_100) == "forest.0.0.0")] <- "Natural"
levels (habitatSR_100) [which (levels (habitatSR_100) == "grassland.0.0.0")] <- "Natural"
levels (habitatSR_100) [which (levels (habitatSR_100) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaFD <- seq (range(deltaFD_100)[1], range(deltaFD_100)[2], 0.05)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit <- lapply (crit_deltaFD, function (limiar) {
 
 (table (habitatSR_100 [which (deltaFD_100 <  limiar)]))/length(deltaFD_100)
 
}
)

# desmanchar a lista = 1 valor para cada limiar
delta_FD_crit <- do.call (rbind , delta_FD_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_100 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_100 >=4)] [which (riqueza_100_4 >=5)])
#levels (habitatSR_100) [which (levels (habitatSR_100) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_100) [which (levels (habitatSR_100) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_100) [which (levels (habitatSR_100) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit_cada_habitat <- lapply (crit_deltaFD, function (limiar) {
 
 table (cada_habitatSR_100 [which (deltaFD_100 <  limiar)])/length(deltaFD_100)
 
}
)

# desmanchar a lista
delta_FD_crit_cada_habitat <- do.call (rbind, delta_FD_crit_cada_habitat)
# remover ultima coluna
delta_FD_crit_cada_habitat <- delta_FD_crit_cada_habitat [,-which(colnames(delta_FD_crit_cada_habitat) == "edge.0.0.1")]

# add plot
plot(NA, xlim=c(range(crit_deltaFD)[1], range(crit_deltaFD) [2]),ylim = c (0,1),
    xlab = expression (paste ("Critical values of SES ",Delta,"FD",sep="")),
    ylab="",
    main=expression (paste ("SES ",Delta,"FD",sep="")))

# natural
points (crit_deltaFD, delta_FD_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaFD, delta_FD_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaFD, delta_FD_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaFD, delta_FD_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaFD, delta_FD_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaFD, delta_FD_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaFD, delta_FD_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_100),lwd=2,lty=2)

legend("topleft", legend= c(
 "Major types",
 "Natural", "Human-modified", 
 "",
 "Categories",
 "Forest","Grassland",
 "Crop field","Clear-cut",
 "Grassland edge", "Forest edge","Tree plantation"
 
), 
pch =  c (NA,19,21,NA,NA,19,19,21,21,21,21,21
),
col=c(NA,"black","black",NA,NA,
     rgb(red = 0, green = 1, blue = 0, alpha = 0.5),
     rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5),
     rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5),
     rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5),
     rgb(red = 1, green = 0, blue = 0, alpha = 0.5),
     rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)
),
bty="n",cex=0.7)

dev.off()

}
)

##################################################################
##################################################################
################       500 traps     #############################
##################################################################
##################################################################

lapply (seq(1,length(pools)), function (i) {
  
## considerando somente deltaSR, por enquanto
habitatSR_500 <- (tabela_resultados_500$habitat [which(tabela_resultados_500$pool== pools [i])] [which (riqueza_500 >=4)] [which (riqueza_500_4 >=5)])
levels (habitatSR_500) [which (levels (habitatSR_500) == "forest.0.0.0")] <- "Natural"
levels (habitatSR_500) [which (levels (habitatSR_500) == "grassland.0.0.0")] <- "Natural"
levels (habitatSR_500) [which (levels (habitatSR_500) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaSR <- seq (range(deltaSR_500)[1], range(deltaSR_500)[2], 0.5)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit <- lapply (crit_deltaSR, function (limiar) {
  
  (table (habitatSR_500 [which (deltaSR_500 <  limiar)]))/length(deltaSR_500)
  
}
)

# desmanchar a lista = 1 valor para cada limiar
delta_SR_crit <- do.call (rbind , delta_SR_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_500 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_500 >=4)] [which (riqueza_500_4 >=5)])
#levels (habitatSR_500) [which (levels (habitatSR_500) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_500) [which (levels (habitatSR_500) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_500) [which (levels (habitatSR_500) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit_cada_habitat <- lapply (crit_deltaSR, function (limiar) {
  
  (table (cada_habitatSR_500 [which (deltaSR_500 <  limiar)]) )/length(deltaSR_500)
  
}
)

# desmanchar a lista
delta_SR_crit_cada_habitat <- do.call (rbind, delta_SR_crit_cada_habitat)
# remover ultima coluna
delta_SR_crit_cada_habitat <- delta_SR_crit_cada_habitat [,-which(colnames(delta_SR_crit_cada_habitat) == "edge.0.0.1")]
delta_SR_crit_cada_habitat <- as.data.frame(delta_SR_crit_cada_habitat)


png(file=paste ("criticalVal500_",pools[i],".png"), width=22, height=12, units="cm", res=600, family="serif")
par(mfrow=c(1,2))
# add plot
plot(NA, xlim=c(range(crit_deltaSR)[1], range(crit_deltaSR) [2]),ylim = c (0,1),
     xlab = expression (paste ("Critical values of ",Delta,"SR",sep="")),
     ylab="Proportion of communities",
     main=expression (paste (Delta,"SR",sep="")))

# natural
points (crit_deltaSR, delta_SR_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaSR, delta_SR_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaSR, delta_SR_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaSR, delta_SR_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaSR, delta_SR_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaSR, delta_SR_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaSR, delta_SR_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_500),lwd=2,lty=2)

### considerando  SES delta FD

## considerando somente deltaSR, por enquanto

habitatSR_500 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools[i])] [which (riqueza_500 >=4)] [which (riqueza_500_4 >=5)])
levels (habitatSR_500) [which (levels (habitatSR_500) == "forest.0.0.0")] <- "Natural"
levels (habitatSR_500) [which (levels (habitatSR_500) == "grassland.0.0.0")] <- "Natural"
levels (habitatSR_500) [which (levels (habitatSR_500) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaFD <- seq (range(deltaFD_500)[1], range(deltaFD_500)[2], 0.05)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit <- lapply (crit_deltaFD, function (limiar) {
  
  (table (habitatSR_500 [which (deltaFD_500 <  limiar)]))/length(deltaFD_500)
  
}
)

# desmanchar a lista = 1 valor para cada limiar
delta_FD_crit <- do.call (rbind , delta_FD_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_500 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_500 >=4)] [which (riqueza_500_4 >=5)])
#levels (habitatSR_500) [which (levels (habitatSR_500) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_500) [which (levels (habitatSR_500) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_500) [which (levels (habitatSR_500) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit_cada_habitat <- lapply (crit_deltaFD, function (limiar) {
  
  table (cada_habitatSR_500 [which (deltaFD_500 <  limiar)])/length(deltaFD_500)
  
}
)

# desmanchar a lista
delta_FD_crit_cada_habitat <- do.call (rbind, delta_FD_crit_cada_habitat)
# remover ultima coluna
delta_FD_crit_cada_habitat <- delta_FD_crit_cada_habitat [,-which(colnames(delta_FD_crit_cada_habitat) == "edge.0.0.1")]

# add plot
plot(NA, xlim=c(range(crit_deltaFD)[1], range(crit_deltaFD) [2]),ylim = c (0,1),
     xlab = expression (paste ("Critical values of SES ",Delta,"FD",sep="")),
     ylab="",
     main=expression (paste ("SES ",Delta,"FD",sep="")))

# natural
points (crit_deltaFD, delta_FD_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaFD, delta_FD_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaFD, delta_FD_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaFD, delta_FD_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaFD, delta_FD_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaFD, delta_FD_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaFD, delta_FD_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_500),lwd=2,lty=2)

legend("topleft", legend= c(
  "Major types",
  "Natural", "Human-modified", 
  "",
  "Categories",
  "Forest","Grassland",
  "Crop field","Clear-cut",
  "Grassland edge", "Forest edge","Tree plantation"
  
), 
pch =  c (NA,19,21,NA,NA,19,19,21,21,21,21,21
),
col=c(NA,"black","black",NA,NA,
      rgb(red = 0, green = 1, blue = 0, alpha = 0.5),
      rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5),
      rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5),
      rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
      rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5),
      rgb(red = 1, green = 0, blue = 0, alpha = 0.5),
      rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)
),
bty="n",cex=0.7)

dev.off()

}
)

##################################################################
##################################################################
################       1,000 traps     #############################
##################################################################
##################################################################

lapply (seq(1,length(pools)), function (i) {
## considerando somente deltaSR, por enquanto
habitatSR_1000 <- (tabela_resultados_1000$habitat [which(tabela_resultados_1000$pool== pools [i])] [which (riqueza_1000 >=4)] [which (riqueza_1000_4 >=5)])
levels (habitatSR_1000) [which (levels (habitatSR_1000) == "forest.0.0.0")] <- "Natural"
levels (habitatSR_1000) [which (levels (habitatSR_1000) == "grassland.0.0.0")] <- "Natural"
levels (habitatSR_1000) [which (levels (habitatSR_1000) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaSR <- seq (range(deltaSR_1000)[1], range(deltaSR_1000)[2], 0.5)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit <- lapply (crit_deltaSR, function (limiar) {
  
  (table (habitatSR_1000 [which (deltaSR_1000 <  limiar)]))/length (deltaSR_1000)
  
}
)

# desmanchar a lista = 1 valor para cada limiar
delta_SR_crit <- do.call (rbind , delta_SR_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_1000 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_1000 >=4)] [which (riqueza_1000_4 >=5)])
#levels (habitatSR_1000) [which (levels (habitatSR_1000) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_1000) [which (levels (habitatSR_1000) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_1000) [which (levels (habitatSR_1000) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_SR_crit_cada_habitat <- lapply (crit_deltaSR, function (limiar) {
  
  (table (cada_habitatSR_1000 [which (deltaSR_1000 <  limiar)]) )/length (deltaSR_1000)
  
}
)

# desmanchar a lista
delta_SR_crit_cada_habitat <- do.call (rbind, delta_SR_crit_cada_habitat)
# remover ultima coluna
delta_SR_crit_cada_habitat <- delta_SR_crit_cada_habitat [,-which(colnames(delta_SR_crit_cada_habitat) == "edge.0.0.1")]
delta_SR_crit_cada_habitat <- as.data.frame(delta_SR_crit_cada_habitat)


png(file=paste ("criticalVal1000_",pools[i],".png"), width=22, height=12, units="cm", res=600, family="serif")
par(mfrow=c(1,2))
# add plot
plot(NA, xlim=c(range(crit_deltaSR)[1], range(crit_deltaSR) [2]),ylim = c (0,1),
     xlab = expression (paste ("Critical values of ",Delta,"SR",sep="")),
     ylab="Proportion of communities",
     main=expression (paste (Delta,"SR",sep="")))

# natural
points (crit_deltaSR, delta_SR_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaSR, delta_SR_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaSR, delta_SR_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaSR, delta_SR_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaSR, delta_SR_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaSR, delta_SR_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaSR, delta_SR_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaSR, delta_SR_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_1000),lwd=2,lty=2)

### considerando  SES delta FD

## considerando somente deltaSR, por enquanto

habitatSR_1000 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools[i])] [which (riqueza_1000 >=4)] [which (riqueza_1000_4 >=5)])
levels (habitatSR_1000) [which (levels (habitatSR_1000) == "forest.0.0.0")] <- "Natural"
levels (habitatSR_1000) [which (levels (habitatSR_1000) == "grassland.0.0.0")] <- "Natural"
levels (habitatSR_1000) [which (levels (habitatSR_1000) != "Natural")] <- "Human-modified"

## valores criticos de deltaSR
crit_deltaFD <- seq (range(deltaFD_1000)[1], range(deltaFD_1000)[2], 0.05)

## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit <- lapply (crit_deltaFD, function (limiar) {
  
  (table (habitatSR_1000 [which (deltaFD_1000 <  limiar)]))/length (deltaFD_1000)
  
}
)

# desmanchar a lista = 1 valor para cada limiar
delta_FD_crit <- do.call (rbind , delta_FD_crit)

## mesma logica, so que agora considerando cada habitat em particular
cada_habitatSR_1000 <- (tabela_resultados$habitat [which(tabela_resultados$pool== pools [i])] [which (riqueza_1000 >=4)] [which (riqueza_1000_4 >=5)])
#levels (habitatSR_1000) [which (levels (habitatSR_1000) == "forest.0.0.0")] <- "Natural"
#levels (habitatSR_1000) [which (levels (habitatSR_1000) == "grassland.0.0.0")] <- "Natural"
#levels (habitatSR_1000) [which (levels (habitatSR_1000) != "Natural")] <- "Human-modified"


## obter o numero de comunidades de natural e modificado, de acordo com cada limiar
delta_FD_crit_cada_habitat <- lapply (crit_deltaFD, function (limiar) {
  
  table (cada_habitatSR_1000 [which (deltaFD_1000 <  limiar)]) /length (deltaFD_1000)
  
}
)

# desmanchar a lista
delta_FD_crit_cada_habitat <- do.call (rbind, delta_FD_crit_cada_habitat)
# remover ultima coluna
delta_FD_crit_cada_habitat <- delta_FD_crit_cada_habitat [,-which(colnames(delta_FD_crit_cada_habitat) == "edge.0.0.1")]

# add plot
plot(NA, xlim=c(range(crit_deltaFD)[1], range(crit_deltaFD) [2]),ylim = c (0,1),
     xlab = expression (paste ("Critical values of SES ",Delta,"FD",sep="")),
     ylab="",
     main=expression (paste ("SES ",Delta,"FD",sep="")))

# natural
points (crit_deltaFD, delta_FD_crit_cada_habitat[,2],pch=19,col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) # forest
points (crit_deltaFD, delta_FD_crit_cada_habitat[,3],pch=19,col=rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5)) # grassland

## human-modified
points (crit_deltaFD, delta_FD_crit_cada_habitat[,4],pch=21,col=rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5)) # crop field
points (crit_deltaFD, delta_FD_crit_cada_habitat[,6],pch=21,col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # clear-cut
points (crit_deltaFD, delta_FD_crit_cada_habitat[,1],pch=21,col=rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5)) # grassland edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,7],pch=21,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5)) # forest edge
points (crit_deltaFD, delta_FD_crit_cada_habitat[,5],pch=21,col=rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)) # tree plantation

## botar os pontos para o grupo maior de habitats
points (crit_deltaFD, delta_FD_crit[,1],pch=21,col="black",cex=1.2)
points (crit_deltaFD, delta_FD_crit[,2],pch=19,col="black",cex=1.2)
#abline (v=mean(deltaSR_1000),lwd=2,lty=2)

legend("topleft", legend= c(
  "Major types",
  "Natural", "Human-modified", 
  "",
  "Categories",
  "Forest","Grassland",
  "Crop field","Clear-cut",
  "Grassland edge", "Forest edge","Tree plantation"
  
), 
pch =  c (NA,19,21,NA,NA,19,19,21,21,21,21,21
),
col=c(NA,"black","black",NA,NA,
      rgb(red = 0, green = 1, blue = 0, alpha = 0.5),
      rgb(red = 0.8, green = 0.6, blue = 0.15, alpha = 0.5),
      rgb(red = 0.5, green = 0.6, blue = 0.6, alpha = 0.5),
      rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
      rgb(red = 0.7, green = 0.5, blue = 0.7, alpha = 0.5),
      rgb(red = 1, green = 0, blue = 0, alpha = 0.5),
      rgb(red = 0, green = 1, blue = 0.3, alpha = 0.5)
),
bty="n",cex=0.7)

dev.off()

}
)

