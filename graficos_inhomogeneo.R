#######################################################
#######################################################
#####                                             #####
#####       GRAFICOS MODELO INHOMOGENEO           #####
#####                                             #####
#######################################################
#######################################################


# Se cargan los datos
load("datos_TFG_PMRS.Rdata")

# Para leer las funciones
source("modelo_inhomogeneo.R")


# Se cargan los paquetes ecespa y spatstat
library(ecespa)
library (spatstat)

# Selecciona parcela

P1 <- plantas2

# se definen simulaciones (para test de asociaciones)
nsim <- 19


#-------
# 1) Se ajusta inghomogeneidad y saca funciones Kcross inhomogenas y envueltas


# Vector de sigmas a probar
# OJO sigmas en las unidades de la parcela. Las unidades son metros
sigmas <- seq(30, 270, by=10)*0.01
sigmas


prueba <- all_cross_auto2(P1, nmin=15, sigmas=sigmas, nsim.sigmas=9, nsim.test=nsim, plotenv=FALSE)
names(prueba)
# alternativamente se pueden probar varias sigmas nm?s, hasta conseguir la reweighted stationarity



# par(mfrow=c(9,9)); par(mar=c(1,1,1,1)); for(i in 1:9) for(j in 1:9) plot(prueba$wIJ.test[[i]][[j]], sqrt(./pi)-r~r, legend=F, main="")

#-------
#2) Se saca promedio ponderado de envueltas ij y ji para cada par de especies
sp.ok <- (names(prueba$wIJ.test))
lista.test <- prueba$wIJ.test # (obtenido en  paso 1)

# numero de individuos por especies (necesario para la media ponderada)
nsp <- table(P1$marks)[match(sp.ok, names(table(P1$marks)))]

# envueltas ponderadas
ejemplo.wEnvelopelist <- wEnvelopelist(lista.test, sp.ok, nsp)

# representacion
dev.new()
lista.env <- ejemplo.wEnvelopelist$wenvelopelist   
len <- length(lista.env)
setpanels(len^2)
par(mar=c(0,0,0,0)); 
for(i in 1:len) for(j in 1:len) plot(lista.env[[i]][[j]], sqrt(./pi)-r~r, legend=F, main="", xaxt="n",yaxt="n")




#-------
#3)  Grafica porcentajes de asociaciones

# lista conj las envueltas ponderadas 
lista.env <- ejemplo.wEnvelopelist$wenvelopelist
# Figura 3 (porcentaje de asiociaciones significativas/no significativas) en funcion de r
figura <- prop.asoc(lista.env, alpha=0.05)  # alpha es el umbral critico del testpointwise

figura


#-------
#4)  En funcion del grafico anterior, elije radio(s) y calcula la(s) matriz (es) de asociacion

R <- 0.25

Aij <- as.dist(sacaA(lista.env, R=R)) # OJO: puede hacerse con las funciones centradas (lista.ev.centradas)

round( Aij,2)

# se representa histograma de porcentajes de valores de asociacion
perchist <- function(x,main=NULL,...){
  if(is.null(main)) title<- deparse(substitute(x))
  cosa <- hist(x, plot=FALSE)
  cosa$counts <- cosa$counts*100/sum(cosa$counts)
  plot(cosa, ..., main=title)
}
dev.new()
perchist(Aij, 
         main = "Lij", 
         xlab = "Fuerza de asociación espacial (Lij) para r=0,25m", 
         ylab = "Porcentaje pares de especies", 
         col = "#4CA6FF")
