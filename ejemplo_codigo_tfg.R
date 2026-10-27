#######################################################
#######################################################
#####                                             #####
#####      SEMEJANZA FUNCIONAL Y ASOCIACIÓN       #####
#####      ESPACIAL EN UN MATORRAL GIPSÓFILO      #####
#####                                             #####
#######################################################
#######################################################


##### 1.) Automatización del cáculo de envueltas cruzadas  #####
################################################################

# Función para simular a partir de un modelo nulo

simu2 <- function(pp2, i, lambdai){
  u <- split(pp2)
  u[[i]] <- rpoispp(lambdai)[u[[i]]$window ]
  split(pp2) <- u
  return(pp2)
}

# Se genera tabla con abundancias de cada especie ordenadas de más a menos

abundancias <- table(plantas2$marks)[rev(order(table(plantas2$marks)))]
abundancias

# Nombres de las especies con más de 50 individuos (van a ser en total 19 especies a considerar)

especies <- names(abundancias)[abundancias >=50]
especies

# Número de especies

n <- length(especies)
n

# Se prepara una lista vacía para almacenar los test cruzados(envueltas) para cada especie

envueltas <- vector (mode="list", length=n)
names(envueltas) <- especies

# se definen el número de simulaciones (elegimos 99, que son las que hace por defecto, pero se pueden aumentar o disminuir con esta definición

nsim <- 99

# Se generan varios bucles anidados:
# a) se instala y se carga el paquete necesario para la estimación de de la densidad espacial de cada especie
# a) primero se selecciona la especie focal, y luego se cruzan todas las demás con ella
# b) cuando se acaba con una especie, comenzamos con la siguiente
# c) se repite con el resto de especies hasta la última

install.packages("spatstat")
library(spatstat)

for(spfocal in especies){ # se elige eligiendo especie a especie:
  # "spfocal" es cada uno de los nombres de las especies
  # se ajusta el modelo nulo para la especie focal
  
  spfocal.lambda <- density.ppp(plantas2[plantas2$marks== spfocal])
  
  # se genera una sublista para cada especie (se almacenan las envueltas cruzadas para cada especie focal
  
  envueltas [[spfocal]] <- vector (mode="list", length=n)
  names(envueltas [[spfocal]] ) <- especies
  
  # comienza el bucle para calcular las envueltas de la función pcf cruzada
  # entre la especie focal y cada una de las otras especies (se van a llamar "sp")
  
  for (sp in especies){
    
  # se pone un condicional para que sólo haga envueltas si "sp" está entre las especies diferente a "spfocal"
    if(sp != spfocal){
      
    # se guardan el resultado de las envueltas en la posición correspondiente
    # (la especie "sp") de la sublista
      
      envueltas [[spfocal]] [[sp]] <- envelope (plantas2, pcfcross, i=spfocal, j=sp,
      simulate=expression(simu2(pp2=plantas2, i=spfocal, lambdai=spfocal.lambda)),
      nsim=nsim, savefuns=TRUE)
   }# acaba el condicional
  }# acaba el bucle para cada especie focal
}# acaba el bucle total
# se guardan los resultados en un workspcace para no tener que volverlos a calcularlos

save(envueltas, file="envueltas.Rdata")


##### 2.) Cálculo matriz Aij                               #####
################################################################

milog <- function(x) ifelse(x>0, log(x),0)

# Se calcula Aij(r=0,2 m) y Aij(r=2 m):
##  a) Se definen los valores de r para los que hacer el sumatorio (integral) en ambos casos.
## r es la misma en todas las envueltas, así que el cálculo lo hacemos a partir de una envuelta cualquiera.

rok2 <-envueltas [[1]] [[2]]$r <= 0.2
rok1 <-envueltas [[1]] [[2]]$r <= 2

## b) Define matriz Aij R=0,2 Y R= 2 (en principio las llenamos del valor NA.)
## número de especies

n <- length(especies)
Aij1 <- Aij2 <- matrix (NA, nrow=n, ncol=n)

## c) se ponen nombres a las filas y columnas de ambas matrices (los mismos nombres en ambas matrices)

dimnames(Aij2) <- dimnames(Aij1) <- list(especies, especies)

## c) se genera el bucle para el cálculo de Aij para r=0,2 y r=2

for(spfocal in especies){
  for (sp in especies){ # sólo hazlo para especies diferentes
    if(sp != spfocal){
      Aij2[spfocal, sp] <- sum(milog(envueltas[[spfocal]][[sp]]$obs[rok2][-1]))
      Aij1[spfocal, sp] <- sum(milog(envueltas[[spfocal]][[sp]]$obs[rok1][-1]))
    }
  }
}



##### 3.) Creación de los histogramas con representación del PORCENTAJE de PARES para r=0,2 y r=2 #####
#######################################################################################################

total <- (n^2)-n
histAij2 <- hist(Aij2)
histAij2$counts <- (histAij2$counts/total)*100
plot(histAij2, main="", las=1, ylab="Proporción pares de especies (%)",
     xlab="Fuerza de asociación espacial (Aij) para r=0,2 m",
     col=c("#4CA6FF", "#4CA6FF", "#4CA6FF", "#4CA6FF", "#4CA6FF"))




total <- (n^2) - n
histAij1 <- hist(Aij1)
histAij1$counts <- (histAij1$counts / total) * 100
plot(histAij1, main = "", las = 1, 
     ylab = "Proporción pares de especies (%)",
     xlab = "Fuerza de asociación espacial (Aij) para r=2 m",
     col = "#4CA6FF")
            
            

##### 4.) Proporción de asociaciones negativas, positivas y neutras (o aleatorias) estimadas   #####
##### mediante la función g_ij (r) y en base al modelo nulo de independencia Poisson homogéneo #####
####################################################################################################

# a) Primero se ponen todas las envueltas en una lista única (envueltas es una lista de listas)
total <- (n^2)-n
envueltas2 <- vector (mode="list", length=total)
# iniciamos el contador "i"
i <- 0
for(spfocal in especies){
  for (sp in especies){
    # sólo hazlo para especies diferentes
    if(sp != spfocal){
      i <- i+1
      envueltas2[[i]] <- envueltas[[spfocal]][[sp]]
    }
  }
}
# despues usamos la función proportion.idar

install.packages("idar")
library(idar)
#para r=2
cosa <-proportion.idar(envueltas2)
# reordena los datos

nr <- length(cosa$r)
datos <- data.frame(r=rep(cosa$r, 3),
                    value = c(cosa$percentage$p.accumulators,
                              cosa$percentage$p.repellers,
                              cosa$percentage$p.neutrals),
                    group = rep(c("Asociaciones positivas", "Asociaciones negativas", "Asociaciones neutras"), each=nr))
datos$group <- factor(datos$group, levels=c("Asociaciones positivas", "Asociaciones neutras", "Asociaciones negativas"))

# Representación del gráfico.
library(ggplot2)
colores <- c("blue", "lightgrey", "orange")
figura <-  ggplot(datos, aes(x=r, y=value, fill=group)) +
  geom_area ()+
  scale_fill_manual (values= colores )+
  xlab ("Radio de vecindad espacial (m)")+
  ylab ("Proporcion pares de especies (%)")+
  theme_bw() + 
  theme(panel.grid = element_blank())
figura



                               