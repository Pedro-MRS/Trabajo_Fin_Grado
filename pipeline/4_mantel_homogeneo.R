#######################################################
#######################################################
#####                                             #####
#####             TEST DE MANTEL                  #####
#####                                             #####
#######################################################
#######################################################


# Se carga el paquete FD
library (FD)


#--------------------------------------------------------------
# 1) Test de Mantel y su representacion grafica

# A) Los rasgos funcionales ("traits") de las 19 especies que se estan estudiando.
# Calculo de las correlaciones entre las matrices de asociacion espacial y las matrices de 
# disimilitud de rasgos


# Se carga el paquete ecodist
library(ecodist)

# Se cargan los datos
load("traits.Rdata")

head(traits) # muestra los primeros rasgos
names (traits) # vector con el nombre de los rasgos

# Se calculan las distancias de Gower

# distancia "multitrait" (todos los traits a la vez)
multitrait.d <- gowdis(traits)

# distancia individual para cada trait
SLA.d <- gowdis(traits[,"SLA", drop=FALSE])
LDMC.d <- gowdis(traits[,"LDMC", drop=FALSE])
Seedmass.d <- gowdis(traits[,"seedmass", drop=FALSE])
LT.d <- gowdis(traits[,"LT", drop=FALSE])
h.d <- gowdis(traits[,"h", drop=FALSE])

# B) Para poder hacer los test de Mantel se necesita que las matrices de asociación:
# Para obtener los resultados del test de Mantel:
  # mantelr: para el valor del coeficiente "r" de correlación de Mantel
  # hay que fijarse en el valor pval3: el p valor para la hipótesis nula para "r = 0"
  # Los intervalos de confianza para el coeficiente de correlación van a ser ilim.2.5% y ulim.97.5%.
  
# tienen que tener formato "matriz de distancia":
Aij2.d <- as.dist(Aij2) # para r=0,2 m
Aij1.d <- as.dist(Aij1) # para r=2 m

# Test de mantel individual para cada trait en cada modelo
# Para r=0,25 m
mantel(Aij1.d~ multitrait.d)
mantel(Aij1.d~ SLA.d)
mantel(Aij1.d~ LDMC.d)
mantel(Aij1.d~ Seedmass.d)
mantel(Aij1.d~ LT.d)
mantel(Aij1.d~ h.d)


# para r= 2 m
mantel(Aij2.d~ multitrait.d)
mantel(Aij2.d~ SLA.d)
mantel(Aij2.d~ LDMC.d)
mantel(Aij2.d~ Seedmass.d)
mantel(Aij2.d~ LT.d)
mantel(Aij2.d~ h.d)



# 2) Gráficos Test de Mantel para ambos modelos 
# Se almacenan todos los test de Mantel (con todos los rasgos)  en una data.frame ("resultados"),
# con una columna que se corresponda con el nombre del rasgo
# Esto facilita su representacion grafica con ggplot2

# se agrupan los resultados en una data.frame llamada resultados (para Aij1) para el modelo homogeno (r=2 m). (para Aij1), r =2m
resultados1 <- data.frame(rbind(mantel (Aij1.d~ multitrait.d),
                                mantel(Aij1.d~ SLA.d), mantel(Aij1.d~ LDMC.d), mantel(Aij1.d~ Seedmass.d), mantel(Aij1.d~ LT.d),
                                mantel(Aij1.d~ h.d)))

# se coloca una columna con el nombre de los rasgos de cada resultado:
resultados1$trait <- c("multitrait","SLA","LDMC","Seedmass", "LT", "h")

# se representa gráficamente con ggplot2:
library(ggplot2)
ggplot(resultados1, aes(y=trait,x=mantelr)) +
  geom_point() + xlim(-1,1)+
  geom_errorbarh(aes(xmin=llim.2.5., xmax=ulim.97.5. ), height=0)+
  geom_vline(xintercept = 0, lty=2)+
  xlab("Mantel r") + ylab ("Rasgo")+
  theme_light()+
  theme(panel.grid = element_blank()
  )

# se agrupan los resultados en una data.frame llamada resultados (para Aij2) para el modelo inhomogeneo (r=0,2 m).
resultados2 <- data.frame(rbind(mantel (Aij2.d~ multitrait.d),
mantel(Aij2.d~ SLA.d), mantel(Aij2.d~ LDMC.d), mantel(Aij2.d~ Seedmass.d),
mantel(Aij2.d~ LT.d), mantel(Aij2.d~ h.d)))
  
# se coloca una columna con el nombre de los rasgos de cada resultados.
resultados2$trait <- c("multitrait","SLA","LDMC","Seedmass", "LT", "h")
  
# se representa gráficamente con ggplot2:
  
ggplot(resultados2, aes(y=trait,x=mantelr)) +
    geom_point() + xlim(-1,1)+
    geom_errorbarh(aes(xmin=llim.2.5., xmax=ulim.97.5. ), height=0)+
    geom_vline(xintercept = 0, lty=2)+
    xlab("Mantel r") + ylab ("Rasgo")+
    theme_light()+          
    theme(panel.grid = element_blank()
    )
  
