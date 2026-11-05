#######################################################
#######################################################
#####                                             #####
#####             MODELO INHOMOGENEO              #####
#####                                             #####
#######################################################
#######################################################


#-------------------------------------------------------------------------- 
# 1) Se calcula lambda()
#--------------------------------------------------------------------------
# funciOn para estimar el ancho de banda (sigma)  que mejor estima la heterogeneidad
# de un patrOn ( a partir de un vector de posibles valores)
# Argumentos:
# Patrón: ppp univariado
# sigma.min: MInimo valor a considerar para el ancho de banda
# sigma.step: intervalos para definir el vector de sigmas
# nsim: numero de simulaciones para calcular las envueltas usadas en el LF.gof
# rmin: limite inferior de integracion para estimar LF.gof.
# rmin: limite superior de integracion para estimar LF.gof.

calcula.lambda <- function(patron, sigmas= NULL, sigma.min, sigma.max, 
                           sigma.step, nsim=19, rmin=NULL, rmax=NULL){
  require(ecespa) # para LF.gof
  
  # Si no se da el vector de sigmas se calcula
  # se mete calculo automatico sino no se da nada (min, max, step)
  
  if(is.null(sigmas)) sigmas <- seq(from = sigma.min, to = sigma.max, 
                                    by = sigma.step)
  u2 <-NULL
  
  for (i in sigmas){
    lambda <- density.ppp(patron, sigma=i, positive=TRUE)
    a <- envelope(Y=patron,fun= Kinhom, sigma=i,
                  nsim= nsim,savefuns=TRUE, simulate =
                    expression(rpoispp(lambda)))
    
    u2 <- c(u2, LF.gof(a,rmin=rmin, rmax=rmax)$u) 
  }
  
  sigma.ok <- sigmas[which.min(u2)]
  lambda <- density.ppp(patron, sigma=sigma.ok, positive=TRUE)
  attributes(lambda)$extra <- list(pp= patron, nsim=nsim,
                                   sigmas=sigmas, rmin=rmin, rmax=rmax)
  return (lambda)
}
#----------------------------------------------


#--------------------------------------------------------------------------
# 2) Elección sigma()
#--------------------------------------------------------------------------

# ALternativa a la anterior.  Consigue funciones inhomogeneas horizontales
# en modelos poisson inhomogeneo
# Función para elegir ancho de banda (sigma) para 
#    a) corregir Kinhom y 
#    b) estimar lambda(u)  de un Poisson inhomogeneo que se simular
# Se basa en que si la estimación  es correcta,
#   - la funcion Kinhom(r) del modelo Poisson inhom = pir^2
#   - por lo tanto L_mmean(r) =0 (horizontal)

# Argumentos
# pp: pp univariado
# sigmas: secuencia de valores de sigma a elegir
# nsim: numero de simulaciones del Poisson inhomogeneo para estimar L_mmena(r)
# seed: semilla para resultados repetuibles
# ...: argumentos adicionales para Kcross.inhom
# plotenv: Si debe plotear los diferentes test para cada sigma
# title: titulo a poner en el marco de la ventana gráfica (devnew(title=title)

#-----------------------------------------------------
#-----------------------------------------------------
# Funcion para elegir sigma

elige.sigma <- function(pp, sigmas, nsim=19,seed=1, ..., plotenv=FALSE, title=""){
  if(plotenv){
    setpanels(length(sigmas), title=title)
    par(mar=c(1,1,1,1))
  }
  
  # incializa matrices
  
  mmeanL <- NULL
  obsL <- NULL
  
  # envueltas para cada sigma
  
  set.seed(seed)
  for( s in sigmas){
    lambda <- density.ppp(pp, positive=TRUE, sigma=s)
    cosa <- envelope(pp, Kinhom, sigma=s,  nsim=nsim,
                     simulate=expression(rpoispp(lambda)),
                     positive =TRUE, verbose = FALSE,...)
    mmeanL <- cbind(mmeanL, sqrt(cosa$mmean/pi)-cosa$r)
    obsL <- cbind(obsL, sqrt(cosa$obs/pi)-cosa$r)
    if(plotenv) {
      plot(cosa,  sqrt(./pi)-r~r, legend=F, main="")
      
      # dibujar el poligono de desviaciones
      
      x<- cosa$r
      y1 <- 0*cosa$r^2
      y2 <- sqrt(cosa$mmean/pi)- cosa$r
      polygon(c(x,rev(x)),c(y2,rev(y1)),col=transparenta("orange",0.3), border=NA)
      lines(x, y1, col="darkgreen", lty=2, lwd=2)
      lines(x, y2, col=2, lty=2, lwd=2)
      
      # abline(h=0, lty=3, col=4)
    }
  }
  
  # Elige el mejor  sigma
  # valor medio de L(r) medio (línea roja)
  mmean.m <- apply(mmeanL,2,mean)
  # desviacion estandard de L(r) medio (linea roja)
  mmean.sd <- apply(mmeanL,2,sd)
  # Se elige la que produce una menor combinacion de media y sd
  # sigma.ok <- which.min(abs(mmean.sd)*abs(mmean.m))
  
  # Otra opcion seria la que tiene una desviacion menor RSS):
  RSS <- apply(mmeanL, 2, function(x) sum(x^2))
  sigma.ok <- which.min(RSS)
  
  # representacion de la desviacion del sigma elegido y de las otras sigmas
  # descarta las que el ultimo valor esta por debajo de 0 (sobreajuste)
  sobreajuste <- obsL[nrow(obsL),] < 0
  # descarta las que el ultimo valor esta por debajo del valor esperado
  sobreajuste2 <- obsL[nrow(obsL),] < mmeanL[nrow(mmeanL),]
  # "pendiente" = difernecia entre el valor para r=0 y rmax
  # para descartar las que tienen un valor excesivo
  pediente <- mmeanL[1,]-mmeanL[nrow(mmeanL),]
  # descarta combinado todo lo anterior
  descarta <- sobreajuste |sobreajuste2 | abs(pediente) > 0.001
  mantiene <- !descarta
  diferencias <- colSums(abs(mmeanL))
  
  #dev.new(title="desviaciones")
  #plot(sigmas[mantiene], diferencias[mantiene], type="o", ylab="desviación")
  #abline(h=0, lty=2)
  #points(sigmas[sigma.ok], diferencias[sigma.ok], pch=19, col=2)
  
  return(sigmas[sigma.ok])
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# 3) Funcion to compute envelopes for the Kcross.inhom con modelos nulos 
# "sin dejar fijadas las lambdas"
#  basados en kernel estimated intensities
# *** deja fijo el patron i  y simula el j ***de acuerdo con un Poisson inhomogeneo 
#  basado en una superficie de intensidad estimada con un kernel de ancho de banda sigmaJ
# Argumentos: Los mismos que en la funcion Kcross inhomogenea, y ademas:
# sigmaI: ancho de banda para estimar la inhomogeneidad del patron i
# sigmaJ: ancho de banda para estimar la inhomogeneidad del patron J y paras generar el
#                Poisson inhomogeneo que se simula.
# verbose: numera las simulaciones que se van realizando para cada test particular

Kci.env <- function(X, i, j, sigmaI=NULL, sigmaJ=NULL,correction="iso",..., nsim=19, nrank = 1, verbose=FALSE){
  library(spatstat.utils)
  # Kobs
  lambdaI <- density.ppp(X[X$marks==i], sigma=sigmaI, at="points", positive = TRUE)
  lambdaJ <- density.ppp(X[X$marks==j], sigma=sigmaJ, at="points", positive = TRUE)
  lambdaI2 <- density.ppp(X[X$marks==i], sigma=sigmaI, positive = TRUE) # im para simular
  lambdaJ2 <- density.ppp(X[X$marks==j], sigma=sigmaJ, positive = TRUE) # im para simular
  
  #      Kobs <- Kcross.inhom(X=X, i=i, j=j, lambdaI= lambdaI, lambdaJ=lambdaJ,
  #                                     correction=correction, ...)
  #  Se hacen unas envueltas de pega  para tener el frame para almacenar las nuestras
  Kenv <- envelope(X, Kcross.inhom, i=i, j=j, lambdaI= lambdaI2, 
                   lambdaJ=lambdaJ2, correction=correction, nsim=3, verbose=FALSE, 
                   savefuns=TRUE,...)
  # Calcular Kcross.inhom con at points y añadir a  la anterior
  Kobs <- Kcross.inhom(X, i=i, j=j, lambdaI= lambdaI, lambdaJ=lambdaJ, 
                       correction=correction,...)
  Kenv$obs <-  data.frame(Kobs)[,3]
  
  # Se guardan los atributos de las simfuns para modificarlos
  atsimfuns <- attributes(attributes(Kenv)$simfuns)
  atKenv <- attributes(Kenv)
  
  simvals <- NULL
  sstate <- list() # copiado de envelopeEngine
  for(s in 1:nsim){
    if(verbose) sstate <- progressreport(s, nsim, state = sstate)
    # simulate new sp pair: split, simulate y recompo?n el X
    # divide
    u <- split(X)
    # randomize (Poisson inhomogeneo) the point pattern of the label
    newj <- rpoispp(lambdaJ2)
    u[[j]] <- newj
    # recompose back the splited as a multivariate pp
    split(X) <- u
    # estimate inhomogenous intensity of the new simulated j 
    #lambdaJsim <- density.ppp(newj, sigma=sigmaJ, positive=TRUE, at="points") 
    # la frase anterior falla en ventanas irregulares (por los puntos que quedan fuera)
    # Mejor la siguiente, y por el mismo motivo, recalcula lambdaI (por algún motivo puede dejar fuera algun punto original)
    lambdaJsim <- density.ppp(X[X$marks==j], sigma=sigmaJ, positive=TRUE, at="points")
    lambdaI <- density.ppp(X[X$marks==i], sigma=sigmaI, at="points", positive = TRUE)
    # Calcula y almacena la funcion KCi simulada
    simvals <- cbind(simvals,
                     data.frame(Kcross.inhom(X=X, i=i, j=j, lambdaI= lambdaI, 
                                             lambdaJ=lambdaJsim,    correction=correction, ...))[,3])
  }
  # se extrae las envueltas al estilo spatstat
  
  nsim <- ncol(simvals)
  if (nrank == 1L) {
    lohi <- apply(simvals, 1L, range)
  } else {
    lohi <- apply(simvals, 1L, orderstats, k = c(nrank, 
                                                 nsim - nrank + 1L))
  }
  Kenv$lo <- lohi[1L, ]
  Kenv$hi <- lohi[2L, ]
  names(Kenv)[names(Kenv)=="theo"] <- "mmean"
  Kenv$mmean <- apply(simvals, 1L, mean, na.rm=TRUE)
  
  # Se crea el fv de funciones simuladas y se asigna sus atributos
  simvals <- data.frame(r = Kenv$r, simvals)
  class(simvals) <- c(class(simvals), "fv")
  
  atsimfuns$dotnames <- paste("sim",1:nsim, sep="")
  atsimfuns$desc <- c("distance argument r", paste("Simulation", 1:nsim))
  atsimfuns$labl <- c("r", atsimfuns$dotnames)
  atsimfuns$names <-   atsimfuns$labl 
  
  attributes(simvals) <- atsimfuns
  
  # Se actualiza los atributos del objeto envelope creado
  atKenv$names[3] <- "mmean"
  atKenv$labl[3] <- "{%s[%s]^{mmean}}(r)" 
  atKenv$desc[3] <- "mean  value of %s from simulations" 
  atKenv$dotnames <- c("obs" , "mmean", "hi",  "lo" )
  atKenv$einfo$nrank <- nrank
  atKenv$einfo$nsim <- nsim
  atKenv$einfo$Nsim <- nsim
  atKenv$einfo$nsim2 <- nsim
  atKenv$simfuns <- simvals
  atKenv$einfo$csr <- FALSE
  atKenv$einfo$csr.theo <- FALSE
  atKenv$einfo$use.theory <- FALSE
  atKenv$einfo$simtype <- "rpois"
  attributes(Kenv) <- atKenv
  
  return(Kenv)
}




#--------------------------------------------------------------------------
# 4) wKcross.inhom()
#--------------------------------------------------------------------------
# funcion para calcular la funcion Kcross inhomogenea ponderada 
#  (el pormedio ponderado entre Kij y Kji)

# Argumentos: los mismos que Kcross.inhom

wKcross.inhom <- function (X, i, j, lambdaI = NULL, lambdaJ = NULL, ...,
                           r = NULL, breaks = NULL, correction = c("border", "isotropic", "Ripley", 
                                                                   "translate"), sigma = NULL, varcov = NULL, lambdaIJ = NULL, 
                           lambdaX = NULL, update = TRUE, leaveoneout = TRUE) {
  
  ij <- Kcross.inhom(X=X, i=i, j=j, lambdaI = lambdaI, 
                     lambdaJ = lambdaJ, ..., r = r, breaks = breaks, 
                     correction = correction, sigma = sigma, varcov = varcov,
                     lambdaIJ = lambdaIJ , lambdaX = lambdaX, update = update,
                     leaveoneout = leaveoneout)
  
  ji <- Kcross.inhom(X=X, i=j, j=i, lambdaI = lambdaJ, 
                     lambdaJ = lambdaI, ..., r = r, breaks = breaks, 
                     correction = correction, sigma = sigma, varcov = varcov,
                     lambdaIJ = lambdaIJ , lambdaX = lambdaX, update = update,
                     leaveoneout = leaveoneout)	
  
  ni <- (X[X$marks==i])$n
  nj <- (X[X$marks==j])$n
  
  wKij <- eval.fv(((ni*ij) + (nj*ji))/(ni+nj))
  
  return(wKij)
}


#--------------------------------------------------------------------------
# 5) binhom() , biPCI () etc.
#--------------------------------------------------------------------------
# funciones para simular patrones multivariados 

# En un ppp multivariado (pp), se dejan todos fijos menos uno (identificado por "label")
# que se simulan mediante un Poisson inhomogeneo de intensidad "lamnbda"

binhom <- function(pp, label, lambda, nsim=1){
  lista.simu <- list()
  for(i in 1:nsim){
    # split the multivariate pp
    u <- split(pp)
    # randomize (Poisson inhomogeneo) the point pattern of the label
    u[[label]] <- rpoispp(lambda)
    # recompose back the splited as a multivariate pp
    split(pp) <- u
    lista.simu[[i]] <- pp
  }
  if(length(lista.simu)>1) return(lista.simu) else return(lista.simu[[1]])   
  
} 


# En un ppp multivariado (pp), se dejan todos fijos menos uno (identificado por "label")
# que se simulan mediante un Poisson inhomogeneo estimado del original con ancho de banda "sigma"

binhomS <- function(pp, label, sigma, nsim=1){
  lista.simu <- list()
  for(i in 1:nsim){
    # split the multivariate pp
    u <- split(pp)
    # estimate inhomogenous intensity
    lambda <- density.ppp(u[[label]] , sigma=sigma)
    # randomize (Poisson inhomogeneo) the point pattern of the label
    u[[label]] <- rpoispp(lambda)
    # recompose back the splited as a multivariate pp
    split(pp) <- u
    lista.simu[[i]] <- pp
  }
  if(length(lista.simu)>1) return(lista.simu) else return(lista.simu[[1]])   
  
} 




# En un ppp multivariado (pp), se dejan todos fijos menos uno (identificado por "label")
# que se simulan mediante un Poisson cluster inhomogeneo "model" (ajusatdo con ipc.estK)
biPCI <- function(pp, label,model, nsim=1){
  lista.simu <- list()
  for(i in 1:nsim){
    # split the multivariate pp
    u <- split(pp)
    # randomize (Poisson inhomogeneo) the point pattern of the label
    u[[label]] <- rIPCP(model)
    # recompose back the splited as a multivariate pp
    split(pp) <- u
    lista.simu[[i]] <- pp
  }
  if(length(lista.simu)>1) return(lista.simu) else return(lista.simu[[1]])   
  
} 



#--------------------------------------------------------------------------
# 6)  Prepara datos all_cross_auto2
# 
#--------------------------------------------------------------------------


# Funcion para automatizar el calculo de los test de independencia (envueltas bivariadas) 
#  entre todos los pares de especies (seleccionadas) de un patron multivariado
# Argumentos :
# pp: patron multivariado. Las marks son exclusivamente un vector de nombres de sp.
# nmin= numero minimo de individuos para considerar una especie
# lambdas? o models?: para ajustar los modelso nulos y/o intensidades de las sp...
#
# nsim.sigmas: numero de simulaciones para los LF.gofs en la seleccion de las sigmas
# nsim.test: numero de simulaciones para las envueltas de lso test de independencia
# plotenv= Logical (plot or not thje test for choosing sigmas
# seed: random seed for reproducible results
# correction: correccion del efecto borde
# nrank: rango para las envueltas 
# ... : argumentos paralas Kcross etc.

all_cross_auto2 <- function (pp, nmin, sigmas, nsim.sigmas, nsim.test, 
                             correction= "iso", plotenv, seed=1, nrank=1, ...){
  
  sp.ok <- names(which(table(pp$marks) >= nmin))
  pp.ok <- pp[pp$marks %in% sp.ok]
  pp.ok$marks <- factor(pp.ok$marks)
  
  pp <- pp.ok # se dejan solo las especies que cumplen las condiciones
  
  # select sigmas   for each species CON SELECT.SIGMAS
  
  lista.sigmas <- list()
  for (i in 1:length(sp.ok)) {
    cat("calculando sigma para ", sp.ok[i], "\n")
    pp.sp <- unmark(pp[pp$marks==sp.ok[i]])
    lista.sigmas[[i]] <-  elige.sigma(pp.sp, sigmas, plotenv=plotenv,
                                      seed=seed, nsim=nsim.sigmas, title=sp.ok[i])
  }
  names(lista.sigmas) <- sp.ok
  
  # Test de independencia entre todos los pares de especies   CON KCI.ENV
  
  # Se realizan los test entre todos los pares de especies
  wIJ.test <- list()
  for (i in 1:length(sp.ok)) {
    wIJ.test[[i]] <- list()
    for (j in 1:length(sp.ok)) {
      cat("calculando test de independencia entre ", sp.ok[i], "y ",  sp.ok[j],"\n")
      if (i !=j){
        wIJ.test[[i]][[j]] <- Kci.env(pp.ok, i=sp.ok[i], j=sp.ok[j],
                                      sigmaI=lista.sigmas[[sp.ok[i]]], sigmaJ=lista.sigmas[[sp.ok[i]]] ,
                                      correction=correction, nsim=nsim.test, nrank = nrank,...)
      }
      if (i ==j){
        wIJ.test[[i]][[j]] <- envelope (pp.ok[pp.ok$marks==sp.ok[i]], Kinhom, 
                                        sigma=lista.sigmas[[sp.ok[i]]],   correction=correction, nsim=nsim.test, 
                                        nrank = nrank, verbose=FALSE, savefuns=TRUE,...)
      }
    }
    names(wIJ.test[[i]]) <- sp.ok
  }
  names(wIJ.test) <- sp.ok
  
  return(list(lista.sigmas= lista.sigmas,
              wIJ.test = wIJ.test))
  
}




#--------------------------------------------------------------------------
# 7)  Preparación de datos (version previa de la anterior): otra forma de seleccionar lambdas
#--------------------------------------------------------------------------


# Funcion para automatizar el calculo de los test de independencia (envueltas bivariadas) 
#  entre todos los pares de especies (seleccionadas) de un patron multivariado
# Argumentos :
# pp: patron multivariado. Las marks son exclusivamente un vector de nombres de sp.
# nmin= numero minimo de individuos para considerar una especie
# lambdas? o models?: para ajustar los modelos nulos y/o intensidades de las sp...
#
# nsim.sigmas: numero de simulaciones para los LF.gofs en la seleccion de las sigmas
# nsim.test: numero de simulaciones para las envueltas de los test de independencia
# ... : argumentos paralas Kcross etc.
all_cross_auto <- function (pp, nmin, sigmas, nsim.sigmas, nsim.test, ...){
  sp.ok <- names(which(table(pp$marks) >= nmin))
  pp.ok <- pp[pp$marks %in% sp.ok]
  pp.ok$marks <- factor(pp.ok$marks)
  
  pp <- pp.ok # se dejan solo las especies que cumplen las condiciones
  
  # select sigmas  and  estimate lambdas for each species
  # Alternative: select best model with select.spm **********
  lista.lambdas <- list()
  for (i in 1:length(sp.ok)) {
    pp.sp <- unmark(pp[pp$marks==sp.ok[i]])
    lista.lambdas[[i]] <- calcula.lambda(pp.sp, sigmas=sigmas, 
                                         nsim=nsim.sigmas)
  }
  names(lista.lambdas) <- sp.ok
  
  # test de independencia entre todos los pares de especies
  
  # Primero,  se generan listas de patrones multivariados simulados
  # para cada especie focal
  options(warn = 1)  # Para evitar que los warnings se conviertan en error
  lista.simulados <- list()
  for (i in 1:length(sp.ok)) {
    lista.simulados[[i]] <- binhom(pp.ok, label=sp.ok[i], 
                                   lambda=lista.lambdas[[sp.ok[i]]], nsim=nsim.test)
  }
  names(lista.simulados) <- sp.ok
  
  # Ahora se hacen los test entre todos los pares de especies
  wIJ.test <- list()
  for (i in 1:length(sp.ok)) {
    wIJ.test[[i]] <- list()
    for (j in 1:length(sp.ok)) {
      wIJ.test[[i]][[j]] <- envelope(pp.ok, wKcross.inhom,
                                     i=sp.ok[i], j=sp.ok[j], lambdaI = lista.lambdas[[sp.ok[i]]], 
                                     lambdaJ = lista.lambdas[[sp.ok[j]]], ...,
                                     #lambdaJ = lista.lambdas[[sp.ok[j]]], 
                                     simulate= lista.simulados[[sp.ok[i]]],
                                     nsim=nsim.test, savefuns=TRUE)
    }
    names(wIJ.test[[i]]) <- sp.ok
  }
  names(wIJ.test) <- sp.ok
  
  return(list(lista.simulados = lista.simulados,
              lista.lambdas= lista.lambdas,
              wIJ.test = wIJ.test))
  
}


#----------------------------------------------------------------------------------
#7) Funciones para copmprobar que  comprueba que las funciones observadas-ponderadas  en todos los pares reciprocos  IJ y JI son las mismas.
# ------------------------------------------------------------------------

# Saca el mfrow para un numero indeterminado de paneles
setpanels <- function(npanels, title=""){
  sq <- sqrt(npanels)
  if(abs(sq-round(sq))>0) {
    nr <- min(round(sq), npanels/round(sq))
    nc <- max(round(sq), npanels/round(sq))
  } else nr <- nc <- sq
  dev.new(title=title)
  par(mfrow= c(nr, nc))
}

# ------------------------------------------------------------------------
# Funcion para ponderar funciones observadas y envueltas
# Se puede usar para promediar funciones ponderadas o no ponderadas
# Argumentos:
#sim.i: envelopes (con simfuns) con la especie i simulada
#sim.j:  envelopes (con simfuns) con la especie j simulada
# ni: numero de puntos de la espcie i
# nj: numero de  puntos de la especie j

wenvelope <- function(sim.i, sim.j, ni, nj){
  # 1) envueltas y funcion promedio esperada
  promedio <- (data.frame(sim.i)[,-1]*ni+
                 data.frame(sim.j)[,-1]*nj)/
    (ni+nj)
  
  # fv. con envueltas y funcion promedio ponderadas 
  wsim <- sim.i
  wsim[,2:5] <- promedio   # DEVOLVER EL OBJETO "WSIM"
  
  # 2) Rango respecto a las envueltas
  all.i <- cbind(data.frame(sim.i)$obs, 
                 data.frame(attributes(sim.i)$simfuns)[,-1])
  all.j <- cbind(data.frame(sim.j)$obs,
                 data.frame(attributes(sim.j)$simfuns)[,-1])
  rank.i <- apply(	all.i, 1, function(x) rank(x)[1])
  rank.j <- apply(	all.j, 1, function(x) rank(x)[1])
  wrank <- (rank.i*ni  + rank.j*nj) /(ni+nj)
  
  
  
  return(list(wenvelope = wsim, wrank = wrank))
  # el numero de simulaciones esta en  attributes(wenvelope)$einfo$nsim
}
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# se procesan las envueltas  para cada par de especies  en una parcela
# (saca el promedio de las  envueltas de  los mismos pares (calculados con simulaciones de cada focal)

# Argumentos:
# lista.test: lista de listas con los test de cada especie focal con el resto
# sp.ok: nombres de las especies testadas (dan nombre a los elementos de lista.test
# nsp : vector (con names) o tabla con el numero de individuos de las especies


wEnvelopelist <- function(lista.test, sp.ok, nsp){
  # copia el frame. Se reescriben dos veces cada ponderada
  # para evitar problemas computacionales (realmente se necesita el equivalente 
  # a una semimatriz
  
  result.env <- result.rank <- lista.test 
  for(i in sp.ok){
    for(j in sp.ok){
      # se hacen los calculos solo si no se han hecho ya
      if(!inherits(result.env[[i]][[j]], "ponderada")){
        ni <- nsp[i] 
        nj <- nsp[j] 
        # estas dos son la misma funcion ponderada
        # pero con envueltas simulando i y con envueltas simulando j
        sim.i <- lista.test[[i]][[j]]
        sim.j <- lista.test[[j]][[i]]
        cosa <- wenvelope(sim.i, sim.j, ni, nj)
        result.env[[i]][[j]] <- result.env[[j]][[i]] <- cosa$wenvelope
        result.rank[[i]][[j]] <- result.rank[[j]][[i]] <- cosa$wrank
        class(result.env[[i]][[j]]) <- class(result.env[[j]][[i]]) <- 
          c(class(result.env[[i]][[j]]),"ponderada")
      }
    }
  }
  return(list(wenvelopelist= result.env, wranklist = result.rank))
}
# ------------------------------------------------------------------------

#--------------------------------------------------------------
# Se generan matrices de asociacion (la matriz A de He y Biswas, pero basada en la funcion L) 
# para un radio determinado
# Argumentos:
# R: valor del radio concreto
# lista.env: lista de listas con las funciones Kcross (o envueltas)  de cada par de especies


sacaA <- function(lista.env, R){
  nsp <- length(lista.env)
  A <- matrix(NA, nr=nsp, nc= nsp)
  dimnames(A) <- list(names(lista.env), names(lista.env))
  
  # Que valor de la funcion? El que se corresponde con el R elegido
  r <- data.frame(lista.env[[1]][[1]])$r
  Rp <- match(R, r)
  if(is.na(Rp)) Rp <- sum(r<= R)
  # R real (para la transformacion de K en L)
  R <- r[Rp]
  
  for( spr in 1: nsp){
    for( spc in 1: nsp){
      A[spr, spc] <- data.frame(lista.env[[spr]][[spc]])[Rp,"obs"]
    }
  }
  # se transforma en  L(R)
  A <- sqrt(A/pi)-R
  return(A)
}

# ------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# 8) Representación grafica porcentajes de asociaciones

# funcion para coger una lista de listas (lista de especies focales
# y,  para cada especie focal, lista de los test con todas sus contrincantes)
# extrater la lista de test *inter-especificos unicos
# y devolver el grafico de porcentaje de asociaciones +, - y 0
# Argumentos:
# lista.env: lista de listas con el contenido mencionado anteriormente
# con los elementos de las listas con nombre (nombre de la especie correspondiente)
# alpha: nivel de significacion para las desviaciones de independencia

prop.asoc <- function(lista.env, alpha=.5, returnenv=FALSE){
  
  #  Primero se ponen todas las envueltas en una unica lista
  #  descartando repeticiones en el caso de envueltas ponderadas simetricas
  especies <- names(lista.env)
  nsp <- length(lista.env)
  total <- ((nsp^2)-nsp)/2
  envueltas2 <-  vector (mode="list", length=total)
  
  # se inicia el contador "i"
  i <- 0
  for(spfocal in especies){
    for (sp in especies){
      # solo se hacen para especies diferentes
      if(sp != spfocal){
        i <- i+1
        envueltas2[[i]] <- lista.env[[spfocal]][[sp]]
        names(envueltas2)[[i]] <- paste(spfocal, sp)
      }
    }
    especies <- especies[-1]
  }
  
  # despues se usa la funcion proportion.idar
  library(idar)
  cosa <-proportion.idar(envueltas2, alfa=alpha)
  
  # se reordenan los datos
  nr <- length(cosa$r)
  datos <- data.frame(r=rep(cosa$r, 3), 
                      value = c(cosa$percentage$p.accumulators, 
                                cosa$percentage$p.repellers,
                                cosa$percentage$p.neutrals),
                      group = rep(c("Asociaciones positivas", 
                                    "Asociaciones negativas","Asociaciones neutras"),
                                  each=nr))
  datos$group <- factor(datos$group, 
                        levels=c("Asociaciones positivas", 
                                 "Asociaciones neutras",
                                 "Asociaciones negativas"))
  
  # Representacion del grafico. Puede buscar maneras de mejorar el aspecto (otros colores, titulos, etc)
  
  library(ggplot2)
  colores <- c("blue", "lightgrey", "orange")
  figura <-  ggplot(datos, aes(x=r, y=value, fill=group)) +
    geom_area ()+
    scale_fill_manual (values= colores, name="grupo")+
    xlab ("Radio de vecindad espacial (m)")+
    ylab ("Proporcion pares de especies (%)")+
    theme_bw() + 
    theme(panel.grid = element_blank())
  # return(figura)
  # comprobación que se ha realizado correctamente
  if(returnenv) return(list(figura=figura, envueltas2 = envueltas2))
  return(figura)
  
}
# ------------------------------------------------------------------------

#--------------------------------------------------------------
# 9) Calculo de distancias de rasgos funcionales

# Funcion para calcular de golpe todas las matrices de distancias de traits
calculaFD <- function(traits, fd="gower"){
  # calcularemos distancias para todos los rasgos y una distancia multitrait
  ndis <- ncol(traits)
  result <- vector(mode="list", length=ndis+1)
  names(result) <- c(names(traits),"multitrait")
  if(fd=="gower"){
    # distancia de Gower
    library(FD)
    # distancias para los rasgos individuales
    for ( i in 1:ndis){
      result[[i]] <- gowdis(traits[,i, drop=FALSE])
    }
    # añade la distancia multitrait
    result[[ndis+1]] <- gowdis(traits)
  }
  
  if(fd!="gower"){
    # distancia euclidea
    # distancias para los rasgos individuales
    for ( i in 1:ndis){
      result[[i]] <- dist(traits[,i, drop=FALSE])
    }
    # añade la distancia multitrait
    result[[ndis+1]] <- dist(traits)
  }
  return(result)
}

# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# Función para centrar funciones K, envueltas y simfuns.
# Esta funcion se hace para "centrar" envueltas etc en donde la función L(r) no quedaba horizontal.
# Sin embargo no deberia emplearse y lo que habria que hacer es ajustar correctamente las
# superficies de inhomogeneidad (eligiendo el sigma adecuado) para que  la L_mean(r) sea horizontal 
# L_mmena(r) = 0
# Argumentos:
# x: fv con funciones o envueltas


centraK<- function(x,...){
  # saca r
  r <- x$r
  
  # se guardan los atributos para centrarlos y devolverselos
  att <- attributes(x)
  simfuns <- att$simfuns
  
  
  # se transforma a L
  x <- eval.fv(sqrt(x/pi)-r)
  simfuns <- eval.fv(sqrt(simfuns/pi)-r)
  
  xmean <- x$mmean
  
  # se centra la función L
  x <- eval.fv(x-xmean)
  simfuns <- eval.fv(simfuns-xmean)
  
  # se vuelve a convertirla en la función K
  x<- eval.fv(pi*(x+r)^2)
  simfuns <- eval.fv(pi*(simfuns+r)^2)
  
  att$simfuns <- simfuns
  attributes(x) <- att
  return(x)
}



#-------------------------------------------------------------------------  
# Función para haacer transparentes los colores
transparenta <- function(x, alpha=0.5) apply(col2rgb(x)/255,2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))


