



### PARALLEL QUASSE  routine
### for a set of simulated phylogenies with simulated traits
### used to estimate TYPE I and TYPE II ERRORS

### TX ... trait vector
### tre ... phylogeny
### stromecky, znakyx ... datasets with phylogenies and traits
###                       (later designated as tre, TX)

rm (list = ls())

ptm  <- proc.time ()


### SET THE PARALLEL COMPUTATION
install.packages (c("foreach", "multicore", "doMC"), repos = "http://cran.us.r-project.org")
library (multicore)
library (foreach)
library (doMC)

### REQUIRED NUMBER OF CORES
registerDoMC(5) ### specify the number of requested cores (here, it is 5)



### START THE MAIN CALCULATION
install.packages (c("geiger","diversitree"), repos = "http://cran.us.r-project.org")
library (geiger)
library (diversitree)

### INPUT DATA AND PARAMETERS

### Load the input data
load ("smaz_stromecky")
load ("smaz_znakyx")

### specify the number of iterations
ITER <- length(stromecky)

### specify the strength of density dependence (values 1-6)
##DENS_STR <- 1


### check the number of cores that are running
cat("\n","\n"," >>>> Number of cores currently running:", getDoParWorkers() ,"\n","\n")


### START THE PARALLEL FOREACH LOOP and save its results to 'vztah_all'
vztah_all <- foreach (i = 1:ITER, .combine = rbind) %dopar% {
  
  # ultr <- FALSE
  # while (ultr == FALSE){
  # rand_num <- ceiling (runif (n = 1, min = 10, max = 24990))
  # stromx <- siml_saveFRC [[ rand_num ]] ### randomly select one of the trees
  # strom  <- stromx [[DENS_STR]] ### specifies the strength of density dependence: DENS_STR values between 1 and 6
  # tre    <- strom$phylogeny
  # TX     <- as.numeric (na.omit (strom$trait))
  # ultr   <- is.ultrametric(tre)
  # }
  
  
  
  ### vyhozeni nulovych koncovych vetvi a odpovidajicich TX hodnot
  #if (sum (which (tre$edge.length == 0) > 0)){
  #   vetve <- cbind (tre$edge, tre$edge.length)
  #   vyhod <- vetve [which (tre$edge.length == 0), 2]
  #
  #   tre   <- drop.tip (tre, vyhod)
  #   TX    <- TX[-vyhod]
  #   }
  
  tre <- stromecky[[i]]
  TX <- znakyx [[i]]
  
  tre$tip.label <- 1:length(tre$tip.label)
  
  states <- TX
  names (states) <- tre$tip.label
  states.sd <- sd (states)
  
  ### (1) SPECIATION LINEARLY DEPENDENT ON TRAIT VALUES
  linearx <- make.linear.x (min (states), max(states))
  lik     <- make.quasse(tre, states, states.sd,
                         lambda = linearx,         ### speciation rate function
                         mu = constant.x,          ### extinction rate function
                         sampling.f = 1)           ### proportion of included species
  
  lik.nodrift <- constrain(lik, drift ~ 0)
  p <- starting.point.quasse(tre, states)
  p.start.lin <- c(0.2, 0.7, 0.1, 5) ### original: c(p[1], 0, p[2:3])
  
  ### search limits
  lower <- c(
    0, ### lower limit in search for minimal speciation rate
    0, ### lower limit in search for maximal speciation rate
    0, ### lower limit in search for minimal extinction rate
    0  ### lower limit for the diffusion coefficient
  )
  upper <- c(
    1, ### upper limit in search for minimal speciation rate
    1, ### upper limit in search for maximum speciation rate
    1, ### upper limit in search for maximal extinction rate
    50  ### upper limit for the diffusion coefficient
  )
  
  ### previously used values (the search range here is too wide)
  ### original:  lower <- c(-Inf, -Inf, 0, 0)
  ### original: upper <- c( Inf,  Inf, Inf, Inf)
  
  cat ("Control point A", "\n")
  
  fit.linear <- find.mle(lik.nodrift, p.start.lin,
                         lower = lower, upper = upper, verbose = 0)
  
  cat ("Control point B", "\n")
  
  ### (2) SPECIATION INDEPENDENT OF TRAIT  VALUES
  lik.nodrift.const <- constrain (lik.nodrift, l.m ~ 0)
  p.start.const  <- c(0.5, 0.1, 5)  ### original: c(p[1],p[2:3])
  
  lower <- c(
    0, ### min constant speciation
    0, ### min constant extinction
    0  ### min diffusion
  )
  upper <- c(
    1, ### max constant speciation
    1, ### max constant extinction
    50 ### max diffusion
  )
  ### original: lower <- c(-Inf, 0, 0)
  ### original: upper <- c( Inf, Inf, Inf)
  fit.const <- find.mle (lik.nodrift.const, p.start.const,
                         upper = upper, lower = lower, verbose = 0)
  
  cat ("Control point C", "\n")
  
  
  ### COMPARE AND SAVE THE RESULTS
  
  anovax <- anova(fit.const, fit.linear)
  
  ### save the immediate results to the 'vztah' dataframe
  vztah <- c(NA, NA)
  if (anovax$AIC[1] > anovax$AIC[2]) { vztah[1] <- "LIN" } else { vztah[1] <- "KONST" }
  if ( anovax$Pr[2] > 0.05)          { vztah[2] <- "N-SIG"} else { vztah[2] <- "SIG"   }
  
  
  
  cat ("Control point D", "\n")
  
  
  ### save important file to the directory
  save (anovax, file = paste ("QUASSE_anova", i, sep = ""))
  write.table (vztah, file=paste ("vztah", i, sep = ""))
  
  #cat ("Control point E", "\n")
  
  print (paste(">>> Run", i , "out of", ITER,">>> Result:", vztah[1], vztah[2], sep = " "))
  
  ### save the looping result to vztah_all through the "foreach" routine
  c (i, NA, vztah[1], vztah[2])
  
}

write.table (vztah_all, file = "vztah_all")

proc.time () -  ptm



### SUMMARIZE RESULTS
LIN_SIG <- LIN_NSIG <- KONST_SIG <- KONST_NSIG <- 0
for (j in 1:length(stromecky)){
  if (vztah_all[j,3] == "LIN" & vztah_all[j,4] == "SIG")    { LIN_SIG   <- LIN_SIG  + 1}
  if (vztah_all[j,3] == "LIN" & vztah_all[j,4] == "N-SIG")  { LIN_NSIG  <- LIN_NSIG + 1}
  if (vztah_all[j,3] == "KONST" & vztah_all[j,4] == "SIG")      { KONST_SIG    <- KONST_SIG  + 1}
  if (vztah_all[j,3] == "KONST" & vztah_all[j,4] == "N-SIG")    { KONST_NSIG   <- KONST_NSIG + 1}
}
cat ("\n","RESULT SUMMARY","\n")
cat ("LIN-SIG:",    LIN_SIG,"\n")
cat ("LIN-NSIG:",   LIN_NSIG,"\n")
cat ("KONST-SIG:",  KONST_SIG,"\n")
cat ("KONST-NSIG:", KONST_NSIG,"\n")










