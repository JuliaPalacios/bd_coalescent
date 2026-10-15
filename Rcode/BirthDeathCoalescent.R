library("PhaseTypeR")
##---------------------------------------------
## Define forward-time birth-death rate matrix
##---------------------------------------------
BDRateMfct <- function(Nmx,Brate,Drate){
  BDRateM <- matrix(0,nrow=Nmx,ncol=Nmx)
  ## Births
  for (i in 1:(Nmx-1)){
    BDRateM[i,i+1] <- i*Brate
    BDRateM[(i+1),i] <- (i+1)*Drate
    BDRateM[i,i] <- -sum(BDRateM[i,])
  }
  BDRateM[Nmx,(Nmx-1)] <- Nmx*Brate
  BDRateM[Nmx,Nmx] <- -Nmx*Brate-Nmx*Drate
  BDRateM <- round(BDRateM,digits=10)
  return(BDRateM)
}
##-----------------------------------------------------------
## Sampling of labelled individuals from birth-death events
## Store the sampling in list called history
##-----------------------------------------------------------
rHistory <- function(events){
  alive <- 1    ## Current individual alive
  next_id <- 2  ## Next individual to be born
  history <- list()
  history[[1]] <- alive
  for ( i in 1:length(events) ) {
    event <- events[i]
    if (event == "birth") {
      alive <- c(alive, next_id) ## Update list of individuals alive
      next_id <- next_id + 1     ## Next individual to be born
    } else if (event == "death") {
      j <- sample(1:length(alive), size = 1) ## Sample individual to die
      alive <- alive[-j]                      ## Remove dead individual from list
    }
    history[[i+1]] <- alive       ## Update history of individuals
  }
  return(history)
}
##----------------------------------------------------------------
## Identify the coalescent time for a randomly sampled pair
## from the final population
##----------------------------------------------------------------
rPairCoalTime <- function(history,event.tm0){
  pair <- sample(history[[length(history)]], 2)
  ## Find all times when both individuals were alive
  times <- which(sapply(history, function(x) {all(pair %in% x) }))
  ## Find minimum time
  first_time <- min(times)
  coal.tm <- event.tm0[first_time]
  out <- list()
  out$coal.tm <- coal.tm
  out$pair <- pair
  return(out)
}
##-----------------------------------------------------
## Illustration of functions
##-----------------------------------------------------
## Define the birth-death process as PhaseType object
Nmx <- 6 ; Brate <- 0.6 ; Drate <- 0.3
BDRateM <- BDRateMfct(Nmx,Brate,Drate)
bgn.vec <- c( 1,rep(0,nrow(BDRateM)-1) )
PHobj <- PH( BDRateM,init_probs=bgn.vec ) 
## Full simulation from the PhaseType object
## (includes states and dwell-times before absorption)
sim <- rFullPH(PHobj) ## Full simulation 
## Corresponding event times
event.tm <- cumsum( sim$time )
event.tm0 <- c(0,event.tm)
## Identify if transition is a birth or death
bd.indx <- pmax(diff(sim$state),0)+1 ## 1=death, 2=birth
events <- c("death","birth")[bd.indx]
## Plot the event times and number of individuals
plot(event.tm0,c(sim$state,Nmx),type="s",
     xlim=range(event.tm0),ylim=c(0,Nmx+1),
     xlab="Forward in time",ylab="Number of individuals alive")
## Simulate a history and add to plot 
history <- rHistory(events)
event.tm.mid <- (event.tm0[-length(event.tm0)] + event.tm0[-1]) / 2
## Add history to plot
check.even <- function(x) {x==floor(x/2)*2}
for (i in 1:length(history)){
  text(event.tm.mid[i],sim$state[i],
       paste(history[[i]],collapse=","),pos=ifelse(!check.even(i),1,3),cex=0.5)
}
## Add coal.time and pair to plot
pair.tm <- rPairCoalTime(history,event.tm0)
abline(v=pair.tm$coal.tm,col="red",lwd=2)
text(pair.tm$coal.tm,2,paste(pair.tm$pair,collapse=","),pos=4,col="red")
## Add coalescent time to plot
tmEnd <- event.tm0[length(event.tm0)]
arrows(x0=tmEnd,y0=1,x1=pair.tm$coal.tm,y1=1,lwd=2,col="red",length=0.1)
##---------------------------------------------------------------
## Finally simulate pairwise coalescent times
##---------------------------------------------------------------
nSim <- 10000
tmV <- rep(0,nSim)
for (i in 1:nSim){
  sim <- rFullPH(PHobj)  
  event.tm <- cumsum( sim$time )
  event.tm0 <- c(0,event.tm)
  bd.indx <- pmax(diff(sim$state),0)+1 ## 1=death, 2=birth
  events <- c("death","birth")[bd.indx]
  history <- rHistory(events)
  pair.tm <- rPairCoalTime(history,event.tm0)
  tmEnd <- event.tm0[length(event.tm0)]
  tmV[i] <- tmEnd-pair.tm$coal.tm
}
brks <- seq(0,max(tmV),len=100)
hist( tmV,breaks=brks ,probability=TRUE)
## Fit the geometric distribution (assuming rate is constant for all times)
mean(tmV)
x.val <- seq(0,max(tmV),len=1000)
dns.exp <- dexp(x.val,rate=1/mean(tmV))
points(x.val,dns.exp,type="l",col="red",lwd=2)
## Fit an Erlang distribution
mn <- mean(tmV) ; vr <- var(tmV) ; k <- mn^2/vr
dns.Erl <- dgamma(x.val, shape = k, rate = k/mn) 
points(x.val,dns.Erl,type="l",col="blue",lwd=2)
##---------------------------------------
## ODE approximation for coalescent rate (NEED TO WORK ON THIS)
##---------------------------------------
CoalRate <- (2*Brate/Nmx)*exp( -(Drate-Brate*brks) )
dns.ode <- CoalRate*exp(-cumsum(CoalRate)*brks[2])*brks[2]
points(brks,dns.ode,type="l",col="blue",lwd=2)

##---------------------------------------
## Coalescent rate with PH reversal (eq. in 0.2)
##---------------------------------------
#N(0)=5 but B is in terms of N+1=6, fix later 
B <- function(i, N, rho){
  sum(rho^seq(0,N-i+1))/sum(rho^seq(0,N-i))
  #1
}

tildeS <- function(N,mu, eta) {
  stopifnot(length(mu) == 1, length(eta) == 1,
            is.finite(mu), is.finite(eta), mu != 0)
  
  rho <- eta / mu
  
  S <- matrix(0, nrow = N, ncol = N)
  
  S[1,1] <- -mu
  S[1,2] <- eta * (1 / B(2,N,rho))
  
  S[2,1] <- 2 * mu * B(2,N,rho)
  S[2,2] <- -2 * (mu + eta)
  S[2,3] <- 2 * eta * (1/B(3,N,rho))
  
  S[3,2] <- 3 * mu * B(3,N,rho)
  S[3,3] <- -3 * (mu + eta)
  S[3,4] <- 3 * eta * (1/B(4,N,rho))
  
  S[4,3] <- 4 * mu * B(4,N,rho)
  S[4,4] <- -4 * (mu + eta)
  S[4,5] <- 4 * eta *(1/B(5,N,rho))
  
  S[5,4] <- 5 * mu * B(5,N,rho)
  S[5,5] <- -5 * (mu + eta)
  
  S
}

N<-6
mu<-.6
eta<-.3
rho<-eta/mu

library(expm)
trans_c<-function(t,N,mu,eta){
  Qrev<-tildeS(N,mu,eta)
  expm(Qrev*t)[N,]
}

rate<-function(t,n,N,mu,eta){
  trans<-trans_c(t,N,mu,eta)
  rho<-eta/mu
  lambda<-n*(n-1)*mu*sum(trans[n:N]*(1/((n-1):(N-1)))*sapply(n:N, B, N = N, rho = rho))
  return(lambda)
}
library("phylodyn")
coal_time<-rep(0,Nsim)
for (j in 1:Nsim){
  coal_time<-coalsim(c(0),2,traj=rate,n=2,N=)
}
