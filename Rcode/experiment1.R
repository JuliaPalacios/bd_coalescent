#Generate S

Smat<-function(N,mu,eta){
  stopifnot(length(mu) == 1, length(eta) == 1,
            is.finite(mu), is.finite(eta), mu != 0)
  
  S <- matrix(0, nrow = N, ncol = N)
  
  S[1,1] <- -mu
  S[1,2] <- mu
  
  S[2,1] <- 2 * eta
  S[2,2] <- -2 * (mu + eta)
  S[2,3] <- 2 *mu
  
  S[3,2] <- 3 * eta
  S[3,3] <- -3 * (mu + eta)
  S[3,4] <- 3 * mu
  
  S[4,3] <- 4 * eta
  S[4,4] <- -4 * (mu + eta)
  S[4,5] <- 4 * mu
  
  S[5,4] <- 5 * eta
  S[5,5] <- -5 * (mu + eta)
  
  S
}

U<-solve(-Smat(N,mu,eta))
u<-U[1,]

Stilde<-solve(diag(u))%*%t(Smat(N,mu,eta))%*%diag(u)
#N is population size at sampling time (abs. 6 , N=5)
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

N<-5
mu<-.6
eta<-.3
rho<-eta/mu

library(expm)
trans_c<-function(t,N,mu,eta){
  Qrev<-tildeS(N,mu,eta)
  expm(Qrev*t)[N,]
}

trans_c2<-function(t,N,mu,eta){
  U<-solve(-Smat(N,mu,eta))
  u<-U[1,]
  Qrev<-solve(diag(u))%*%t(Smat(N,mu,eta))%*%diag(u)
  expm(Qrev*t)[N,]
}



exp_val<-function(t,n,N,mu,eta){
  trans<-trans_c(t,N,mu,eta)
  rho<-eta/mu
  lambda<-n*(n-1)*mu*sum(trans[n:N]*(1/((n-1):(N-1)))*sapply(n:N, B, N = N, rho = rho))
  return(lambda)
}

exp_val2<-function(t,n,N,mu,eta){
  trans<-trans_c2(t,N,mu,eta)
  rho<-eta/mu
  lambda<-n*(n-1)*mu*sum(trans[n:N]*(1/((n-1):(N-1)))*sapply(n:N, B, N = N, rho = rho))
  return(lambda)
}


lambda_ode <- function(t, n, N, mu, eta) {
  Nt<-max(2,N*exp(-(mu-eta)*t))
  2 * choose(n,2)*mu * (1/(Nt-1))
}

# lambda_ode <- function(t, n, N, mu, eta) {
#   2 * choose(n,2)*mu * exp((mu - eta) * t) / (N)
# }
x <- seq(0, 2, by = 0.01)

y1 <- sapply(x, exp_val, n = 2, N = 5, 0.6, 0.3)
y2 <- sapply(x, exp_val, n = 2, N = 5, 0.6, 0.5)
#y3 <- sapply(x, exp_val, n = 2, N = 5, 0.5, 0.6)
y4 <- sapply(x,lambda_ode,n=2, N=5, 0.6,0.3)
y5 <- sapply(x,lambda_ode,n=2, N=5, 0.6,0.5)
#y6 <- sapply(x,lambda_ode,n=2, N=5, 0.5,0.6)
y7 <- sapply(x, exp_val2, n = 2, N = 5, 0.6, 0.3)
y8 <- sapply(x, exp_val2, n = 2, N = 5, 0.6, 0.5)
#y9 <- sapply(x, exp_val2, n = 2, N = 5, 0.5, 0.6)

pdf("exp1.pdf", width = 10, height = 4)

plot(x, y1, type = "l", lwd = 2, col = "gray",
     xlab = "x", ylab = "lambda(x)",ylim=c(0,1),
     main = "lambda(x), n=2, N=6")

lines(x, y2, lwd = 2, col = "red")
#lines(x, y3, lwd = 2, col = "blue")
lines(x, y4, lwd = 2, col = "gray",lty=2)
lines(x, y5, lwd = 2, col = "red",lty=2)
#lines(x, y6, lwd = 2, col = "blue",lty=2)
lines(x, y7, lwd = 2, col = "gray",lty=3)
lines(x, y8, lwd = 2, col = "red",lty=3)
#lines(x, y9, lwd = 2, col = "darkblue",lty=5)

legend("topright",
       legend = c("(mu=0.6, eta=0.3)", "(0.6, 0.5)", "(0.5, 0.6)"),
       col    = c("black", "red", "blue"),
       lwd    = 2, lty = 1, bty = "n")
dev.off()


#######Simulate coalescent 
## ------------------------------------------------------------
## Birth–death simulation with individual labels (1..5),
## then sample 2 extant individuals at the first time N hits 5,
## then compute coalescence time measured BACKWARDS from present:
## ------------------------------------------------------------

simulate_bd_tree_to_N <- function(mu, eta, Nmax = 5) {
  #if(!is.null(seed)) set.seed(seed)
  stopifnot(Nmax >= 1)
  
  t <- 0
  n <- 1L
  alive <- 1L
  next_id <- 2L
  
  parent <- rep(NA_integer_, Nmax); parent[1] <- 0L
  birth_time <- rep(NA_real_, Nmax); birth_time[1] <- 0
  
  log <- data.frame(time=numeric(0), event=character(0),
                    id=integer(0), parent=integer(0), n=integer(0),
                    stringsAsFactors = FALSE)
  
  while(n > 0L && n < Nmax) {
    rate_birth <- n * mu
    rate_death <- n * eta
    if (n==1){rate_death<-0} #Cannot go extinct 
    rate_tot   <- rate_birth + rate_death
    
    # next event time
    t <- t + rexp(1, rate = rate_tot)
    
    if(runif(1) < rate_birth / rate_tot) {
      # Birth: choose a parent among alive, create a new ID
      p <- sample(alive, 1)
      id_new <- next_id
      next_id <- next_id + 1L
      
      parent[id_new] <- p
      birth_time[id_new] <- t
      
      alive <- c(alive, id_new)
      n <- n + 1L
      
      log <- rbind(log, data.frame(time=t, event="birth",
                                   id=id_new, parent=p, n=n))
    } else {
      # Death: remove a uniformly chosen alive individual
      v <- sample(alive, 1)
      alive <- setdiff(alive, v)
      n <- n - 1L
      
      log <- rbind(log, data.frame(time=t, event="death",
                                   id=v, parent=parent[v], n=n))
    }
  }
  
  if(n == 0L) stop("Extinction before reaching Nmax.")
  list(mu=mu, eta=eta, Nmax=Nmax,
       T_present=t, alive=alive, parent=parent,
       birth_time=birth_time, log=log)
}

get_lineage <- function(parent, tip) {
  lin <- integer(0)
  x <- tip
  while(x != 0L) {
    lin <- c(lin, x)
    x <- parent[x]
  }
  lin  # tip -> ... -> root
}

coalescence_time_backwards_two <- function(parent, birth_time, a, b, T_present) {
  la <- get_lineage(parent, a)
  lb <- get_lineage(parent, b)
  
  # MRCA id = common ancestor closest to the tips
  common <- intersect(la, lb)
  ra <- match(common, la)  # 1 = tip
  rb <- match(common, lb)
  mrca <- common[which.min(pmax(ra, rb))]
  
  # node just below MRCA on each path (toward each tip)
  ka <- match(mrca, la)
  kb <- match(mrca, lb)
  
  # if one tip IS the MRCA of the other, there is no branching event "giving rise to both"
  if(ka == 1L || kb == 1L) {
    return(list(
      a=a, b=b, mrca_id=mrca,
      t_branch=NA_real_, tau_backwards=NA_real_,
      note="One sampled individual is ancestor of the other (no split event)."
    ))
  }
  
  child_a <- la[ka - 1L]
  child_b <- lb[kb -  1L]
  
  # Branching event time (forward) that produced the split:
  # in this birth model, the split occurs when the *later* of the two child lineages was born
  t_branch <- min(birth_time[child_a], birth_time[child_b])
  
  # Time from present back to that branching event:
  tau <- T_present - t_branch
  
  list(
    a=a, b=b,
    mrca_id = mrca,
    child_on_path_to_a = child_a,
    child_on_path_to_b = child_b,
    t_branch = t_branch,
    T_present = T_present,
    tau_backwards = tau
  )
}

run_once <- function(mu, eta, Nmax = 5) {
  sim <- simulate_bd_tree_to_N(mu, eta, Nmax=Nmax)
  
  
  # sample 2 individuals from those alive at first hitting Nmax
  samp <- sample(sim$alive, 2)
  
  coal <- coalescence_time_backwards_two(
    parent = sim$parent,
    birth_time = sim$birth_time,
    a = samp[1], b = samp[2],
    T_present = sim$T_present
  )
  coal$tau_backwards
  list(sample=samp, coalescence=coal$tau_backwards, sim=sim)
}

## ---------------- Example ----------------
Nsim<-100
coal<-rep(0,Nsim)
for (j in 1:Nsim){
out <- run_once(mu=1.0, eta=0.2, Nmax=5)
out$coalescence
coal[j]<-out$coalescence$tau_backwards                    # event history
}



