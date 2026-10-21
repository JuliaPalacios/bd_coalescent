# Exact likelihood and MLE for all merger times of one genealogy.
# Observation convention: times = c(t_n, ..., t_2), t_(n+1) = 0.
# N is the absorbing forward population threshold; sampling population is N-1.
# Requires expm. Sourcing this file does not run simulations or create plots.

bd_B <- function(N, mu, eta) {
  m <- N - 1L
  rho <- eta / mu
  b <- numeric(m)
  b[m] <- 1 + rho
  if (m > 2L) for (k in seq.int(m - 1L, 2L)) {
    b[k] <- 1 + rho - rho / b[k + 1L]
  }
  b
}

simulate_bd_genealogies <- function(n, N, mu, eta, nsim = 100L) {
  stopifnot(length(n)==1L, n>=2, n==floor(n), length(N)==1L,
            N>n, N==floor(N), length(mu)==1L, is.finite(mu), mu>0,
            length(eta)==1L, is.finite(eta), eta>=0,
            length(nsim)==1L, nsim>=0, nsim==floor(nsim))
  m <- N - 1L
  b <- bd_B(N, mu, eta)
  k <- seq.int(2L, m)
  down <- up <- numeric(m)
  down[k] <- k * mu * b[k]
  if (m > 2L) {
    k <- seq.int(2L, m - 1L)
    up[k] <- k * eta / b[k + 1L]
  }
  times <- matrix(0, nsim, n - 1L,
                  dimnames=list(NULL, paste0('t_', seq.int(n, 2L))))
  for (r in seq_len(nsim)) {
    population <- m; lineages <- n; elapsed <- 0
    while (lineages > 1L) {
      total <- down[population] + up[population]
      elapsed <- elapsed + rexp(1L, total)
      if (runif(1L) < up[population] / total) {
        population <- population + 1L
      } else {
        p <- lineages * (lineages - 1) / (population * (population - 1))
        population <- population - 1L
        if (runif(1L) < p) {
          times[r, n-lineages+1L] <- elapsed
          lineages <- lineages - 1L
        }
      }
    }
  }
  times
}

# Prepare the no-merger generators once for each candidate parameter pair.
# Pure spectral calculations are unsafe for very ill-conditioned diagonal
# similarities. 'spectral' therefore explicitly uses expm for those stages.
# 'spectral_raw' is provided only for numerical diagnostics, not fitting.
bd_prepare <- function(n, N, mu, eta,
                       method=c('expm','spectral','spectral_raw'),
                       condition_limit=1e5) {
  method <- match.arg(method)
  stopifnot(n>=2, n==floor(n), N>n, N==floor(N),
            is.finite(mu), mu>0, is.finite(eta), eta>=0)
  b <- bd_B(N, mu, eta)
  m <- N - 1L
  stages <- vector('list', n-1L)
  for (l in seq.int(n, 2L)) {
    k <- seq.int(l, m); size <- length(k)
    c_rate <- l * (l-1) * mu * b[k] / (k-1)
    diagonal <- -k * (mu+eta)
    K <- diag(diagonal, nrow=size, ncol=size)
    log_r <- numeric(size)
    if (size > 1L) {
      i <- seq_len(size-1L)
      upper <- k[i] * eta / b[k[i]+1L]
      lower <- (k[i]+1L)*mu*b[k[i]+1L] *
        (1-l*(l-1)/(k[i]*(k[i]+1L)))
      K[cbind(i,i+1L)] <- upper
      K[cbind(i+1L,i)] <- lower
      if (eta>0) log_r <- c(0,cumsum(.5*(log(lower)-log(upper))))
    }
    stage <- list(K=K, c=c_rate, spectral=FALSE, fallback=FALSE,
                  condition=if (eta==0 && size>1L) Inf else exp(diff(range(log_r))))
    use_spectral <- method!='expm' && (size==1L || eta>0)
    if (method=='spectral' && stage$condition>condition_limit) use_spectral <- FALSE
    if (use_spectral) {
      H <- diag(diagonal, nrow=size, ncol=size)
      if (size>1L) {
        H[cbind(i,i+1L)] <- H[cbind(i+1L,i)] <- sqrt(upper*lower)
      }
      eig <- eigen(H, symmetric=TRUE)
      stage$values <- eig$values
      stage$vectors <- eig$vectors
      stage$r <- exp(log_r-max(log_r))
      stage$spectral <- TRUE
    }
    stage$fallback <- method!='expm' && !stage$spectral
    stages[[n-l+1L]] <- stage
  }
  list(stages=stages, n=n, N=N, mu=mu, eta=eta, method=method)
}

bd_loglik_prepared <- function(times, prepared, details=FALSE, score=FALSE) {
  n <- prepared$n; m <- prepared$N-1L
  stopifnot(is.numeric(times), length(times)==n-1L,
            all(is.finite(times)), all(diff(c(0,times))>0))
  delta <- diff(c(0,times))
  alpha <- numeric(m-n+1L); alpha[length(alpha)] <- 1
  loglik <- 0; fallback <- 0L
  alpha_dot <- numeric(length(alpha)); score_total <- 0
  for (i in seq_len(n-1L)) {
    stage <- prepared$stages[[i]]
    size <- length(alpha)
    used_fallback <- stage$fallback
    if (stage$spectral) {
      shift <- max(stage$values)
      exponentials <- exp((stage$values-shift)*delta[i])
      coefficients <- drop(crossprod(stage$vectors, alpha*stage$r))
      a <- drop(stage$vectors %*% (coefficients*exponentials)) / stage$r
      if (score) {
        coefficients_dot <- drop(crossprod(stage$vectors,alpha_dot*stage$r))
        a_dot <- drop(stage$vectors %*% (coefficients_dot*exponentials))/stage$r
      }
      # A valid no-merger distribution is nonnegative. Reject large negative
      # entries; tiny negative roundoff is clipped before normalization.
      invalid <- any(!is.finite(a)) || min(a) < -1e-12*max(1,max(abs(a)))
      if (invalid && prepared$method=='spectral_raw') return(if(details) list(loglik=NA_real_,fallback=fallback) else NA_real_)
      if (invalid) used_fallback <- TRUE
    }
    if (!stage$spectral || used_fallback) {
      shift <- -min(stage$c)
      shifted_K <- stage$K
      diag(shifted_K) <- diag(shifted_K)-shift
      E <- expm::expm(shifted_K*delta[i], method='Higham08.b')
      a <- drop(alpha %*% E)
      if (score) a_dot <- drop(alpha_dot %*% E)
    }
    fallback <- fallback+as.integer(used_fallback)
    a <- pmax(a,0)
    g <- sum(a*stage$c)
    if (!is.finite(g) || g<=0) return(if(details) list(loglik=-Inf,fallback=fallback) else -Inf)
    loglik <- loglik+log(g)+shift*delta[i]
    if (score) {
      # K and c scale linearly with total=mu+eta at fixed shape.
      a_dot <- a_dot+drop(a %*% stage$K)*delta[i]
      v_dot <- (a_dot+a)*stage$c
      score_i <- sum(v_dot)/g
      score_total <- score_total+score_i
      if (i<n-1L) alpha_dot <- c((v_dot-a*stage$c*score_i)/g,0)
    }
    if (i<n-1L) alpha <- c(a*stage$c/g,0)
  }
  if (details) list(loglik=loglik,fallback=fallback,stages=n-1L,
                   score_total=score_total) else loglik
}

coalescent_loglik <- function(times, N, mu, eta,
                             method=c('expm','spectral','spectral_raw'),
                             details=FALSE, condition_limit=1e5, score=FALSE) {
  prepared <- bd_prepare(length(times)+1L,N,mu,eta,match.arg(method),condition_limit)
  bd_loglik_prepared(times,prepared,details,score)
}

# Coalescent times alone identify the unordered pair {mu,eta}:
# mu B_k(eta/mu) = eta B_k(mu/eta). Use mu>=eta as a convention.
# Optimize log(total rate) and shape=4*mu*eta/(mu+eta)^2 in [0,1].
# shape=0 (pure birth) and shape=1 (equal rates) are included boundaries.
fit_coalescent_mle <- function(times, N, method=c('expm','spectral'),
                               starts=rbind(c(0,.9),c(0,1)),
                               control=list(maxit=150, factr=1e8, pgtol=1e-5),
                               shape_step=1e-5) {
  method <- match.arg(method)
  stopifnot(is.numeric(times), length(times)>=1L, all(is.finite(times)),
            all(diff(c(0,times))>0), N>length(times)+1L,
            is.matrix(starts),ncol(starts)==2L)
  unpack <- function(par) {
    total <- exp(par[1L]); root <- sqrt(max(0,1-par[2L]))
    c(rate_large=total*(1+root)/2,rate_small=total*(1-root)/2)
  }
  evaluations <- 0L; fallbacks <- 0L; attempted_stages <- 0L
  cache <- new.env(parent=emptyenv())
  evaluate <- function(par) {
    key <- paste(sprintf('%.15g',par),collapse=',')
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    evaluations <<- evaluations+1L
    rates <- unpack(par)
    value <- coalescent_loglik(times,N,rates[1L],rates[2L],method,
                              details=TRUE,score=TRUE)
    fallbacks <<- fallbacks+value$fallback
    attempted_stages <<- attempted_stages+length(times)
    assign(key,value,envir=cache)
    value
  }
  objective <- function(par) {
    value <- evaluate(par)
    if (is.finite(value$loglik)) -value$loglik else 1e6+sum(par^2)
  }
  gradient <- function(par) {
    value <- evaluate(par)
    if (!is.finite(value$loglik)) return(2*par)
    h <- shape_step
    lo <- hi <- par
    lo[2L] <- max(0,par[2L]-h)
    hi[2L] <- min(1,par[2L]+h)
    vlo <- evaluate(lo)$loglik; vhi <- evaluate(hi)$loglik
    if (!is.finite(vlo)) { vlo <- value$loglik; lo <- par }
    if (!is.finite(vhi)) { vhi <- value$loglik; hi <- par }
    shape_score <- if (hi[2L]>lo[2L]) (vhi-vlo)/(hi[2L]-lo[2L]) else 0
    -c(value$score_total,shape_score)
  }
  begin <- proc.time()
  # Scale bounds follow the observed time units. The generous interval can
  # be inspected with scale_boundary; it is not a biological rate constraint.
  rate_scale <- 1/median(diff(c(0,times)))
  log_bounds <- log(rate_scale)+c(-14,14)
  starts[,1L] <- pmin(pmax(starts[,1L],log_bounds[1L]),log_bounds[2L])
  fits <- lapply(seq_len(nrow(starts)),function(i) {
    optim(starts[i,],objective,gr=gradient,method='L-BFGS-B',
          lower=c(log_bounds[1L],0),upper=c(log_bounds[2L],1),control=control)
  })
  best <- fits[[which.min(vapply(fits,function(x)x$value,numeric(1)))]]
  # Validate with a finer shape derivative; curvature is high near equal rates.
  shape_step <- min(shape_step,1e-7)
  projected <- gradient(best$par)
  if(best$par[2L]<1e-7) projected[2L] <- min(projected[2L],0)
  if(best$par[2L]>1-1e-7) projected[2L] <- max(projected[2L],0)
  if(max(abs(projected))>1e-3) {
    refined <- optim(best$par,objective,gr=gradient,method='L-BFGS-B',
                     lower=c(log_bounds[1L],0),upper=c(log_bounds[2L],1),
                     control=list(maxit=300,factr=1e6,pgtol=1e-6))
    fits[[length(fits)+1L]] <- refined
    if(refined$value<=best$value) best <- refined
    projected <- gradient(best$par)
    if(best$par[2L]<1e-7) projected[2L] <- min(projected[2L],0)
    if(best$par[2L]>1-1e-7) projected[2L] <- max(projected[2L],0)
  }
  elapsed <- proc.time()-begin
  rates <- unpack(best$par)
  list(estimates=rates,loglik=-best$value,convergence=best$convergence,
       message=best$message,projected_gradient=max(abs(projected)),
       shape=best$par[2L],total=exp(best$par[1L]),
       shape_boundary=best$par[2L]<1e-6 || best$par[2L]>1-1e-6,
       scale_boundary=any(abs(best$par[1L]-log_bounds)<1e-5),
       method=method,evaluations=evaluations,fallbacks=fallbacks,
       fallback_fraction=if(attempted_stages) fallbacks/attempted_stages else 0,
       cpu_seconds=unname(sum(elapsed[c('user.self','sys.self')])),
       elapsed_seconds=unname(elapsed['elapsed']),optimizations=fits)
}
