##### Birth-Death Coalescent
##### Author: JP +Codex

make_coal_model <- function(n, N, mu, eta) {
  stopifnot(
    length(N) == 1L, is.finite(N), N >= 3, N == floor(N),
    length(n) == 1L, is.finite(n),
    n >= 2, n <= N - 1, n == floor(n),
    length(mu) == 1L, is.finite(mu), mu > 0,
    length(eta) == 1L, is.finite(eta), eta > 0
  )
  
  m <- N - 1L
  rho <- eta / mu
  stopifnot(is.finite(rho))
  
  # b[i] = B(i, N, rho).
  # Recurrence avoids repeatedly summing powers and handles rho = 1.
  # verified
  b <- numeric(m)
  b[m] <- 1 + rho
  if (m > 2L) {
    for (i in seq.int(m - 1L, 2L)) {
      b[i] <- 1 + rho - rho / b[i + 1L]
    }
  }
  
  i <- seq_len(m - 1L)
  upper <- i * eta / b[i + 1L]
  lower <- (i + 1L) * mu * b[i + 1L]
  
  diagonal <- -seq_len(m) * (mu + eta)
  diagonal[1L] <- -mu
  # Q has lower diagonal l values and upper diagonal a
  # Q = D H D^(-1), with H symmetric, D is a diagonal matrix
  # This gives a real, orthogonal eigendecomposition.
  #\[H=D^{-1}QD\], To enforce symmetry we get: \[
  #\frac{d_{i+1}}{d_i}=\sqrt{\frac{\ell_i}{a_i}}.
  #\]
  log_d <- c(0, cumsum(0.5 * (log(lower) - log(upper))))
  d <- exp(log_d - max(log_d))
  stopifnot(all(is.finite(d)), all(d > 0))
  
  H <- diag(diagonal)
  off <- sqrt(upper) * sqrt(lower)
  H[cbind(i, i + 1L)] <- off
  H[cbind(i + 1L, i)] <- off
  
  eig <- eigen(H, symmetric = TRUE)
  stopifnot(all(eig$values < 0))
  
  # lambda(t) = e_m' exp(Q t) w.
  w <- numeric(m)
  j <- seq.int(n, m)
  ##this is ok for now, but needs to go in the simulation function, right now n=2 so it does not matter
  w[j] <- n * (n - 1) * mu * b[j] / (j - 1)
  coefficients <- d[m] * eig$vectors[m, ] *
    drop(crossprod(eig$vectors, w / d))
  
  rates <- eig$values
  
  check_times <- function(t) {
    stopifnot(is.numeric(t), !anyNA(t), all(t >= 0))
  }
  
  lambda <- function(t) {
    check_times(t)
    drop(exp(outer(t, rates)) %*% coefficients)
  }
  
  ##will be used for simulation later
  cumulative_hazard <- function(t) {
    check_times(t)
    drop(expm1(outer(t, rates)) %*% (coefficients / rates))
  }
  
  list(
    lambda = lambda,
    trajectory = function(t) 1 / lambda(t),
    cumulative_hazard = cumulative_hazard,
    total_hazard = sum(-coefficients / rates),
    n = n
  )
}

model03 <- make_coal_model(n = 2, N = 6, mu = 0.6, eta = 0.3)
model05 <- make_coal_model(n = 2, N = 6, mu = 0.6, eta = 0.5)

x <- seq(0, 2, by = 0.01)

y1 <- model03$lambda(x)
y2 <- model05$lambda(x)

model03 <- make_coal_model(n = 2, N = 6, mu = 0.6, eta = 0.3)
model05 <- make_coal_model(n = 2, N = 6, mu = 0.6, eta = 0.5)

x <- seq(0, 2, by = 0.01)

y1 <- model03$lambda(x)
y2 <- model05$lambda(x)


# pmax() makes the ODE approximation work on a vector of times.
lambda_ode_fast <- function(t, n, N, mu, eta) {
  Nt <- pmax(2, (N - 1) * exp(-(mu - eta) * t))
  n * (n - 1) * mu / (Nt - 1)
}

y3 <- lambda_ode_fast(x, 2, 6, 0.6, 0.3)
y4 <- lambda_ode_fast(x, 2, 6, 0.6, 0.5)

pdf("exp1.pdf", width = 10, height = 4)

plot(x, y1, type = "l", lwd = 2, col = "gray",
     xlab = "x", ylab = "lambda(x)",ylim=c(0,1),
     main = "lambda(x), n=2, N=6")

lines(x, y2, lwd = 2, col = "red")
lines(x, y3, lwd = 2, col = "gray",lty=2)
lines(x, y4, lwd = 2, col = "red",lty=3)

legend("topright",
       legend = c("(mu=0.6, eta=0.3)", "(0.6, 0.5)", "(0.5, 0.6)"),
       col    = c("black", "red", "blue"),
       lwd    = 2, lty = 1, bty = "n")
dev.off()

##Now coalescent simulation by transformation
simulate_two <- function(model, nsim = 1000L) {
  stopifnot(
    model$n == 2,
    length(nsim) == 1L, is.finite(nsim),
    nsim >= 0, nsim == floor(nsim)
  )
  
  targets <- rexp(nsim)
  result <- rep(Inf, nsim)
  finite <- which(targets < model$total_hazard)
  
  if (!length(finite)) return(result)
  
  # Find one upper bracket for every finite draw.
  hi <- 1
  largest <- max(targets[finite])
  while (model$cumulative_hazard(hi) < largest) {
    hi <- 2 * hi
    if (!is.finite(hi)) stop("Could not bracket coalescence times.")
  }
  
  result[finite] <- vapply(targets[finite], function(target) {
    uniroot(
      function(t) model$cumulative_hazard(t) - target,
      interval = c(0, hi),
      tol = 1e-9
    )$root
  }, numeric(1))
  
  result
}

set.seed(123)
coal_times <- simulate_two(model03, 1000)

hist(coal_times)
abline(v=mean(coal_times),col="red")

# Probability of no finite coalescence under this model:
exp(-model03$total_hazard)
