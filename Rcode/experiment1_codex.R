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
  
  # K_n tracks population sizes n,...,m before the first sample
  # coalescence. Only non-coalescing downward transitions remain.
  states <- seq.int(n, m)
  size <- length(states)
  c_rate <- n * (n - 1) * mu * b[states] / (states - 1)
  diagonal <- -states * (mu + eta)
  K <- diag(diagonal, nrow = size, ncol = size)

  log_d <- numeric(size)
  H <- K
  if (size > 1L) {
    i <- seq_len(size - 1L)
    k <- states[i]
    upper <- k * eta / b[k + 1L]
    p_coal <- n * (n - 1) / ((k + 1) * k)
    lower <- (k + 1) * mu * b[k + 1L] * (1 - p_coal)
    K[cbind(i, i + 1L)] <- upper
    K[cbind(i + 1L, i)] <- lower

    # K = D H D^(-1), with H symmetric.
    log_d <- c(0, cumsum(0.5 * (log(lower) - log(upper))))
    off <- sqrt(upper) * sqrt(lower)
    H[cbind(i, i + 1L)] <- off
    H[cbind(i + 1L, i)] <- off
  }
  d <- exp(log_d - max(log_d))
  stopifnot(all(is.finite(d)), all(d > 0))

  eig <- eigen(H, symmetric = TRUE)
  rates <- eig$values
  stopifnot(all(rates < 0))

  # S(t) = beta exp(K_n t) 1; f(t) = beta exp(K_n t) c.
  initial <- d[size] * eig$vectors[size, ]
  survival_coefficients <- initial *
    drop(crossprod(eig$vectors, 1 / d))
  density_coefficients <- initial *
    drop(crossprod(eig$vectors, c_rate / d))

  check_times <- function(t) {
    stopifnot(is.numeric(t), !anyNA(t), all(t >= 0))
  }

  # Factor out the slowest exponential to keep the hazard and
  # log survival accurate even when survival itself underflows.
  leading_rate <- max(rates)
  log_survival <- function(t) {
    check_times(t)
    result <- rep(-Inf, length(t))
    finite <- is.finite(t)
    scaled <- drop(exp(outer(t[finite], rates - leading_rate)) %*%
                     survival_coefficients)
    stopifnot(all(scaled > 0))
    result[finite] <- leading_rate * t[finite] + log(scaled)
    result[t == 0] <- 0
    result
  }

  lambda <- function(t) {
    check_times(t)
    weights <- exp(outer(t, rates - leading_rate))
    # At infinity only the leading eigenmode survives.
    weights[, rates == leading_rate] <- 1
    drop(weights %*% density_coefficients) /
      drop(weights %*% survival_coefficients)
  }

  list(
    lambda = lambda,
    trajectory = function(t) 1 / lambda(t),
    survival = function(t) exp(log_survival(t)),
    density = function(t) lambda(t) * exp(log_survival(t)),
    cumulative_hazard = function(t) -log_survival(t),
    total_hazard = Inf,
    K = K,
    states = states,
    absorption_rates = c_rate,
    n = n
  )
}

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
     xlab = "x", ylab = "lambda(x)",ylim=range(c(0, y1, y2, y3, y4)),
     main = "lambda(x), n=2, N=6")

lines(x, y2, lwd = 2, col = "red")
lines(x, y3, lwd = 2, col = "gray",lty=2)
lines(x, y4, lwd = 2, col = "red",lty=2)

legend("topright",
       legend = c("Conditional PH (mu=0.6, eta=0.3)",
                  "Conditional PH (mu=0.6, eta=0.5)",
                  "ODE (mu=0.6, eta=0.3)",
                  "ODE (mu=0.6, eta=0.5)"),
       col = c("gray", "red", "gray", "red"),
       lwd = 2, lty = c(1, 1, 2, 2), bty = "n")
dev.off()

##Now coalescent simulation by transformation
simulate_two <- function(model, nsim = 1000L) {
  stopifnot(
    model$n == 2,
    length(nsim) == 1L, is.finite(nsim),
    nsim >= 0, nsim == floor(nsim)
  )
  
  targets <- rexp(nsim)
  if (!length(targets)) return(numeric(0))

  # H(t) = -log S(t) tends to infinity, so all draws are finite.
  hi <- 1
  largest <- max(targets)
  while (model$cumulative_hazard(hi) < largest) {
    hi <- 2 * hi
    if (!is.finite(hi)) stop("Could not bracket coalescence times.")
  }
  
  vapply(targets, function(target) {
    uniroot(
      function(t) model$cumulative_hazard(t) - target,
      interval = c(0, hi),
      tol = 1e-9
    )$root
  }, numeric(1))
}

set.seed(123)
coal_times <- simulate_two(model03, 1000)

hist(coal_times)
abline(v=mean(coal_times),col="red")

# Eventual coalescence has probability one:
exp(-model03$total_hazard)

