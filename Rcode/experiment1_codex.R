##### Birth-Death Coalescent
##### Author: JP +Codex

# Hazard of the next merger, conditional on n distinct lineages and
# initial_population individuals at the start of this waiting interval.
# Initially initial_population = N - 1 (sampling just before hitting N).
make_coal_model <- function(n, N, mu, eta, initial_population = N - 1L) {
  stopifnot(
    length(N) == 1L, is.finite(N), N >= 3, N == floor(N),
    length(n) == 1L, is.finite(n),
    n >= 2, n <= N - 1, n == floor(n),
    length(initial_population) == 1L, is.finite(initial_population),
    initial_population == floor(initial_population),
    initial_population >= n, initial_population <= N - 1,
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

  # beta is concentrated at initial_population, which need not be m
  # after a previous merger. S(t) = beta exp(K_n t) 1 and
  # f(t) = beta exp(K_n t) c, so lambda(t) = f(t) / S(t).
  initial_index <- initial_population - n + 1L
  initial <- d[initial_index] * eig$vectors[initial_index, ]
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
    n = n,
    N = N,
    mu = mu,
    eta = eta,
    B = b,
    initial_population = initial_population
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

# Simulate all n-1 mergers by following the joint population/lineage CTMC.
# This samples the same waiting laws as the conditional PH hazards above,
# while retaining the population size needed for each subsequent interval.
# No independent restart of the population process occurs at a merger.
simulate_coalescent <- function(model, nsim = 1000L) {
  stopifnot(
    length(nsim) == 1L, is.finite(nsim),
    nsim >= 0, nsim == floor(nsim)
  )

  n <- model$n
  m <- model$N - 1L
  states <- seq.int(2L, m)
  down <- up <- numeric(m)
  down[states] <- states * model$mu * model$B[states]
  if (m > 2L) {
    k <- seq.int(2L, m - 1L)
    up[k] <- k * model$eta / model$B[k + 1L]
  }

  labels <- paste0(seq.int(n, 2L), "_to_", seq.int(n - 1L, 1L))
  times <- matrix(0, nrow = nsim, ncol = n - 1L,
                  dimnames = list(NULL, labels))
  population_sizes <- matrix(NA_integer_, nrow = nsim, ncol = n - 1L,
                             dimnames = list(NULL, labels))
  for (replicate in seq_len(nsim)) {
    population <- model$initial_population
    lineages <- n
    elapsed <- 0
    while (lineages > 1L) {
      total_rate <- up[population] + down[population]
      elapsed <- elapsed + rexp(1L, rate = total_rate)
      if (runif(1L) < up[population] / total_rate) {
        population <- population + 1L
      } else {
        # A population merger coalesces a uniformly chosen pair.
        p_coal <- lineages * (lineages - 1) /
          (population * (population - 1))
        population <- population - 1L
        if (runif(1L) < p_coal) {
          column <- n - lineages + 1L
          times[replicate, column] <- elapsed
          population_sizes[replicate, column] <- population
          lineages <- lineages - 1L
        }
      }
    }
  }

  waiting_times <- times
  if (n > 2L) {
    waiting_times[, 2:(n - 1L)] <-
      times[, 2:(n - 1L), drop = FALSE] -
      times[, 1:(n - 2L), drop = FALSE]
  }
  list(
    # Absolute times backward from sampling, in increasing order.
    times = times,
    # Durations between consecutive mergers (first measured from sampling).
    waiting_times = waiting_times,
    # Population sizes immediately after each merger.
    population_sizes = population_sizes
  )
}

set.seed(123)
coal_times <- simulate_two(model03, 1000)

hist(coal_times)
abline(v=mean(coal_times),col="red")

# Eventual coalescence has probability one:
#exp(-model03$total_hazard)


# General sample: each row has n-1 ordered times backward from sampling.
model4 <- make_coal_model(n = 4, N = 6, mu = 0.6, eta = 0.3)
coal4 <- simulate_coalescent(model4, nsim = 1000)
head(coal4$times)
head(coal4$waiting_times)
boxplot(coal4$times)
# Conditional hazard after the first merger in the first replicate:
# elapsed time below is measured since that merger, not since sampling.
after_first <- make_coal_model(
  n = model4$n - 1L, N = model4$N, mu = model4$mu, eta = model4$eta,
  initial_population = coal4$population_sizes[1L, 1L]
)
after_first$lambda(c(0, 0.5, 1))
