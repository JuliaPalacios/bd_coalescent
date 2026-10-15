myfun <- function(i, N, rho){
  (1 - rho^(N - i)) / (1 - rho^(N - i - 1))
}

rhos <- c(0.2, 0.5, 0.8)
cols <- c("cyan","steelblue", "darkorange", "firebrick")

plot_for_N <- function(N){
  i <- 1:(N - 2)
  ylist <- lapply(rhos, function(r) myfun(i, N, r))
  ylim <- range(unlist(ylist), finite = TRUE)
  
  plot(i, ylist[[1]], type = "l", lwd = 2, col = cols[1],
       xlab = "i", ylab = expression(frac(A[i], A[i+1])),
       ylim = ylim, main = paste("N =", N))
  
  for(k in 2:length(rhos)){
    lines(i, ylist[[k]], lwd = 2, col = cols[k])
  }
  
  legend("topright", legend = paste0("rho = ", rhos),
         col = cols, lwd = 2, bty = "n")
}

pdf("Ai_over_Ai1_N20_N50.pdf", width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1) + 0.1)

plot_for_N(20)
plot_for_N(50)

dev.off()

###same but for rho>1
myfun <- function(i, N, rho){
  sum(rho^seq(0,N-i-1))/sum(rho^seq(0,N-i-2))
}

rhos <- c(1.2, 2)
cols <- c("cyan","steelblue", "darkorange", "firebrick")

plot_for_N <- function(N){
  i <- 1:(N - 2)
  ylist <- lapply(rhos, function(r) myfun(i, N, r))
  ylim <- range(unlist(ylist), finite = TRUE)
  
  plot(i, ylist[[1]], type = "l", lwd = 2, col = cols[1],
       xlab = "i", ylab = expression(frac(A[i], A[i+1])),
       ylim = ylim, main = paste("N =", N))
  
  for(k in 2:length(rhos)){
    lines(i, ylist[[k]], lwd = 2, col = cols[k])
  }
  
  legend("topright", legend = paste0("rho = ", rhos),
         col = cols, lwd = 2, bty = "n")
}

pdf("Ai_over_Ai1_N20_N50_2.pdf", width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1) + 0.1)

plot_for_N(20)
plot_for_N(50)

dev.off()
