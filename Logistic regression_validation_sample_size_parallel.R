# simulation-based approach to calculate the sample size required for 
# external validation of a particular prediction model for a binary outcome
# coded in parallel
# Gary Collins (02-December-2020)

library(doParallel)
library(foreach)
library(rms)
rm(list = ls())
dev.off(dev.list()["RStudioGD"])

set.seed(200520)

## needs the rms package to be installed 

hanley_se_auc <- function(x, y){
  ### standard error for the C-statistic based on Hanley's formula
  n1 <- sum(y == 0)
  n2 <- sum(y == 1)
  Q1 <- x / (2 - x)
  Q2 <- 2 * x^2 / (1 + x)
  sqrt((x * (1 - x) + (n1 - 1) * (Q1 - x^2) + (n2 - 1) * (Q2 - x^2)) / (n1 * n2))
}

my_fun <- function(){    
  # generate the linear predictor
  LP <- rnorm(n[k], mu.v, sd.v[j])
  pr <- 1 / (1 + exp(-LP))
  
  # generate the outcome
  y  <- rbinom(n[k], 1, pr)
  
  df <- data.frame(y = y, LP = LP)
  
  # fit calibration model
  fit <- lrm(y~LP, data = df, x = T, y = T)
  
  # get the c-statistic, its standard error and 95% CI width
  c.statistic    <- as.numeric(fit$stats['C'])
  c.statistic.se <- hanley_se_auc(c.statistic, y)
  c.width        <- 2 * qnorm(0.025, lower = F) * c.statistic.se
  
  # get the calibration slope, its standard error and 95% CI width
  slope       <- as.numeric(coef(fit)[2])
  slope.se    <- as.numeric(sqrt(diag(vcov(fit)))[2])
  slope.width <- 2 * qnorm(0.025, lower = F) * slope.se

  events <- table(y)[2]
  
  myList <- list(n = n[k], sd.v = sd.v[j], events = events, c.statistic, c.width, slope, slope.width)
  
  return(myList)
}

E         <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
base.prob <- 0.01
mu.v      <- qlogis(base.prob, lower = T)
sd.v      <- c(0.2, 0.4, 0.6, 0.8, 1.0) # standard deviation of the linear predictor
n         <- E / boot::inv.logit(mu.v)
n.E.scenarios  <- length(E)
n.sd.scenarios <- length(sd.v)

NSIM    <- 500 # number of simulations
n_cores <- detectCores()
cl      <- makeCluster(n_cores - 1)
registerDoParallel(cl)
system.time(tmp <- foreach(k = 1:length(n), .combine = 'cbind') %:% foreach(j = 1:length(sd.v), .combine = 'cbind') %:% foreach(i = 1:NSIM, .packages = 'rms', .combine = 'rbind') %dopar% {rslt <- my_fun()})
stopCluster(cl)
OUT <- matrix(colMeans(matrix(unlist(tmp), ncol = 7 * 55, byrow = F)), ncol = 7,byrow = T)
OUT <- data.frame(OUT)
names(OUT) <- c("n", "sd(LP)", "events", "c.statistic", "c (width)", "slope", "slope (width)")

par(mfrow = c(1, 2))
plot(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51), 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51), 5], type = 'b', xlab = 'Average number of events', ylab = 'Average 95% CI width for the C-statistic', ylim = c(0, max(OUT[, 5])), xlim = c(0, max(OUT[, 3])), pch = 4, cex.lab = 0.75, xaxt = 'n')
axis(side = 1, at = E[-1], cex.axis = 0.75)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 1, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 1, 5], type = 'b', col = 2, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 2, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 2, 5], type = 'b', col = 3, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 3, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 3, 5], type = 'b', col = 4, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 4, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 4, 5], type = 'b', col = 5, cex = 0.5, pch = 4)
grid()
legend("topright", col = 1:6, legend = as.character(sd.v), lty = 1, title = expression(paste('SD(LP) ', sigma)), cex = 0.75)
abline(v = 100, lty = 2)
abline(v = 200, lty = 2)

plot(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51), 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51), 7], type = 'b', xlab = 'Average number of events', ylab = 'Average 95% CI width for the calibration slope', ylim = c(0, max(OUT[, 7])), xlim = c(0, max(OUT[, 3])), pch = 4, cex.lab = 0.75, xaxt = 'n')
axis(side = 1, at = E[-1], cex.axis = 0.75)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 1, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 1, 7], type = 'b', col = 2, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 2, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 2, 7], type = 'b', col = 3, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 3, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 3, 7], type = 'b', col = 4, cex = 0.5, pch = 4)
lines(OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 4, 3], OUT[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51) + 4, 7], type = 'b', col = 5, cex = 0.5, pch = 4)
grid()
legend("topright", col = 1:6, legend = as.character(sd.v), lty = 1, title = expression(paste('SD(LP) ', sigma)), cex = 0.75)
abline(v = 100, lty = 2)
abline(v = 200, lty = 2)

