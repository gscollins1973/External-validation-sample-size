# simulation-based approach to calculate the sample size required for 
# external validation of a particular prediction model for a binary outcome
# Gary Collins (02-December-2020)

set.seed(200520)

## needs the rms package to be installed 

E         <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
base.prob <- 0.05
mu.v      <- qlogis(base.prob, lower = T)
sd.v      <- c(0.2, 0.4, 0.6, 0.8, 1.0) # standard deviation of the linear predictor
n         <- E / boot::inv.logit(mu.v)
NSIM      <- 500 # number of simulations
n.E.scenarios  <- length(E)
n.sd.scenarios <- length(sd.v)

hanley_se_auc <- function(x, y){
  ### standard error for the C-statistic based on Hanley's formula
  n1 <- sum(y == 0)
  n2 <- sum(y == 1)
  Q1 <- x / (2 - x)
  Q2 <- 2 * x^2 / (1 + x)
  sqrt((x*(1-x) + (n1-1)*(Q1-x^2)+(n2-1)*(Q2-x^2)) / (n1*n2))
}

events         <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store the number of events
c.statistic    <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store c-Statistic
c.statistic.se <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store the standard error of the C-statistic (Hanley)
c.width        <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store CI widths of the C-statistic
slope          <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store the calibration slope
slope.se       <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store the standard error of the calibration slope
slope.width    <- array(dim = c(NSIM, n.E.scenarios, n.sd.scenarios)) # store the CI widths of the calibration slope

system.time(for(k in 1:n.sd.scenarios){
  for(j in 1:n.E.scenarios){
    for(i in 1:NSIM){
      
      # generate the linear predictor
      LP <- rnorm(n[j], mu.v, sd.v[k])
      pr <- 1 / (1 + exp(-LP))
      # generate the outcome
      y  <- rbinom(n[j], 1, pr)
    
      df <- data.frame(y = y, LP = LP)
      
      # fit calibration model
      fit <- rms::lrm(y~LP, data = df, x = T, y=T)
      
      # get the calibration slope, its standard error and 95% CI width
      slope[i,j,k]       <- as.numeric(coef(fit)[2])
      slope.se[i,j,k]    <- as.numeric(sqrt(diag(vcov(fit)))[2])
      slope.width[i,j,k] <- 2 * qnorm(0.025, lower = F) * slope.se[i, j, k]
      
      # get the c-statistic, its standard error and 95% CI width
      c.statistic[i,j,k]    <- as.numeric(fit$stats['C'])
      c.statistic.se[i,j,k] <- hanley_se_auc(c.statistic[i,j,k], y)
      c.width[i,j,k]        <- 2 * qnorm(0.025, lower = F) * c.statistic.se[i,j,k]
      
      # store the number outcomes simulated 
      events[i,j,k] <- table(y)[2]
    }
  }
})

# plot the results
par(mfrow = c(1, 2))
plot(colMeans(events[,,1]),  colMeans(c.width[,,1]), cex = 0.5, type = 'b', xlim = c(0, max(colMeans(events[,,5]))), ylim = c(0, max(c.width)), xlab = 'Average number of events', ylab = 'Average 95% CI width for the C-statistic', pch = 4, cex.lab=0.75, xaxt = 'n')
axis(side = 1, at = E[-1], cex.axis=0.75)
lines(colMeans(events[,,2]), colMeans(c.width[,,2]), cex = 0.5, type = 'b', col = 2, pch = 4)
lines(colMeans(events[,,3]), colMeans(c.width[,,3]), cex = 0.5, type = 'b', col = 3, pch = 4)
lines(colMeans(events[,,4]), colMeans(c.width[,,4]), cex = 0.5, type = 'b', col = 4, pch = 4)
lines(colMeans(events[,,5]), colMeans(c.width[,,5]), cex = 0.5, type = 'b', col = 5, pch = 4)
grid()
legend("topright", col = 1:6, legend = as.character(sd.v), lty = 1, title = expression(paste('SD(LP) ', sigma)), cex = 0.75)
abline(v = 100, lty = 2)
abline(v = 200, lty = 2)

plot(colMeans(events[,,1]),  colMeans(slope.width[,,1]), cex = 0.5, type = 'b', xlim = c(0, max(colMeans(events[,,5]))), ylim = c(0, max(slope.width)), xlab = 'Average number of events', ylab = 'Average 95% CI width for the calibration slope', pch = 4, cex.lab=0.75, xaxt = 'n')
axis(side = 1, at = E[-1], cex.axis=0.75)
lines(colMeans(events[,,2]), colMeans(slope.width[,,2]), cex = 0.5, type = 'b', col = 2, pch = 4)
lines(colMeans(events[,,3]), colMeans(slope.width[,,3]), cex = 0.5, type = 'b', col = 3, pch = 4)
lines(colMeans(events[,,4]), colMeans(slope.width[,,4]), cex = 0.5, type = 'b', col = 4, pch = 4)
lines(colMeans(events[,,5]), colMeans(slope.width[,,5]), cex = 0.5, type = 'b', col = 5, pch = 4)
grid()
legend("topright", col = 1:6, legend = as.character(sd.v), lty = 1, title = expression(paste('SD(LP) ', sigma)), cex = 0.75)
abline(v = 100, lty = 2)
abline(v = 200, lty = 2)



