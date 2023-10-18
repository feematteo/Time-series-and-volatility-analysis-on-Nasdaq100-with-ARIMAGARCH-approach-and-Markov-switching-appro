library(MSGARCH)
# 1 ----
data("dem2gbp", package = "MSGARCH")
ms2.garch.n <- CreateSpec(variance.spec = list(model = "sGARCH"),
                          distribution.spec = list(distribution = "norm"),
                           switch.spec = list(K = 2))
fit.ml <- FitML(ms2.garch.n, data = dem2gbp)
summary(fit.ml)


set.seed(1234)
fit.mcmc <- FitMCMC(ms2.garch.n, data = dem2gbp)
summary(fit.mcmc)


pred <- predict(fit.ml, nahead = 5, do.return.draw = TRUE)
pred$vol  #, is a numeric vector of length h, containing the standard
                                     #deviations of the distributions yT +j |IT

pred$draw     #, is a matrix of dimension h × nsim of the simulated MSGARCH process


#Quantitative risk management

# The VaR measures
# the threshold value such that the probability of observing a loss larger or equal to it in a
# given time horizon is equal to α. The ES measures the expected loss below the VaR level.

#Both VaR and ES are a matrix of dimension nahead × length(alpha)
risk <- Risk(fit.ml, alpha = c(0.01, 0.05), nahead = 5)


#the difference between predict() and simulate(). The function predict() is used for 
#forecasting the conditional volatility and the predictive distribution
#(using the argument do.return.draw = TRUE) while the function simulate() aims at generating 
#simulation paths for a given MSGARCH model

# 2 ----
data("SMI", package = "MSGARCH")
#in-sample analysis of the daily log-returns of the major equity index
#for the Swiss market
plot(SMI, xlab = "Date (year)")
# in financial time series, such as volatility clustering and presence of outliers, are evident 
#Furthermore, we also note that large (absolute) returns are more frequent 
#at the start  and at the end than in the middle 
#This suggests that the conditional variance is time-varying according to a regime-switching specification

## Create MS(2)-GJR-std specification
ms2.gjr.s <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                        distribution.spec = list(distribution = "std"),
                        switch.spec = list(K = 2),
                        constraint.spec = list(regime.const = "nu"))

## ML estimation
fit.ml <- FitML(ms2.gjr.s, data = SMI)

## Summary
summary(fit.ml)

# Fitted parameters:
#   Estimate Std. Error  t value  Pr(>|t|)
# alpha1_1   0.0005     0.0088   0.0568 4.774e-01
# alpha2_1   0.2137     0.0619   3.4505 2.798e-04
# beta_1     0.5264     0.0995   5.2921 6.045e-08
# 
# alpha1_2   0.0052     0.0169   0.3050 3.802e-01
# alpha2_2   0.1516     0.0381   3.9771 3.488e-05
# beta_2     0.8716     0.0354  24.6378    <1e-16

#Parameter estimates indicate that the evolution of the volatility process is heterogeneous
#across the two regimes. the two regimes report different unconditional volatility levels:

# the first regime
# is characterized by: (i) low unconditional volatility, (ii) strong volatility reaction to past
# negative returns, and (iii) low persistence of the volatility process

# the second
# regime is characterized by: (i) high unconditional volatility, (ii) weak volatility reaction to
# past negative returns, and (iii) high persistence of the volatility process

## Unconditional vol
set.seed(1234)
sqrt(250) * sapply(ExtractStateFit(fit.ml), UncVol)

## Smoothed probabilities in regime 2 and volatility
smoothed.prob <- State(fit.ml)$SmoothProb[, 1, 2, drop = TRUE]
vol <- sqrt(250) * Volatility(fit.ml)


op <- par(mfrow = c(2,1),
          oma = c(1,1,0,0) + 0.0,
          mar = c(2,2,2,2) + 0.0)
plot(as.vector(SMI), axes = FALSE, ann = FALSE)
par(new = TRUE)
ylabel <- expression(paste("Pr(", s[t], " = 2 | ", hat(psi), ", ", I[t], ")"))
plot(zoo::zoo(smoothed.prob, order.by = zoo::index(SMI)), plot.type = "single")
title(main = "Smoothed probabilities")
plot(zoo::zoo(vol, order.by = zoo::index(SMI)), plot.type = "single")
title(main = "Volatility (%)")
par(op)

#when the smoothed probabilities of regime two are near one, the filtered volatility 
#of the process sharply increases

# note that the Markov chain evolves persistently over time
# 
# Stable probabilities:
#   State 1 State 2 
# 0.5407  0.4593 

#!!!
# ML estimation can be difficult for MSGARCH-type models
# MCMC procedures can be used to explore the joint posterior distribution of
# the model parameters, thus avoiding convergence to local maxima commonly encountered
# via ML estimation

## MCMC estimation
nmcmc <- 12500
nburn <- 5000
nthin <- 5
ctr <- list(nmcmc = nmcmc, nburn = nburn,
            nthin = nthin, par0 = fit.ml$par)
fit.mcmc <- FitMCMC(ms2.gjr.s, data = SMI, ctr = ctr)
summary(fit.mcmc)


#Convergence of the chain
par(mfrow = c(3, 4))
coda::traceplot(fit.mcmc$par)
coda::heidel.diag(fit.mcmc$par)
coda::acfplot(fit.mcmc$par)

# in The Bayesian approach distributions of nonlinear functions of the model parameters 
# can be obtained at low cost by simulating from the joint posterior distribution, 
# and parameter uncertainty can be integrated in the forecasts through the predictive distribution


# the estimation method is a random-walk Metropolis-Hastings algorithm with coerced acceptance rate

draws <- as.matrix(fit.mcmc$par)

sel <- c("alpha1_1", "alpha2_1")
tmp <- draws[, sel]
par.mle <- fit.ml$par[sel]
par.bay <- apply(tmp, 2, mean)
xlim <- range(c(tmp[,1], par.mle[1]))
ylim <- range(c(tmp[,2], par.mle[2]))
par(mfrow = c(1, 1))
plot(tmp, xlim = xlim, ylim = ylim,)
grid()
par(new = TRUE)
points(par.bay[1], par.bay[2],col = 'green', lwd = 2,
       pch = 15, xlim = xlim, ylim = ylim)
points(par.mle[1], par.mle[2], col = "red", lwd = 2,
       pch = 15, xlim = xlim, ylim = ylim)
#The green square reports the posterior mean while the red
#triangle reports the ML estimate


# for each draw in
# the posterior sample we can compute the unconditional volatility in each regime, to get its
# posterior distribution

## This function computes the unconditional volatility
## for a GJR model with symmeetric disturbances
f_ucvol <- function(par) {
  if (is.vector(par)) {
    par <- matrix(data = par, nrow = 1, dimnames = list(1, names(par)))
  }
  ucvol_1 <- sqrt(250) * par[,"alpha0_1"] / (1 - (par[,"alpha1_1"] + 0.5 * par[,"alpha2_1"] + par[,"beta_1"]))
  ucvol_2 <- sqrt(250) * par[,"alpha0_2"] / (1 - (par[,"alpha1_2"] + 0.5 * par[,"alpha2_2"] + par[,"beta_2"]))
  out <- list(ucvol_1 = ucvol_1, ucvol_2 = ucvol_2)
  return(out)
}

## Compute unconditional volatility
ucvol.draws <- f_ucvol(draws)
ucvol.bay   <- lapply(ucvol.draws, mean)
ucvol.mle   <- f_ucvol(fit.ml$par)

#graph -> posterior distributions of the unconditional annualized volatility in each regime

n <- length(ucvol.draws$ucvol_1)
op <- par(mar = c(2, 2, 2, 2),
          mfrow = c(1, 2),
          oma = c(2, 2, 0.2, 0.2))
hist(ucvol.draws$ucvol_1, nclass = round(10 * log(n)), prob = TRUE,
       xlab = "Volatility (%)", main= '')
title(main = "Regime 1")
lines(density(ucvol.draws$ucvol_1), col = "black", lwd = 2)
rug(ucvol.draws$ucvol_1); box()
points(ucvol.bay$ucvol_1, 0, pch = 15, col = "green", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_1, 0, pch = 17, col = "red", lwd = 2, cex = 3)
hist(ucvol.draws$ucvol_2, nclass = round(10 * log(n)), prob = TRUE,
      xlab = "Volatility (%)",
     ylab = "", cex.lab = 1.5, cex.axis = 1.5, main = "")
rug(ucvol.draws$ucvol_2); box()
points(ucvol.bay$ucvol_2, 0, pch = 15, col = "green", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_2, 0, pch = 17, col = "red", lwd = 2, cex = 3)
title(main = "Regime 2", cex.main = 1.5)
lines(density(ucvol.draws$ucvol_2), col = "black", lwd = 2)
par(op)

# In the low-volatility regime, the distribution is centered
# around 8.9% per annum. For the high-volatility regime, the distribution is centered around
# 32.7% per annum

## Quantiles of unconditional volatility
sapply(ucvol.draws, quantile, probs = c(0.025, 0.975))

# The 95% confidence bands given by the Bayesian approach are 
#        ucvol_1   ucvol_2
# 2.5%   7.98409   23.90536
# 97.5%  10.61048  44.25162

###include parameter uncertainty in the one-step ahead predictive density of MSGARCH models

#solving problems
save(fit.ml,file = 'probl_1/a.txt')
rm(fit.ml)

load("probl_1/a.txt")
fit.ml$spec <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                          distribution.spec = list(distribution = "std"),
                          switch.spec = list(K = 2),
                          constraint.spec = list(regime.const = "nu"))


save(fit.mcmc,file = 'probl_1/b.txt')
rm(fit.mcmc)

load("probl_1/b.txt")
fit.mcmc$spec <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                          distribution.spec = list(distribution = "std"),
                          switch.spec = list(K = 2),
                          constraint.spec = list(regime.const = "nu"))
## Impact of paramter uncertainty in pred
nmesh <- 1000
x <- seq(from = -5, to = 0, length.out = nmesh)
pred.mle <- as.vector(PredPdf(fit.ml, x = x, nahead = 1))
pred.bay <- as.vector(PredPdf(fit.mcmc, x = x, nahead = 1))

pred.draws <- matrix(data = NA, nrow = nrow(draws), ncol = nmesh)
for (i in 1:nrow(draws)) {
  
  save(ms2.gjr.s,file = 'probl_1/c.txt')
  rm(ms2.gjr.s)
  
  load("probl_1/c.txt")
  ms2.gjr.s <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                              distribution.spec = list(distribution = "std"),
                              switch.spec = list(K = 2),
                              constraint.spec = list(regime.const = "nu"))
  
  tmp <- PredPdf(ms2.gjr.s, par = draws[i,], x = x, data = SMI, nahead = 1)
  pred.draws[i,] <- as.vector(tmp)
}
#we want to evaluate the one-step ahead predictive density in the range of values from −5 to 0

xlim <- c(-4, -1.2)
ylim <- c(0, 0.1)
par(mfrow = c(1, 1))
matplot(x, t(pred.draws), xlim = xlim, ylim = ylim,
        type = "l",col = "lightsteelblue", xlab = "Return (%)", ylab = "Predictives",
        lty = 1.5, las = 1, cex.axis = 1.5, cex.lab = 1.5)
title(main = "Left-tail forecast of SMI index return", cex.main = 1.5)
lines(x, pred.bay, xlim = xlim, ylim = ylim,
      type = "l", lty = "solid", col = "green", lwd = 3)
lines(x, pred.mle, xlim = xlim, ylim = ylim,
      type = "l", pch = "o", lty = "dashed", col = "red", lwd = 3)
legend("topleft", c("MCMC draws", "Bayesian","ML"),
       col = c("lightsteelblue", "blue", "red"), lwd = 3,
       lty = c(1, 1, 2), bty = "n", cex = 2)
box()

# The Bayesian predictive density (solid green line) is a particular average 
# of the predictive densities that can be formed with individual posterior 
# MCMC draws (thin solid blue lines). It is generally more
# conservative than the predictive density with plugged ML estimates (dashed red line) and
# offers additional flexibility by accounting for all likely scenarios within the model structure









##backtest analysis of the two-regime model against its single-regime counterpart

## Create GJR-std specification for comparison
gjr.s <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                    distribution.spec = list(distribution = "std"),
                    switch.spec = list(K = 1))

models <- list(gjr.s, ms2.gjr.s)

n.ots    <- 1000 # number of out-of-sample evaluation
n.its    <- 1500 # fit sample size
alpha    <- 0.05 # risk Level
k.update <- 100  # estimation frequency

## Initialization
VaR   <- matrix(NA, nrow = n.ots, ncol = length(models))
y.ots <- matrix(NA, nrow = n.ots, ncol = 1)
model.fit <- vector(mode = "list", length = length(models))

# iterate over out-of-sample time
for (i in 1:n.ots) {
  cat("Backtest - Iteration: ", i, "\n")
  y.its    <- SMI[i:(n.its + i - 1)] # in-sample data
  y.ots[i] <- SMI[n.its + i]         # out-of-sample data
  
  # iterate over models
  for (j in 1:length(models)) {
    
    # do we update the model estimation
    if (k.update == 1 || i %% k.update == 1) {
      cat("Model", j, "is reestimated\n")
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its,
                              ctr = list(do.se = FALSE))
    #  models <- list(gjr.s, ms2.gjr.s)
      
      save(model.fit , file = "probl_1/d.txt")
      rm(model.fit )
      load("probl_1/d.txt")
      models <- list( CreateSpec(variance.spec = list(model = "gjrGARCH"),
                                  distribution.spec = list(distribution = "std"),
                                  switch.spec = list(K = 1)),  CreateSpec(variance.spec = list(model = "sGARCH"),
                                                                          distribution.spec = list(distribution = "norm"),
                                                                          switch.spec = list(K = 2))    )
      model.fit[[j]]$spec  <-as.vector( models[[j]])
    }
    else  {
      
      save(model.fit , file = "probl_1/d.txt")
      rm(model.fit )
      load("probl_1/d.txt")
      models <- list( CreateSpec(variance.spec = list(model = "gjrGARCH"),
                                 distribution.spec = list(distribution = "std"),
                                 switch.spec = list(K = 1)),  CreateSpec(variance.spec = list(model = "sGARCH"),
                                                                         distribution.spec = list(distribution = "norm"),
                                                                         switch.spec = list(K = 2))    )
      model.fit[[j]]$spec  <-as.vector( models[[j]])
      
    }
    # save(model.fit , file = "probl_1/e.txt")
    # rm(model.fit )
    # load("probl_1/e.txt")
    # model.fit[[j]]$spec  <- models[[j]]
    # calculate VaR 1-step ahead
    VaR[i,j] <- Risk(model.fit[[j]]$spec, par = model.fit[[j]]$par,
                     data = y.its,
                     n.ahead = 1,
                     alpha   = alpha,
                     do.es   = FALSE,
                     do.its  = FALSE)$VaR
  }
}

library("GAS")
CC.pval <- DQ.pval <- vector("double", length(models))
for (j in 1:length(models)) {
  test <- GAS::BacktestVaR(data  = y.ots,
                           VaR   = VaR[,j],
                           alpha = alpha)
  
  CC.pval[j] <- test$LRcc[2]
  DQ.pval[j] <- test$DQ$pvalue
}
names(CC.pval) <- names(DQ.pval) <- c("GJR-std", "MS2-GJR-std")

print(CC.pval)
print(DQ.pval)

sink()

library("zoo")
time.index <- zoo::index(SMI)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR, order.by = time.index)

par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, lwd = 1, xlab = "Date (year)",
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
lines(VaR[,1], type = 'l', col = "red", lwd = 3, lty = "dashed")
lines(VaR[,2], type = 'l', col = "blue", lwd = 3)
legend("topleft", legend = c("VaR 5% - GJR-std", "VaR 5% - MS2-GJR-std"),
       col = c("red", "blue"), lwd = 3, cex = 1.5, lty = c("dashed", "solid"))
abline(h = 0)
title("Backtesing VaR at 5% risk level", cex.main = 1.5)


