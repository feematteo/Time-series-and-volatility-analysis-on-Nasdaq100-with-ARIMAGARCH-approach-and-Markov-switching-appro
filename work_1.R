library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(zoo)

nasdaq <- Ad(getSymbols("^IXIC",from= "2002-01-01",to= "2021-11-28"))
chartSeries(nasdaq)

par(mfrow=c(2,1))
acf(log(nasdaq),lag.max = 100)
pacf(log(nasdaq),lag.max = 100)
#From the ACF plot, we observe that the plot decays to zero slowly, 
#meaning the shock affects the process permanently. We can conclude that we need 
#to perform time series analysis on the daily return (log return) of the stock prices.
library(fBasics)


return <- CalculateReturns(nasdaq)
return <- return[-c(1),]
chart_Series(return)
#zero mean   high volatility 
#standard stationarity won't work
print(basicStats(return))
#From the basic statistics of the log return of the stock prices, 
#we observe that the mean is 0 and the distribution of log returns has large kurtosis(fat tails)

chart.Histogram(return, method= c('add.density',"add.normal"),
                colorset = c("blue","red","black"))
legend("topright",legend = c("return","kernel","normal_dist")
       ,fill = c("blue","red","black"))

qqnorm(return)
qqline(return, col = 2)
#seems to be more skewed than the normal distribution 
#no normal distribution
#the series has a somewhat normal distribution with fat tails at both ends.

par(mfrow=c(2,1))
#Time plot of square of log return of prices
plot(return^2,type='l', ylab = "square of stock price return", main="Plot of daily nasdaq price squared return")

#Time plot of absolute value of log return of prices
plot(abs(return),type='l', ylab = "square of stock price return", main="Plot of daily nasdaq price absolute return")

#During the years 2002, 2008-2009, 2013 and 2020, there is spike in volatility indicating non-constant conditional volatility. 
#The high volatility doesn't decrease as fast because the negative shocks have a effect on the process.

#ACF plot of log return of prices
par(mfrow=c(2,1))
acf(return)

#ACF plot of square of log return of prices
acf(return^2)

#The statistics showed that the mean was constant and nearly 0. This is further confirmed by the time series plot. 
#The ACF plot further shows that since, the log stock price returns are not correlated, the mean is constant for the time series.
#However, the squared stock price return values have high correlation. 
#Thus,we may conclude that the log returns process has a strong non-linear dependence.

#  theoretical background mean models ----

#synthetic i.i.d. data
# generate Gaussian synthetic return data
library(mvtnorm)
set.seed(400)
N <- 100
T_max <- 1000
mu <- runif(N)
U <- t(rmvnorm(n = round(0.7*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
X <- rmvnorm(n = T_max, mean = mu, sigma = Sigma)

# now loop over subsets of the samples
error_mu_vs_T <- error_Sigma_vs_T <- NULL
T_sweep <- ceiling(seq(1.01*N, T_max, length.out = 20))
for (T_ in T_sweep) {
  X_ <- X[1:T_, ]
  # sample estimates
  mu_sm <- colMeans(X_)
  Sigma_scm <- cov(X_)
  # compute errors
  error_mu_vs_T    <- c(error_mu_vs_T,    norm(mu_sm     - mu, "2"))
  error_Sigma_vs_T <- c(error_Sigma_vs_T, norm(Sigma_scm - Sigma, "F"))
}
names(error_mu_vs_T) <- names(error_Sigma_vs_T) <- paste("T =", T_sweep)

# plots
plot(T_sweep, error_mu_vs_T, type = "b", pch = 20, col = "blue",
     main = "Error in estimation of mu", xlab = "T", ylab = "error")

plot(T_sweep, error_Sigma_vs_T, type = "b", pch = 20, col = "blue",
     main = "Error in estimation of Sigma", xlab = "T", ylab = "error")

#Univariate ARMA model
library(rugarch)

# specify an AR(1) model with given coefficients and parameters
arma_fixed_spec <- arfimaspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
                              fixed.pars = list(mu = 0.01, ar1 = -0.9, sigma = 0.2))
true_params <- unlist(arma_fixed_spec@model$fixed.pars)
# simulate one path
T <- 2000
set.seed(50)
path_arma <- arfimapath(arma_fixed_spec, n.sim = T)
# convert to xts and plot
synth_log_returns <- xts(path_arma@path$seriesSim, order.by = as.Date("2010-01-01") + 0:(T-1))
plot(synth_log_returns, main = "Synthetic log-returns from ARMA model", lwd = 1.5)
synth_log_prices <- xts(diffinv(synth_log_returns)[-1], order.by = index(synth_log_returns))
plot(synth_log_prices, main = "Synthetic log-prices from ARMA model", lwd = 1.5)


#now estimate the parameters
estim_coeffs_vs_T <- error_coeffs_vs_T <- NULL
T_sweep <- ceiling(seq(100, T, length.out = 20))
for (T_ in T_sweep) {
  arma_fit <- arfimafit(spec = arma_spec, data = synth_log_returns[1:T_])
  estim_coeffs_vs_T <- rbind(estim_coeffs_vs_T, coef(arma_fit))
  error_coeffs_vs_T <- rbind(error_coeffs_vs_T, abs(coef(arma_fit) - true_params)/true_params)
}
rownames(error_coeffs_vs_T) <- rownames(estim_coeffs_vs_T) <- paste("T =", T_sweep)

# plots
matplot(T_sweep, estim_coeffs_vs_T, 
        main = "Estimated ARMA coefficients", xlab = "T", ylab = "value",
        type = "b", pch = 20, col = rainbow(3))
legend("center", inset = 0.01, legend = colnames(estim_coeffs_vs_T), pch = 20, col = rainbow(3))


matplot(T_sweep, 100*error_coeffs_vs_T, 
        main = "Relative error in estimated ARMA coefficients", xlab = "T", ylab = "error (%)",
        type = "b", pch = 20, col = rainbow(3))
legend("topright", inset = 0.01, legend = colnames(error_coeffs_vs_T), pch = 20, col = rainbow(3))

true_params
# mu   ar1  sigma 
#0.01 -0.90  0.20 

#The coefficients seem to be well estimated after T=800 sample, Furthermore the estimation of μ isn't noisy.

#In the previous example, we implicitly assumed we knew the order of the ARMA model. 
#In practice, the order is unknown and one has to try different combinations of orders. 
#The higher the order, the better the fit, but this will inevitable produce overfitting. 
#Many methods have been developed to penalize the increase of the order complexity to avoid overfitting (AIC, BIC, SIC).
arma_fit <- autoarfima(data = synth_log_returns, 
                       ar.max = 3, ma.max = 3, include.mean = TRUE, 
                       criterion = "BIC", method = "partial")
# see the ranking of the combinations
arma_fit$rank.matrix
##choose the best
armaOrder <- arma_fit$rank.matrix[1, c("AR","MA")]

# estimate model excluding the out of sample
out_of_sample <- round(T/2)
dates_out_of_sample <- tail(index(synth_log_returns), out_of_sample)
arma_spec = arfimaspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE))
arma_fit <- arfimafit(spec = arma_spec, data = synth_log_returns, out.sample = out_of_sample)

# forecast log-returns along the whole out-of-sample
arma_fore <- arfimaforecast(arma_fit, n.ahead = 1, n.roll = out_of_sample-1)
forecast_log_returns <- xts(arma_fore@forecast$seriesFor[1, ], dates_out_of_sample)

# recover log-prices
prev_log_price <- head(tail(synth_log_prices, out_of_sample+1), out_of_sample)
forecast_log_prices <- xts(prev_log_price + arma_fore@forecast$seriesFor[1, ], dates_out_of_sample)

par(mfrow=c(2,1))
# plot of log-returns
plot(cbind("fitted"   = fitted(arma_fit),
           "forecast" = forecast_log_returns,
           "original" = synth_log_returns), 
     col = c("blue", "red", "black"), lwd = c(0.5, 0.5, 2),
     main = "Forecast of synthetic log-returns", legend.loc = "topleft")
# plot of log-prices
plot(cbind("forecast" = forecast_log_prices,
           "original" = synth_log_prices), 
     col = c("red", "black"), lwd = c(0.5, 0.5, 2),
     main = "Forecast of synthetic log-prices", legend.loc = "topleft")
# --- --- --- --- --- --- --- --- --- --- --- ---
# prepare training and test data ----

T <- nrow(return)
T_trn <- round(0.7*T)
T_tst <- T - T_trn
return_trn <- return[1:T_trn]
return_tst <- return[-c(1:T_trn)]
{ plot(return, main = "Returns", lwd = 1.5)
  addEventLines(xts("training", index(return[T_trn])), srt=90, pos=2, lwd = 2, col = "blue") }
#Now, we use the training data (i.e., for t=1,…,Ttrn) to fit different models 
#(note that out-of-sample data is excluded by indicating out.sample = T_tst).
library(rugarch)

# fit AR(1) model
ar_spec <- arfimaspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE))
ar_fit <- arfimafit(spec = ar_spec, data = return, out.sample = T_tst)
coef(ar_fit)

# fit ARMA(2,2) model
arma_spec <- arfimaspec(mean.model = list(armaOrder = c(2,2), include.mean = TRUE))
arma_fit <- arfimafit(spec = arma_spec, data = return, out.sample = T_tst)
coef(arma_fit)

# fit ARMA(2,1) model
arma_spec_2 <- arfimaspec(mean.model = list(armaOrder = c(2,1), include.mean = TRUE))
arma_fit_2 <- arfimafit(spec = arma_spec_2, data = return, out.sample = T_tst)
coef(arma_fit_2)

# fit MA(1) model
ma_spec <- arfimaspec(mean.model = list(armaOrder = c(0,1), include.mean = TRUE))
ma_fit <- arfimafit(spec = ma_spec, data = return, out.sample = T_tst)
coef(ma_fit)

#use the different models to forecast the log-returns:
dates_out_of_sample <- tail(index(return), T_tst)

# forecast with AR(1) model
ar_fore_return <- xts(arfimaforecast(ar_fit, n.ahead = 1, n.roll = T_tst - 1)@forecast$seriesFor[1, ], 
                          dates_out_of_sample)

# forecast with ARMA(2,2) model
arma_fore_return <- xts(arfimaforecast(arma_fit, n.ahead = 1, n.roll = T_tst - 1)@forecast$seriesFor[1, ], 
                            dates_out_of_sample)

# forecast with ARMA(2,2) model
arma_2_fore_return <- xts(arfimaforecast(arma_fit_2, n.ahead = 1, n.roll = T_tst - 1)@forecast$seriesFor[1, ], 
                        dates_out_of_sample)

# forecast with MA(1) model
ma_fore_return <- xts(arfimaforecast(ma_fit, n.ahead = 1, n.roll = T_tst - 1)@forecast$seriesFor[1, ], 
                      dates_out_of_sample)

#compute the forecast errors
error_var <- rbind( "AR(1)"                  = c(var(return - fitted(ar_fit)),
                                                var(return - ar_fore_return)),
                   "ARMA(2,2)"              = c(var(return - fitted(arma_fit)),
                                                var(return - arma_fore_return)),
                   "ARMA(2,1)"              = c(var(return - fitted(arma_fit_2)),
                                                var(return - arma_2_fore_return)),
                   "MA(1,0)"              = c(var(return - fitted(ma_fit)),
                                                var(return - ma_fore_return)))
colnames(error_var) <- c("in-sample", "out-of-sample")
print(error_var)

#           in-sample    out-of-sample
#AR(1)     0.0001984002  0.0001678625
#ARMA(2,2) 0.0001978884  0.0001698360
#ARMA(2,1) 0.0001980006  0.0001695621
#MA(1,0)   0.0001983303  0.0001677420
#The important quantity is the out-of-sample error: we can see that increasing the model complexity may give few results. 
#It seems that the simplest AR or MA model is good enough in terms of error of forecast returns.

error_return <- cbind(return - ar_fore_return,
                          return - arma_fore_return,
                          return - arma_2_fore_return,
                          return - ma_fore_return)

names(error_return) <- c("AR", "ARMA", "ARMA2", 'MA')
plot(error_return, col = c("red", "green", "blue",'purple'), lwd = c(0.2, 2, 3, 4),
     main = "Out-of-sample error of static return forecast for different models", 
     legend.loc = "bottomleft")
#since we are not refitting the models, the error tends to become worse as the time advances

#Rolling-window comparison ----

#concept of static forecast vs rolling forecast

# model specification of an ARMA(2,2)
specx <- arfimaspec(mean.model = list(armaOrder = c(2,2), include.mean = TRUE))

# static fit and forecast
ar_static_fitx <- arfimafit(spec = specx, data = return, out.sample = T_tst)
ar_static_fore_logreturnsx <- xts(arfimaforecast(ar_static_fitx, n.ahead = 1, n.roll = T_tst - 1)@forecast$seriesFor[1, ],
                                 dates_out_of_sample)

# rolling fit and forecast
modelrollx <- arfimaroll(spec = specx, data = return, n.ahead = 1, 
                        forecast.length = T_tst, refit.every = 50, refit.window = "moving")
ar_rolling_fore_logreturnsx <- xts(modelrollx@forecast$density$Mu, dates_out_of_sample)

par(mfrow=c(2,1))
# plot of forecast
plot(cbind("static forecast"  = ar_static_fore_logreturnsx,
           "rolling forecast" = ar_rolling_fore_logreturnsx),
     col = c("black", "red"), lwd = 2,
     main = "Forecast with ARMA(2,2) model", legend.loc = "topleft")


# plot of forecast error
error_logreturnsx <- cbind(return - ar_static_fore_logreturnsx,
                          return - ar_rolling_fore_logreturnsx)
names(error_logreturnsx) <- c("rolling forecast", "static forecast")
plot(error_logreturnsx, col = c("black", "red"), lwd = 2,
     main = "Forecast error with ARMA(2,2) model", legend.loc = "topleft")
#We can observe the effect of the rolling-window process in keeping track with the time series is very smooth.



# rolling forecast with AR(1) model
ar_rolling_fore_return <- xts(arfimaroll(ar_spec, data = return, n.ahead = 1, forecast.length = T_tst, 
                                             refit.every = 50, refit.window = "moving")@forecast$density$Mu, 
                                  dates_out_of_sample)

# rolling forecast with ARMA(2,2) model
arma_rolling_fore_return <- xts(arfimaroll(arma_spec, data = return, n.ahead = 1, forecast.length = T_tst, 
                                               refit.every = 50, refit.window = "moving")@forecast$density$Mu, 
                                    dates_out_of_sample)

# rolling forecast with ARMA(2,1) model
arma_2_rolling_fore_return <- xts(arfimaroll(arma_spec, data = return, n.ahead = 1, forecast.length = T_tst, 
                                           refit.every = 50, refit.window = "moving")@forecast$density$Mu, 
                                dates_out_of_sample)

# rolling forecast with MA(1) model
ma_rolling_fore_return <- xts(arfimaroll(ma_spec, data = return, n.ahead = 1, forecast.length = T_tst, 
                                         refit.every = 50, refit.window = "moving")@forecast$density$Mu, 
                              dates_out_of_sample) 

rolling_error_var <- rbind("AR(1)"                  = c(var(return - fitted(ar_fit)),
                                                        var(return - ar_rolling_fore_return)),
                           "ARMA(2,2)"              = c(var(return - fitted(arma_fit)),
                                                        var(return - arma_rolling_fore_return)),
                           "ARMA(2,1)"              = c(var(return - fitted(arma_fit)),
                                                        var(return - arma_2_rolling_fore_return)),
                           "MA(1,0)"              = c(var(return - fitted(ma_fit)),
                                                        var(return - ma_rolling_fore_return)))

colnames(rolling_error_var) <- c("in-sample", "out-of-sample")
print(rolling_error_var)

#          in-sample     out-of-sample
#AR(1)     0.0001984002  0.0001675942
#ARMA(2,2) 0.0001978884  0.0001701219
#ARMA(2,1) 0.0001978884  0.0001701219
#MA(1,0)   0.0001983303  0.0001674599


error_return <- cbind(return - ar_rolling_fore_return,
                          return - arma_rolling_fore_return,
                          return - arma_2_rolling_fore_return,
                      return - ma_rolling_fore_return)
names(error_return) <- c("AR", "ARMA", "ARMA2", 'MA')
plot(error_return, col = c("red", "green", "blue",'purple'), lwd = c(0.5, 2, 3,4),
     main = "Error of rolling forecast for different models", legend.loc = "topleft")

#We see that now all the models track the time series. 
#Also, we don’t observe any significant difference among the models

#compare the static vs rolling basis errors:
barplot(rbind(error_var[, "out-of-sample"], rolling_error_var[, "out-of-sample"]), 
        col = c("darkblue", "darkred"), 
        legend = c("static forecast", "rolling forecast"), 
        main = "Out-of-sample forecast error for different models", 
        xlab = "method", ylab = "variance", beside = TRUE)
#a rolling refitting not only does not hurt but in some cases is better
#even though here there is no significant difference

# --- --- --- --- --- --- --- --- --- --- 
# ARIMA fit ----
library(SBAGM) #lower AIC
ARIMAAIC(data = return, p=4, q=4, d=0, season=list(order=c(0,0,0))) 

library(forecast)
auto.arima(y = return)

arima_mod <- arima(return,order=c(1,0,0))
# Call:
#   arima(x = return, order = c(1, 0, 0))
# 
# Coefficients:
#     ar1    intercept
#    -0.1005   5e-04
#s.e. 0.0141   2e-04
# 
# sigma^2 estimated as 0.0001889:  log likelihood = 14372.2,  aic = -28738.4

res_arima=arima_mod$res
par(mfcol=c(3,1))
plot(res_arima,main='Residuals')
acf_res=acf(res_arima,main='ACF
Residuals',lag.max=100,ylim=c(-0.5,1))
pacfres=pacf(res_arima,main='PACF
Residuals',lag.max=100,ylim=c(-0.5,1))
#The residual plot, ACF and PACF do not have any significant lag, indicating ARIMA(1,1,0) is a
#good model to represent the series   
Box.test(res_arima, lag = c(2), type = c( "Ljung-Box"), fitdf = 1)             
Box.test(res_arima, lag = c(5), type = c( "Ljung-Box"), fitdf = 1)                 
Box.test(res_arima, lag = c(10), type = c( "Ljung-Box"), fitdf = 1)   
Box.test(res_arima, lag = c(25), type = c( "Ljung-Box"), fitdf = 1) 

res_arima2=res_arima^2
par(mfcol=c(3,1))
plot(res_arima2,main='Residuals Squared')
acf_res2=acf(res_arima2,main='ACF Squared
Residuals',lag.max=100,ylim=c(-0.5,1))
pacfres2=pacf(res_arima2,main='PACF Squared
Residuals',lag.max=100,ylim=c(-0.5,1)) #Squared residuals plot shows cluster of volatility 
#PACF and ACF cuts off after lag 10 and 20 even though some remaining lags are significants


chart.RollingPerformance(R= return["2002::2021"],width = 22
                         ,FUN = "sd.annualized",scale=252,
                         main = "NASDAQ's monthly volatility")
#needs stochastic model for conditional volatility

#Variance models ----

#theoretical background ----
garch_fixed_spec <- ugarchspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
                               variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                               fixed.pars = list(mu = 0.005, ar1 = -0.9, 
                                                 omega = 0.001, alpha1 = 0.3, beta1 = 0.65))
true_params <- unlist(garch_fixed_spec@model$fixed.pars)
true_params
#  mu    ar1  omega  alpha1  beta1 
#0.005 -0.900  0.001  0.300  0.650 
T <- 2000
set.seed(50)
path_garch <- ugarchpath(garch_fixed_spec, n.sim = T)

synth_log_returns <- xts(path_garch@path$seriesSim, order.by = as.Date("2010-01-01") + 0:(T-1))
synth_volatility <- xts(path_garch@path$sigmaSim, order.by = as.Date("2010-01-01") + 0:(T-1))
{ plot(synth_log_returns, main = "Synthetic log-returns from GARCH model", lwd = 1.5)
  lines(synth_volatility, col = "red", lwd = 2) }

# specify a GARCH model
garch_spec <- ugarchspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
                         variance.model = list(model = "sGARCH", garchOrder = c(1,1)))
estim_coeffs_vs_T <- error_coeffs_vs_T <- NULL
T_sweep <- ceiling(seq(100, T, length.out = 20))
for (T_ in T_sweep) {
  garch_fit <- ugarchfit(spec = garch_spec, data = synth_log_returns[1:T_])
  error_coeffs_vs_T <- rbind(error_coeffs_vs_T, abs((coef(garch_fit) - true_params)/true_params))
  estim_coeffs_vs_T <- rbind(estim_coeffs_vs_T, coef(garch_fit))
}
rownames(error_coeffs_vs_T) <- rownames(estim_coeffs_vs_T) <- paste("T =", T_sweep)

# plots

matplot(T_sweep, estim_coeffs_vs_T, 
        main = "Estimated GARCH coefficients", xlab = "T", ylab = "value",
        type = "b", pch = 20, col = rainbow(5))
legend("center", inset = 0.01, legend = colnames(estim_coeffs_vs_T), pch = 20, col = rainbow(5))



matplot(T_sweep, 100*error_coeffs_vs_T, 
        main = "Relative error in estimated GARCH coefficients", xlab = "T", ylab = "error (%)",
        type = "b", pch = 20, col = rainbow(5), ylim=c(-10,250))
legend("topright", inset = 0.01, legend = colnames(error_coeffs_vs_T), pch = 20, col = rainbow(5))

#the coefficients seem to be well estimated after T=800 samples.

#forecast
# estimate model excluding the out-of-sample
out_of_sample <- round(T/2)
dates_out_of_sample <- tail(index(synth_log_returns), out_of_sample)
garch_spec <- ugarchspec(mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
                         variance.model = list(model = "sGARCH", garchOrder = c(1,1)))
garch_fit <- ugarchfit(spec = garch_spec, data = synth_log_returns, out.sample = out_of_sample)

# forecast log-returns along the whole out-of-sample
garch_fore <- ugarchforecast(garch_fit, n.ahead = 1, n.roll = out_of_sample-1)
forecast_log_returns <- xts(garch_fore@forecast$seriesFor[1, ], dates_out_of_sample)
forecast_volatility <- xts(garch_fore@forecast$sigmaFor[1, ], dates_out_of_sample)

par(mfrow=c(2,1))
# plot of log-returns
plot(cbind("fitted"   = fitted(garch_fit),
           "forecast" = forecast_log_returns,
           "original" = synth_log_returns), 
     col = c("blue", "red", "black"), lwd = c(0.5, 0.5, 2),
     main = "Forecast of synthetic log-returns", legend.loc = "topleft")
# plot of volatility log-returns
plot(cbind("fitted volatility"   = sigma(garch_fit),
           "forecast volatility" = forecast_volatility,
           "log-returns"         = synth_log_returns), 
     col = c("blue", "red", "black"), lwd = c(2, 2, 1),
     main = "Forecast of volatility of synthetic log-returns", legend.loc = "topleft")
#--- --- --- --- --- --- --- --- --- --- --- --- --- 

#differents methods ----
#costant
var_constant <- var(return_trn) 
par(mfrow=c(2,1))
plot(cbind(sqrt(var_constant), return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Constant")
#install.packages('RcppRoll')
library(RcppRoll)

#MA
lookback_var <- 20
var_t <- roll_meanr(return_trn^2, n = lookback_var, fill = NA)
plot(cbind(sqrt(var_t), return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on simple rolling means of squares (lookback=20)")
return_trn_std <- return_trn/sqrt(var_t)
var_ma <- var(return_trn_std, na.rm = TRUE) * tail(var_t, 1)

#EWMA
library(forecast)
fit_ets <- ets(return_trn^2, model = "ANN")
std_t <- as.numeric(sqrt(fit_ets$fitted))
plot(cbind(std_t, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on EWMA of squares")
return_trn_std_2 <- return_trn/std_t
var_ewma <- var(return_trn_std_2, na.rm = TRUE) * tail(std_t, 1)^2

#Multiplicative ETS
fit_ets2 <- ets(1e-6 + return_trn^2, model = "MNN")
std_t2 <- as.numeric(sqrt(fit_ets2$fitted))
plot(cbind(std_t, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on ETS(M,N,N) of squares")
return_trn_std_3 <- return_trn/std_t2
var_ets_mnn <- var(return_trn_std_3, na.rm = TRUE) * tail(std_t2, 1)^2


appgarch(data = return, methods = c('sGARCH', 'eGARCH', 'gjrGARCH'),distributions =  c("norm", "std", "snorm",'ged'), 
         aorder = c(1,0), gorder = c(1, 1), stepahead = 2 ) #eGARCH - snorm    rmse_mean 0.01620853

appgarch(data = return, methods = c('sGARCH', 'eGARCH', 'gjrGARCH'),distributions =  c("norm", "std", "snorm",'ged'), 
         aorder = c(1,0), gorder = c(2, 1), stepahead = 2 ) #eGARCH - snorm    rmse_mean 0.01621036

appgarch(data = return, methods = c('sGARCH', 'eGARCH', 'gjrGARCH'),distributions =  c("norm", "std", "snorm",'ged'), 
         aorder = c(1,0), gorder = c(3, 2), stepahead = 2 ) #eGARCH - snorm    rmse_mean 0.01620979


#ARCH
library(fGarch)
arch_fit <- fGarch::garchFit(formula = ~ garch(5,0), return_trn, trace = FALSE)
std_t4 <- arch_fit@sigma.t
plot(cbind(std_t4, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on ARCH(5)")
var_arch <- tail(std_t4, 1)^2

#GARCH
garch_fit <- fGarch::garchFit(formula = ~ garch(1,1), return_trn, trace = FALSE)
std_t5 <- garch_fit@sigma.t
plot(cbind(std_t5, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on GARCH(1,1)")
var_garch <- tail(std_t5, 1)^2

# garch_fit2 <- fGarch::garchFit(formula = ~ garch(2,1), return_trn, trace = FALSE,cond.dist = 'snorm')
# std_t6 <- garch_fit2@sigma.t
# plot(cbind(std_t6, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
#      main = "Envelope based on GARCH(2,1)S")
# var_garch2 <- tail(std_t6, 1)^2
# 
# garch_fit3 <- fGarch::garchFit(formula = ~ garch(3,2), return_trn, trace = FALSE,cond.dist = 'snorm')
# std_t7 <- garch_fit3@sigma.t
# plot(cbind(std_t7, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
#      main = "Envelope based on GARCH(3,2)S")
# var_garch3 <- tail(std_t7, 1)^2

#SV
#install.packages("stochvol")
library(stochvol)
res_sv <- svsample(return_trn - mean(return_trn))
std_t_sv <- res_sv$summary$sd[, 1]
plot(cbind(std_t_sv, return_trn), col = c("red", "black"), lwd = c(2.5, 1.5),
     main = "Envelope based on stochastic volatility")
var_sv <- tail(std_t_sv, 1)^2

#compare the error in the estimation of the variance by each method for the out-of-sample period 
error_all <- c("MA"         = abs(var_ma      - var(return_tst)),
               "EWMA"       = abs(var_ewma    - var(return_tst)),
               "ETS(M,N,N)" = abs(var_ets_mnn - var(return_tst)),
               "ARCH(5)"    = abs(var_arch    - var(return_tst)),
               "GARCH(1,1)" = abs(var_garch   - var(return_tst)),
               "SV"         = abs(var_sv      - var(return_tst)))

print(error_all)
barplot(error_all, main = "Error in estimation of out-of-sample variance", col = rainbow(8))

#Rolling-window comparison----
error_sv <- error_garch <- error_garch2 <- error_arch <- error_ets_mnn <- error_ewma <- error_ma <- NULL

# rolling window
#install.packages("roll")
#library(roll)
lookback_ <- 200
len_tst <- 40
for (i in seq(lookback_, T-len_tst, by = len_tst)) {
  return_trn <- return[(i-lookback_+1):i]
  return_tst <- return[(i+1):(i+len_tst)]
  var_tst <- var(return_tst)
  
  # MA
  var_t <- roll_meanr(return_trn^2,n =20,fill = NA )
  var_fore <- var(return_trn/sqrt(var_t), na.rm = TRUE) * tail(var_t, 1)
  error_ma <- c(error_ma, abs(var_fore - var_tst))
  
  # EWMA
  fit_ets <- ets(return_trn^2, model = "ANN")
  std_t <- as.numeric(sqrt(fit_ets$fitted))
  var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  error_ewma <- c(error_ewma, abs(var_fore - var_tst))
  
  # ETS(M,N,N)
  fit_ets <- ets(1e-6 + return_trn^2, model = "MNN")
  std_t <- as.numeric(sqrt(fit_ets$fitted))
  var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  error_ets_mnn <- c(error_ets_mnn, abs(var_fore - var_tst))
  
  # ARCH
  arch_fit <- fGarch::garchFit(formula = ~ garch(5,0), return_trn, trace = FALSE)
  std_t <- as.numeric(arch_fit@sigma.t)
  var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  error_arch <- c(error_arch, abs(var_fore - var_tst))
  
  # GARCH
  garch_fit <- fGarch::garchFit(formula = ~ garch(1,1), return_trn, trace = FALSE)
  std_t <- as.numeric(garch_fit@sigma.t)
  var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  error_garch <- c(error_garch, abs(var_fore - var_tst))
  
  # GARCH 3,2 S
  #garch_fit2 <- fGarch::garchFit(formula = ~ garch(3,2), return_trn,cond.dist = 'snorm', trace = FALSE)
  #std_t <- as.numeric(garch_fit2@sigma.t)
  #var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  #error_garch2 <- c(error_garch2, abs(var_fore - var_tst))
  
  # SV
  res <- svsample(return_trn - mean(return_trn))
  std_t <- res$summary$sd[, 1]
  var_fore <- var(return_trn/std_t, na.rm = TRUE) * tail(std_t, 1)^2
  error_sv <- c(error_sv, abs(var_fore - var_tst))
}
error_all <- c("MA"         = mean(error_ma),
               "EWMA"       = mean(error_ewma),
               "ETS(M,N,N)" = mean(error_ets_mnn),
               "ARCH(5)"    = mean(error_arch),
               "GARCH(1,1)" = mean(error_garch),
               "SV"         = mean(error_sv))
print(error_all)
barplot(error_all, main = "Error in estimation of variance", col = rainbow(8))


#Multivariate GARCH ----
library(quantmod)
stock_namelist <- c("QQQ", "SPY", "AAPL")

# download data from YahooFinance
prices <- xts()
for (stock_index in 1:length(stock_namelist))
  prices <- cbind(prices, Ad(getSymbols(stock_namelist[stock_index], 
                                        from= "2002-07-01",to= "2021-11-28", auto.assign = FALSE)))
colnames(prices) <- stock_namelist
indexClass(prices) <- "Date"
logreturns <- diff(log(prices))[-1]

# plot the four series of log-prices
plot(log(prices), col = c("black", "blue", "green",'pink'),
     main = "Log-prices of the two ETFs + AAPL", legend.loc = "bottomright")

##define the model
library(rmgarch)

# specify i.i.d. model for the univariate time series
ugarch_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), 
              include.mean = FALSE), variance.model = list(model = "sGARCH", 
                                          garchOrder = c(1,1)))

# specify DCC model
dcc_spec <- dccspec(uspec = multispec(replicate(ugarch_spec, n = 3)),
                    VAR = TRUE, lag = 3,
                    model = "DCC", dccOrder = c(1,1))
# estimate model
garchdcc_fit <- dccfit(dcc_spec, data = logreturns, solver = "nlminb")

# extract time-varying covariance and correlation matrix
dcc_cor <- rcor(garchdcc_fit)
dim(dcc_cor) #3    3 4887
#plot the time-varying correlations
corr_t <- xts(cbind(dcc_cor[1, 2, ], dcc_cor[1, 3, ], dcc_cor[2, 3, ]), order.by = index(logreturns))
colnames(corr_t) <- c("QQQ vs SPY", "QQQ vs AAPL", "SPY vs AAPL")
plot(corr_t, col = c("black", "red", "blue"),
     main = "Time-varying correlations", legend.loc = "bottomright")

#We see the correlation between the two ETFs is extremely high and quite stable. 
#The correlation between AAPL and the ETFs is cyclical, shifting between 0.8 and 0.4


#model specification ----
model_specify <- ugarchspec(mean.model = list(armaOrder=c(1,0), include.mean = TRUE),
                            variance.model = list(model="sGARCH",
                                                  garchOrder=c(1,1)),
                            distribution.model = "snorm")
model_fit <- ugarchfit(data = return,spec = model_specify)
par(mfrow = c(3, 4))
for (i in 1:12) {
  plot(model_fit, which = i)
}


model_specify2 <- ugarchspec(mean.model = list(armaOrder=c(1,0), include.mean = TRUE),
                            variance.model = list(model="eGARCH",
                                                  garchOrder=c(1,1)),
                            distribution.model = "snorm")
model_fit2 <- ugarchfit(data = return,spec = model_specify2)
par(mfrow = c(3, 4))
for (i in 1:12) {
  plot(model_fit2, which = i)
}

model_specify3 <- ugarchspec(mean.model = list(armaOrder=c(1,0), include.mean = TRUE),
                            variance.model = list(model="eGARCH",
                                                  garchOrder=c(1,1)),
                            distribution.model = "sstd")
model_fit3 <- ugarchfit(data = return,spec = model_specify3)
par(mfrow = c(3, 4))
for (i in 1:12) {
  plot(model_fit3, which = i)
}


#MSGARCH ----

library(MSGARCH)
appmsgarch(data =return, methods = c( "sGARCH", "eGARCH", "gjrGARCH", "tGARCH"),
           distributions = c("norm", "snorm", "std",'ged'),
           stepahead = 2 )  #package SBAGM
# tGARCH-tGARCH   ged-ged , rmse and mae

mrs_spec <-CreateSpec(variance.spec = list(model=c("tGARCH")),distribution.spec = list(distribution=c("ged")),
                      switch.spec = list(do.mix = FALSE, K = 2),
                      constraint.spec = list(regime.const = "nu") )
summary(mrs_spec)

fitml <- FitML(spec = mrs_spec,data = return) #Maximum Likelihood estimation.

## Unconditional vol
sqrt(250) * sapply(ExtractStateFit(fitml), UncVol)
# 0.3001438 0.4693295


## Smoothed probabilities in regime 2 and volatility
smoothed.prob <- State(fitml)$SmoothProb[, 1, 2, drop = TRUE]
vol <- sqrt(250) * Volatility(fitml)


op <- par(mfrow = c(2,1),
          oma = c(1,1,0,0) + 0.0,
          mar = c(2,2,2,2) + 0.0)
plot(as.vector(return), axes = FALSE, ann = FALSE)
par(new = TRUE)
ylabel <- expression(paste("Pr(", s[t], " = 2 | ", hat(psi), ", ", I[t], ")"))
plot(zoo::zoo(smoothed.prob, order.by = zoo::index(return)), plot.type = "single")
title(main = "Smoothed probabilities")
plot(zoo::zoo(vol, order.by = zoo::index(return)), plot.type = "single")
title(main = "Volatility (%)")
par(op)

## MCMC estimation
nmcmc <- 12500
nburn <- 5000
nthin <- 5
ctr <- list(nmcmc = nmcmc, nburn = nburn,
            nthin = nthin, par0 = fitml$par)
fit.mcmc <- FitMCMC(mrs_spec, data = return, ctr = ctr)# Markov Chain Monte Carlo / Bayesian estimation
summary(fit.mcmc)



draws <- as.matrix(fit.mcmc$par)

sel <- c("alpha1_1", "alpha2_1")
tmp <- draws[, sel]
par.mle <- fitml$par[sel]
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
# reports the ML estimate


## This function computes the unconditional volatility
##
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
ucvol.mle   <- f_ucvol(fitml$par)

#graph -> posterior distributions of the unconditional annualized volatility in each regime

n <- length(ucvol.draws$ucvol_1)
op <- par(mar = c(2, 2, 2, 2),
          mfrow = c(1, 2),
          oma = c(2, 2, 0.2, 0.2))
hist(ucvol.draws$ucvol_1, nclass = round(10 * log(n)), prob = TRUE,
     xlab = "Volatility (%)", main= '',ylim = range(0,0.08),xlim = range(-30,30),breaks = 1000)
title(main = "Regime 1")
lines(density(ucvol.draws$ucvol_1), col = "black", lwd = 2)
rug(ucvol.draws$ucvol_1); box()
points(ucvol.bay$ucvol_1, 0, pch = 15, col = "green", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_1, 0, pch = 17, col = "red", lwd = 2, cex = 3)
hist(ucvol.draws$ucvol_2, nclass = round(10 * log(n)), prob = TRUE,
     xlab = "Volatility (%)",
     ylab = "", cex.lab = 1.5, cex.axis = 1.5, main = "",ylim = range(0,10),xlim = range(-1.5,1),breaks = 50)
rug(ucvol.draws$ucvol_2); box()
points(ucvol.bay$ucvol_2, 0, pch = 15, col = "green", lwd = 2, cex = 3)
points(ucvol.mle$ucvol_2, 0, pch = 17, col = "red", lwd = 2, cex = 3)
title(main = "Regime 2", cex.main = 1.5)
lines(density(ucvol.draws$ucvol_2), col = "black", lwd = 2)
par(op)




## Quantiles of unconditional volatility
sapply(ucvol.draws, quantile, probs = c(0.025, 0.975))
#           ucvol_1   ucvol_2
# 2.5%  -0.6255203   -8121.995
# 97.5%  3.9397738   11331.303


###include parameter uncertainty in the one-step ahead predictive density of MSGARCH models

#solving problems
save(fitml,file = 'probl_1/aa.txt')
rm(fitml)

load("probl_1/aa.txt")
fitml$spec <- CreateSpec(variance.spec = list(model=c("tGARCH")),distribution.spec = list(distribution=c("ged")),
                         switch.spec = list(do.mix = FALSE, K = 2),
                         constraint.spec = list(regime.const = "nu") )


save(fit.mcmc,file = 'probl_1/bb.txt')
rm(fit.mcmc)

load("probl_1/bb.txt")
fit.mcmc$spec <- CreateSpec(variance.spec = list(model=c("tGARCH")),distribution.spec = list(distribution=c("ged")),
                            switch.spec = list(do.mix = FALSE, K = 2),
                            constraint.spec = list(regime.const = "nu") )
## Impact of paramter uncertainty in pred
nmesh <- 1000
x <- seq(from = -5, to = 0, length.out = nmesh)
pred.mle <- as.vector(PredPdf(fitml, x = x, nahead = 1))
pred.bay <- as.vector(PredPdf(fit.mcmc, x = x, nahead = 1))

pred.draws <- matrix(data = NA, nrow = nrow(draws), ncol = nmesh)
for (i in 1:nrow(draws)) {
  
  save(mrs_spec,file = 'probl_1/cc.txt')
  rm(mrs_spec)
  
  load("probl_1/cc.txt")
  mrs_spec <- CreateSpec(variance.spec = list(model=c("tGARCH")),distribution.spec = list(distribution=c("ged")),
                         switch.spec = list(do.mix = FALSE, K = 2),
                         constraint.spec = list(regime.const = "nu") )
  
  tmp <- PredPdf(mrs_spec, par = draws[i,], x = x, data = return, nahead = 1)
  pred.draws[i,] <- as.vector(tmp)
}

#we want to evaluate the one-step ahead predictive density in the range of values from −5 to 0

xlim <- c(-0.4, 0.05)
ylim <- c(0, 0.3)
par(mfrow = c(1, 1))
matplot(x, t(pred.draws), xlim = xlim, ylim = ylim,
        type = "l",col = "lightsteelblue", xlab = "Return (%)", ylab = "Predictives",
        lty = 1.5, las = 1, cex.axis = 1.5, cex.lab = 1.5)
title(main = "Left-tail forecast of Nasdaq return", cex.main = 1.5)
lines(x, pred.bay, xlim = xlim, ylim = ylim,
      type = "l", lty = "solid", col = "green", lwd = 3)
lines(x, pred.mle, xlim = xlim, ylim = ylim,
      type = "l", pch = "o", lty = "dashed", col = "red", lwd = 3)
legend("topleft", c("MCMC draws", "Bayesian","ML"),
       col = c("lightsteelblue", "green", "red"), lwd = 3,
       lty = c(1, 1, 2), bty = "n", cex = 2)
box()

##backtest analysis of the two-regime model against its single-regime counterpart

## Create tgrach-ged specification for comparison
mrs_spec2 <- CreateSpec(variance.spec = list(model = "tGARCH"),
                    distribution.spec = list(distribution = "ged"),
                    switch.spec = list(K = 1))

models <- list(mrs_spec2, mrs_spec)

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
  y.its    <- return[i:(n.its + i - 1)] # in-sample data
  y.ots[i] <- return[n.its + i]         # out-of-sample data
  
  # iterate over models
  for (j in 1:length(models)) {
    
    # do we update the model estimation
    if (k.update == 1 || i %% k.update == 1) {
      cat("Model", j, "is reestimated\n")
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its,
                              ctr = list(do.se = FALSE))
      
      save(model.fit , file = "probl_1/dd.txt")
      rm(model.fit )
      load("probl_1/dd.txt")
      models <- list( CreateSpec(variance.spec = list(model = "tGARCH"),
                                 distribution.spec = list(distribution = "ged"),
                                 switch.spec = list(K = 1)),  CreateSpec(variance.spec = list(model = "tGARCH"),
                                                                         distribution.spec = list(distribution = "ged"),
                                                                         switch.spec = list(K = 2))    )
      model.fit[[j]]$spec  <-as.vector( models[[j]])
    }
    else  {
      
      save(model.fit , file = "probl_1/dd.txt")
      rm(model.fit )
      load("probl_1/dd.txt")
      models <- list( CreateSpec(variance.spec = list(model = "tGARCH"),
                                 distribution.spec = list(distribution = "ged"),
                                 switch.spec = list(K = 1)),  CreateSpec(variance.spec = list(model = "tGARCH"),
                                                                         distribution.spec = list(distribution = "ged"),
                                                                         switch.spec = list(K = 2))    )
      model.fit[[j]]$spec  <-as.vector( models[[j]])
      
    }
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
names(CC.pval) <- names(DQ.pval) <- c("Tgarch-ged", "MS2-Tgarch-ged")

print(CC.pval)
print(DQ.pval)

#sink()


library("zoo")
time.index <- zoo::index(return)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR, order.by = time.index)

par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, lwd = 1, xlab = "Date (year)",
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
lines(VaR[,1], type = 'l', col = "red", lwd = 3, lty = "dashed")
lines(VaR[,2], type = 'l', col = "blue", lwd = 3)
legend("topleft", legend = c("VaR 5% - Tgarch-ged", "VaR 5% - MS2-Tgarch-ged"),
       col = c("red", "blue"), lwd = 3, cex = 1.5, lty = c("dashed", "solid"))
abline(h = 0)
title("Backtesing VaR at 5% risk level", cex.main = 1.5)







#####

pred <- predict(firml, nahead = 5, do.return.draw = TRUE)
pred$vol
risk <- Risk(firml, alpha = c(0.01, 0.05), nahead = 5)
#the five-step ahead VaR and ES
#measures at the 1% and 5% risk levels are computed

smooth_prob <- State(firml)$SmoothProb[, 1, 2, drop = TRUE]


plot(smooth_prob,type='l')
vol <- sqrt(500) * Volatility(firml)
plot(vol)

# xx <- seq(from = -5, to = 0, length.out = 1000)
# pred.ml <- as.vector(PredPdf(object = firml,x = xx,nahead = 1))
# 
# #install.packages('tidyverse')
# library(tidyverse)
# ggplot(try2)
### ### #

#Pairs Trading ----

library(egcm)  
egcm.set.default.pvalue(0.01)
SPY_prices <- Ad(getSymbols("SPY", from = "2013-01-01", to = "2013-12-31", auto.assign = FALSE))
IVV_prices <- Ad(getSymbols("IVV", from = "2013-01-01", to = "2013-12-31", auto.assign = FALSE))
plot(cbind(SPY_prices, IVV_prices), legend.loc = "topleft", main = "ETF prices")

resx <- egcm(SPY_prices, IVV_prices)
summary(resx)
plot(resx)



