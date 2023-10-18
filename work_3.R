library(forecast)
#install.packages("expsmooth")
library(expsmooth)


fit1 <- ets(bonds)
etsfit <- ets(usnetelec)
fit3 <- ets(ukcars)
fit4 <- ets(visitors)

par(mfrow=c(2,2))
plot(forecast(fit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")
plot(forecast(etsfit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")
plot(forecast(fit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")
plot(forecast(fit4),xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")

# Four time series showing point forecasts and 80% & 95% prediction intervals 
#obtained using exponential smoothing state space models
print(etsfit)
accuracy(etsfit)

# Fit an ETS model to the first 45 months
fit <- ets(usnetelec[1:45])
# Apply the same model to the last 10 months without re-estimating the parameters
test <- ets(usnetelec[46:55],model=fit)
# Look at the in-sample accuracy
accuracy(test)
# Look at the out-of-sample accuracy
accuracy(forecast(fit,10), usnetelec[46:55])

#ARIMA
aafit1 <- auto.arima(bonds,max.P=0,max.Q=0,D=0,approximation=FALSE)
arimafit <- auto.arima(usnetelec)
aafit3 <- auto.arima(ukcars,approximation=FALSE)
aafit4 <- auto.arima(log(visitors),approximation=FALSE)


# Figure 2
par(mfrow=c(2,2))
plot(forecast(aafit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")
plot(forecast(arimafit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")
plot(forecast(aafit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")
plot(forecast(aafit4),lambda=0,xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")

# COMPARISONS
summary(forecast(etsfit))
summary(forecast(arimafit))



