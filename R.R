#CO2 emissions from the industrial sector in the USA (Millions of Tm)

#Source: US Energy Information Administration 
#https://www.eia.gov/totalenergy/data/monthly/ Environment/ Industrial Sector

#ARIMA MODEL
serie=ts(read.table("co2IndUSA.dat"),start=1990,freq=12)
plot(serie,main="CO2 emissions from the industrial sector in the USA", ylab="millions of Tm")
abline(v=1985:2020,lty=3,col=4)
#It seems to be a stationary pattern, as the series has been decreased on some periods. 
 
#1. IDENTIFICATION
#a) Determine the transformations that make the series stationary. Justify with graphs and numerical results the selected transformations.
#IS THE VARIANCE CONSTANT?
#Plot of means-variances
m=apply(matrix(serie,nr=12),2,mean)
v=apply(matrix(serie,nr=12),2,var)
plot(m,v,xlab="Anual means",ylab="Anual variances",main="serie")
abline(lm(v~m),col=2,lty=3,lwd=2)
#There is not a clear increase of the variances as the level of the series (the means) increases.

w1 <- window(serie,start=1990,end=2000)
w2 <- window(serie,start=2000,end=2010)
w3 <- window(serie,start=2010,end=2019)

par(mfrow=c(1,3))
boxplot(w1,ylim = c(min(serie),max(serie)))
boxplot(w2,ylim = c(min(serie),max(serie)))
boxplot(w3,ylim = c(min(serie),max(serie)))
#Variance appears to be constant: it seems not to change with the level of the series

#Boxplot per years (measure of dispersion)
boxplot(serie~floor(time(serie)))
#We can see that the widths of the boxplots don't change with the level of the series, as we can appreciate wide boxes both on the bottom and on the upper part of the Y axis. As it's not completely clear whether the variance is constant or not, we will better not apply a log-transformation to our data. Thus, we will work with the original series.

#IS THERE SEASONALITY PRESENT?
#Decomposition in basic components (trend+seasonality+residual)
plot(decompose(serie))
#In the trend plot, we can see that the evolution of the series through the years (without the seasonal component) has been irregular. We can also perceive a seasonal pattern if we look at the seasonal plot. And the noise is changing (random plot). The observed plot is obtained by adding the random, seasonal and trend plots.

#Analysis of the seasonal indices (mensual plot)
monthplot(serie)
#We can see seasonal patterns. That is to say, a pattern repeats itself every s=12 months. In January, August and December, there were more C02 emissions. Thus, we will apply a seasonal difference.

#Seasonal difference of order 12 to remove the seasonal component
d12serie=diff(serie,12)

#IS THE MEAN CONSTANT?
plot(d12serie)
abline(h=0)
#The mean doesn't look constant, so we will apply the necessary regular differences so as to transform the series into stationary.

#Regular differenciation to make the mean constant
d1d12serie=diff(d12serie,1)
plot(d1d12serie)
abline(h=0)

mean(d1d12serie)
#Now, the mean looks constant and centered at zero. However, we can appreciate some outliers. Let's see if other regular differentiation is needed.

#New regular differenciation for a possible non-stationarity
d1d1d12serie=diff(d1d12serie,1)
plot(d1d1d12serie)
abline(h=0)
#We can see that the series didn't change much. As follows, we will check if this last regular differentiation was needed through the comparison of the variances.

#ARE WE OVERDIFFERENCIATING?
#Variances comparison
var(serie)
var(d12serie) #seasonal difference
var(d1d12serie) #seasonal difference + regular difference
var(d1d1d12serie) #seasonal difference + 2 regular differences
#We can see that the second regular difference was not needed, as the variance is artificially increased.

par(mfrow=c(1,2))
acf(serie,ylim=c(-1,1),col=c(2,rep(1,11)),lwd=2,lag.max=72)
acf(d1d12serie,ylim=c(-1,1),col=c(rep(1,11),2),lwd=2,lag.max=72)
par(mfrow=c(1,1))
#In conclusion, the series d1d12serie is stationary, as the mean and variance are constant, there is no seasonality present, and the ACF decays rapidly towards 0. 

#b) Analyze the ACF and PACF of the series to identify at least two possible models. Explain in which features from the graphics you rely to identify these models.
par(mfrow=c(1,2))
acf(d1d12serie,ylim=c(-1,1),col=c(2,rep(1,11)),lwd=2,lag.max=72)
pacf(d1d12serie,ylim=c(-1,1),col=c(rep(1,11),2),lwd=2,lag.max=72)
par(mfrow=c(1,1))
#For the regular component (we will only look at the first lags):**
#In the PACF, we can identify 3 significant lags (different from zero). Also, we can appreciate an alternate (positive lag-negative lag) exponential decreasing pattern in the ACF. Thus, we can propose an AR(3) model for the regular part. We can identify a clear MA(1) or a MA(5) for the regular part by looking at the ACF's first lags. Also, we can clearly see an exponential decreasing pattern in the PACF. We can propose an ARMA (1,1) model for the regular part as well. However, we will start fitting the simplest models, AR(3) and MA(1).

#FOR THE SEASONAL COMPONENT (WE WILL ONLY LOOK AT THE RED LAGS)
#In the ACF, we can appreciate one clearly significant red lag. As for the PACF, we can appreciate three significant red lags. However, we can only identify an exponential decreasing pattern in the PACF. Thus, we will propose an MA(1).

#2. ESTIMATION
#a) Use R to estimate the identified models.
#ARIMA(3,0,0)(0,0,1)12
(mod1=arima(d1d12serie,order=c(3,0,0),seasonal=list(order=c(0,0,1),period=12)))
#Let's check the significance of the estimated parameters.

abs(mod1$coef/sqrt(diag(mod1$var.coef)))
#Both the parameter ar3 and the intercept (mean) are not significant with a t-ratio |t| < 2. Let's start removing the intercept from the model by applying both a regular and a seasonal difference.

#ARIMA(3,1,0)(0,1,1)12
(mod2=arima(serie,order=c(3,1,0),seasonal=list(order=c(0,1,1),period=12)))

abs(mod2$coef/sqrt(diag(mod2$var.coef)))
#We can see that all the parameters but the ar3 are significant and the AIC has slightly decreased. Let us remove the non significant parameter from the model by fitting an $ARIMA(2,1,0)(0,1,1)_{12}$.

#ARIMA(2,1,0)(0,1,1)12
(mod3=arima(serie,order=c(2,1,0),seasonal=list(order=c(0,1,1),period=12)))

abs(mod3$coef/sqrt(diag(mod3$var.coef)))
#Now, all the parameters are significant, but the AIC has slightly increased. Hoever, we will remove this parameter from the model as the increase in the AIC is very little and this model is simple. Let's fit the other of the proposed models, an $ARIMA(0,0,1)(0,0,1)_{12}$.

#ARIMA(0,0,1)(0,0,1)12
(mod4=arima(d1d12serie,order=c(0,0,1),seasonal=list(order=c(0,0,1),period=12)))

abs(mod4$coef/sqrt(diag(mod4$var.coef)))
#In this case, the intercept is not significant. Thus, we will remove it from the model by applying a regular and a seasonal difference.

#ARIMA(0,1,1)(0,1,1)12
(mod5=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12)))

abs(mod5$coef/sqrt(diag(mod5$var.coef)))
#We can see that the two parameters are significant and the AIC decreased in contrast with the last model. 

#3. VALIDATION
################# Validation #################################
validation=function(model,dades){
  s=frequency(get(model$series))
  resid=model$residuals
  par(mfrow=c(2,2),mar=c(3,3,3,3))
  #Residuals plot
  plot(resid,main="Residuals")
  abline(h=0)
  abline(h=c(-3*sd(resid),3*sd(resid)),lty=3,col=4)
  #Square Root of absolute values of residuals (Homocedasticity)
  scatter.smooth(sqrt(abs(resid)),main="Square Root of Absolute residuals",
                 lpars=list(col=2))
  
  #Normal plot of residuals
  qqnorm(resid)
  qqline(resid,col=2,lwd=2)
  
  ##Histogram of residuals with normal curve
  hist(resid,breaks=10,freq=F)
  curve(dnorm(x,mean=mean(resid),sd=sd(resid)),col=2,add=T)
  
  
  #ACF & PACF of residuals
  par(mfrow=c(1,2))
  acf(resid,ylim=c(-1,1),lag.max=60,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid,ylim=c(-1,1),lag.max=60,col=c(rep(1,s-1),2),lwd=1)
  par(mfrow=c(1,1))
  
  #ACF & PACF of square residuals 
  par(mfrow=c(1,2))
  acf(resid^2,ylim=c(-1,1),lag.max=60,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid^2,ylim=c(-1,1),lag.max=60,col=c(rep(1,s-1),2),lwd=1)
  par(mfrow=c(1,1))
  
  #Ljung-Box p-values
  par(mar=c(2,2,1,1))
  tsdiag(model,gof.lag=7*s)
  cat("\n--------------------------------------------------------------------\n")
  
  ##Shapiro-Wilks Normality test
  print(shapiro.test(resid(model)))
  cat("\nNormality Tests\n")
  cat("\n--------------------\n")
 
  ##Shapiro-Wilks Normality test
  print(shapiro.test(resid(model)))

  suppressMessages(require(nortest,quietly=TRUE,warn.conflicts=FALSE))
  ##Anderson-Darling test
  print(ad.test(resid(model)))
  
  suppressMessages(require(tseries,quietly=TRUE,warn.conflicts=FALSE))
  ##Jarque-Bera test
  print(jarque.bera.test(resid(model)))
  
  cat("\nHomoscedasticity Test\n")
  cat("\n--------------------\n")
  suppressMessages(require(lmtest,quietly=TRUE,warn.conflicts=FALSE))
  ##Breusch-Pagan test
  obs=get(model$series)
  print(bptest(resid(model)~I(obs-resid(model))))
  
  cat("\nIndependence Tests\n")
  cat("\n--------------------\n")
  
  ##Durbin-Watson test
  print(dwtest(resid(model)~I(1:length(resid(model)))))
  
  ##Ljung-Box test
  cat("\nLjung-Box test\n")
  print(t(apply(matrix(c(1:4,(1:4)*s)),1,function(el) {
    te=Box.test(resid(model),type="Ljung-Box",lag=el)
    c(lag=(te$parameter),statistic=te$statistic[[1]],p.value=te$p.value)})))
  
  #Sample ACF vs. Teoric ACF
  par(mfrow=c(2,2),mar=c(3,3,3,3))
  acf(dades, ylim=c(-1,1) ,lag.max=36,main="Sample ACF")
  
  plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36),ylim=c(-1,1), 
       type="h",xlab="Lag",  ylab="", main="ACF Teoric")
  abline(h=0)
  
  #Sample PACF vs. Teoric PACF
  pacf(dades, ylim=c(-1,1) ,lag.max=36,main="Sample PACF")
  
  plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
       type="h", xlab="Lag", ylab="", main="PACF Teoric")
  abline(h=0)
  par(mfrow=c(1,1))
}
################# End of Validation #################################

dades=d1d12serie
model=mod3
validation(model,dades)

#1. ARIMA(2,1,0)(0,1,1)12
#a) Perform the complete analysis of residuals, justifying the premises from the corresponding graphical results.
#NORMALITY ASSUMPTION
#Looking at the Q-Q plot, we can accept the normality assumption, as the residuals are following the straight line. The histogram of the residuals also confirms the normality assumption.
#The Shapiro-Wilk, Anderson-Darling and Jarque Bera normality tests also confirm the normality of the residuals, with a p-value higher than 0.05.

#CONSTANT VARIANCE ASSUMPTION
#In the plot of the residuals, we can appreciate a constant variance, but we can also see some residuals. As for the plot of square root of absolute residuals, we can confirm the constant variance of the residuals, as the red line is almost horizontal.
#In the ACF and PACF of the square of the residuals, we can see three significant lags, that could be explained by the presence of outliers or a volatility problem.
#Finally, the studentized Breusch-Pagan test rejects the homoscedasticity of the residuals, with a p-value smaller than 0.05.

#INDEPENDENCE OF THE RESIDUALS ASSUMPTION
#The ACF and PACF of the residuals are compatible with a random noise because there are no significant lags. Thus, we can accept the independence of the residuals. The Durbin-Watson test also accepts the independence assumption, with a p-value higher than 0.05.
#As for the plot of the p-values for the Ljung-Box statistic, we can see that all the lags are non-significant, thus we can consider that the observed autocorrelation structure is well explained by the model. That is to say, all the lags can be considered white noise, something that we can also check in the Ljung-Box test. Thus, we can confirm the independence of the residuals.

#CONCORDANCE TEORIC AND SAMPLE ACF/PACF
#We can appreciate a clear agreement between the teoric (the adjusted model) and sample (the original series) ACF and PACF.

#CONCLUSION
#Although we cannot completely validate this model because of the non-constant variance, instead of proposing a new model, we will first treat for outliers and calendar effects to see whether the variance can be considered constant or not.

#b) Include data from the expressions of the AR and MA infinite models, discuss if they are stationary and / or invertible.
#Stationary and Invertible
  cat("\nModul of AR Characteristic polynomial Roots: ", 
      Mod(polyroot(c(1,-model$model$phi))),"\n")
  cat("\nModul of MA Characteristic polynomial Roots: ",
      Mod(polyroot(c(1,model$model$theta))),"\n")
  
  #Model expressed as an MA infinity (psi-weights)
  psis=ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36)
  names(psis)=paste("psi",1:36)
  cat("\nPsi-weights (MA(inf))\n")
  cat("\n--------------------\n")
  print(psis[1:20])
  
  plot(psis,type="h",main="Psis Weights - MA infinite")
  
  #Model expressed as an AR infinity (pi-weights)
  pis=-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36)
  names(pis)=paste("pi",1:36)
  cat("\nPi-weights (AR(inf))\n")
  cat("\n--------------------\n")
  print(pis[1:20])
  
  plot(pis,type="h",main="Pis Weights - AR infinite")
#We can see that the moduls of AR characteristic polynomial roots are all higher than 1. It means that all the complex roots of the characteristic polynomial of the AR part are outside the unit circle, so we can conclude that the model is causal. Remember than all the AR models are invertible.
#As for the moduls of MA characteristic polynomial roots, they are also higher than 1. It means that all the complex roots of the characteristic polynomial of the MA part lie outside the unit circle, so we can conclude that the model is invertible. Remember that all the MA models are causal.
#Looking at the $\pi$ and $\psi$ weights, they drop down to zero, thus we can confirm that the transformed series is causal and invertible. Thus, we can use this model to make predictions.

#2. ARIMA(0,1,1)(0,1,1)12**
dades=d1d12serie
model=mod5
validation(model,dades)

#a) Perform the complete analysis of residuals, justifying the premises from the corresponding graphical results.
#NORMALITY ASSUMPTION
#Looking at the Q-Q plot, we can accept the normality assumption, as the residuals are following the straight line. The histogram of the residuals also confirms the normality assumption.
#The Shapiro-Wilk, Anderson-Darling and Jarque Bera normality tests also confirm the normality of the residuals, with a p-value higher than 0.05.

#CONSTANT VARIANCE ASSUMPTION
#In the plot of the residuals, we can appreciate a constant variance, but we can also see some residuals. As for the plot of square root of absolute residuals, we can confirm the constant variance of the residuals, as the red line is almost horizontal.
#In the ACF and PACF of the square of the residuals, we can see three significant lags, that could be explained by the presence of outliers or a volatility problem.
#Finally, the studentized Breusch-Pagan test rejects the homoscedasticity of the residuals, with a p-value smaller than 0.05.

#INDEPENDENCE OF THE RESIDUALS ASSUMPTION
#The ACF and PACF of the residuals are not totally compatible with a random noise because there are three significant lags. However, the Durbin-Watson test accepts the independence assumption, with a p-value higher than 0.05.
#As for the plot of the p-values for the Ljung-Box statistic, we can see that all the lags up to k = 26 are non-significant, thus we can consider that the observed autocorrelation structure is not completely well explained by the model. Thus, we cannot totally confirm the independence of the residuals.

**Concordance teoric and sample ACF/PACF**:
#CONCORDANCE TEORIC AND SAMPLE ACF/PACF
#We can appreciate a clear agreement between the teoric (the adjusted model) and sample (the original series) ACF and PACF.

#CONCLUSION
#Although we cannot completely validate this model because of the non-constant variance, instead of proposing a new model, we will first treat for outliers and calendar effects to see whether the variance can be considered constant or not.

#b) Include data from the expressions of the AR and MA infinite models, discuss if they are stationary and / or invertible.
#Stationary and Invertible
  cat("\nModul of AR Characteristic polynomial Roots: ", 
      Mod(polyroot(c(1,-model$model$phi))),"\n")
  cat("\nModul of MA Characteristic polynomial Roots: ",
      Mod(polyroot(c(1,model$model$theta))),"\n")
  
  #Model expressed as an MA infinity (psi-weights)
  psis=ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36)
  names(psis)=paste("psi",1:36)
  cat("\nPsi-weights (MA(inf))\n")
  cat("\n--------------------\n")
  print(psis[1:20])
  
  plot(psis,type="h",main="Psis Weights - MA infinite")
  
  #Model expressed as an AR infinity (pi-weights)
  pis=-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36)
  names(pis)=paste("pi",1:36)
  cat("\nPi-weights (AR(inf))\n")
  cat("\n--------------------\n")
  print(pis[1:20])
  
  plot(pis,type="h",main="Pis Weights - AR infinite")
#We can see that we don't have moduls of AR characteristic polynomial for the AR part, but it is well known that all the MA models are causal.
#As for the moduls of MA characteristic polynomial roots, they are also higher than 1. It means that all the complex roots of the characteristic polynomial of the MA part lie outside the unit circle, so we can conclude that the model is invertible.
#Looking at the $\pi$ and $\psi$ weights, they drop down to zero, thus we can confirm that the transformed series is causal and invertible. Thus, we can use this model to make predictions.

#4. FORECAST
#1- ARIMA(2,1,0)(0,1,1)12
#a) Check the stability of the model and evaluates their foresight, reserving the last 12 observations.
ultim=c(2018,12)
pdq=c(2,1,0)
PDQ=c(0,1,1)

serie2=window(serie,end=ultim) #Incomplete series: leave last year out!
serie1=window(serie,end=ultim+c(1,0)) #Complete series

(modA=arima(serie1,order=pdq,seasonal=list(order=PDQ,period=12)))

(modB=arima(serie2,order=pdq,seasonal=list(order=PDQ,period=12)))
#We can conclude that the model is stable, as both the complete and the incomplete series share similar values (in sign, magnitude and signification). In practice, this means that the correlation structure has not changed in the last year, and that the use of the complete series for making predictions is reliable.

#b) Get long term forecasts for the twelve months following the last observation and include confidence intervals.
#OUT-OF-SAMPLE PREDICTIONS
#How close are the predicted values from the observed?
pred=predict(modB,n.ahead=12)

pr<-ts(c(tail(serie2,1),pred$pred),start=ultim,freq=12) #point predictions
se<-ts(c(0,pred$se),start=ultim,freq=12) #standard errors for point predictions

#Prediction Intervals (back transformed to original scale using exp-function)
tl<-ts(pr-1.96*se,start=ultim,freq=12)
tu<-ts(pr+1.96*se,start=ultim,freq=12)
pr<-ts(pr,start=ultim,freq=12) #predictions in original scale

#Plot of the original series and out-of-sample predictions (only time window 2015-2019 shown)
ts.plot(serie,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-3,+2),type="o",main="Model ARIMA(2,1,0)(0,1,1)12")
abline(v=(ultim[1]-3):(ultim[1]+2),lty=3,col=4)
#The point predictions (in red) are all located between the confidence limits (in blue) and they are very close to the observations (in black), so we can consider that the model makes good predictions

#Point and Interval predictions for year 2019, sing the series without the last 12 observations, and prediction errors
(previs=window(cbind(tl,pr,tu,serie,error=round(serie-pr,3)),start=ultim))

obs=window(serie,start=ultim)
#Measures of acuracy: how exact are the predictions with respect to the observations
(mod.EQM1=sqrt(sum(((obs-pr)/obs)^2)/12)) #RMSPE (Root Mean Square Percentage Error)
(mod.EAM1=sum(abs(obs-pr)/obs)/12) #MAPE (Mean Absolute Percentage Error)
#The RMSPE and the MAPE measures indicate that there is a 2.12% and a 1.60% of error in the predictions. As both percentages are smaller than the rule of thumb value of 10%, we can consider that the current model has a reasonably good prediction ability and that we can proceed to perform long term predictions.

#LONG-TERM PREDICTIONS
#Predictions for the year 2020 with the complete series
pred=predict(modA,n.ahead=12)
pr<-ts(c(tail(serie,1),pred$pred),start=ultim+c(1,0),freq=12)
se<-ts(c(0,pred$se),start=ultim+c(1,0),freq=12)

#Intervals
tl1<-ts((pr-1.96*se),start=ultim+c(1,0),freq=12)
tu1<-ts((pr+1.96*se),start=ultim+c(1,0),freq=12)
pr1<-ts((pr),start=ultim+c(1,0),freq=12)

ts.plot(serie,tl1,tu1,pr1,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(ultim[1]-2,ultim[1]+3),type="o",main="Model ARIMA(2,1,0)(0,1,1)12")
abline(v=(ultim[1]-2):(ultim[1]+3),lty=3,col=4)

#We apply an ARIMA model to this series, so as to find the predictions for the year 2020
(previs1=window(cbind(tl1,pr1,tu1),start=ultim+c(1,0)))

#1. ARIMA(0,1,1)(0,1,1)12
#a) Check the stability of the model and evaluates their foresight, reserving the last 12 observations.
ultim=c(2018,12)
pdq=c(0,1,1)
PDQ=c(0,1,1)

serie2=window(serie,end=ultim) #Incomplete series: leave last year out!
serie1=window(serie,end=ultim+c(1,0)) #Complete series

(modA=arima(serie1,order=pdq,seasonal=list(order=PDQ,period=12)))

(modB=arima(serie2,order=pdq,seasonal=list(order=PDQ,period=12)))
#We can conclude that the model is stable, as both the complete and the incomplete series share similar values (in sign, magnitude and signification). In practice, this means that the correlation structure has not changed in the last year, and that the use of the complete series for making predictions is reliable.

#b) Get long term forecasts for the twelve months following the last observation and include confidence intervals.
#OUT-OF-SAMPLE PREDICTIONS
#How close are the predicted values from the observed?
pred=predict(modB,n.ahead=12)

pr<-ts(c(tail(serie2,1),pred$pred),start=ultim,freq=12) #point predictions
se<-ts(c(0,pred$se),start=ultim,freq=12) #standard errors for point predictions

#Prediction Intervals (back transformed to original scale using exp-function)
tl<-ts(pr-1.96*se,start=ultim,freq=12)
tu<-ts(pr+1.96*se,start=ultim,freq=12)
pr<-ts(pr,start=ultim,freq=12) #predictions in original scale

#Plot of the original series and out-of-sample predictions (only time window 2015-2019 shown)
ts.plot(serie,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-3,+2),type="o",main="Model ARIMA(0,1,1)(0,1,1)12")
abline(v=(ultim[1]-3):(ultim[1]+2),lty=3,col=4)
#The point predictions (in red) are all located between the confidence limits (in blue) and they are very close to the observations (in black), so we can consider that the model makes good predictions

#Point and Interval predictions for year 2019, sing the series without the last 12 observations, and prediction errors
(previs=window(cbind(tl,pr,tu,serie,error=round(serie-pr,3)),start=ultim))

obs=window(serie,start=ultim)
#Measures of acuracy: how exact are the predictions with respect to the observations
(mod.EQM2=sqrt(sum(((obs-pr)/obs)^2)/12)) #RMSPE (Root Mean Square Percentage Error)
(mod.EAM2=sum(abs(obs-pr)/obs)/12) #MAPE (Mean Absolute Percentage Error)
#The RMSPE and the MAPE measures indicate that there is a 2.18% and a 1.71% of error in the predictions. As both percentages are smaller than the rule of thumb value of 10%, we can consider that the current model has a reasonaly good prediction ability and that we can proceed to perform long term predictions.

#LONG-TERM PREDICTIONS
#Predictions for the year 2020 with the complete series
pred=predict(modA,n.ahead=12)
pr<-ts(c(tail(serie,1),pred$pred),start=ultim+c(1,0),freq=12)
se<-ts(c(0,pred$se),start=ultim+c(1,0),freq=12)

#Intervals
tl1<-ts((pr-1.96*se),start=ultim+c(1,0),freq=12)
tu1<-ts((pr+1.96*se),start=ultim+c(1,0),freq=12)
pr1<-ts((pr),start=ultim+c(1,0),freq=12)

ts.plot(serie,tl1,tu1,pr1,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(ultim[1]-2,ultim[1]+3),type="o",main="Model ARIMA(0,1,1)(0,1,1)12")
abline(v=(ultim[1]-2):(ultim[1]+3),lty=3,col=4)

#We apply an ARIMA model to this series, so as to find the predictions for the year 2020
(previs2=window(cbind(tl1,pr1,tu1),start=ultim+c(1,0)))

#SELECTION OF THE BEST MODEL FOR FORECASTING
#3 predictions obtained by the Box-Jenkins methodology
par(mfrow=c(2,2),mar=c(3,3,1,1))
ts.plot(serie,previs1,lty=c(1,2,1,2),col=c(1,4,2,4),xlim=c(2016,2021),type="o",main="Model ARIMA(2,1,0)(0,1,1)12")
abline(v=2016:2021,lty=3,col=4,ylim=c(15,280))
ts.plot(serie,previs2,lty=c(1,2,1,2),col=c(1,4,2,4),xlim=c(2016,2021),type="o",main="Model ARIMA(0,1,1)(0,1,1)12")
abline(v=2016:2021,lty=3,col=4,ylim=c(15,280))
par(mfrow=c(1,1))

resul=data.frame(
  par=c(length(coef(mod3)),length(coef(mod5))),
  Sigma2Z=c(mod3$sigma2,mod5$sigma2),
  AIC=c(AIC(mod3),AIC(mod5)),
  BIC=c(BIC(mod3),BIC(mod5)),
  RMSPE=c(mod.EQM1,mod.EQM2),
  MAPE=c(mod.EAM1,mod.EAM2),
  meanLength=c(sum(previs1[,3]-previs1[,1]),sum(previs2[,3]-previs2[,1]))/12)
row.names(resul)=c("ARIMA(2,1,0)(0,1,1)12","ARIMA(0,1,1)(0,1,1)12")

resul
#We can conclude that the best model for forecasting is the ARIMA(0,1,1)(0,1,1)12 model, as it provides similar adequacy measures as the ARIMA(2,1,0)(0,1,1)12 and it is simpler.

#5. OUTLIER TREATMENT
#a) Analyze whether the Calendar Effects are significant for this series.
source("CalendarEffects.r")

data=c(start(serie)[1],start(serie)[2], length(serie)) #start year[1], month[2]

#We will next fit two auxiliary series, one for the Easter effect and other for the trading days effect. Then, we will use these two series as exogenous variables to fit an ARIMAX model, an ARIMA model with eXogenous variables.
#Auxiliar series Trading days effect
(wTradDays=Wtrad(data))
#A positive value means that there is no balance between the trading days and the weekend days (we expect to find a balance of 5 trading days/2 weekend days). For example, in May 1990, we have 3 extra trading days per each 2 weekend days.On the other hand, the a negative value means that there is an excess of weekend days. In December 1990, we can find 4 trading days less. 
#-> All the months but February will need a trading days correction (in February, we will always find 20 trading days/8 weekend days).

#Auxiliar series Easter effects
(wEast=Weaster(data))
#In April 1990, we can find a +0.50, which means that the Easter was celebrated in April that year. On the contrary, we can find a -0.50 in March. The opposite pattern can be found the next year, 1991.

#We will now build a third auxiliary series, in which we will indicate the moment in which the CO2 emissions started an exponential decrease that still remains. The 1 means that the "intervention" is active. 
#Auxiliar series Intervention analysis
recession=ts(rep(0,length(serie)),start=start(serie),freq=frequency(serie))
recession[229:length(serie)]=1
recession
#In 2009, the CO2 emissions show an important decline. The causes (the "intervention") was the impact of the recession on industrial production and overall energy use.

(mod1=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))) #Original series 

(mod1TD=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=wTradDays)) #Trading days

abs(mod1TD$coef/sqrt(diag(mod1TD$var.coef)))
#The trading days effect is significant, and the AIC and the $\sigma^2_z$ are smaller than in the case of the originals series. Thus, we can correct the original series by trading days. The coefficient estimate indicates that the trading days affect the series by increasing the CO2 emissions.

(mod1Ea=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=wEast)) #Easter

abs(mod1Ea$coef/sqrt(diag(mod1Ea$var.coef)))
#On the contrary, the easter effect is not significant. Thus, it won't be necessary to correct the series for Easter effects, as they are not affecting the behavior of the original series. Also, the AIC has increased with respect to the previous model.

(mod1EC=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,wEast))) #Trading days + Easter

abs(mod1EC$coef/sqrt(diag(mod1EC$var.coef)))
#We can see that the trading days are still significant, whereas the Easter effects are not. The AIC decreased with reespect to the last model.

(mod1ECI=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,wEast,recession))) #Trading days + Easter + Intervention analysis 

abs(mod1ECI$coef/sqrt(diag(mod1ECI$var.coef)))
#We can see that just the trading days effect and the intervention are significant. Also, the AIC is smaller than in the previous model. Thus, we will correct the original series by Trading Days and Intervention effects. The 2009 recession has been negative in the sense that it has affected the series by decreasing the CO2 emissions (negative coefficient estimate).

(mod1TDI=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,recession))) #Trading days + Intervention analysis

abs(mod1TDI$coef/sqrt(diag(mod1TDI$var.coef)))
#We can see that both the Trading Days and the intervention effects are significant. This is the model with the best acuracy. Thus, we will use the trading days and the intervention analysis auxiliary series to correct the original one.

#CORRECTED SERIES BY TRADING DAYS

#Significant coefficients
EfecTD=coef(mod1TDI)["wTradDays"]*wTradDays #Trading days calendar effect
EfecIN=coef(mod1TDI)["recession"]*recession #Intervention analysis

#Series corrected by trading days calendar effects
serieTDI=serie-EfecTD-EfecIN #we have to substract the trading days effect to the original series

#We have to apply the regular and seasonal differences to the series corrected by trading days, as we have previously done to the original series
d1d12serieTDI=diff(diff(serieTDI,12))

#REIDENTIFICATION OF THE MODEL: ACF AND PACF
par(mfrow=c(1,2))
acf(d1d12serieTDI,ylim=c(-1,1),lag.max=72,col=c(2,rep(1,11)),lwd=2)
pacf(d1d12serieTDI,ylim=c(-1,1),lag.max=72,col=c(rep(1,11),2),lwd=2)
par(mfrow=c(1,1))
#We can identify a clear MA(1) both for the regular and the seasonal part. We can appreciate an exponential decrease both in the PACF both for the regular (first lags) and the seasonal (red lags) components.

#FINAL ARIMA MODEL (0,1,1)(0,1,1)12 (CORRECTED BY CALENDAR EFFECTS)

(mod1TDI=arima(serie,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,recession)))

abs(mod1TDI$coef/sqrt(diag(mod1TDI$var.coef)))

EfecTD=coef(mod1TDI)["wTradDays"]*wTradDays
EfecIN=coef(mod1TDI)["recession"]*recession
plot(EfecTD+EfecIN) #Plot of the trading days and intervention effect 
#We can see that, when the recession occurred, the CO2 emissions dropped.

#Linearized time series/teoretic (in red)
plot(serie)
lines(serieTDI,col=2)
abline(v=1992:2018,lty=3,col=4)
#We can see that the CO2 emissions will have been higher if the recession hadn't taken place.

#Differences between the original series and the corrected by trading days and intervention analysis (exogenous variables)
window(cbind(wTradDays,recession,serie,serieTDI,resid(mod1TDI)),start=2005,end=2007)

#VALIDATION
dades=d1d12serieTDI
model=mod1TDI
validation(model,dades)
#We can see that the variance is still non constant, but this model is better than the first one. We will later apply the outliers treatment to our series to finally solve the heterocedasticity problem.

#Stationary and Invertible
  cat("\nModul of AR Characteristic polynomial Roots: ", 
      Mod(polyroot(c(1,-model$model$phi))),"\n")
  cat("\nModul of MA Characteristic polynomial Roots: ",
      Mod(polyroot(c(1,model$model$theta))),"\n")
  
  #Model expressed as an MA infinity (psi-weights)
  psis=ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36)
  names(psis)=paste("psi",1:36)
  cat("\nPsi-weights (MA(inf))\n")
  cat("\n--------------------\n")
  print(psis[1:20])
  
  plot(psis,type="h",main="Psis Weights - MA infinite")
  
  #Model expressed as an AR infinity (pi-weights)
  pis=-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36)
  names(pis)=paste("pi",1:36)
  cat("\nPi-weights (AR(inf))\n")
  cat("\n--------------------\n")
  print(pis[1:20])
  
  plot(pis,type="h",main="Pis Weights - AR infinite")

#FORECASTING
#IS THE MODEL STABLE?
ultim=c(2018,12)

serie1=window(serie,end=ultim+c(1,0))
serie2=window(serie,end=ultim)
wTradDays2=window(wTradDays,end=ultim)
recession2=window(recession,end=ultim)

(modTDI=arima(serie1,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,recession)))

(modTDI2=arima(serie2,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays2,recession2)))
#We can conclude that the model is stable.

#OUT-OF-SAMPLE PREDICTIONS

pred=predict(modTDI2,n.ahead=12,newxreg=window(cbind(wTradDays,recession),start=c(ultim[1]+1,1)))
predic=pred$pr
pr<-ts(c(tail(serie2,1),predic),start=ultim,freq=12)
se<-ts(c(0,pred$se),start=ultim,freq=12)

#Intervals
tl<-ts(pr-1.96*se,start=ultim,freq=12)
tu<-ts(pr+1.96*se,start=ultim,freq=12)
pr<-ts(pr,start=ultim,freq=12)

ts.plot(serie,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-3,2),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA")
abline(v=(ultim[1]-3):(ultim[1]+2),lty=3,col=4)

(previs.lin=window(cbind(tl,pr,tu,serie,error=round(serie-pr,3)),start=ultim))

obs=window(serie,start=ultim)
#RMSPE and MAPE
(mod.EQM2=sqrt(sum(((obs-pr)/obs)^2)/12))
(mod.EAM2=sum(abs(obs-pr)/obs)/12)
#We had slightly more exact predictions with the first model.

#LONG-TERM PREDICTIONS

data2=c(ultim[1]+2, 1, 12)

wTradDays2=Wtrad(data2)
recession2=ts(rep(1,12),start=c(ultim[1]+2),freq=12)

pred=predict(modTDI,n.ahead=12,newxreg=data.frame(wTradDays2,recession2))
predic=pred$pr
pr<-ts(c(serie[length(serie)],predic),start=ultim+c(1,0),freq=12)
se<-ts(c(0,pred$se),start=ultim+c(1,0),freq=12)

#Intervals
tl2<-ts(pr-1.96*se,start=ultim+c(1,0),freq=12)
tu2<-ts(pr+1.96*se,start=ultim+c(1,0),freq=12)
pr2<-ts(pr,start=ultim+c(1,0),freq=12)

ts.plot(serie,tl2,tu2,pr2,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-2,3),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA")
abline(v=(ultim[1]-2):(ultim[1]+3),lty=3,col=4)

mod1TDI=modTDI
(previs3=window(cbind(tl2,pr2,tu2),start=ultim+c(1,0)))

#b) For the final model selected, apply the automatic detection of outliers. Try the interpretation of detected outliers.

#ARIMA MODEL WITH OUTLIER TREATMENT (AND CALENDAR EFFECTS)

#Outliers detection
source("atipics2.r")

#automatic outlier detection
modTDI.atip=outdetec(mod1TDI,dif=c(1,12),crit=2.8,LS=F) 
#We have one regular difference (d = 1) and one seasonal difference of order 12 (D = 12)
#crit = 2.8 (this threshold indicates the number of outliers we want: more if we low the crit value)
#LS = F (we don't want to find level shift, because we don't want the changes in the trend to be interpreted as changes in the level because it will lead to a worse detection of outliers)

#Table of outliers found and the date
atipics=modTDI.atip$atip[order(modTDI.atip$atip[,1]),]
months=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

#Estimated residual variance after outliers detection and treatment
modTDI.atip$sigma2

#Table with detected outliers, their types, magnitud, statistic values, cronology and perc.Obs
data.frame(atipics,Fecha=paste(months[(atipics[,1]-1)%%12+1],start(serie)[1]+((atipics[,1]-1)%/%12)),perc.Obs=atipics[,3])
#We have found 7 outliers, 4 Transitory Changes (TC) and 3 Additive Outliers (AO). W_coeff is the $\beta$ related to the outlier and ABS_L_Ratio os the significance of this $\beta$. We can also appreciate the date in which the outlier was produced. 

#c) With the linearized series, get forecasts for the original series by using the linearized model and compare them with those obtained previously.

#LINEARIZED SERIES, HAVING REMOVED THE OUTLIERS EFFECT
serie.lin=lineal(serie,modTDI.atip$atip)
#An ARIMA model is identified and fitted to this linearized series. The found model, once validated, can be used for performing forecasting. It would be expected that one gets better predictions; particularly more precise ones. 

#FINAL MODEL: ARIMA(0,1,1)(0,1,1)12 having removed the outliers and applied the calendar effects corrections

```{r}
#Estimation of the identified ARIMA model to be linearized series
(mod1TDI.lin=arima(serie.lin,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,recession)))

abs(mod1TDI$coef/sqrt(diag(mod1TDI$var.coef)))
#We can see that all the parameters of this model are significant. 

#Plot of the outliers effect
EfecTD=coef(mod1TDI.lin)["wTradDays"]*wTradDays
EfecIN=coef(mod1TDI.lin)["recession"]*recession
serieTDI.lin=serie.lin-EfecTD-EfecIN
plot(serie-serieTDI.lin)

#In the previous plot, we can find the different types of outliers that have been found.

#Plot together (in original scale) the observed and the theoretical/linearized series (in red, without outliers)
plot(serieTDI.lin,col=2) #in red, the linearized series (without outliers) but in original scale
lines(serie) #in black, the original series

#VALIDATION
dades=diff(diff(serieTDI.lin,12))
model=mod1TDI.lin
validation(model,dades)
#Now, the outliers have disappeared and the normality, independence and homocedasticity of the residuals are clear. Thus, this model is valid.

#Stationary and Invertible
  cat("\nModul of AR Characteristic polynomial Roots: ", 
      Mod(polyroot(c(1,-model$model$phi))),"\n")
  cat("\nModul of MA Characteristic polynomial Roots: ",
      Mod(polyroot(c(1,model$model$theta))),"\n")
  
  #Model expressed as an MA infinity (psi-weights)
  psis=ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36)
  names(psis)=paste("psi",1:36)
  cat("\nPsi-weights (MA(inf))\n")
  cat("\n--------------------\n")
  print(psis[1:20])
  
  plot(psis,type="h",main="Psis Weights - MA infinite")
  
  #Model expressed as an AR infinity (pi-weights)
  pis=-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36)
  names(pis)=paste("pi",1:36)
  cat("\nPi-weights (AR(inf))\n")
  cat("\n--------------------\n")
  print(pis[1:20])
  
  plot(pis,type="h",main="Pis Weights - AR infinite")

#FORECASTING
#IS THE MODEL STABLE?
ultim=c(2018,12)

serie.lin1=window(serie.lin,end=ultim+c(1,0))
serie.lin2=window(serie.lin,end=ultim)
wTradDays2=window(wTradDays,end=ultim)
recession2=window(recession,end=ultim)

(modTDI.lin=arima(serie.lin1,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays,recession)))

(modTDI.lin2=arima(serie.lin2,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12),xreg=data.frame(wTradDays2,recession2)))
#The model is stable.

#OUT-OF-SAMPLE PREDICTIONS
pred=predict(modTDI.lin2,n.ahead=12,newxreg=window(cbind(wTradDays,recession),start=c(ultim[1]+1,1)))
predic=pred$pr
wLS=sum(modTD.atip$atip[modTDI.atip$atip$type_detected=="LS" & modTD.atip$atip$Obs<=length(serie)-12,3])
pr<-ts(c(tail(serie2,1),predic+wLS),start=ultim,freq=12)

se<-ts(c(0,pred$se),start=ultim,freq=12)

#Intervals
tl<-ts(pr-1.96*se,start=ultim,freq=12)
tu<-ts(pr+1.96*se,start=ultim,freq=12)
pr<-ts(pr,start=ultim,freq=12)

ts.plot(serie,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-3,2),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA+Outliers")
abline(v=(ultim[1]-3):(ultim[1]+2),lty=3,col=4)

(previs.lin=window(cbind(tl,pr,tu,serie,error=round(serie-pr,3)),start=ultim))

obs=window(serie,start=ultim)
#RMSPE and MAPE
(mod.EQM3=sqrt(sum(((obs-pr)/obs)^2)/12))
(mod.EAM3=sum(abs(obs-pr)/obs)/12)
#The predictions are now more exact.

#LONG-TERM PREDICTIONS

data3=c(ultim[1]+2, 1, 12)

wTradDays3=Wtrad(data3)
wEast3=Weaster(data3)
pred=predict(modTDI.lin,n.ahead=12,newxreg=data.frame(wTradDays3,wEast3))
predic=pred$pr
wLS=sum(modTDI.atip$atip[modTDI.atip$atip$type_detected=="LS",3])
pr<-ts(c(serie[length(serie)],predic+wLS),start=ultim+c(1,0),freq=12)
se<-ts(c(0,pred$se),start=ultim+c(1,0),freq=12)

#Intervals
tl3<-ts(pr-1.96*se,start=ultim+c(1,0),freq=12)
tu3<-ts(pr+1.96*se,start=ultim+c(1,0),freq=12)
pr3<-ts(pr,start=ultim+c(1,0),freq=12)

ts.plot(serie,tl3,tu3,pr3,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=ultim[1]+c(-2,3),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA+Outliers")
abline(v=(ultim[1]-2):(ultim[1]+3),lty=3,col=4)

mod1TDI.lin=modTDI.lin
(previs4=window(cbind(tl3,pr3,tu3),start=ultim+c(1,0)))

#MODEL COMPARISON
#3 predictions obtained by the Box-Jenkins methodology
par(mfrow=c(2,2),mar=c(3,3,1,1))
ts.plot(serie,previs2,lty=c(1,2,1,2),col=c(1,4,2,4),xlim=c(2016,2021),type="o",main="Model ARIMA(0,1,1)(0,1,1)12")
abline(v=2016:2021,lty=3,col=4,ylim=c(15,280))
ts.plot(serie,previs3,lty=c(1,2,1,2),col=c(1,4,2,4),xlim=c(2016,2021),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA")
abline(v=2016:2021,lty=3,col=4,ylim=c(15,280))
ts.plot(serie,previs4,lty=c(1,2,1,2),col=c(1,4,2,4),xlim=c(2016,2021),type="o",main="Model ARIMA(0,1,1)(0,1,1)12+TD+IA+Outliers")
abline(v=2016:2021,lty=3,col=4,ylim=c(15,280))
par(mfrow=c(1,1))

resul=data.frame(
  par=c(length(coef(mod1)),length(coef(mod1TDI)),length(coef(mod1TDI.lin))+nrow(modTDI.atip$atip)),
  Sigma2Z=c(mod1$sigma2,mod1TDI$sigma2,mod1TDI.lin$sigma2),
  AIC=c(AIC(mod1),AIC(mod1TDI),AIC(mod1TDI.lin)+2*nrow(modTDI.atip$atip)),
  BIC=c(BIC(mod1),BIC(mod1TDI),BIC(mod1TDI.lin)+log(length(serie)-13)*nrow(modTDI.atip$atip)),
  RMSPE=c(mod.EQM1,mod.EQM2,mod.EQM3),
  MAPE=c(mod.EAM1,mod.EAM2,mod.EAM3),
  meanLength=c(sum(previs1[,3]-previs1[,1]),sum(previs2[,3]-previs2[,1]),sum(previs3[,3]-previs3[,1]))/12)
row.names(resul)=c("ARIMA(0,1,1)(0,1,1)12","ARIMA(0,1,1)(0,1,1)12+TD+IA","ARIMA(0,1,1)(0,1,1)12+TD+IA+Outliers")
resul

#We will choose the last model to make the CO2 emissions forecasting for the year 2020.

#Let's finally plot together predicted values using observed series with (blue) and series without outliers (green).

ts.plot(serie,tl1,tu1,pr1,tl2,tu2,pr2,lty=c(1,2,2,1,2,2,1),col=c(1,4,4,2,3,3,6),xlim=ultim[1]+c(1,3),type="o",main="CO2 emissions in USA")
legend("topleft",c("ARIMA(0,1,1)(0,1,1)_12","ARIMA(0,1,1)(0,1,1)_12 with outlier treatment"),col=c(4,3),lty=1,lwd=2)
abline(v=ultim[1]+1:3,lty=3,col=4)


