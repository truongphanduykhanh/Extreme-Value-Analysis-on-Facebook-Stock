# ====================================================================================
#                               Khanh Truong | Paul Dublanche
#                              Extreme Value Theory - Project
#                                          12/2018
# ====================================================================================

# Import library and data ------------------------------------------------------------
rm(list=ls())
library('evd')
library('quantmod')

getSymbols("FB",from='2012-06-01',to='2018-12-01')
head(FB) 
tail(FB)
(plot.price <- plot.xts(FB$FB.Adjusted,main='Facebook Stock Price'))

# Prepare data -----------------------------------------------------------------------
FB$nlogreturn <- -diff(log(FB$FB.Adjusted)) # calculate negative log-return
(plot.return <- plot.xts(FB$nlogreturn,main='Facebook Negative Log-return'))

FB <- FB[!is.na(FB$nlogreturn),] # remove first row
FB <- data.frame(date=index(FB), coredata(FB)) # transform to dataframe 

FB$Year <- as.numeric(substr(FB$date,1,4)) # extract year
FB$Month <- as.numeric(substr(FB$date,6,7)) # extract month
FB$Weekday <- weekdays(FB$date) # extract week-day

month.max <- aggregate(nlogreturn ~ Year+Month, data=FB, FUN=max) # monthly max
summary(month.max$nlogreturn)


### ==================================================================================
### BLOCK MAXIMA =====================================================================
### ==================================================================================

# Fit models =========================================================================
library(MASS) # library for function 'fitdistr'
(normal <- fitdistr(month.max$nlogreturn, densfun="normal")) # Gaussian
(gev <- fgev(month.max$nlogreturn)) # GEV

# Gaussian versus GEV
hist(month.max$nlogreturn, pch=20, breaks=25, prob=TRUE, main=NULL,
     xlab='Negative log-return')
curve(dnorm(x, normal$estimate['mean'], normal$estimate['sd']),col="blue",lwd=2,add=T)
curve(dgev(x, gev$estimate['loc'], gev$estimate['scale'], gev$estimate['shape']),
      col="red", lwd=2, add=T)
legend("topright", legend=c("GEV fit","Gaussian fit"), col=c("red","blue"), lwd=2,
       cex=0.8, box.lty = 0)

# Model checking =====================================================================

# Significance of parameters ---------------------------------------------------------
confint(gev) # Symmetric confidence intervel quick way

# Symmetric confidence intervel manual way 
c(gev$estimate['loc']-qnorm(0.975)*gev$std.err['loc'],
  gev$estimate['loc']+qnorm(0.975)*gev$std.err['loc'])
c(gev$estimate['scale']-qnorm(0.975)*gev$std.err['scale'],
  gev$estimate['scale']+qnorm(0.975)*gev$std.err['scale'])
c(gev$estimate['shape']-qnorm(0.975)*gev$std.err['shape'],
  gev$estimate['shape']+qnorm(0.975)*gev$std.err['shape']) # contain 0

# Parameter 'shape' ------------------------------------------------------------------
# Wald test
(p_value <- 2*pnorm(abs(gev$estimate['shape'])/gev$std.err['shape'],lower.tail=FALSE))
# likelihood ratio test, H0: 2 models are the same
gumble <- fgev(month.max$nlogreturn, shape=0)
anova(gev,gumble) 
# likelihood ratio test by hand
pchisq((abs(gev$deviance-gumble$deviance)), df=1, lower.tail=FALSE)
# => 'shape' is significantly different from 0 => keep gev

# Plot GEV and Gumbel 
hist(month.max$nlogreturn, breaks=25, prob=TRUE, main=NULL,
     xlab = 'Negative log-return')
curve(dgev(x, gev$estimate['loc'], gev$estimate['scale'],
           gev$estimate['shape']), col="red", lwd=2, add=T)
curve(dgev(x, gumble$estimate['loc'], gumble$estimate['scale'], 0),
      col="blue", lwd=2, add=T)
legend("topright", legend=c("GEV fit","Gumbel fit"), col=c("red","blue"),
       lwd=2, cex=0.8, box.lty=0)
# GEV fits more, especially at the peak

par(mfrow=c(1,2))
qq(gev, main="QQ Plot of GEV fit")
qq(gumble, main="QQ Plot of Gumbel fit")
# Understandingly pointwise CI of Gumbel is smaller than GEV
# It is because one parameter is fixed
# But Gumble is not fit with large values (underestimate the large values)

qq(gev, xlim=c(0.01,0.03), ylim=c(0,0.06), main="QQ Plot of GEV fit")
qq(gumble, xlim=c(0.01,0.03), ylim=c(0,0.06), main="QQ Plot of Gumble fit")
# zoom in. very alike

# Choose GEV since:
# 1.'shape' is significant different from 0 (wald 0.05, likelihood 0.01)
# 2. As a risk analyst, we are more risk-adverse,
# we choose model which fits more to the worst case (large values)


# Return period - Return level =======================================================
gev <- fgev(month.max$nlogreturn) # recall GEV to compare

gev.1year <- fgev(month.max$nlogreturn, prob = 1/12) # 1-year return level
gev.0.5year <- fgev(month.max$nlogreturn, prob = 1/6) # 0.5-year return level
gev.3year <- fgev(month.max$nlogreturn, prob = 1/36, std.err = FALSE) # 3-year rl
# standard error 3-year rl can't be estimated because information matrix is singular

# Check re-parameterization in model gev ---------------------------------------------
loc <- gev$estimate['loc']
scale <- gev$estimate['scale']
shape <- gev$estimate['shape']

# check re-parameterization in three models
loc - (scale/shape)*(1-(-log(1-1/12))^(-shape)) # 1 year
gev.1year$estimate['quantile'] # similar to 'quantile' of re-parameterized

loc - (scale/shape)*(1-(-log(1-1/6))^(-shape))  # 0.5 year
gev.0.5year$estimate['quantile'] # similar to 'quantile' of re-parameterized

loc - (scale/shape)*(1-(-log(1-1/36))^(-shape)) # 3 years
gev.3year$estimate['quantile'] # similar to 'quantile' of re-parameterized

rm(loc,scale,shape)

# Return period - Return level Plot ==================================================
rl(gev,main = NULL) # return level plot

# Check return level curve -----------------------------------------------------------
rl.1year <- qgev(1/12, gev$estimate["loc"], gev$estimate["scale"],
  gev$estimate["shape"], lower.tail=FALSE) # estimate 1-year return level
period <- function(T){
  return(-1/log(1- 1/(T)))} # x-axis of rl(gev) is not T but period(T) 
abline(v=period(12), h=rl.1year, col='red') # Check return level curve

# Check empirical points -------------------------------------------------------------
largest <- max(month.max$nlogreturn)
second.largest <- sort(month.max$nlogreturn)[nrow(month.max)-1]

abline(v=period(nrow(month.max)+1), h=largest, col='blue') # check max
abline(v=period((nrow(month.max)+1)/2), h=second.largest, col='blue') # second max

# Intepretation: 0.07358201  is expected to be exceeded every 1 years.
# Equivalently: There is probability 1/12 (0.083) that
# 0.07358201 is exceeded each month  


# Confidence intervel of Return level ================================================
plot(profile(gev.1year, "quantile"), main='Log-likelihood of 1-Year Return Level')
# locator(2)
abline(v=c(0.06125737,0.09350722),col='red') # get C.I from locator(2)
text(0.065,190,'0.061') # add text on the plot
text(0.097,190,'0.094') # add text on the plot
confint(gev.1year)['quantile',] # compare to symmetric C.I

# For appendix: gev.0.5year models obtains more data
# So log-likelihood function would be smoothier.
plot(profile(gev.0.5year, "quantile"), main='Log-likelihood of 6-Month Return Level')


### ==================================================================================
### PEAK OVER THRESHOLD ==============================================================
### ==================================================================================
plot.return # stationary

# Choose threshold ===================================================================
sum(FB$nlogreturn>0.1,na.rm = TRUE)
sum(FB$nlogreturn>0.1,na.rm = TRUE)/nrow(FB) # very small sample

sum(FB$nlogreturn>0.08,na.rm = TRUE)
sum(FB$nlogreturn>0.08,na.rm = TRUE)/nrow(FB) # small sample

# Plot to select threshold
# All the fixed numbers below are just for nice plot purpose
par(mfrow=c(1,2))
tcplot(FB$nlogreturn,tlim=c(0,0.08),ylims = c(-0.195,0.2),which = 1)
abline(h=0,col='red',xpd=TRUE)

tcplot(FB$nlogreturn,tlim=c(0,0.08),which = 2)
abline(h=0.36,col='red',xpd=TRUE)
abline(v=0.0467,col='red')

threshold <- 0.0467 # select threshold


# Fit GPD distribution ===============================================================
(gpd <- fpot(FB$nlogreturn, threshold, npp = 252)) # Generalized Pareto
(pp <- fpot(FB$nlogreturn, threshold, npp = 252, model = "pp")) # model point process

# Check relationship of parameters between 2 models
pp$estimate['scale']+pp$estimate['shape']*(threshold-pp$estimate['loc'])
gpd$estimate['scale'] # very similar to method Generalized Pareto

par(mfrow=c(1,1))
qq(gpd) # modeling checking

par(mfrow=c(1,1)) # plot histogram and the estimator
exceedances <- FB[FB$nlogreturn>=threshold,'nlogreturn']
hist(exceedances, breaks=10, prob=TRUE,
     main=NULL,xlab = 'Negative log-return',ylim = c(0,40))
curve(dgpd(x,pp$estimate['loc'],pp$estimate['scale'],pp$estimate['shape']),
      col="blue", lwd=2, add=T)
legend("topright", legend=c("GPD fit"), col=c("blue"),
       lwd=2, cex=0.8, box.lty = 0)

# Generalized Pareto 1-year
(gpd.1year <- fpot(FB$nlogreturn, threshold, npp = 252, mper = 1))  

# Profile log-likelihood
plot(profile(gpd),'shape')
# locator(2)
abline(v=c(0.05151051,0.94064803),col='red') # get C.I from locator(2)
text(0.15,100,'0.052') # add text on the plot
text(1.05,100,'0.941') # add text on the plot

confint(gpd)['shape',]

(p_value <- 2*(1-pnorm(abs(gpd$estimate['shape'])/gpd$std.err['shape']))) # Wald test


# Return period ======================================================================
par(mfrow=c(1,1))

(gpd <- fpot(FB$nlogreturn, threshold, npp = 252))
(gpd.0.5year <- fpot(FB$nlogreturn, threshold, mper = 0.5, npp = 252))
(gpd.1year <- fpot(FB$nlogreturn, threshold, mper = 1, npp = 252))
(gpd.3year <- fpot(FB$nlogreturn, threshold, mper = 3, npp = 252))

plot(profile(gpd.1year),"rlevel")
# locator(2)
abline(v=c(0.06475264,0.09296453),col='red',lty=2) # get C.I from locator(2)
text(0.07,102.5,'0.065') # add text on the plot
text(0.098,102.5,'0.093') # add text on the plot

confint(gpd.1year)['rlevel',]

rl(gpd)


### ==================================================================================
### COMPARE GEV AND GPD ==============================================================
### ==================================================================================
# Compare estimates
gev.0.5year$estimate['quantile']
gpd.0.5year$estimate['rlevel']

gev.1year$estimate['quantile']
gpd.1year$estimate['rlevel']

gev.3year$estimate['quantile']
gpd.3year$estimate['rlevel']

# Compare CI
confint(gev.0.5year)['quantile',]
confint(gpd.0.5year)['rlevel',]

confint(gev.1year)['quantile',]
confint(gpd.1year)['rlevel',]

confint(gev.3year)['quantile',]
confint(gpd.3year)['rlevel',]