library(lintr)
library(devtools) 
library(RDataTracker)
library(MASS)
library(multcomp)
library(lattice)
library(DHARMa)
library(astsa)
library(vcd)

# automated format check and debug tool
lint("Phenology_model.R")
ddg.run("Phenology_model.R")
ddg.display()

#------------------------
# Read in and format data
#------------------------
# Count data
df <- read.csv("ByCatch.csv")
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
# get rid of incomplete records
df <- na.omit(df)

# subset by species code
Code <- factor("AMBBIS", levels = levels(df$Species))
cap.dates <- subset(df, df$Species == Code & df$Direction.of.Travel == "IN")
# get counts
cap.dates <- hist(cap.dates[, 1], breaks = 1000, freq = T)
cap.dates <- data.frame(as.Date(cap.dates$mids, origin = "1970-01-01"),
                        cap.dates$counts)
colnames(cap.dates) <- c("Date", "Count")
cap.dates <- cap.dates[cap.dates$Count > 0, ]

# get all dates
all.dates <- unique(df[, 1])
all0.dates <- all.dates[ !all.dates %in% cap.dates$Date ]
all0.dates <- data.frame(all0.dates, rep(0, length(all0.dates)))
colnames(all0.dates) <- c("Date", "Count")

dat <- rbind(cap.dates, all0.dates)


# covariate data
dd1 <- read.csv("10-11_All.csv")
dd2 <- read.csv("11-12_All.csv")
dd3 <- read.csv("12-13_All.csv")
dd4 <- read.csv("13-14_All.csv")
dd5 <- read.csv("14-15_All.csv")
dd6 <- read.csv("15-16_All.csv")

dd1$Date <- as.Date(dd1$Date, format = "%m/%d/%Y") - 1
dd2$Date <- as.Date(dd2$Date, format = "%m/%d/%Y") - 1
dd3$Date <- as.Date(dd3$Date, format = "%m/%d/%Y") - 1
dd4$Date <- as.Date(dd4$Date, format = "%m/%d/%Y") - 1
dd5$Date <- as.Date(dd5$Date, format = "%m/%d/%Y") - 1
dd6$Date <- as.Date(dd6$Date, format = "%m/%d/%Y") - 1

covs <- rbind(dd1, dd2, dd3, dd4, dd5, dd6)
# check that all cov dates are included in count data 
subset(dat$Date, !dat$Date %in% covs$Date)

# merge based on date
all.data <- merge(dat, covs)

# append year column
tick <- as.Date(c("2010-06-15",
                 "2011-06-15",
                 "2012-06-15",
                 "2013-06-15",
                 "2014-06-15",
                 "2015-06-15",
                 "2016-06-15"))

year <- vector()
for (i in 1:(length(tick) - 1)) {
  x <- length(subset(all.data$Date, all.data$Date > tick[i] &
                                    all.data$Date < tick[i + 1]))
  year <- append(year, rep(i, x))
}

all.data <- cbind(all.data, year)
all.data$year <- as.factor(all.data$year)


all.data <- cbind(all.data, jds = as.numeric(format(all.data$Date, "%j")))

for (i in 1:length(all.data$jds)) {
  if (all.data$jds[i] <= 180) {
    all.data$jds[i] <- all.data$jds[i] + 366
  }
}

all.data$jds <- all.data$jds - 255
min(all.data$jds)


# proportion of individuals histogram
Nsum <- by(all.data$Count, INDICES = all.data$year, sum)
prop <- all.data$Count / Nsum[all.data$year]

plot(all.data$Date, prop, type = "h",
     xlab = "Date", ylab = "Proportion of Individuals", lwd = 2)

#------------
#Analysis
#------------

# Null Model
fit0 <- glm.nb(Count ~ 1, all.data)

summary(fit0)

1 - pchisq(summary(fit0)$deviance,
           summary(fit0)$df.residual) # GOF test really bad!


# Poisson Model
fit1 <- glm(Count ~ ppt + tmin + ppt * tmin + jds + year, all.data,
            family = "poisson")

summary(fit1)

1 - pchisq(summary(fit1)$deviance,
           summary(fit1)$df.residual) # GOF test really bad!

poipred <- cbind(all.data,
      Mean <- predict(fit1, newdata = all.data, type = "response"),
      SE <- predict(fit1, newdata = all.data, type = "response", se.fit = T)$se.fit)


# Negative binomial model
fit2 <- glm.nb(Count ~ ppt + tmin + ppt * tmin + year + jds, all.data)
summary(fit2)

# see if year as a whole is significant (it is)
m2 <- update(fit2, . ~ . - jds)
anova(fit2, m2)

# pairwise comparison of specific years 
summary(glht(fit2, mcp(year = "Tukey")))

1 - pchisq(summary(fit2)$deviance,
           summary(fit2)$df.residual) # GOF test really good!

anova(fit2, fit0) # significantly better than null model

# predictions for plotting
nbpred <- cbind(all.data,
               Mean = predict(fit2, newdata = all.data, type = "response"),
               SE = predict(fit2, newdata = all.data, type = "response", se.fit = T)$se.fit)


xy <- data.frame(ppt = sort(rep(seq(1, 150, length.out = 100), 100)),
                 tmin = rep(seq(-5, 25, length.out = 100), 100))
ys <- sort(rep(1:6, 10000))
ds <- sort(rep(all.data$jds, length.out = 60000))

xy <- cbind(rbind(xy, xy, xy, xy, xy, xy), ys, ds)
colnames(xy)  <- c("ppt", "tmin", "year", "jds")
xy[, 3] <- as.factor(xy[, 3])

cpred <- predict(fit2, newdata = xy, type = "response")
cpred[cpred > 80] <- 80
cplot <- cbind(xy, cpred)

contourplot(cpred ~ ppt * tmin, data = cplot, cuts = 40, labels = F,
            xlab = "Precipitation (mm)", ylab = "Min Temperature")

# diagnostic plot from ver hoef and boveng
res <- (nbpred$Count - nbpred$Mean) ^ 2
plot(nbpred$Mean, nbpred$SE ^ 2)
plot(nbpred$Mean, res)

plot(nbpred$Count, nbpred$SE ^ 2, xlim = c(0, 20))


# Polynomial NB model
plys <- poly(all.data$ppt, 3)
plys1 <- poly(all.data$tmin, 3)
plys2 <- poly(all.data$jds, 3)

fit3 <- glm.nb(Count ~ ppt + plys1 + ppt * tmin + plys2 + year, all.data)
summary(fit3)

anova(fit2, fit3) # polynomial model not better!

exp(0.003) / (1 + exp(0.003))
y <- function(x) x ^ 3 + x
plot(1:100, y(1:100))

# Simulated NB vs real data
n <- 867
data.nb <- rnbinom(n, size = fit3$theta, mu = exp(coef(fit3)[1]))



fits <- hist(data.nb, breaks = 20)
fits$counts <- sqrt(fits$counts) #sqrt frequencies for plotting


sum(data.nb)
sum(fits$counts)
sum(all.data$Count)
raw <- hist(all.data$Count, breaks = 70)
raw$counts <- sqrt(raw$counts)   #sqrt frequencies for plotting


# plot histograms side by side
par(mfrow = c(1, 2))
plot(raw,  xlim = c(0, 40), ylim = c(0, 30), ylab = expression(sqrt(Frequency)),
     xlab = "Captures", main = "Obsevered Data")
plot(fits, xlim = c(0, 40), ylim = c(0, 30), ylab = expression(sqrt(Frequency)),
     xlab = "Captures", main = "Predicted Data")

x <- fit2$model$Count
# rootogram to assess model fit
rootogram(~x, dfun = function(x)
          dnbinom(x, size = fit2$theta, mu = exp(coef(fit2)[1])),
                  probability = F, transormation = sqrt)
rootogram(fit1)
mean(poipred$SE)
mean(nbpred$SE)
# Here you see the 'danger' of ignoring overdispersion in the Poisson model. 
# The SE estimates are lower for the Poisson model than for the negative binomial model, 
# which increases the likelihood of incorrectly detecting a significant effect in the Poisson model.


# chi sq test to see how good the model (fit2) is (it's good!)
raw$counts ^ 2
obs <- c(763,  31,  22,  11,  9,   3,   28)
sum(obs)

fits$counts ^ 2
exs <- c(669,  40,  23,  22,  16,   13, 84)
sum(exs)
cbind(obs, exs)

chisq.test(cbind(obs, exs))


# check residual plots
simulationOutput <- simulateResiduals(fittedModel = fit2, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput)


testZeroInflation(simulationOutput = simulationOutput)

plotResiduals(all.data$ppt,  simulationOutput$scaledResiduals)
plotResiduals(all.data$tmin,  simulationOutput$scaledResiduals)






#---------------------------------------------------------
# cross correlation lag analysis (need to be every day!!)
#---------------------------------------------------------

# Count data
df <- read.csv("ByCatch.csv")
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
# get rid of incomplete records
df <- na.omit(df)

# subset by species code
Code <- factor("AMBBIS", levels = levels(df$Species))
cap.dates <- subset(df, df$Species == Code & df$Direction.of.Travel == "IN")
# get counts
cap.dates <- hist(cap.dates[, 1], breaks = 1000, freq = T)
cap.dates <- data.frame(as.Date(cap.dates$mids, origin = "1970-01-01"),
                        cap.dates$counts)
colnames(cap.dates) <- c("Date", "Count")
cap.dates <- cap.dates[cap.dates$Count > 0, ]

# get all dates
all.dates <- unique(df[, 1])
all0.dates <- all.dates[ !all.dates %in% cap.dates$Date ]
all0.dates <- data.frame(all0.dates, rep(0, length(all0.dates)))
colnames(all0.dates) <- c("Date", "Count")

dat <- rbind(cap.dates, all0.dates)

# all dates without sampling
ds <- seq(as.Date("2010-10-30"), as.Date("2016-06-04"), by = "days")
nads <- ds[ !ds %in% dat$Date]
nads <- data.frame(as.Date(nads), rep(NA, length(nads)))
names(nads) <- c("Date", "Count")
dat <- rbind(dat, nads)

# covariate data
dd1 <- read.csv("10-11_All.csv")
dd2 <- read.csv("11-12_All.csv")
dd3 <- read.csv("12-13_All.csv")
dd4 <- read.csv("13-14_All.csv")
dd5 <- read.csv("14-15_All.csv")
dd6 <- read.csv("15-16_All.csv")

dd1$Date <- as.Date(dd1$Date, format = "%m/%d/%Y") - 1
dd2$Date <- as.Date(dd2$Date, format = "%m/%d/%Y") - 1
dd3$Date <- as.Date(dd3$Date, format = "%m/%d/%Y") - 1
dd4$Date <- as.Date(dd4$Date, format = "%m/%d/%Y") - 1
dd5$Date <- as.Date(dd5$Date, format = "%m/%d/%Y") - 1
dd6$Date <- as.Date(dd6$Date, format = "%m/%d/%Y") - 1

covs <- rbind(dd1, dd2, dd3, dd4, dd5, dd6)
# check that all cov dates are included in count data 
subset(dat$Date, !dat$Date %in% covs$Date)

# merge based on date
all.data <- merge(dat, covs, by = "Date", all = T)

par(mfrow = c(1, 2))
ccf(x = all.data$ppt, y = all.data$Count, na.action = na.pass)
ccf(x = all.data$tmin, y = all.data$Count, na.action = na.pass)
 
 lag2.plot(all.data$tmax, all.data$Count, 10)
 

 
 set.seed(123)
 x <- arima.sim(model = list(0.2, 0, 0.5), n = 100)
 y <- arima.sim(model = list(0.4, 0, 0.4), n = 100)
 plot(x, y)
 ccf(x, y, type = "correlation")
 
 
 
 
 
 
 # test diagnostics and simulations
 
 
 data(quine)
 fit <- goodfit(quine$Days)
 summary(fit)
 rootogram(fit)
 
 
fit <- goodfit(all.data$Count, type = "nbinomial")
summary(fit)
rootogram(fit)

fit <- goodfit(all.data$Count, type = "poisson")
summary(fit)
rootogram(fit)


Ord_plot(all.data$Count)
distplot(all.data$Count, type = "nbinomial")
distplot(all.data$Count, type = "poisson")
influencePlot(fit2)
influencePlot(fit1)

x <- rnbinom(n = 100000, size = .1, mu = 2)
fit <- goodfit(x, type = "nbinomial")
summary(fit)
rootogram(fit)

fit <- goodfit(x, type = "poisson")
summary(fit)
rootogram(fit)