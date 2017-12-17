library(lintr)
library(devtools) 
library(RDataTracker)
library(MASS)
library(multcomp)
library(lattice)
library(DHARMa)
library(astsa)
library(vcd)
library(plotly)
library(visreg)
library(AICcmodavg)
library(RColorBrewer)

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
cap.dates <- subset(df, df$Species == Code)

#& df$Direction.of.Travel == "OUT" & df$Direction.of.Travel == "IN"


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

all.data$jds <- all.data$jds - 254
min(all.data$jds)

all.data[all.data$ppt > 200,]
all.data = all.data[-526,]

max(all.data$tmin)
min(all.data$tmin)
mean(all.data$tmin)
max(all.data$ppt)
min(all.data$ppt)
mean(all.data$ppt)

x = all.data[all.data$Count > 0, ]
max(x$tmin)
min(x$tmin)
mean(x$tmin)
max(x$ppt)
min(x$ppt)
mean(x$ppt)


# proportion of individuals histogram
Nsum <- by(all.data$Count, INDICES = all.data$year, sum)
prop <- all.data$Count / Nsum[all.data$year]

plot(all.data$Date, prop, type = "h",
     xlab = "Date", ylab = "Proportion of Individuals", lwd = 2)

#------------
#Analysis
#------------

varnames <- c("ppt", "tmin", "ppt * tmin", "jds", "year", "I(jds^2)", "I(tmin^2)", "I(ppt^2)")
rhs <- unlist( sapply(1:length(varnames),function(k) apply(combn(varnames,k),2,paste,collapse=" + ") ) )
formulae <- as.formula( quote( paste("Count ~", rhs) ) )

aics = vector()
for(i in 1:length(formulae)){
  fit = glm.nb(formulae[i], all.data)
  aics[i] = fit$aic
}


aicsP = vector()
for(i in 1:length(formulae)){
  fit = glm(formulae[i], all.data, family = "poisson")
  aicsP[i] = fit$aic
}

# Poisson Model
fit1 <- glm(formulae[252], all.data, family = "poisson")

# Negative binomial model
fit2 <- glm.nb(formulae[252], all.data)
summary(fit2)


fit3 <- glm.nb(Count ~ ppt + tmin + jds + year + I(jds^2) + I(tmin^2) + I(ppt^2) +jds*ppt, all.data)
summary(fit3)



# AIC weights
aics = unique(aics)
daics = aics - aics[162]
sort(daics)
ws = exp(-.5 * daics) / sum(exp(-.5 * daics)) 
sort(ws)


# get R to print decimals not exponetials
options("scipen"=100, "digits"=4)
round(ws, 3)


# pairwise comparison of specific years 
summary(glht(fit2, mcp(year = "Tukey")))


# predictions for plotting
xy <- data.frame(ppt = sort(rep(seq(1, 150, length.out = 100), 100)),
                 tmin = rep(seq(-8, 25, length.out = 100), 100))
ys <- sort(rep(1:6, 10000))
ds <- sort(rep(all.data$jds, length.out = 60000))

xy <- cbind(rbind(xy, xy, xy, xy, xy, xy), ys, ds)
colnames(xy)  <- c("ppt", "tmin", "year", "jds")
xy[, 3] <- as.factor(xy[, 3])

cpred <- predict(fit2, newdata = xy, type = "response")
max(cpred)
cpred[cpred > 80] <- 80
cplot <- cbind(xy, cpred)

contourplot(cpred ~ ppt * tmin, data = cplot, cuts = 230, labels = T,
            xlab = "Precipitation (mm)", ylab = "Min Temperature")

plot_ly(x = cplot$ppt, y = cplot$tmin, z = cplot$cpred, type = "contour", contours = list(coloring = "Greens"))

plot_ly(x = cplot$ppt, y = cplot$tmin, z = cplot$cpred, 
        type = "contour", colors = "Greys",
        contours = list(start = 0, end = 10, size = 1)) %>% 
  colorbar(title = "Nightly Captures", titleside = "right") %>% 
  layout(xaxis = list(title = "Precipitation (mm)"), 
         yaxis = list(title = "Min Temperature"))


# diagnostic plot from ver hoef and boveng
res <- (nbpred$Count - nbpred$Mean) ^ 2
plot(nbpred$Mean, nbpred$SE ^ 2)
plot(nbpred$Mean, res)

plot(nbpred$Count, nbpred$SE ^ 2, xlim = c(0, 20))



cpred <- predict(fit3, type = "response")
cpred[cpred > 80] <- 80
cplot <- cbind(all.data, cpred)

contourplot(cpred ~ ppt * tmin, data = cplot, cuts = 2, labels = T,
            xlab = "Precipitation (mm)", ylab = "Min Temperature")

plot_ly(x = cplot$ppt, y = cplot$tmin, z = cplot$cpred, 
        type = "contour",autocontour = TRUE, line = list(smoothing = 0.85),
        contours = list(coloring = "heatmap", showlabels = TRUE, start = 1, size = 3)) %>% 
  colorbar(title = "Nightly Captures")


#----------------------------
# Simulated NB vs real data
raw <- hist(all.data$Count, breaks = -.5:89.5)
raw$counts <- sqrt(raw$counts)   #sqrt frequencies for plotting
n <- length(all.data$Count)

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(raw, xlim = c(0, 100), ylim = c(0,30), col = "grey70", freq = T,
      ylab = expression(sqrt(Frequency)), xlab = " ", main = " ", 
      axes = FALSE)
ps = seq(0, 90, by = 5)
points(ps, sqrt(dnbinom(ps, size = fit2$theta, mu = exp(sum(coef(fit2)[c(1:4, 10:12)]))) * n), 
       col = "black", pch = 16, cex = 0.8) 
lines(0:90, sqrt(dnbinom(0:90, size = fit2$theta, mu = exp(sum(coef(fit2)[c(1:4, 10:12)]))) * n), 
      col = "black", lty = "dashed", lwd = 2)

axis(1)
axis(2)
mtext("Nightly Captures", side = 1, line = 2.8, cex = 1.5)
legend(30,25, c('observed', 'predicted'), pch = c(15,16), lty = c(NA,"dashed"), 
       col=c('grey70',"black"), pt.cex=1.5, lwd = 2, cex=1.5, bty='n') 




1 - pchisq(summary(fit2)$deviance,
           summary(fit2)$df.residual)

# 
# max(data.nb)
# 
fits <- hist(data.nb, breaks = 100)
# fits$counts <- sqrt(fits$counts) #sqrt frequencies for plotting
# 
# 
# sum(data.nb)
# sum(fits$counts)
# sum(all.data$Count)
# raw <- hist(all.data$Count, breaks = -0.5:80.5)
# raw$counts <- sqrt(raw$counts)   #sqrt frequencies for plotting
# 
# 
# # plot histograms side by side
# par(mfrow = c(1, 2))
# plot(raw, ylim = c(0, 30), ylab = expression(sqrt(Frequency)),
#      xlab = "Captures", main = "Obsevered Data")
# plot(fits, xlim = c(0, 40), ylim = c(0, 30), ylab = expression(sqrt(Frequency)),
#      xlab = "Captures", main = "Predicted Data")



x <- fit8$model$Count
# rootogram to assess model fit
rootogram(~x, dfun = function(x)
          dnbinom(x, size = fit8$theta, mu = exp(coef(fit8)[1])),
                  probability = F, transormation = sqrt)
rootogram(fit8)
mean(poipred$SE)
mean(nbpred$SE)
# Here you see the 'danger' of ignoring overdispersion in the Poisson model. 
# The SE estimates are lower for the Poisson model than for the negative binomial model, 
# which increases the likelihood of incorrectly detecting a significant effect in the Poisson model.




# check residual plots
simulationOutput <- simulateResiduals(fittedModel = fit8, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput)


testZeroInflation(simulationOutput = simulationOutput)

plotResiduals(all.data$ppt,  simulationOutput$scaledResiduals)
plotResiduals(all.data$tmin,  simulationOutput$scaledResiduals)



# visreg(fit2, "ppt", by = "tmin", scale = "response")
# 
# visreg(fit2, "tmin", type = "conditional", scale = "response", gg = T, partial = T)
# visreg(fit2, "ppt", type = "conditional", scale = "response", by = "tmin", gg = T, partial = T)
# summary(fit2)
# 
# visreg2d(fit3, "ppt", "tmin", xlab = "Precipitation (mm)", ylab = "Minimum Temperature", main = "")
# 
# visreg(fit8, "jds", partial = T, by = "tmin")
# visreg2d(fit2, "ppt", "jds")
# visreg2d(fit3, "ppt", "tmin")

visreg(fit8, "jds", axes = F, xlab = "Day of Year", ylab = "Relative Count")
visreg2d(fit8, "ppt", "tmin", scale = "response")
visreg2d(fit2, "ppt", "tmin", scale = "response", key.axes = axis(2, labels = F), 
         main = "", xlab = "", ylab = "")
mtext("Precipitation / mm", side = 1, line = 3, cex = 1.5, adj = 0.2)
mtext(expression(paste("Minimum Temperature / ",degree,"C")), side = 2, line = 3, cex = 1.5, las = 0)


visreg(fit2, scale = "response")

visreg2d(fit2, "ppt", "tmin", scale = "response", key.axes = axis(2, labels = T), 
         main = "", xlab = "", ylab = "", type = "contrast", nn = 100, 
         col = colorRampPalette(brewer.pal(9,"Greys"))(28))

#col = gray.colors(20)

# visreg(fit3, "ppt", partial = F, by = "tmin", scale = "response")
# visreg2d(fit3, "ppt", "tmin", scale = "response")
# visreg(fit3, scale = "response")
# summary(fit3)
# 
# visreg(fit3, "ppt", scale = "response")
# 
# 
# visreg2d(fit3, "jds", "ppt", scale = "response")
# visreg2d(fit3, "jds", "tmin", scale = "response")
# visreg2d(fit3, "ppt", "year", scale = "response")
# visreg2d(fit3, "tmin", "year", scale = "response")

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