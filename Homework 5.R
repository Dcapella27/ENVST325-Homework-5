library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)
library(forecast)
library(jtools)
library(lubridate)
library(patchwork)

# Total inches of water evapo that occurred in one month
ETdat <- read.csv("/cloud/project/activity06/ETdata.csv")

# Greenhouse gas emssisions from anthropogenic land use
ghg <- read.csv("/cloud/project/activity05/Deemer_GHG_Data.csv")

# ETall <- ETdata %>%
#   group_by(date, crop) %>%
#   summarise( ET.in = mean(Ensemble.ET))
# 
# Pistachios <- ETall %>%
#   filter(crop == "Pistachios")
# 
# Pistachios_ts <- na.omit(ts(Pistachios$ET.in,
#                  start = c(2016, 1),
#                  frequency = 12))
# Pistachios_ts
# Pistachios_dec <- decompose(Pistachios_ts)
# plot(Pistachios_dec)
# Pistachios_dec
# 
# acf(Pistachios_ts, lag.max = 24)
# pacf.plot <- pacf(na.omit(Pistachios_ts))
# pacf.plot
# 
# model1 <- arima(Pistachios_ts, 
#                 order = c(1, 0, 0))
# # model1
# model3 <- arima(Pistachios_ts, 
#                 order = (3, 0, 0))
# model3
# 
# AR_fit1 <- Pistachios_ts - residuals(model1) 
# AR_fit4 <- Pistachios_ts - residuals(model4)
# #plot data
# plot(almond_y)
# # plot fit
# points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
# points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
# legend("topleft", c("data","AR1","AR4"),
#        lty=c(1,2,2), lwd=c(1,2,2), 
#        col=c("black", "tomato3","darkgoldenrod4"),
#        bty="n")

# Question 1
ghg$transformedCo2 <- 1/(ghg$co2 + 1000)
# log transform carbon dioxoide fluxes
ghg$log.mean.depth <- log(ghg$mean.depth)
ghg$log.age <- log(ghg$age)
ghg$log.chlorophyll.a <- log(ghg$chlorophyll.a)
ghg$log.surface.area <- log(ghg$surface.area)
ghg$log.volume <- log(ghg$volume + 1)
ghg$log.DIP <- log(ghg$DIP + 1)
ghg$log.precipitation <- log(ghg$precipitation)
ghg$log.runoff <- log(ghg$runoff + 1)
ghg$Region <- as.factor(ghg$Region)
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(transformedCo2 ~ log.mean.depth + log.age + 
                 log.chlorophyll.a + surface.area + 
                 log.volume + DIP + runoff + log.precipitation + 
                 airTemp + Region, data=ghg) 
summary(mod.full)
res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)
# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)
# shapiro-wilks test
shapiro.test(res.full)
plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)
reg.data <- data.frame(ghg$transformedCo2, ghg$log.mean.depth, ghg$log.age,
                       ghg$log.chlorophyll.a, ghg$surface.area,
                       ghg$log.volume, ghg$DIP, ghg$runoff, ghg$log.precipitation, 
                       ghg$airTemp)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)
# run stepwise
summ(mod.full)

pistachios <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(pistachios, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

pistachios_ts <- ts(pistachios$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit
# decompose almond ET time series
pistachios_dec <- decompose(pistachios_ts)
# plot decomposition
plot(pistachios_dec)

corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(corn, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

corn_ts <- ts(corn$ET.in, # data
                    start = c(2016,1), #start year 2016, month 1
                    #first number is unit of time and second is observations within a unit
                    frequency= 12) # frequency of observations in a unit
# decompose almond ET time series
corn_dec <- decompose(corn_ts)
# plot decomposition
plot(corn_dec)

fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(fallow, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

fallow_ts <- ts(fallow$ET.in, # data
                    start = c(2016,1), #start year 2016, month 1
                    #first number is unit of time and second is observations within a unit
                    frequency= 12) # frequency of observations in a unit
# decompose almond ET time series
fallow_dec <- decompose(fallow_ts)
# plot decomposition
plot(fallow_dec)

acf(na.omit(pistachios_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
acf(na.omit(corn_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
acf(na.omit(fallow_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(pistachios_ts))
pacf.plot <- pacf(na.omit(corn_ts))
pacf.plot <- pacf(na.omit(fallow_ts))

pistachios_y <- na.omit(pistachios_ts)
pistachios.mod <- arima(pistachios_y , # data 
                order = c(3,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
pistachios.mod

corn_y <- na.omit(corn_ts)
corn.mod <- arima(corn_y , # data 
                        order = c(2,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
corn.mod

fallow_y <- na.omit(fallow_ts)
fallow.mod <- arima(fallow_y , # data 
                        order = c(3,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
fallow.mod


pistachios_fit <- pistachios_y - residuals(pistachios.mod) 
#plot data
plot(pistachios_y)
# plot fit
points(pistachios_fit, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR3"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

fallow_fit <- fallow_y - residuals(fallow.mod) 
#plot data
plot(fallow_y)
# plot fit
points(fallow_fit, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR3"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

corn_fit <- corn_y - residuals(corn.mod) 
#plot data
plot(corn_y)
# plot fit
points(corn_fit, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR2"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newPistachios <- forecast(pistachios.mod)
newPistachios
newCorn <- forecast(corn.mod)
newCorn
newFallow <- forecast(fallow.mod)
newFallow


newPistachiosF <- data.frame(newPistachios)
newFallowF <- data.frame(newFallow)
newCornF <- data.frame(newCorn)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newFallowF$dateF <- ymd(paste(years,"/",month,"/",1))
newCornF$dateF <- ymd(paste(years,"/",month,"/",1))
newPistachiosF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
pistachiosP <-ggplot() +
  geom_line(data = pistachios, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachios$date[1]),newPistachiosF$dateF[24])+  # Plotting original data
  geom_line(data = newPistachiosF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newPistachiosF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title = "Pistachios")

# make a plot with data and predictions including a prediction interval
cornP <- ggplot() +
  geom_line(data = corn, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(corn$date[1]),newCornF$dateF[24])+  # Plotting original data
  geom_line(data = newCornF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newCornF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title = "Corn" )

# make a plot with data and predictions including a prediction interval
fallowP <- ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newFallowF$dateF[24])+  # Plotting original data
  geom_line(data = newFallowF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newFallowF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title = "Fallow")


pistachiosP
cornP
fallowP
