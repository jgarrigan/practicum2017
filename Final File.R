if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("plot3D")) {
  install.packages("plot3D", dependencies = TRUE)
  library(plot3D)
}
if (!require("xlsx")) {
  install.packages("xlsx", dependencies = TRUE)
  library(xlsx)
}
if (!require("zoo")) {
  install.packages("zoo", dependencies = TRUE)
  library(zoo)
}
if (!require("lubridate")) {
  install.packages("lubridate", dependencies = TRUE)
  library(lubridate)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}
if (!require("mtsdi")) {
  install.packages("mtsdi", dependencies = TRUE)
  library(mtsdi)
}
if (!require("VIM")) {
  install.packages("VIM", dependencies = TRUE)
  library(VIM)
}
if (!require("MissMech")) {
  install.packages("MissMech", dependencies = TRUE)
  library(MissMech)
}
if (!require("mice")) {
  install.packages("mice", dependencies = TRUE)
  library(mice)
}
if (!require("xts")) {
  install.packages("xts", dependencies = TRUE)
  library(xts)
}
#for calculating the cross-corr matrix
if (!require("MTS")) {                   
  install.packages("MTS", dependencies = TRUE)
  library(MTS)
}
#MODWT transform function
if (!require("wavelets")) {
  install.packages("wavelets", dependencies = TRUE)
  library(wavelets)
}
#visualising missing data
if (!require("Amelia")) {
  install.packages("Amelia", dependencies = TRUE)
  library(Amelia)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("GGally")) {
  install.packages("GGally", dependencies = TRUE)
  library(GGally)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("readxl")) {
  install.packages("readxl", dependencies = TRUE)
  library(readxl)
}
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}
if (!require("W2CWM2C")) {
  install.packages("W2CWM2C", dependencies = TRUE)
  library(W2CWM2C)
}
if (!require("colorspace")) {
  install.packages("colorspace", dependencies = TRUE)
  library(colorspace)
}
if (!require("waveslim")) {
  install.packages("waveslim", dependencies = TRUE)
  library(waveslim)
}
if (!require("FGN")) {
  install.packages("FGN", dependencies = TRUE)
  library(FGN)
}
if (!require("wavemulcor")) {
  install.packages("wavemulcor", dependencies = TRUE)
  library(wavemulcor)
}
if (!require("ggseas")) {
  install.packages("ggseas", dependencies = TRUE)
  library(ggseas)
}

setwd("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data")

library(readxl)
LX0088_Raw <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/LX0088 Raw.xlsx", 
                         col_types = c("skip", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric", 
                                       "numeric"), na = "NA")

LX0088_Raw <- as.data.frame(LX0088_Raw)

LX0088.Neighbour.Hrly.RTWP <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Hourly Final with NA.xlsx", 
                         col_types = c("skip", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric", 
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric",
                                       "numeric", "numeric")
                         , na = "NA")


LX0088.Neighbour.Raw.RTWP <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Raw Final with NA.xlsx",
                                        col_types = c("skip", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric", "numeric",
                                                      "numeric", "numeric")
                                        , na = "NA")

RNC20_30Hourly <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/RNC20_30Hourly.xlsx")

DO0182U09A3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09A3_24hrs.xlsx", 
                               header = TRUE,
                               sheetIndex = 1)

DO0182U09B3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09B3_24hrs.xlsx",
                               header = TRUE,
                               sheetIndex = 1)

DO0182U09C3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09C3_24hrs.xlsx",
                               header = TRUE,
                               sheetIndex = 1)

Raw_Normalised_RTWP_Finally_Got_There_020817 <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Raw Normalised RTWP Finally Got There 020817.xlsx",sheet = "Normalised")


#visualisation of the missing values for LX0088.Neighbour.Raw.RTWP
par(mar=c(7,5,5,5))
missmap(LX0088.Neighbour.Raw.RTWP)

#Boxplot for LX0088 

boxplot(LX0088_Raw, 
        horizontal = TRUE,
        par(mar=c(4,10,2,1)),
        las=2,
        xlab = "RTWP Levels (dBm)",
        main = "Boxplot of RTWP values for LX0088",
        col = rainbow(9),
        cex.main = 1.5,
        cex.lab = 1.5,
        cex.axis = 1.5)

#Correlation plot between variables for LX0088
M <- cor(LX0088_Raw)
corrplot(M, 
         method = "number", 
         type = "lower",
         title = "Correlation plot between variables for LX0088",
         mar=c(0,0,1,0))

#cross-correlation plot for LX0088
par(mfrow=c(3,1),mar=c(5,4,4,4),oma=c(0,0,0,0))
ccfU09A3_U09B3 <- ccf(LX0088_Raw$LX0088U09A3,LX0088_Raw$LX0088U09B3,
                      type = "correlation",
                      plot = TRUE,
                      lag.max = 1200,
                      main = "CCF for LX0088U09A3 & LX0088U09B3",
                      ylab = "CCF",
                      xlab = "Lag (15 min intervals)")

ccfU09A3_U09C3 <- ccf(LX0088_Raw$LX0088U09A3,LX0088_Raw$LX0088U09C3,
                      type = "correlation",
                      plot = TRUE,
                      lag.max = 1200,
                      main = "CCF for LX0088U09A3 & LX0088U09C3",
                      ylab = "CCF",
                      xlab = "Lag (15 min intervals)")

ccfU09A3_U09C3 <- ccf(LX0088_Raw$LX0088U09B3,LX0088_Raw$LX0088U09C3, 
                      type = "correlation",
                      plot = TRUE, 
                      lag.max = 1200,
                      main = "CCF for LX0088U09B3 & LX0088U09C3",
                      ylab = "CCF",
                      xlab = "Lag (15 min intervals)")

par(mfrow=c(1,1))
my.date <- seq(as.POSIXlt(strptime('20/07/2017 18:15',"%d/%m/%Y %H:%M"),tz="GMT"), as.POSIXlt(strptime('23/07/2017 21:15',"%d/%m/%Y %H:%M"),tz="GMT"),length.out = 7)
my.date <- format(my.date, "%a %H:%M")
#RTWP vs. Time plot for LX0088
plot(LX0088_Raw[500:800,1],
     ylab="RTWP (dBm)",
     xlab="Time",
     type="l",
     col = "#e41a1c",
     lwd = 3,
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5,
     xaxt='n',
     ylim = c(-107, -85),
      main = "LX0088 RTWP vs. Time")
lines(LX0088_Raw[500:800,2] , col = "#377eb8",lwd=3)
lines(LX0088_Raw[500:800,3] , col = "#4daf4a",lwd=3)
legend(0,-85,
       ncol = 1,
       cex = 1.5,
       legend = c("Sec A - U900MHz",
                  "Sec B - U900MHz",
                  "Sec C - U900MHz"),
                  #"F3 - U2100MHz"),
       fill = c("#e41a1c","#377eb8","#4daf4a"),
       title = "Sectors of LX0088")
axis(1, at=c(0,50,100,150,200,250,300),labels = my.date, cex.axis=1.5)

#autocorrelation plot to determine if process is stationary

LX0088U09A3.diff1 <- diff(LX0088_Raw$LX0088U09A3,differences = 1)
LX0088U21A1.diff1 <- diff(LX0088_Raw$LX0088U21A1,differences = 1)

par(mar=c(5,5,5,3), oma = c(0,0,0,0),mfrow = c(2,1))
acf(LX0088U09A3.diff1, lag.max = 96,
    main = "Autocorrelation Function Plot of LX0088U09A3",
    xlab = "Lag index (15 mins)",
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 1.5)
legend(90,1.1,ncol = 1, cex = 1.5,legend = "(a)",box.lwd = 0,box.col = "white", bg = "transparent")


acf(LX0088U21A1.diff1, lag.max = 96,
    main = "Autocorrelation Function Plot of LX0088U21A1",
    xlab = "Lag index (15 mins)",
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 1.5)
legend(90,1.1,ncol = 1, cex = 1.5,legend = "(b)",box.lwd = 0,box.col = "white", bg = "transparent")

par(mar=c(5,5,5,3), oma = c(0,0,0,0),mfrow = c(2,1))
pacf(LX0088U09A3.diff1, lag.max = 100,
     main = "Partial Autocorrelation Function Plot of LX0088U09A3",
     xlab = "Lag index (15 mins)",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5)
legend(90,-0.05,ncol = 1, cex = 1.5,legend = "(a)",box.lwd = 0,box.col = "white", bg = "transparent")

pacf(LX0088U21A1.diff1, lag.max = 100,
     main = "Partial Autocorrelation Function Plot of LX0088U21A1",
     xlab = "Lag index (15 mins)",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5)
legend(90,-0.10,ncol = 1, cex = 1.5,legend = "(b)",box.lwd = 0,box.col = "white", bg = "transparent")


# acf(LX0088_Raw$LX0088U21A1, lag.max = 96,
#     main = "Autocorrelation Function Plot of LX0088U21A1",
#     xlab = "Lag index (15 mins)")

# #plot of the partial autocorrelation function
# pacf(LX0088_Raw$LX0088U09A3, lag.max = 672,
#      main = "Partial Autocorrelation Function Plot of LX0088U09A3",
#      xlab = "Lag index (15 mins)")


# date.range <- parse_date_time(c("15/07/2017 13:45","29/07/2017 13:15"), "%d%m%y HM")
# my.dates <- seq(date.range[1], date.range[2], "15 min")

# LX0088U09A3.datetime <- data.frame(
#   x = as.numeric(1:1325),
#   y = as.numeric(LX0088_Raw$LX0088U09A3)) 
# 
# colnames(LX0088U09A3.datetime) <- c("Index","RTWP")
# 
# LX0088U09A3.datetime %>%
# ggsdc(aes(x=LX0088U09A3.datetime[,1],y=LX0088U09A3.datetime[,2],colour="red"),
#       frequency = 96,method = "stl",
#       facet.titles = c("Original Time Series", "Underlying Trend",
#                        "Regular Seasonal Impacts", "Residual Randomness"))
  

#season decompostion of LX0088
LX0088U09A3ts <- ts(as.numeric(LX0088_Raw$LX0088U09A3), frequency = 96)
LX0088U09A3.stl <- stl(LX0088U09A3ts, s.window = "periodic", robust = TRUE)
plot(LX0088U09A3.stl, 
     main = "Seasonal Decompostion for LX0088U09A3",
     col.range = "red",
     set.pars = list(cex = 10,
                     cex.axis = 1,
                     cex.lab = 1, 
                     cex.main = 1, 
                     mar = c(0, 6, 0, 6), 
                     oma = c(6, 0, 4, 0), 
                     tck = -0.02,
                     lwd=2, 
                     mfrow = c(4,1)))

LX0088U09B3ts <- ts(as.numeric(LX0088_Raw$LX0088U09B3), frequency = 96)
LX0088U09B3.stl <- stl(LX0088U09B3ts, s.window = "periodic", robust = TRUE)
plot(LX0088U09B3.stl, main = "Seasonal Decompostion for LX0088U09B3")

LX0088U09C3ts <- ts(as.numeric(LX0088_Raw$LX0088U09C3), frequency = 96)
LX0088U09C3.stl <- stl(LX0088U09C3ts, s.window = "periodic", robust = TRUE)
plot(LX0088U09C3.stl, main = "Seasonal Decompostion for LX0088U09C3")

#imputation of the hourly RTWP for the neighbours of LX0088

LX0088.Neigh.Hrly.RTWP.withNA.ts <- ts(LX0088.Neighbour.Hrly.RTWP, frequency = 24)
#imputation of the raw RTWP for the neighbours of LX0088
LX0088.Neigh.Raw.RTWP.withNA.ts <- ts(LX0088.Neighbour.Raw.RTWP, frequency = 96)

# #read in the test data containing the original and modified copy.
# Imputation_Test <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Imputation Test.xlsx", 
#                               na = "NA")
# 
# Imputation_Test <- ts(Imputation_Test,frequency = 96)
# 
# variable2bImputed <- ~KD0027U09A3+KD0027U09C3+KD0046U09A3+KD0046U09C3+KD0103U09B3+KD0103U09C3+KD0166U09C3+KD0177U09C3+LX0017U09B3+LX0022U09A3+LX0022U09B3+LX0022U09C3+LX0023U09A3+LX0023U09B3+LX0060U09A3+LX0060U09B3+LX0060U09C3+LX0073U09A3+LX0074U09A3+LX0074U09B3+LX0074U09C3+LX0079U09A3+LX0080U09A3+LX0080U09B3+LX0080U09C3+LX0088U09A3+LX0088U09B3+LX0088U09C3+OF0028U09B3+OF0061U09B3+KD0056U09A3+KK0125U09A3+LX0017U09A3+LX0017U09C3+LX0025U09A3+LX0096U09A3+LX0097U09A3+Modified
# 
# imputed.spline <- mnimput(variable2bImputed,
#                           Imputation_Test,
#                           eps=1e-3,
#                           ts=TRUE,
#                           method="spline")
# 
# arcontrol<-list(order=cbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0)), period=NULL)
# imputed.arima <- mnimput(variable2bImputed,
#                           Imputation_Test,
#                           eps=1e-3,
#                           ts=TRUE,
#                           method="arima",
#                        ar.control = arcontrol)
# 
# imputed.gam <- mnimput(variable2bImputed,
#                           Imputation_Test,
#                           eps=1e-3,
#                           ts=TRUE,
#                           method="gam",
#                           ga.control=list())
# 
# rmse <- function(error)
# {
#   sqrt(mean(error^2))
# }
# 
# # Function that returns Mean Absolute Error
# mae <- function(error)
# {
#   mean(abs(error))
# }
# 
# predicted.data.spline <- imputed.spline$filled.dataset$Modified
# original.data <-  Imputation_Test[,"LX0096U09A3"] 
#   
# # Calculate error
# prediction.error.spline <- original.data - predicted.data.spline
# 
# rmse(prediction.error.spline)
# mae(prediction.error.spline)


missing.colnames.Hourly <- ~KD0027U09A3+KD0027U09C3+KD0046U09A3+KD0046U09C3+KD0103U09B3+KD0103U09C3+KD0166U09C3+KD0177U09C3+LX0017U09B3+LX0022U09A3+LX0022U09B3+LX0022U09C3+LX0023U09A3+LX0023U09B3+LX0060U09A3+LX0060U09B3+LX0060U09C3+LX0073U09A3+LX0074U09A3+LX0074U09B3+LX0074U09C3+LX0079U09A3+LX0080U09A3+LX0080U09B3+LX0080U09C3+LX0088U09A3+LX0088U09B3+LX0088U09C3+OF0028U09B3+OF0061U09B3
missing.colnames.raw <- ~KD0027U09A3+KD0027U09C3+KD0046U09A3+KD0046U09C3+KD0103U09B3+KD0103U09C3+KD0166U09C3+KD0177U09C3+LX0017U09B3+LX0022U09A3+LX0022U09B3+LX0022U09C3+LX0023U09A3+LX0023U09B3+LX0060U09A3+LX0060U09B3+LX0060U09C3+LX0073U09A3+LX0074U09A3+LX0074U09B3+LX0074U09C3+LX0079U09A3+LX0080U09A3+LX0080U09B3+LX0080U09C3+LX0088U09A3+LX0088U09B3+LX0088U09C3+OF0028U09B3+OF0061U09B3+KD0056U09A3+KK0125U09A3+LX0017U09A3+LX0017U09C3+LX0025U09A3+LX0096U09A3+LX0097U09A3

imputed.Hourly <- mnimput(missing.colnames.Hourly,
                          LX0088.Neigh.Hrly.RTWP.withNA.ts,
                          eps=1e-3,
                          ts=TRUE,
                          method="spline")

imputed.Raw <- mnimput(missing.colnames.raw,
                       LX0088.Neigh.Raw.RTWP.withNA.ts,
                       eps=1e-3,
                       ts=TRUE,
                       method="spline")

#summary(imputed.Hourly)
summary(imputed.Raw)


imputed.Hourly.Vars <- imputed.Hourly$filled.dataset
imputed.Raw.Vars <- imputed.Raw$filled.dataset

original.hrly.comp.vars <- LX0088.Neighbour.Hrly.RTWP[,c("KD0056U09A3",
                                                          "KK0125U09A3",
                                                          "LX0017U09A3",
                                                          "LX0017U09C3",
                                                          "LX0025U09A3",
                                                          "LX0096U09A3",
                                                          "LX0097U09A3")]

original.hrly.comp.vars.ts <- ts(original.hrly.comp.vars, frequency = 24)

hourly.data.NoNA <- ts.union(original.hrly.comp.vars.ts,imputed.Hourly.Vars, dframe = TRUE)

col.names.original <- colnames(original.hrly.comp.vars.ts)
col.names.imputed <- colnames(imputed.Hourly.Vars)

colnames(hourly.data.NoNA)[1:7] <- col.names.original
colnames(hourly.data.NoNA)[8:37] <- col.names.imputed

hourly.data.NoNA.ts <- ts(hourly.data.NoNA, frequency = 24)
raw.data.NoNA.ts <- ts(imputed.Raw.Vars, frequency = 96)

#distribution of Hurst Exponent for RNC20_30
RNC20_30Hourly <- RNC20_30Hourly[-1]

RNC20_30TS <- ts(RNC20_30Hourly, frequency = 25)

HurstExp <- apply(RNC20_30TS, 2, HurstK)

HurstExp <- as.data.frame(HurstExp)

par(mfrow=c(1,1), las=0, mar=c(4,4,3,2)+0.2)

tmp <- density(HurstExp[,1], na.rm = T)
hist(HurstExp[,1],
     #prob = FALSE,
     prob = TRUE, 
     ylim = c(0, max(tmp$y)),
     main = "Distribution of Hurst Exponent Estimates for RNC20_30",
     xlab = "Hurst Exponent Estimate Value",
     freq = FALSE,
     #freq = TRUE,
     col="lightgreen",
     cex.main = 1.5,
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(tmp, col = "darkblue", lwd = 2, cex = 1.5)

library(e1071)
test.skewness <- apply(RNC20_30Hourly,2,skewness,na.rm=TRUE)
test.kurtosis <- apply(RNC20_30Hourly,2,kurtosis,na.rm=TRUE)
par(mfrow=c(2,1), oma=c(0,5,0,1))
hist(test.skewness,
     main = "Distribution of Skewness - Measure of Symmetry",
     xlab = "Skewness Value",
     col = rainbow(1),
     cex.main = 1.5,
     cex.lab = 1.5)
legend(11,2500,ncol = 1,cex = 1.5, legend = "(a)", box.lwd = 0, box.col = "white", bg = "transparent")


hist(test.kurtosis,
     main = "Distribution of Kurtosis - Measure of Tail Shape",
     xlab = "Kurtosis Value",
     col = rainbow(1),
     cex.main = 1.5,
     cex.lab = 1.5)
legend(150,2500,ncol = 1,cex = 1.5, legend = "(b)", box.lwd = 0, box.col = "white", bg = "transparent")

#intra-day correlation plots for DO0182U09A3
# DO0182U09A3_24hrs_Mat <- as.matrix(DO0182U09A3_24hrs)
# DO0182U09B3_24hrs_Mat <- as.matrix(DO0182U09B3_24hrs)
# DO0182U09C3_24hrs_Mat <- as.matrix(DO0182U09C3_24hrs)
# 
# DO0182U09A3_24hrCor <- cor(DO0182U09A3_24hrs_Mat)
# DO0182U09B3_24hrCor <- cor(DO0182U09B3_24hrs_Mat)
# DO0182U09C3_24hrCor <- cor(DO0182U09C3_24hrs_Mat)
# 
LX0088U09A3_Daily_Data <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/LX0088U09A3 Daily Data.xlsx")
LX0088U09A3_Daily_Data <- as.matrix(LX0088U09A3_Daily_Data)
LX0088U09A3_24hrCor <- cor(LX0088U09A3_Daily_Data)

LX0088U21A1_Daily_Data <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/LX0088U21A1 Daily Data.xlsx")
LX0088U21A1_Daily_Data <- as.matrix(LX0088U21A1_Daily_Data)
LX0088U21A1_24hrCor <- cor(LX0088U21A1_Daily_Data)

par( ps=16,mfrow=c(1,2))
corrplot(LX0088U09A3_24hrCor, 
         method = "number", 
         type = "lower",
         mfrow=c(1,1),
         oma=c(0.1,0.1,0.1,0.1),
         mar=c(0,0,1,0),
         title = "Correlation plot between days for LX0088U09A3")

par( ps=16)
corrplot(LX0088U21A1_24hrCor, 
         method = "number", 
         type = "lower",
         mfrow=c(1,1),
         oma=c(0.1,0.1,0.1,0.1),
         mar=c(0,0,1,0),
         title = "Correlation plot between days for LX0088U21A1")

# my.date.corplot <- seq(as.POSIXlt(strptime('15/07/2017 13:45',"%d/%m/%Y %H:%M"),tz="GMT"), as.POSIXlt(strptime('29/07/2017 13:15',"%d/%m/%Y %H:%M"),tz="GMT"),length.out = 14)
# my.date.corplot <- format(my.date.corplot, "%a %H:%M")
my.new.labels <- c("Day 1","Day 2","Day 3","Day 4","Day 5","Day 6","Day 7","Day 8","Day 9","Day 10","Day 11","Day 12","Day 13","Day 14")
par(mfrow=c(1,1),mar=c(10,5,3,5))
plot(LX0088_Raw[,"LX0088U09A3"],
     type="l",
     col="#e41a1c",
     xlab="",
     ylab="RTWP (dBm)",
     main="RTWP Plot for LX0088U21A1 and LX0088U09A3",
     cex.lab=1.5,
     xaxt='n',
     cex.axis=1.5,
     cex.main=1.5,
     lwd=2)
lines(LX0088_Raw[,"LX0088U21A1"],
     type="l",
     col="#377eb8",
     lwd=2)
axis(1, at=c(0,96,192,288,384,480,576,672,768,864,960,1056,1152,1248),labels = my.new.labels, cex.axis=1.5,las=2)
mtext("Days",side = 1, line = 8,cex = 1.5)
legend(600,-87.5,
       ncol = 2,
       cex = 1.5,
       legend = c("Sec A - U900MHz",
                  "Sec A - U2100MHz"),
       fill = c("#e41a1c","#377eb8"),
       title = "Sector A for LX0088")

# par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(1,3,2,1))
# 
# corrplot(DO0182U09B3_24hrCor, 
#          method = "number", 
#          type = "lower",
#          mfrow=c(1,1),
#          oma=c(0,0,0,0),
#          mar=c(0,0,1.5,0), 
#          title = "Correlation plot between days for DO0192U09B3")
# 
# corrplot(DO0182U09C3_24hrCor, 
#          method = "number", 
#          type = "lower",
#          mfrow=c(1,1),
#          oma=c(0,0,0,0),
#          mar=c(0,0,1.5,0), 
#          title = "Correlation plot between days for DO0192U09C3")

#compute the multiple wavelet correlation between variables
#hourly.data.NoNA.ts.df <- diff(hourly.data.NoNA.ts, differences = 1)

LX0088U09A3 <- raw.data.NoNA.ts[,"LX0088U09A3"]
LX0080U09C3 <- raw.data.NoNA.ts[,"LX0080U09C3"]
LX0022U09A3 <- raw.data.NoNA.ts[,"LX0022U09A3"]
LX0017U09B3 <- raw.data.NoNA.ts[,"LX0017U09B3"]
KD0177U09C3 <- raw.data.NoNA.ts[,"KD0177U09C3"]
KD0046U09C3 <- raw.data.NoNA.ts[,"KD0046U09C3"]
OF0028U09B3 <- raw.data.NoNA.ts[,"OF0028U09B3"]
OF0061U09B3 <- raw.data.NoNA.ts[,"OF0061U09B3"]
KD0027U09A3 <- raw.data.NoNA.ts[,"KD0027U09A3"]
LX0080U09A3 <- raw.data.NoNA.ts[,"LX0080U09A3"]

wf <- "la8"
J <- 7

LX0088U09A3.modwt <- waveslim::modwt(LX0088U09A3, wf, J,boundary = "periodic")
LX0088U09A3.modwt.bw <- brick.wall(LX0088U09A3.modwt, wf)
LX0080U09C3.modwt <- waveslim::modwt(LX0080U09C3, wf, J,boundary = "periodic")
LX0080U09C3.modwt.bw <- brick.wall(LX0080U09C3.modwt, wf)
LX0022U09A3.modwt <- waveslim::modwt(LX0022U09A3, wf, J,boundary = "periodic")
LX0022U09A3.modwt.bw <- brick.wall(LX0022U09A3.modwt, wf)
LX0017U09B3.modwt <- waveslim::modwt(LX0017U09B3, wf, J,boundary = "periodic")
LX0017U09B3.modwt.bw <- brick.wall(LX0017U09B3.modwt, wf)
KD0177U09C3.modwt <- waveslim::modwt(KD0177U09C3, wf, J,boundary = "periodic")
KD0177U09C3.modwt.bw <- brick.wall(KD0177U09C3.modwt, wf)
KD0046U09C3.modwt <- waveslim::modwt(KD0046U09C3, wf, J,boundary = "periodic")
KD0046U09C3.modwt.bw <- brick.wall(KD0046U09C3.modwt, wf)
OF0028U09B3.modwt <- waveslim::modwt(OF0028U09B3, wf, J,boundary = "periodic")
OF0028U09B3.modwt.bw <- brick.wall(OF0028U09B3.modwt, wf)
OF0061U09B3.modwt <- waveslim::modwt(OF0061U09B3, wf, J,boundary = "periodic")
OF0061U09B3.modwt.bw <- brick.wall(OF0061U09B3.modwt, wf)
KD0027U09A3.modwt <- waveslim::modwt(KD0027U09A3, wf, J,boundary = "periodic")
KD0027U09A3.modwt.bw <- brick.wall(KD0027U09A3.modwt, wf)
LX0080U09A3.modwt <- waveslim::modwt(LX0080U09A3, wf, J,boundary = "periodic")
LX0080U09A3.modwt.bw <- brick.wall(LX0080U09A3.modwt, wf)


xx <- list(LX0088U09A3.modwt.bw,
           LX0080U09C3.modwt.bw,
           LX0022U09A3.modwt.bw,
           LX0017U09B3.modwt.bw,
           KD0177U09C3.modwt.bw,
           KD0046U09C3.modwt.bw,
           OF0028U09B3.modwt.bw,
           OF0061U09B3.modwt.bw,
           KD0027U09A3.modwt.bw,
           LX0080U09A3.modwt.bw)

Lst <- wave.multiple.correlation(xx, N = length(xx[[1]][[1]]))#, ymaxr = 1)

LX0088.modwt.cor <- Lst$xy.mulcor[1:J,]
YmaxR <- Lst$YmaxR

cell.names <- c("Cell 1",
                "Cell 2", 
                "Cell 3", 
                "Cell 4", 
                "Cell 5", 
                "Cell 6",
                "Cell 7",
                "Cell 8",
                "Cell 9",
                "Cell 10")

par(mfrow=c(1,1), oma=c(1,1,3,0), las=0, mar=c(5,4,4,2)+.1)
matplot(2^(0:(J-1)), LX0088.modwt.cor[-(J+1),], type="b",
        log="x", pch="*LU", xaxt="n", lty=1, col=c(1,4,4),
        xlab="Wavelet Scale",
        ylab = "Wavelet Multiple Correlation",
        main="Estimate of Wavelet Multiscale Multiple Correlation",
        cex = 1.5,
        cex.axis=1.5,
        cex.lab=1.5,
        cex.main=1.5)
axis(side=1, at=2^(0:7),cex.axis=1.5)
abline(h=0)
text(2^(0:7), min(LX0088.modwt.cor[-(J+1),])+0.03,
     labels=cell.names[YmaxR], adj=0.5, cex=1.5)
#axis(side=3, at=2^(0:7),cex.axis=1.5)
#mtext("Wavelet Scale", ylab="Wavelet Multiple Correlation", side = 3, outer = TRUE,line =3)

#wavelet local multiple correlation

N <- length(LX0080U09C3)
b <-trunc(N/3)
t1<- 1:b
t2<- (b+1):(2*b)
t3<- (2*b+1):N

wf<- "la8"
M <- N/2^3 #defines the length of the rolling window
window <- "uniform" #type of weight function 

local.J <- 7#trunc(log2(N))-3 #sets the decomposition level

cor1 <- cor(LX0088U09A3[t1],KD0177U09C3[t1])
cor2 <- cor(LX0088U09A3[t2],KD0177U09C3[t2])
cor3 <- cor(LX0088U09A3[t3],KD0177U09C3[t3])

LX0088U09A3.modwt <- modwt(LX0088U09A3, wf, local.J)
LX0088U09A3.modwt.bw <- brick.wall(LX0088U09A3.modwt, wf)

KD0177U09C3.modwt <- modwt(KD0177U09C3, wf, local.J)
KD0177U09C3.modwt.bw <- brick.wall(KD0177U09C3.modwt, wf)

localxx <- list(LX0088U09A3.modwt.bw,KD0177U09C3.modwt.bw)

xy.mulcor <- wave.local.multiple.correlation(localxx, M, window=window)

val <- as.matrix(xy.mulcor$val)
lo <- as.matrix(xy.mulcor$lo)
up <- as.matrix(xy.mulcor$up)
YmaxR <- as.matrix(xy.mulcor$YmaxR)

#old.par <- par()

scale.names <- paste0("(",c("2-4",
                            "4-8",
                            "8-16",
                            "16-32",
                            "32-64",
                            "64-128",
                            "128-256",
                            "256-512",
                            "512-1024",
                            "1024-2048"),"]")
scale.names <- c(scale.names[1:local.J],"smooth")

title <- paste("Wavelet Local Multiple Correlation")
sub <- paste("first",b,"obs:",round(100*cor1,1),"% correlation;","middle",b,"obs:",
             round(100*cor2,1),"%","rest:",round(100*cor3,1),"%")
xlab <- "time index (15 mins)"
ylab <- "periods"
image2D(z=val, x=1:nrow(val), y=1:ncol(val),
        main=title, #sub=sub,
        xlab=xlab, ylab=ylab, axes=FALSE, clab = expression(varphi),
        rasterImage = TRUE, contour = list(lwd = 2, col = jet.col(11)))
axis(side=1, at=seq(10,nrow(val),by=10), cex.axis=0.75)
axis(side=2, at=1:ncol(val),labels=scale.names, las=1,cex.axis=0.75)

colnames(val)[1:local.J] <- paste0("level",1:local.J)
par(mfrow=c(4,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
for(i in local.J:1) {
  matplot(1:N,val[,i], type="l", lty=1, ylim=c(-1,1), xaxt="n",
          xlab="", ylab="", main=colnames(val)[i])
  if(i<3) {axis(side=1, at=seq(10,N,by=10))}
  #axis(side=2, at=c(-.2, 0, .5, 1))
  lines(lo[,i], lty=1, col=2) ##Add Connected Line Segments to a Plot
  lines(up[,i], lty=1, col=2)
  abline(h=0) ##Add Straight horiz and vert Lines to a Plot
}
par(las=0)
mtext('time', side=1, outer=TRUE, adj=0.5)
mtext('Wavelet Local Multiple Correlation', side=2, outer=TRUE, adj=0.5)

par(old.par)

#wavelet multiple cross-correlation
lmax <- 20 #max lags
n <- length(LX0088U09A3)
lags <- length(-lmax:lmax)

Lst <- wave.multiple.cross.correlation(localxx, lmax)

LX0088.RTWP.cross.cor <- as.matrix(Lst$xy.mulcor[1:local.J,])
YmaxR <- Lst$YmaxR #index number of var whose correlation is calculated against 
#a linear combination of the rest 

cell.names <- c("LX0088U09A3",
                "LX0080U09C3",
                "LX0022U09A3",
                "LX0017U09B3",
                "KD0177U09C3",
                "KD0046U09C3",
                "OF0028U09B3",
                "OF0061U09B3",
                "KD0027U09A3",
                "LX0080U09A3")

rownames(LX0088.RTWP.cross.cor) <- rownames(LX0088.RTWP.cross.cor, do.NULL = FALSE, prefix = "Level ")


#calculates Fishers transformation applied to the sample wavelet multiple 
#correlation to compute the confidence intervals
lower.ci <- tanh(atanh(LX0088.RTWP.cross.cor) - qnorm(0.975) / 
                   sqrt(matrix(trunc(n/2^(1:local.J)), nrow=local.J, ncol=lags)- 3))
upper.ci <- tanh(atanh(LX0088.RTWP.cross.cor) + qnorm(0.975) /
                   sqrt(matrix(trunc(n/2^(1:local.J)), nrow=local.J, ncol=lags)- 3))

par(mfrow=c(4,2), las=1, pty="m", mar=c(4,5,2,1)+.1, oma=c(1.2,1.2,0,0))
for(i in local.J:1) {
  matplot((1:(2*lmax+1)),LX0088.RTWP.cross.cor[i,], type="l", lty=1, ylim=c(-1,1),
          xaxt="n", xlab="", ylab="", main=rownames(LX0088.RTWP.cross.cor)[[i]][1],cex.axis=1.5,cex.main=1.5)
  if(i<3) {axis(side=1, at=seq(1, 2*lmax+1, by=10),
                labels=seq(-lmax, lmax, by=10), cex.axis=1.5)}
  #axis(side=2, at=c(-.2, 0, .5, 1))
  lines(lower.ci[i,], lty=1, col=2) ##Add Connected Line Segments to a Plot
  lines(upper.ci[i,], lty=1, col=2)
  abline(h=0,v=lmax+1) ##Add Straight horiz and vert Lines to a Plot
  text(1,1, labels=cell.names[YmaxR[i]], adj=c(-0,1), cex=1.5)
}
par(las=0)
mtext('Lag Index (15 mins)', side=1, outer=TRUE, adj=0.5, cex = 1.5, line = -1 )
mtext('Wavelet Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5, cex = 1.5, line = -1)

#calcualtes the eigenvalue decompostion of the normalised RTWP data for LX0088 and Neighbours

Raw_Normalised_RTWP_Finally_Got_There_020817 <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Most definitely the last attempt.xlsx",sheet = "Normalised")

RTWP.LX0088.Neighbours <- Raw_Normalised_RTWP_Finally_Got_There_020817[,c("LX0088U09A3",
                                                                          "LX0080U09C3",
                                                                          "LX0022U09A3",
                                                                          "LX0017U09B3",
                                                                          "KD0177U09C3",
                                                                          "KD0046U09C3",
                                                                          "OF0028U09B3",
                                                                          "OF0061U09B3",
                                                                          "KD0027U09A3",
                                                                          "LX0080U09A3")]
RTWP.LX0088.Neighbours.TS <- ts(RTWP.LX0088.Neighbours, frequency = 96)

Cor <- function(x) {
  corr <- cor(x)
  out <- as.data.frame.table(corr)
  with(out, setNames(Freq, paste(Var1, Var2)))
}

#defimes the size of the window to use in the sliding correlation window
windowSize10 <- 10
windowSize20 <- 20
windowSize50 <- 50
windowSize100 <- 100

#calculates the correlation coeff. for all variables for each time period
#by changing the window to less than 10 there will be standard deviation errors
LX0088.Cor.Matrix.window10 <- rollapplyr(RTWP.LX0088.Neighbours.TS,
                                         windowSize10,
                                         Cor,
                                         by.column = FALSE)

LX0088.Cor.Matrix.window20 <- rollapplyr(RTWP.LX0088.Neighbours.TS,
                                         windowSize20,
                                         Cor,
                                         by.column = FALSE)

LX0088.Cor.Matrix.window50 <- rollapplyr(RTWP.LX0088.Neighbours.TS,
                                         windowSize50,
                                         Cor,
                                         by.column = FALSE)

LX0088.Cor.Matrix.window100 <- rollapplyr(RTWP.LX0088.Neighbours.TS,
                                          windowSize100,
                                          Cor,
                                          by.column = FALSE)


#calculate the eigenvalues of the corr matrix at each time interval
#using sliding window
eigen.window10 <- apply(LX0088.Cor.Matrix.window10, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(LX0088.Cor.Matrix.window10))))$values[1:10])

eigen.window20 <- apply(LX0088.Cor.Matrix.window20, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(LX0088.Cor.Matrix.window20))))$values[1:10])

eigen.window50 <- apply(LX0088.Cor.Matrix.window50, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(LX0088.Cor.Matrix.window50))))$values[1:10])

eigen.window100 <- apply(LX0088.Cor.Matrix.window100, 
                         1, 
                         function(x) 
                           eigen(matrix(x,
                                        nrow = sqrt(ncol(LX0088.Cor.Matrix.window100))))$values[1:10])

eigen.window10 <- t(eigen.window10)
eigen.window20 <- t(eigen.window20)
eigen.window50 <- t(eigen.window50)
eigen.window100 <- t(eigen.window100)

#write.xlsx(eigen.window20, "C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/Eigenvalues.GodHelpMe.xlsx")

par(mfrow=c(2,2),oma=c(0,2,0,0),mar=c(5,5,2,1))
plot(eigen.window10[,1], 
     type="l", 
     col="#d7191c",
     lwd=1.5,
     main=expression(paste("(a) ", lambda[1]," 90 min window")),
     xlab="Time Index (15 min intervals)",
     ylab=expression(paste(lambda[1]," magnitude")),
     cex.main = 1.5,
     cex.axis =1.5,
     cex.lab = 1.5)
legend(1200,9,ncol = 1, cex = 1.5,legend = "(a)",box.lwd = 0,box.col = "white", bg = "transparent")
plot(eigen.window20[,1],
     type = "l",
     col="#fdae61",
     lwd=1.5,
     main=expression(paste("(b) ",lambda[1]," 300 min window")),
     xlab="Time Index (15 min intervals)",
     ylab=expression(paste(lambda[1]," magnitude")),
     cex.main = 1.5,
     cex.axis =1.5,
     cex.lab = 1.5)
legend(1200,9.2,ncol = 1, cex = 1.5,legend = "(b)",box.lwd = 0,box.col = "white", bg = "transparent")
plot(eigen.window50[,1],
     type = "l",
     col="#008837",
     lwd=1.5,
     main=expression(paste("(c) ",lambda[1]," 750 min window")),
     xlab="Time Index (15 min intervals)",
     ylab=expression(paste(lambda[1]," magnitude")),
     cex.main = 1.5,
     cex.axis =1.5,
     cex.lab = 1.5)
legend(1150,9,ncol = 1, cex = 1.5,legend = "(c)",box.lwd = 0,box.col = "white", bg = "transparent")
legend(560,7.5,ncol = 1, cex = 2,legend = "*",box.lwd = 0,text.col = "red",box.col = "white", bg = "transparent")
legend(650,9,ncol = 1, cex = 2,legend = "*",box.lwd = 0,text.col = "red",box.col = "white", bg = "transparent")
legend(750,9,ncol = 1, cex = 2,legend = "*",box.lwd = 0,text.col = "red",box.col = "white", bg = "transparent")
plot(eigen.window100[,1],
     type = "l",
     col="#2c7bb6",
     lwd=1.5,
     main=expression(paste("(d) ",lambda[1]," 1500 min window")),
     xlab="Time Index (15 min intervals)",
     ylab=expression(paste(lambda[1]," magnitude")),
     cex.main = 1.5,
     cex.axis =1.5,
     cex.lab = 1.5)
legend(1100,7.8,ncol = 1, cex = 1.5,legend = "(d)",box.lwd = 0,box.col = "white", bg = "transparent")




# #computes the intra day correlation
# #####################Importing 24 hour data for DO0182U09A3#####################
# DO0182U09A3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09A3_24hrs.xlsx", 
#                                header = TRUE,
#                                sheetIndex = 1)
# 
# DO0182U09B3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09B3_24hrs.xlsx",
#                                header = TRUE,
#                                sheetIndex = 1)
# 
# DO0182U09C3_24hrs <- read.xlsx("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/wideRawDF_DU0182U09C3_24hrs.xlsx",
#                                header = TRUE,
#                                sheetIndex = 1)
# 
# DO0182U09A3_24hrs_Mat <- as.matrix(DO0182U09A3_24hrs)
# DO0182U09B3_24hrs_Mat <- as.matrix(DO0182U09B3_24hrs)
# DO0182U09C3_24hrs_Mat <- as.matrix(DO0182U09C3_24hrs)
# 
# #convert the matrix to a correlation matrix
# DO0182U09A3_24hrCor <- cor(DO0182U09A3_24hrs_Mat)
# DO0182U09B3_24hrCor <- cor(DO0182U09B3_24hrs_Mat)
# DO0182U09C3_24hrCor <- cor(DO0182U09C3_24hrs_Mat)
# 
# corrplot(DO0182U09A3_24hrCor, 
#          method = "number", 
#          type = "lower",
#          mfrow=c(1,1),
#          oma=c(0,0,0,0),
#          mar=c(0,0,1.5,0), 
#          title = "Correlation plot between days for DO0192U09A3")
# 
# #par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(1,3,2,1))
# 
# corrplot(DO0182U09B3_24hrCor, 
#          method = "number", 
#          type = "lower",
#          mfrow=c(1,1),
#          oma=c(0,0,0,0),
#          mar=c(0,0,1.5,0), 
#          title = "Correlation plot between days for DO0192U09B3")
# 
# corrplot(DO0182U09C3_24hrCor, 
#          method = "number", 
#          type = "lower",
#          mfrow=c(1,1),
#          oma=c(0,0,0,0),
#          mar=c(0,0,1.5,0), 
#          title = "Correlation plot between days for DO0192U09C3")

#wavelet cross-correlation heatmap
RTWP.LX0088.Neighbours <- Raw_Normalised_RTWP_Finally_Got_There_020817[,c("LX0088U09A3",
                                                                          "LX0080U09C3",
                                                                          "LX0022U09A3",
                                                                          "LX0017U09B3",
                                                                          "KD0177U09C3",
                                                                          "KD0046U09C3",
                                                                          "OF0028U09B3",
                                                                          "OF0061U09B3",
                                                                          "KD0027U09A3",
                                                                          "LX0080U09A3")]


LX0088U09A3 <- raw.data.NoNA.ts[,"LX0088U09A3"]
LX0080U09C3 <- raw.data.NoNA.ts[,"LX0080U09C3"]
LX0022U09A3 <- raw.data.NoNA.ts[,"LX0022U09A3"]
LX0017U09B3 <- raw.data.NoNA.ts[,"LX0017U09B3"]
KD0177U09C3 <- raw.data.NoNA.ts[,"KD0177U09C3"]
KD0046U09C3 <- raw.data.NoNA.ts[,"KD0046U09C3"]
OF0028U09B3 <- raw.data.NoNA.ts[,"OF0028U09B3"]
OF0061U09B3 <- raw.data.NoNA.ts[,"OF0061U09B3"]
KD0027U09A3 <- raw.data.NoNA.ts[,"KD0027U09A3"]
LX0080U09A3 <- raw.data.NoNA.ts[,"LX0080U09A3"]

tenCells <- cbind(LX0088U09A3,
                  LX0080U09C3,
                  LX0022U09A3,
                  LX0017U09B3,
                  KD0177U09C3,
                  KD0046U09C3,
                  OF0028U09B3,
                  OF0061U09B3,
                  KD0027U09A3,
                  LX0080U09A3)

# RTWP.LX0088.Neighbours.TS <- ts(RTWP.LX0088.Neighbours, frequency = 96)
RTWP.LX0088.Neighbours.TS <- ts(tenCells, frequency = 96)

RTWP.LX0088.Neighbours.TS.MODWT <- apply(RTWP.LX0088.Neighbours.TS,2,waveslim::modwt,wf="la8",n.levels=7,boundary="periodic")
# 
# LX0088U09A3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$LX0088U09A3,wf="la8")
# LX0080U09C3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$LX0080U09C3,wf="la8")
# LX0022U09A3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$LX0022U09A3,wf="la8")
# LX0017U09B3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$LX0017U09B3,wf="la8")
# KD0177U09C3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$KD0177U09C3,wf="la8")
# KD0046U09C3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$KD0046U09C3,wf="la8")
# OF0028U09B3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$OF0028U09B3,wf="la8")
# OF0061U09B3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$OF0061U09B3,wf="la8")
# KD0027U09A3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$KD0027U09A3,wf="la8")
# LX0080U09A3.modwt.bw <- brick.wall(RTWP.LX0088.Neighbours.TS.MODWT$LX0080U09A3,wf="la8")

rollingWindow <- 20

modwt2 <- function(...) unlist(head(brick.wall(modwt(...), wf = "la8"), rollingWindow))

rollr <- lapply(RTWP.LX0088.Neighbours.TS, function(x) rollapplyr(x, rollingWindow, FUN = modwt2, wf = "la8", 
                                          n.levels = 4, boundary = "periodic"))

L <- lapply(rollr, function(x) lapply(1:nrow(x), function(i) matrix(x[i,], , 4)))

res <- lapply(L, function(y) lapply(y, function(x) as.list(as.data.frame(x))))

create_4mat <- function(w) {
  # create four 3*3 correlation matrices (one for each level) for window w
  M <- replicate(4, matrix(0, nrow = 10, ncol = 10), simplify = FALSE)
  for (k in 1:4) {
    for (i in 1:10) {
      for (j in (i:10)[-1]) {
        M[[k]][i, j] = wave.correlation(res[[i]][[w]], res[[j]][[w]], N=rollingWindow)[k, 1]
      }
    }
    M[[k]] <- M[[k]] + t(M[[k]]) + diag(1, 10, 10)
  }
  M
}

output <- lapply(1:((length(RTWP.LX0088.Neighbours.TS[,1])-rollingWindow)+1), create_4mat)

eigenvalues1 <- lapply(output, function(x) eigen(x[[1]], symmetric = TRUE, 
                                                 only.values = TRUE)$values)
eigenvalues1.df <- data.frame(matrix(unlist(eigenvalues1), nrow=1325, byrow=T))

eigenvalues2 <- lapply(output, function(x) eigen(x[[2]], symmetric = TRUE, 
                                                 only.values = TRUE)$values)
eigenvalues2.df <- data.frame(matrix(unlist(eigenvalues2), nrow=1325, byrow=T))

eigenvalues3 <- lapply(output, function(x) eigen(x[[3]], symmetric = TRUE, 
                                                 only.values = TRUE)$values)
eigenvalues3.df <- data.frame(matrix(unlist(eigenvalues3), nrow=1325, byrow=T))

eigenvalues4 <- lapply(output, function(x) eigen(x[[4]], symmetric = TRUE, 
                                                 only.values = TRUE)$values)
eigenvalues4.df <- data.frame(matrix(unlist(eigenvalues4), nrow=1325, byrow=T))

par(mfrow=c(4,1))
plot(eigenvalues1.df[,1],type="l")
plot(eigenvalues2.df[,1],type="l")
plot(eigenvalues3.df[,1],type="l")
plot(eigenvalues4.df[,1],type="l")

modwt2 <- function(...) unlist(head(brick.wall(modwt(...),"la8"), 5))
roller.LX0088U09A3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,1],rollingWindow,FUN=modwt2,wf="la8",n.levels=5,boundary="periodic")
roller.LX0080U09C3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,2],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.LX0022U09A3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,3],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.LX0017U09B3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,4],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.KD0177U09C3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,5],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.KD0046U09C3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,6],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.OF0028U09B3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,7],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.OF0061U09B3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,8],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.KD0027U09A3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,9],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")
roller.LX0080U09A3 <- rollapplyr(RTWP.LX0088.Neighbours.TS[,10],rollingWindow,FUN=modwt2,wf="la8",n.levels=7,boundary="periodic")



output.LX0088U09A3 <- apply(roller.LX0088U09A3, 1, function(x) list(matrix(x, 50)))
output.LX0080U09C3 <- apply(roller.LX0080U09C3, 1, function(x) list(matrix(x, 50)))
output.LX0022U09A3 <- apply(roller.LX0022U09A3, 1, function(x) list(matrix(x, 50)))
output.LX0017U09B3 <- apply(roller.LX0017U09B3, 1, function(x) list(matrix(x, 50)))
output.KD0177U09C3 <- apply(roller.KD0177U09C3, 1, function(x) list(matrix(x, 50)))
output.KD0046U09C3 <- apply(roller.KD0046U09C3, 1, function(x) list(matrix(x, 50)))
output.OF0028U09B3 <- apply(roller.OF0028U09B3, 1, function(x) list(matrix(x, 50)))
output.OF0061U09B3 <- apply(roller.OF0061U09B3, 1, function(x) list(matrix(x, 50)))
output.KD0027U09A3 <- apply(roller.KD0027U09A3, 1, function(x) list(matrix(x, 50)))
output.LX0080U09A3 <- apply(roller.LX0080U09A3, 1, function(x) list(matrix(x, 50)))

write.xlsx(output.LX0088U09A3,file = "output.LX0088U09A3.xlsx")



wavelet.coeff.s1 <- cbind(LX0088U09A3.modwt.bw$d1,
                          LX0080U09C3.modwt.bw$d1,
                          LX0022U09A3.modwt.bw$d1,
                          LX0017U09B3.modwt.bw$d1,
                          KD0177U09C3.modwt.bw$d1,
                          KD0046U09C3.modwt.bw$d1,
                          OF0028U09B3.modwt.bw$d1,
                          OF0061U09B3.modwt.bw$d1,
                          KD0027U09A3.modwt.bw$d1,
                          LX0080U09A3.modwt.bw$d1)

#wavelet.coeff.s1 <- na.omit(wavelet.coeff.s1)

wavelet.coeff.s2 <- cbind(LX0088U09A3.modwt.bw$d2,
                          LX0080U09C3.modwt.bw$d2,
                          LX0022U09A3.modwt.bw$d2,
                          LX0017U09B3.modwt.bw$d2,
                          KD0177U09C3.modwt.bw$d2,
                          KD0046U09C3.modwt.bw$d2,
                          OF0028U09B3.modwt.bw$d2,
                          OF0061U09B3.modwt.bw$d2,
                          KD0027U09A3.modwt.bw$d2,
                          LX0080U09A3.modwt.bw$d2)

#wavelet.coeff.s2 <- na.omit(wavelet.coeff.s2)


wavelet.coeff.s3 <- cbind(LX0088U09A3.modwt.bw$d3,
                          LX0080U09C3.modwt.bw$d3,
                          LX0022U09A3.modwt.bw$d3,
                          LX0017U09B3.modwt.bw$d3,
                          KD0177U09C3.modwt.bw$d3,
                          KD0046U09C3.modwt.bw$d3,
                          OF0028U09B3.modwt.bw$d3,
                          OF0061U09B3.modwt.bw$d3,
                          KD0027U09A3.modwt.bw$d3,
                          LX0080U09A3.modwt.bw$d3)

#wavelet.coeff.s3 <- na.omit(wavelet.coeff.s3)

wavelet.coeff.s4 <- cbind(LX0088U09A3.modwt.bw$d4,
                          LX0080U09C3.modwt.bw$d4,
                          LX0022U09A3.modwt.bw$d4,
                          LX0017U09B3.modwt.bw$d4,
                          KD0177U09C3.modwt.bw$d4,
                          KD0046U09C3.modwt.bw$d4,
                          OF0028U09B3.modwt.bw$d4,
                          OF0061U09B3.modwt.bw$d4,
                          KD0027U09A3.modwt.bw$d4,
                          LX0080U09A3.modwt.bw$d4)

#wavelet.coeff.s4 <- na.omit(wavelet.coeff.s4)

wavelet.coeff.s5 <- cbind(LX0088U09A3.modwt.bw$d5,
                          LX0080U09C3.modwt.bw$d5,
                          LX0022U09A3.modwt.bw$d5,
                          LX0017U09B3.modwt.bw$d5,
                          KD0177U09C3.modwt.bw$d5,
                          KD0046U09C3.modwt.bw$d5,
                          OF0028U09B3.modwt.bw$d5,
                          OF0061U09B3.modwt.bw$d5,
                          KD0027U09A3.modwt.bw$d5,
                          LX0080U09A3.modwt.bw$d5)

#wavelet.coeff.s5 <- na.omit(wavelet.coeff.s5)

wavelet.coeff.s6 <- cbind(LX0088U09A3.modwt.bw$d6,
                          LX0080U09C3.modwt.bw$d6,
                          LX0022U09A3.modwt.bw$d6,
                          LX0017U09B3.modwt.bw$d6,
                          KD0177U09C3.modwt.bw$d6,
                          KD0046U09C3.modwt.bw$d6,
                          OF0028U09B3.modwt.bw$d6,
                          OF0061U09B3.modwt.bw$d6,
                          KD0027U09A3.modwt.bw$d6,
                          LX0080U09A3.modwt.bw$d6)

#wavelet.coeff.s6 <- na.omit(wavelet.coeff.s6)

wavelet.coeff.s7 <- cbind(LX0088U09A3.modwt.bw$d7,
                          LX0080U09C3.modwt.bw$d7,
                          LX0022U09A3.modwt.bw$d7,
                          LX0017U09B3.modwt.bw$d7,
                          KD0177U09C3.modwt.bw$d7,
                          KD0046U09C3.modwt.bw$d7,
                          OF0028U09B3.modwt.bw$d7,
                          OF0061U09B3.modwt.bw$d7,
                          KD0027U09A3.modwt.bw$d7,
                          LX0080U09A3.modwt.bw$d7)

#wavelet.coeff.s7 <- na.omit(wavelet.coeff.s7)


write.xlsx(wavelet.coeff.s1, "wavelet.coeff.s1.xlsx")
write.xlsx(wavelet.coeff.s2, "wavelet.coeff.s2.xlsx")
write.xlsx(wavelet.coeff.s3, "wavelet.coeff.s3.xlsx")
write.xlsx(wavelet.coeff.s4, "wavelet.coeff.s4.xlsx")
write.xlsx(wavelet.coeff.s5, "wavelet.coeff.s5.xlsx")
write.xlsx(wavelet.coeff.s6, "wavelet.coeff.s6.xlsx")
write.xlsx(wavelet.coeff.s7, "wavelet.coeff.s7.xlsx")

wavelet_coeff_s1.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s1.result.xlsx", sheet = "Result")
wavelet_coeff_s2.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s2.result.xlsx", sheet = "Result")
wavelet_coeff_s3.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s3.result.xlsx", sheet = "Result")
wavelet_coeff_s4.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s4.result.xlsx", sheet = "Result")
wavelet_coeff_s5.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s5.result.xlsx", sheet = "Result")
wavelet_coeff_s6.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s6.result.xlsx", sheet = "Result")
wavelet_coeff_s7.result <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/wavelet.coeff.s7.result.xlsx", sheet = "Result")

eigen.coeff.s1 <- apply(wavelet_coeff_s1.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s1.result))))$values[1:10])
eigen.coeff.s1 <- t(eigen.coeff.s1)

eigen.coeff.s2 <- apply(wavelet_coeff_s2.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s2.result))))$values[1:10])

eigen.coeff.s2 <- t(eigen.coeff.s2)

eigen.coeff.s3 <- apply(wavelet_coeff_s3.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s3.result))))$values[1:10])

eigen.coeff.s3 <- t(eigen.coeff.s3)

eigen.coeff.s4 <- apply(wavelet_coeff_s4.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s4.result))))$values[1:10])

eigen.coeff.s4 <- t(eigen.coeff.s4)

eigen.coeff.s5 <- apply(wavelet_coeff_s5.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s5.result))))$values[1:10])

eigen.coeff.s5 <- t(eigen.coeff.s5)

eigen.coeff.s6 <- apply(wavelet_coeff_s6.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s6.result))))$values[1:10])

eigen.coeff.s6 <- t(eigen.coeff.s6)

eigen.coeff.s7 <- apply(wavelet_coeff_s7.result, 
                        1, 
                        function(x) 
                          eigen(matrix(x,
                                       nrow = sqrt(ncol(wavelet_coeff_s7.result))))$values[1:10])

eigen.coeff.s7 <- t(eigen.coeff.s7)

eigenvals.comb <- t(cbind(eigen.coeff.s1[,1],
                        eigen.coeff.s2[,1],
                        eigen.coeff.s3[,1],
                        eigen.coeff.s4[,1],
                        eigen.coeff.s5[,1],
                        eigen.coeff.s6[,1],
                        eigen.coeff.s7[,1]))

write.xlsx(eigenvals.comb, "eigenvals.comb.xlsx")

eigenvals_percentages <- read_excel("C:/Users/John/Google Drive/Masters/Practicum/RTWP Datasets/Complete Datasets/Absolute Last Version/RTWP/Final Data/eigenvals.percentages.xlsx", 
                                    col_names = FALSE)

eigenvals_percentages <- as.matrix(eigenvals_percentages)

clnames <- rep("",ncol(eigenvals_percentages))
sq <- seq(96,1322,96)
clnames[sq] <- sq
colnames(eigenvals_percentages) <- clnames   

library(gplots)
source("customHeatmap.r")
par(mfrow=c(1,1),oma=c(3,3,0,6))
myheatmap.2(eigenvals_percentages,
            trace = "none",
            dendrogram = "none",
            Rowv = NULL,
            Colv = NULL,
            density.info = "none",
            margin = c(5,7),
            main = "",
            xlab = "",
            ylab = "",
            lmat = rbind(c(2,3,6),c(4,1,5)),
            lwid = c(0.8, 4, 0.5),
            lhei = c(0.5, 4),
            key = TRUE,
            key.title="",
            key.xlab = "Eigenvalue\n Percentage",
            col = rev(rainbow(20*10, start = 0/6, end = 4/6)),
            cexCol=1.2,
            cexRow = 2,
            cex.axis=2)

title(main=expression(paste("Heatmap of Largest Eigenvalues ",
                            lambda[1],  " Across 7 Wavelet Scales")),
      cex.main=2)
mtext("Time Index (15 min intervals)",side = 1,line = 6, cex = 2)
mtext("Wavelet Scales",side = 4,line = 0, cex = 2)


par(mfrow=c(4,2),oma=c(0,2,0,2))
plot(eigen.coeff.s1[,1],type="l", col="blue",main="Wavelet Scale 1 (30-60 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s2[,1],type="l", col="blue",main="Wavelet Scale 2 (60-120 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s3[,1],type="l", col="blue",main="Wavelet Scale 3 (120-240 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s4[,1],type="l", col="blue",main="Wavelet Scale 4 (240-480 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s5[,1],type="l", col="blue",main="Wavelet Scale 5 (480-960 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s6[,1],type="l", col="blue",main="Wavelet Scale 6 (960-1920 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
plot(eigen.coeff.s7[,1],type="l", col="blue",main="Wavelet Scale 7 (1920-3840 min)", cex.main=2.0,xlab="Time Index (15 mins)",ylab="",cex.axis=1.5,cex.lab=1.5)
mtext(expression(lambda[1]),side = 2,line = 2,cex = 1.5)
