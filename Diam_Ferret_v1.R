########################
### Loading packages ###
########################
library(psych)
library(ggplot2)

###################
### Import data ###
###################

data1 <- read.csv("~/Dropbox/UFRN/LaSiD/SLP/Alunos Iniciacao Cientifica/Bartolomeu/IC/TCC/20132/XBM_corn_5-ASA_1-1.csv", sep=";")
data2 <- read.csv("~/Dropbox/UFRN/LaSiD/SLP/Alunos Iniciacao Cientifica/Bartolomeu/IC/TCC/20132/XBM_corn_5-ASA_1-1.csv", sep=";")
data3 <- read.csv("~/Dropbox/UFRN/LaSiD/SLP/Alunos Iniciacao Cientifica/Bartolomeu/IC/TCC/20132/XBM_corn_5-ASA_1-1.csv", sep=";")
d1 <- data1
d2 <- data2
d3 <- data3

conversionFactor1 = 1.6
conversionFactor2 = 1.6
conversionFactor3 = 1.6
conversionFactor_T = 1.6
sampleCount1 = 500
sampleCount2 = 500
sampleCount3 = 500
sampleCount_T = sum(sampleCount1,sampleCount2, sampleCount3)

#The standards values to conversionFactor are:
# 40x = 1.6
# 10x = 3.6

####################
### Fitting data ###
####################

d1$Micrometer <- d1$Size * conversionFactor1
d1["Micrometer"]

d2$Micrometer <- d2$Size * conversionFactor2
d2["Micrometer"]

d3$Micrometer <- d3$Size * conversionFactor3
d3["Micrometer"]

#Looking at data after processing
d1
d2
d3

###################################
### Processing data to analysis ###
###################################
# For now on, the analysis will be done separetelly

##################
#### To data1 ####
##################

# To remove the first line where the frequence is equal to zero.
i = 1
mustContinue = TRUE
while (mustContinue) {
  if (d1$Freq[i] != 0) {
    mustContinue = FALSE
  } else {
    i = i + 1
  }
}

#Now, d1p is the data after the processing
d1p = d1[i:length(d1$Size),1:3]

#Amplitude

#High value
major_1 <- max(d1p$Micrometer)
major_1

#Low valor
minor_1 <- min(d1p$Micrometer)
minor_1

Amplitude <- major_1 - minor_1
Classes_1 <- 1+3.333*log(sampleCount1)

#Intervals
Interval1 <- Amplitude/Classes_1
Interval_1 <- round(Interval1) 
#the function signif was used to floor the number to a specific number of significative digits

#Classes
nClasses_1 <- floor(length(d1p$Size)/Interval_1)
myClasses_1 <- data.frame(cbind(cinf = numeric(0), csup = numeric(0), cavg = numeric(0), cfreq = numeric(0), cfreqcum = numeric(0), cfreqper = numeric(0), cfreqcumper = numeric(0)))

for (i in 1:nClasses_1) {
  myClasses_1[i,1] = d1p$Micrometer[(i-1)*Interval_1+1]
  myClasses_1[i,2] = d1p$Micrometer[i*Interval_1]
  myClasses_1[i,3] = (myClasses_1[i,1] + myClasses_1[i,2])/2
  myClasses_1[i,4] = 0
  for (j in 1:Interval_1) {
    myClasses_1[i,4] = myClasses_1[i,4] + d1p$Freq[(i-1)*Interval_1+j]
  }
  if (i > 1) {
    myClasses_1[i,5] = myClasses_1[i-1,5] + myClasses_1[i,4]
  } else {
    myClasses_1[i,5] = myClasses_1[i,4]
  }
  myClasses_1[i,6] = (myClasses_1[i,4]/sampleCount1)*100
  myClasses_1[i,7] = (myClasses_1[i,5]/sampleCount1)*100
}

#Extrapolating the last interval of the last class
if (nClasses_1*Interval_1 < length(d1p$Size)) {
  myClasses_1[nClasses_1+1,1] = d1p$Micrometer[nClasses_1*Interval_1+1]
  myClasses_1[nClasses_1+1,2] = myClasses_1[nClasses_2+1,1] + (Interval_1-1) * conversionFactor1
  myClasses_1[nClasses_1+1,3] = (myClasses_1[nClasses_2+1,1] + myClasses_1[nClasses_1+1,2])/2
  myClasses_1[nClasses_1+1,4] = 0
  
  for (i in (nClasses_1*Interval_1+1):length(d1p$Size)) {
    print(i)
    myClasses_1[nClasses_1+1,4] = myClasses_1[nClasses_1+1,4] + d1p$Freq[i]
  }
  myClasses_1[nClasses_1+1,5] = myClasses_1[nClasses_1,5] + myClasses_1[nClasses_1+1,4]
  myClasses_1[nClasses_1+1,6] = (myClasses_1[nClasses_1+1,4]/sampleCount1)*100
  myClasses_1[nClasses_1+1,7] = (myClasses_1[nClasses_1+1,5]/sampleCount1)*100
}

#Now, it is time to evaluate the general data
samples_1 = numeric(0)
for (i in 1:(length(myClasses_1$cfreq))) {
  samples_1 = c(samples_1, rep(myClasses_1[i,3], myClasses_1[i,4]))
}
describe(samples_1)

#Histogram
xAxis_1 = paste(myClasses_1$cinf,"-",myClasses_1$csup)
barplot(myClasses_1$cfreqper, width=1, names.arg=(xAxis_1), xlab="Particle size", ylab="Frequency (%)", ylim=c(0,50))
library(moments)
kurtosis(myClasses_1$cfreqper)

####################
### Fitting data ###
####################
equation_1 <- lm(myClasses_1$cfreqcumper ~ log10(myClasses_1$cavg))
equation_1
plot(equation_1)
anova(equation_1)
summary(equation_1)

#Scatter plot
logplot <- plot(myClasses_1$cavg, myClasses_1$cfreqcumper, xlab="Mean size (micrometer)", ylab="Culmulative frequency (%)", main="Particle Size Distribution", ylim= c(0,100), pch=20, xlim=c(1,500), log="x")
lines(myClasses_1$cavg,myClasses_1$cfreqcumper, type="l")
grid(logplot, lwd=2 , equilogs=FALSE)

#Add abline to plot
abline(equation_1, col="blue")

###############################
### Finding D10,D50 and D90 ###
###############################

D10_1 <- 10^((10-equation_1$coefficients[[1]])/equation_1$coefficients[[2]])
D50_1 <- 10^((50-equation_1$coefficients[[1]])/equation_1$coefficients[[2]])
D90_1 <- 10^((90-equation_1$coefficients[[1]])/equation_1$coefficients[[2]])

##################
### Span Index ###
##################

Span_Index_1 <- (D90_1-D10_1)/D50_1
Span_Index_1

#######################################################################################
#######################################################################################
#######################################################################################

##################
#### To data2 ####
##################

# To remove the first line where the frequence is equal to zero.
i = 1
mustContinue = TRUE
while (mustContinue) {
  if (d2$Freq[i] != 0) {
    mustContinue = FALSE
  } else {
    i = i + 1
  }
}

#Now, d2p is the data after the processing
d2p = d2[i:length(d2$Size),1:3]

#Amplitude

#High value
major_2 <- max(d2p$Micrometer)
major_2

#Low valor
minor_2 <- min(d2p$Micrometer)
minor_2

Amplitude <- major_2 - minor_2
Classes_2 <- 1+3.333*log(sampleCount2)

#Intervals
Interval2 <- Amplitude/Classes_2
Interval_2 <- round(Interval2) 
#the function signif was used to floor the number to a specific number of significative digits

#Classes
nClasses_2 <- floor(length(d2p$Size)/Interval_2)
myClasses_2 <- data.frame(cbind(cinf = numeric(0), csup = numeric(0), cavg = numeric(0), cfreq = numeric(0), cfreqcum = numeric(0), cfreqper = numeric(0), cfreqcumper = numeric(0)))

for (i in 1:nClasses_2) {
  myClasses_2[i,1] = d2p$Micrometer[(i-1)*Interval_2+1]
  myClasses_2[i,2] = d2p$Micrometer[i*Interval_2]
  myClasses_2[i,3] = (myClasses_2[i,1] + myClasses_2[i,2])/2
  myClasses_2[i,4] = 0
  for (j in 1:Interval_2) {
    myClasses_2[i,4] = myClasses_2[i,4] + d2p$Freq[(i-1)*Interval_2+j]
  }
  if (i > 1) {
    myClasses_2[i,5] = myClasses_2[i-1,5] + myClasses_2[i,4]
  } else {
    myClasses_2[i,5] = myClasses_2[i,4]
  }
  myClasses_2[i,6] = (myClasses_2[i,4]/sampleCount2)*100
  myClasses_2[i,7] = (myClasses_2[i,5]/sampleCount2)*100
}

#Extrapolating the last interval of the last class
if (nClasses_2*Interval_2 < length(d2p$Size)) {
  myClasses_2[nClasses_2+1,1] = d2p$Micrometer[nClasses_2*Interval_2+1]
  myClasses_2[nClasses_2+1,2] = myClasses_2[nClasses_2+1,1] + (Interval_2-1) * conversionFactor2
  myClasses_2[nClasses_2+1,3] = (myClasses_2[nClasses_2+1,1] + myClasses_2[nClasses_2+1,2])/2
  myClasses_2[nClasses_2+1,4] = 0
  
  for (i in (nClasses_2*Interval_2+1):length(d2p$Size)) {
    print(i)
    myClasses_2[nClasses_2+1,4] = myClasses_2[nClasses_2+1,4] + d2p$Freq[i]
  }
  myClasses_2[nClasses_2+1,5] = myClasses_2[nClasses_2,5] + myClasses_2[nClasses_2+1,4]
  myClasses_2[nClasses_2+1,6] = (myClasses_2[nClasses_2+1,4]/sampleCount2)*100
  myClasses_2[nClasses_2+1,7] = (myClasses_2[nClasses_2+1,5]/sampleCount2)*100
}

#Now, it is time to evaluate the general data
samples_2 = numeric(0)
for (i in 1:(length(myClasses_2$cfreq))) {
  samples_2 = c(samples_2, rep(myClasses_2[i,3], myClasses_2[i,4]))
}
describe(samples_2)

#Histogram
xAxis_2 = paste(myClasses_2$cinf,"-",myClasses_2$csup)
barplot(myClasses_2$cfreqper, width=1, names.arg=(xAxis_2), xlab="Particle size", ylab="Frequency (%)", ylim=c(0,50))
library(moments)
kurtosis(myClasses_2$cfreqper)

####################
### Fitting data ###
####################
equation_2 <- lm(myClasses_2$cfreqcumper ~ log10(myClasses_2$cavg))
equation_2
plot(equation_2)
anova(equation_2)
summary(equation_2)

#Scatter plot
logplot <- plot(myClasses_2$cavg, myClasses_2$cfreqcumper, xlab="Mean size (micrometer)", ylab="Culmulative frequency (%)", main="Particle Size Distribution", ylim= c(0,100), pch=20, xlim=c(1,500), log="x")
lines(myClasses_2$cavg,myClasses_2$cfreqcumper, type="l")
grid(logplot, lwd=2 , equilogs=FALSE)

#Add abline to plot
abline(equation_2, col="blue")

###############################
### Finding D10,D50 and D90 ###
###############################

D10_2 <- 10^((10-equation_2$coefficients[[1]])/equation_2$coefficients[[2]])
D50_2 <- 10^((50-equation_2$coefficients[[1]])/equation_2$coefficients[[2]])
D90_2 <- 10^((90-equation_2$coefficients[[1]])/equation_2$coefficients[[2]])

##################
### Span Index ###
##################

Span_Index_2 <- (D90_2-D10_2)/D50_2
Span_Index_2

#######################################################################################
#######################################################################################
#######################################################################################

##################
#### To data3 ####
##################

# To remove the first line where the frequence is equal to zero.
i = 1
mustContinue = TRUE
while (mustContinue) {
  if (d3$Freq[i] != 0) {
    mustContinue = FALSE
  } else {
    i = i + 1
  }
}

#Now, d3p is the data after the processing
d3p = d3[i:length(d3$Size),1:3]

#Amplitude

#High value
major_3 <- max(d3p$Micrometer)
major_3

#Low valor
minor_3 <- min(d3p$Micrometer)
minor_3

Amplitude <- major_3 - minor_3
Classes_3 <- 1+3.333*log(sampleCount3)

#Intervals
Interval3 <- Amplitude/Classes_3
Interval_3 <- round(Interval3) 
#the function signif was used to floor the number to a specific number of significative digits

#Classes
nClasses_3 <- floor(length(d3p$Size)/Interval_3)
myClasses_3 <- data.frame(cbind(cinf = numeric(0), csup = numeric(0), cavg = numeric(0), cfreq = numeric(0), cfreqcum = numeric(0), cfreqper = numeric(0), cfreqcumper = numeric(0)))

for (i in 1:nClasses_3) {
  myClasses_3[i,1] = d3p$Micrometer[(i-1)*Interval_3+1]
  myClasses_3[i,2] = d3p$Micrometer[i*Interval_3]
  myClasses_3[i,3] = (myClasses_3[i,1] + myClasses_3[i,2])/2
  myClasses_3[i,4] = 0
  for (j in 1:Interval_3) {
    myClasses_3[i,4] = myClasses_3[i,4] + d3p$Freq[(i-1)*Interval_3+j]
  }
  if (i > 1) {
    myClasses_3[i,5] = myClasses_3[i-1,5] + myClasses_3[i,4]
  } else {
    myClasses_3[i,5] = myClasses_3[i,4]
  }
  myClasses_3[i,6] = (myClasses_3[i,4]/sampleCount3)*100
  myClasses_3[i,7] = (myClasses_3[i,5]/sampleCount3)*100
}

#Extrapolating the last interval of the last class
if (nClasses_3*Interval_3 < length(d3p$Size)) {
  myClasses_3[nClasses_3+1,1] = d3p$Micrometer[nClasses_3*Interval_3+1]
  myClasses_3[nClasses_3+1,2] = myClasses_3[nClasses_3+1,1] + (Interval_3-1) * conversionFactor3
  myClasses_3[nClasses_3+1,3] = (myClasses_3[nClasses_3+1,1] + myClasses_3[nClasses_3+1,2])/2
  myClasses_3[nClasses_3+1,4] = 0
  
  for (i in (nClasses_3*Interval_3+1):length(d3p$Size)) {
    print(i)
    myClasses_3[nClasses_3+1,4] = myClasses_3[nClasses_3+1,4] + d3p$Freq[i]
  }
  myClasses_3[nClasses_3+1,5] = myClasses_3[nClasses_3,5] + myClasses_3[nClasses_3+1,4]
  myClasses_3[nClasses_3+1,6] = (myClasses_3[nClasses_3+1,4]/sampleCount3)*100
  myClasses_3[nClasses_3+1,7] = (myClasses_3[nClasses_3+1,5]/sampleCount3)*100
}

#Now, it is time to evaluate the general data
samples_3 = numeric(0)
for (i in 1:(length(myClasses_3$cfreq))) {
  samples_3 = c(samples_3, rep(myClasses_3[i,3], myClasses_3[i,4]))
}
describe(samples_3)

#Histogram
xAxis_3 = paste(myClasses_3$cinf,"-",myClasses_3$csup)
barplot(myClasses_3$cfreqper, width=1, names.arg=(xAxis_3), xlab="Particle size", ylab="Frequency (%)", ylim=c(0,50))
library(moments)
kurtosis(myClasses_3$cfreqper)

####################
### Fitting data ###
####################
equation_3 <- lm(myClasses_3$cfreqcumper ~ log10(myClasses_3$cavg))
equation_3
plot(equation_3)
anova(equation_3)
summary(equation_3)

#Scatter plot
logplot <- plot(myClasses_3$cavg, myClasses_3$cfreqcumper, xlab="Mean size (micrometer)", ylab="Culmulative frequency (%)", main="Particle Size Distribution", ylim= c(0,100), pch=20, xlim=c(1,500), log="x")
lines(myClasses_3$cavg,myClasses_3$cfreqcumper, type="l")
grid(logplot, lwd=2 , equilogs=FALSE)

#Add abline to plot
abline(equation_3, col="blue")

###############################
### Finding D10,D50 and D90 ###
###############################

D10_3 <- 10^((10-equation_3$coefficients[[1]])/equation_3$coefficients[[2]])
D50_3 <- 10^((50-equation_3$coefficients[[1]])/equation_3$coefficients[[2]])
D90_3 <- 10^((90-equation_3$coefficients[[1]])/equation_3$coefficients[[2]])

##################
### Span Index ###
##################

Span_Index_3 <- (D90_3-D10_3)/D50_3
Span_Index_3

#You have finished, Congrats!#