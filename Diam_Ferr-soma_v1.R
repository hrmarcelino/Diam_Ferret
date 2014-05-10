###############################################################
### Carregar packages para realizar os calculos necessarios ###
###############################################################
library(psych)
library(ggplot2)

########################
### Importando dados ###
########################

data1 <- MPX.F15_c1_2
data2 <- MPX.F15_c2_2
data3 <- MPX.F15_c3_2
d1 <- data1
d2 <- data2
d3 <- data3

#Convertendo os valores observados para unidades metricas (micrometro)
conversionFactor1 = 1.6
conversionFactor2 = 1.6
conversionFactor3 = 1.6
conversionFactor_T = 1.6

#Total de particulas amostradas
sampleCount1 = 500
sampleCount2 = 500
sampleCount3 = 500
sampleCount_T = sum(sampleCount1,sampleCount2, sampleCount3)

#Os valores de conversao utilizados para o microscopio da Farmacognosia (Leica) sao:
# 40x = 1.6
# 10x = 3.6

####################
### Fitting data ###
####################
#Criar coluna com valores convertidos de unidade adimensional para micrometro
d1$Micrometer <- d1$Size * conversionFactor1
d1["Micrometer"]

d2$Micrometer <- d2$Size * conversionFactor2
d2["Micrometer"]

d3$Micrometer <- d3$Size * conversionFactor3
d3["Micrometer"]

#Revelar dados apos a criacao da nova coluna
d1
d2
d3

###################################
### Getting all data togheter ###
###################################

#Obtencao do vetor de maior tamanho para que dTp, tabela que engloba todos os dados, consiga englobar a menor e a maior particula possivel.
dTSize <- max(length(d1$Size),length(d2$Size),length(d3$Size))

dTp <- data.frame()
for (i in 1:dTSize){
  Size = 0
  Freq = 0
  Micrometer = 0 
  if (i <= length(d1[,1])){
    Size = d1[i,1]
    Micrometer = d1[,3]
    Freq = Freq + d1[i,2]
  }
  if (i <= length(d2[,1])){
    Size = d2[i,1]
    Micrometer = d2[i,3]
    Freq = Freq + d2[i,2]
  }
  if (i <= length(d3[,1])){
    Size = d3[i,1]
    Micrometer = d3[i,3]
    Freq = Freq + d3[i,2]
  }
  dTp <- rbind(dTp, cbind(Size,Micrometer,Freq))
}

print(dTp)

###################################
### Processing data to analysis ###
###################################

#Para remover as primeiras linhas, em que a frequencia e igual a zero. Isto e feito para que elas nao atrapalhem o calculo da quantidade de classes.
i = 1
mustContinue = TRUE
while (mustContinue) {
  if (dTp$Freq[i] != 0) {
    mustContinue = FALSE
  } else {
    i = i + 1
  }
}

dTp = dTp[i:length(dTp$Size),1:3]

#Amplitude

#Valor mais alto
major_a <- max(dTp$Micrometer)
major_a

#Valor mais baixo
minor_a <- min(dTp$Micrometer)
minor_a

#Calculo da amplitude e do total do numero de classes que sera utilizado para analise dos dados obtidos.
Amplitude <- major_a - minor_a
Classes_a <- 1+3.333*log(sampleCount_T)

#Determinacao do intervalo de classes.
Intervala <- Amplitude/Classes_a
Interval_a <- round(Intervala) 
#A funcao "signif" foi utilizada para arredondar o numero para um numero significativo.

#A partir de agora a frequencia de cada classe sera calculada, e em seguida serao calculcadas as frequencias acumuladas nominal e percentagens
#Classes
nClasses_a <- floor(length(dTp$Size)/Interval_a)
myClasses_a <- data.frame(cbind(cinf = numeric(0), csup = numeric(0), cavg = numeric(0), cfreq = numeric(0), cfreqcum = numeric(0), cfreqper = numeric(0), cfreqcumper = numeric(0)))

for (i in 1:nClasses_a) {
  myClasses_a[i,1] = dTp$Micrometer[(i-1)*Interval_a+1]
  myClasses_a[i,2] = dTp$Micrometer[i*Interval_a]
  myClasses_a[i,3] = (myClasses_a[i,1] + myClasses_a[i,2])/2
  myClasses_a[i,4] = 0
  for (j in 1:Interval_a) {
    myClasses_a[i,4] = myClasses_a[i,4] + dTp$Freq[(i-1)*Interval_a+j]
  }
  if (i > 1) {
    myClasses_a[i,5] = myClasses_a[i-1,5] + myClasses_a[i,4]
  } else {
    myClasses_a[i,5] = myClasses_a[i,4]
  }
  myClasses_a[i,6] = (myClasses_a[i,4]/sampleCount_T)*100
  myClasses_a[i,7] = (myClasses_a[i,5]/sampleCount_T)*100
}

#Os valores para a ultima classe sao extrapolados para o total do sampleCount. E para tal tambem sao calculados todos os parametros ja calculados anteriormente para as demais classes.
if (nClasses_a*Interval_a < length(dTp$Size)) {
  myClasses_a[nClasses_a+1,1] = dTp$Micrometer[nClasses_a*Interval_a+1]
  myClasses_a[nClasses_a+1,2] = myClasses_a[nClasses_a+1,1] + (Interval_a-1) * conversionFactor3
  myClasses_a[nClasses_a+1,3] = (myClasses_a[nClasses_a+1,1] + myClasses_a[nClasses_a+1,2])/2
  myClasses_a[nClasses_a+1,4] = 0
  
  for (i in (nClasses_a*Interval_a+1):length(dTp$Size)) {
    print(i)
    myClasses_a[nClasses_a+1,4] = myClasses_a[nClasses_a+1,4] + dTp$Freq[i]
  }
  myClasses_a[nClasses_a+1,5] = myClasses_a[nClasses_a,5] + myClasses_a[nClasses_a+1,4]
  myClasses_a[nClasses_a+1,6] = (myClasses_a[nClasses_a+1,4]/sampleCount_T)*100
  myClasses_a[nClasses_a+1,7] = (myClasses_a[nClasses_a+1,5]/sampleCount_T)*100
}

#Na funcao a seguir serao obtidos dados de estatistica descritiva atraves do package(psych).
samples_a = numeric(0)
for (i in 1:(length(myClasses_a$cfreq))) {
  samples_a = c(samples_a, rep(myClasses_a[i,3], myClasses_a[i,4]))
}
describe(samples_a)

# Obtencao de histograma utilizando o lattice.
xAxis_a = paste(myClasses_a$cinf,"-",myClasses_a$csup)
barplot(myClasses_a$cfreqper, width=1, names.arg=(xAxis_a), xlab="Particle size", ylab="Frequency (%)", ylim=c(0,50))

#Para as proximas versoes
#1) Utilizar o package(ggplot2) para obtencao dos graficos;
#2) Dimensionar os graficos de modo que as legendas de cada eixo fiquem legiveis.

####################
### Fitting data ###
####################

#Sera realizada a regressao linear para obtencao da relacao entre tamanho e frequencia.
equation_a <- lm(myClasses_a$cfreqcumper ~ log10(myClasses_a$cavg))
equation_a
plot(equation_a)
anova(equation_a)
summary(equation_a)

#Plotagem do grafico de dispersao para obtencao da imagem da relacao entre frequencia e tamanho
logplot <- plot(myClasses_a$cavg, myClasses_a$cfreqcumper, xlab="Mean size (micrometer)", ylab="Culmulative frequency (%)", main="Particle Size Distribution", ylim= c(0,100), pch=20, xlim=c(1,500), log="x")
lines(myClasses_a$cavg,myClasses_a$cfreqcumper, type="l")
grid(logplot, lwd=2 , equilogs=FALSE)

#Adicionar linha de tendencia, obtida atraves da regressao linear, ao grafico de dispersao
abline(equation_a, col="blue")

###############################
### Finding D10,D50 and D90 ###
###############################

# Obtencao dos percentis para calculo do Span Index (Reynald, 2011)
D10_a <- 10^((10-equation_a$coefficients[[1]])/equation_a$coefficients[[2]])
D50_a <- 10^((50-equation_a$coefficients[[1]])/equation_a$coefficients[[2]])
D90_a <- 10^((90-equation_a$coefficients[[1]])/equation_a$coefficients[[2]])

##################
### Span Index ###
##################

Span_Index_a <- (D90_a-D10_a)/D50_a
Span_Index_a

###########
### PdI ###
###########

#Calculo do indece de polispersao, segundo manual NanoComposix Guidelines for DLS Measurements and Analysis.

x <- describe(samples_a)
media <- x[,3]
desvio_padrao <- x[,4]
PdI <- ((x[,4])/(x[,3]))^2
PdI
