##analisi di sopravvivenza popolazione canadese, valutazione modello con miglior bontà di adattamento e proiezione di mortalità
library(demography)
library(StMoMo)
##range di età scelto per la popolazione canadese: 55:100
##periodo: 1950:2016
##carico dataset popolazione canadese, ricavato dall'human mortality database
canada <- read.demogdata("CanadaMx_1x1.txt","CanadaExposures_1x1.txt", type="mortality", label="Canada")
##estraggo range età
canada.=extract.ages(canada, 55:100)
##estraggo periodo
.canada.=extract.years(canada., 1950:2016)
##converto l'oggetto demogdata in stmomodata
canadast=StMoMoData(.canada., series="total")
##inizializzo i dati
canadaini=central2initial(canadast)
##carico le funzioni dei vari modelli che verranno applicati
##LC
LC <- lc(link = "logit")
##CBD
CBD <- cbd()
## APC, RH, M7
RH <- rh(link = "logit", cohortAgeFun = "1")
APC <- apc(link = "logit")
M7 <- m7()
##PLAT model
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1",
               constFun = constPlat)
##costruiamo la matrice pesata 0-1
ages.fit=55:100
wxt <- genWeightMat(ages = ages.fit, years = canadaini$years, clip = 3)
##modello Lee-Carter
LCfit <- fit(LC, data = canadaini, ages.fit = ages.fit, wxt = wxt)
##per ottenere il grafico dei parametri
plot(LCfit)
##analisi dei residui
LCres=residuals(LCfit)
##grafici per i residui
plot(LCres) ##scatterplot
plot(LCres, type="colourmap") ##mappa di calore
##modello APC
APCfit <- fit(APC, data = canadaini, ages.fit = ages.fit, wxt = wxt)
##grafico parametri APCfit
plot(APCfit)
##analisi residui APC
APCres=residuals(APCfit)
##grafici residui APC
plot(APCres) ##scatter plot
plot(APCres, type="colourmap") ##mappa di calore
##modello CBD
CBDfit <- fit(CBD, data = canadaini, ages.fit = ages.fit, wxt = wxt)
##analisi dei residui CBD
CBDres=residuals(CBDfit)
CBDres
##grafici residui CBD
plot(CBDres) ##scatter plot
plot(CBDres, type="colourmap") ##mappa di calore
##M7 model
M7fit <- fit(M7, data = canadaini, ages.fit = ages.fit, wxt = wxt)
##grafici dei parametri del modello M7
plot(M7fit)
##analisi dei residui M7
M7res=residuals(M7fit)
##grafici dei residui M7
plot(M7res) ##scatter plot
plot(M7res, type="colourmap") ##mappa di calore
##modello RH
RHfit <- fit(RH, data = canadaini, ages.fit = ages.fit, wxt = wxt, start.ax = LCfit$ax, start.bx = LCfit$bx, start.kt = LCfit$kt)
##grafici parametri modello RH
plot(RHfit)
##analisi dei residui RH
RHres=residuals(RHfit)
##grafici dei residui RH
plot(RHres) ##scatter plot
plot(RHres, type="colourmap") ##mappa di calore
##modello PLAT
PLATfit <- fit(PLAT, data = canadaini, ages.fit = ages.fit, wxt = wxt)
##grafici dei parametri modello PLAT
plot(PLATfit)
##analisi dei residui modello PLAT
PLATres=residuals(PLATfit)
##grafici residui modello PLAT
plot(PLATres) ##scatter plot
plot(PLATres, type="colourmap") ##mappa di calore
##oltre all'analisi dei residui, usiamo altri due criteri per valutare la bontà di adattamento dei modelli, AIC e BIC.
##AIC dei vari fit
AIC(LCfit)
[1] 41271.37
AIC(APCfit)
[1] 40481.57
AIC(CBDfit)
[1] 58085.05
AIC(RHfit)
[1] 34056.6
AIC(M7fit)
[1] 36443.1
AIC(PLATfit)
[1] 35107.43
##BIC dei vari modelli
BIC(LCfit)
[1] 42217.99
BIC(APCfit)
[1] 41783.93
BIC(CBDfit)
[1] 58892.99
BIC(RHfit)
[1] 35636.31
BIC(M7fit)
[1] 38276.04
BIC(PLATfit)
[1] 36801.7
##in base ai valori AIC e BIC, bisogna preferire il modello con i valori più bassi. Nel caso della popolazione canadese il modello meglio adattato è quello RH.
##adesso eseguiamo una proiezione di 35 anni sui dati canadesi utilizzando il fitting del modello RH.
RHfor <- forecast(RHfit, h = 35, gc.order = c(1, 1, 0))
RHfor
Stochastic Mortality Model forecast
Call: forecast.fitStMoMo(object = RHfit, h = 35, gc.order = c(1, 1,  
Call:     0))

Binomial model with predictor: logit q[x,t] = a[x] + b1[x] k1[t] + g[t-x]

kt model: mrwd
gc model: ARIMA(1,1,0) with drift
Jump-off method: fit
Data:  Canada
Series:  total
Years in forecast: 2017 - 2051
Ages in forecast: 55 - 100 
##Grafico della proiezione
plot(RHfor, only.kt = TRUE)
##fine.
