##scegliamo il modello meglio adattato sulla provincia di Alberta
library(demography)
library(StMoMo)
library(lifecontingencies)
alberta <- read.demogdata("AlbertaMx_1x1.txt",
  "AlbertaExposures_1x1.txt", type="mortality", label="Alberta")
alberta.=extract.years(alberta, 1950:2015) 
albertast=StMoMoData(alberta., series="total")
albertaini=central2initial(albertast)
##richiamo le funzioni dei vari modelli
LC <- lc(link = "logit")
CBD <- cbd()
RH <- rh(link = "logit", cohortAgeFun = "1")
APC <- apc(link = "logit")
M7 <- m7()
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
ages.fit <- 0:100
wxt <- genWeightMat(ages = ages.fit, years = albertaini$years, clip = 3)
#modello LC
LCfit <- fit(LC, data = albertaini, ages.fit = ages.fit, wxt = wxt)
##residui
LCres=residuals(LCfit)
plot(LCres, type="colourmap") ##mappa di calore
##modello APC
APCfit <- fit(APC, data = albertaini, ages.fit = ages.fit, wxt = wxt)
APCres=residuals(APCfit) ##residui
plot=(APCres, type="colourmap") ##mappa di calore
##modello CBD
CBDfit <- fit(CBD, data = albertaini, ages.fit = ages.fit, wxt = wxt)
CBDres=residuals(CBDfit) ##residui
plot(CBDres, type="colourmap") ##mappa di calore
##modello M7
M7fit <- fit(M7, data = albertaini, ages.fit = ages.fit, wxt = wxt)
M7res=residuals(M7fit)
plot(M7res, type="colourmap") ##mappa di calore
##modello RH
RHfit <- fit(RH, data = albertaini, ages.fit = ages.fit, wxt = wxt, start.ax = LCfit$ax, start.bx = LCfit$bx, start.kt = LCfit$kt)
RHres=residuals(RHfit) ##residui
plot(RHres, type="colourmap") ##mappa di calore
##modello PLAT
PLATfit <- fit(PLAT, data = albertaini, ages.fit = ages.fit, wxt = wxt)
PLATres=residuals(PLATfit) ##residui
plot(PLATres, type="colourmap") ##mappa di calore
##adesso, tramite il criterio AIC, scegliamo il modello meglio adattato alla popolazione della provincia di Alberta.
AIC(LCfit) ##49013.71
AIC(APCfit) ##50784.65
AIC(CBDfit) ##329198.8
AIC(M7fit) ##143712.3
AIC(RHfit) ##48231.87
AIC(PLATfit) ##50510.05
##Il modello RH, con AIC più basso, si riconferma il modello meglio adattato. 
##adesso proietto per i prossimi 50 anni l'RHfit
RHfor=forecast(RHfit, h=50)
##unisco tassi di mortalità storici e proiettati per ricavare la lifetable, scelgo il 1950 come coorte di nascita.
rh_historical_rates <- extractCohort(fitted(RHfit, type = "rates"), cohort = 1950) ##tassi storici
rh_forecasted_rates <- extractCohort(RHfor$rates, cohort = 1950) ##tassi proiettati
rh_rates_1950 <- c(rh_historical_rates,rh_forecasted_rates)
##trasformo i tassi di mortalità in probabilità di morte
rh_qx_1950<-mx2qx(rh_rates_1950)
##ricavo la lifetable
rh_lifetable_1950<-probs2lifetable(probs=rh_qx_1950,type = "qx",name = paste("RH","1950","lt",sep="_"))
##verifico la veridicità del modello provando ad estrarre l'aspettativa di vita di una persona di età x=65.
exn(rh_lifetable_1950,x=65) ##7,87, poco realistico, cambio modello, scegliendo il LC, secondo con AIC più basso.
##modello LC
##proietto i dati per 50 anni
LCfor=forecast(LCfit, h=50)
lc_historical_rates <- extractCohort(fitted(LCfit, type = "rates"),  cohort = 1950)
lc_forecasted_rates <- extractCohort(LCfor$rates, cohort = 1950)
lc_rates_1950 <- c(lc_historical_rates,lc_forecasted_rates)
lc_qx_1950<-mx2qx(lc_rates_1950)
lc_lifetable_1950<-probs2lifetable(probs=lc_qx_1950,type = "qx", name = paste("LC","1950","lt",sep="_"))
##verifichiamo l'aspettativa di vita di una persona di età x=65 
exn(lc_lifetable_1950,x=65) ##21.30, risultato più realistico, accetto il modello LC
##ricavo tavola attuariale, tasso tecnico del 2,5%
lc_acttbl_1950<-new("actuarialtable",x=lc_lifetable_1950@x,lx=lc_lifetable_1950@lx, interest=0.025,name="LC ActTbl")
##calcoliamo l'APV, actuarial present value, potizziamo una coorte di studio con stessa età di partenza per il contratto pensionistico (x=25), età di pensionamento (67).
pensionamento=c(68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101)
renditapensionistica=axn(actuarialtable = lc_acttbl_1950,x=pensionamento)
renditapensionistica
 [1] 15.366041 14.909958 14.451521 13.990589 13.529557 13.063575 12.596295 12.125262 11.649508 11.166185 10.688365 10.203428  9.730719  9.268193  8.810324  8.360852
[17]  7.909910  7.456077  7.029958  6.609334  6.191765  5.791729  5.412014  5.045143  4.710232  4.376123  4.057964  3.713897  3.372700  3.024631  2.605808  2.221539
[33]  1.691487  1.000000
##per calcolare i contributi che in età lavorativa gli impiegati devono pagare al fondo annualmente, per m=42 anni fino ad m=0,  partendo dall'eta x=25 fino ad arrivare ad x=67, con n=1 pagamento annuo.
##creiamo il vettore degli anni lavorativi
work=c(42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0) 
##creiamo il vettore delle età
workage=c(25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67)
pagamenti=axn(actuarialtable=lc_acttbl_1950, x=workage, m=work) / + axn(actuarialtable=lc_acttbl_1950, x=workage, n=1) 
pagamenti
 [1]  4.858318  4.985292  5.115031  5.248275  5.384895  5.525111  5.669491  5.817324  5.969352  6.125205  6.285546  6.450615  6.620214  6.794407  6.973923  7.158753
[17]  7.349132  7.545441  7.747805  7.956737  8.172803  8.396094  8.626750  8.865638  9.111830  9.367990  9.632850  9.908307 10.193851 10.490829 10.799024 11.120624
[33] 11.454523 11.805283 12.173336 12.553545 12.952961 13.371975 13.810396 14.271810 14.759424 15.274611 15.823003
##grafici andamenti prestazioni
plot(pagamenti)
plot(renditapensionistica)


