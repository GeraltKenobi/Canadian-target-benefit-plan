library(demography)
library(StMoMo)
library(lifecontingencies)
##Provincia dell'Ontario
Ontario <- read.demogdata("OntarioMx_1x1.txt",
  "OntarioExposures_1x1.txt", type="mortality", label="Ontario")
ontario.=extract.years(Ontario, 1950:2015)
ontariost=StMoMoData(ontario., series="total")
ontarioini=central2initial(ontariost)
##faccio il fitting dei dati con i vari modelli 
##richiamo le funzioni
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
##età massima 100
ages.fit=0:100
wxt <- genWeightMat(ages = ages.fit, years = ontarioini$years, clip = 3)
LCfit <- fit(LC, data = ontarioini, ages.fit = ages.fit, wxt = wxt) 
APCfit <- fit(APC, data = ontarioini, ages.fit = ages.fit, wxt = wxt)
CBDfit <- fit(CBD, data = ontarioini, ages.fit = ages.fit, wxt = wxt) 
 M7fit <- fit(M7, data = ontarioini, ages.fit = ages.fit, wxt = wxt) 
 PLATfit <- fit(PLAT, data = ontarioini, ages.fit = ages.fit, wxt = wxt)
 RHfit <- fit(RH, data = ontarioini, ages.fit = ages.fit, wxt = wxt,  start.ax = LCfit$ax, start.bx = LCfit$bx, start.kt = LCfit$kt)
 ##valuto il modello meglio adattato in base agli AIC
 AIC(LCfit) ##66780.38
 AIC(APCfit) ##68528.08
 AIC(CBDfit) ##1058332
 AIC(M7fit) ##416620
 AIC(RHfit) ##59925.76
 AIC(PLATfit) ##67966.7
 ##scelgo il modello RH.
 ##proietto i dati per 50 anni, ricavo i tassi di mortalità storici e proiettati e li trasformo in probabilità di morte, coorte di nascita=1950.
 RHfor=forecast(RHfit, h=50)
rh_historical_rates <- extractCohort(fitted(RHfit, type = "rates"), cohort = 1950)
rh_forecasted_rates <- extractCohort(RHfor$rates, cohort = 1950)
rh_rates_1950 <- c(rh_historical_rates,rh_forecasted_rates)
rh_qx_1950<-mx2qx(rh_rates_1950)
##ricavo la lifetable
rh_lifetable_1950<-probs2lifetable(probs=rh_qx_1950,type = "qx",name = paste("RH","1950","lt",sep="_"))
##valuto la veridicità della lifetable estraendo l'aspettativa di vita di una testa d'età x=65.
exn(rh_lifetable_1950,x=65) ##25.25, risultato realistico, accetto il modello RH.
##estraggo tavola attuariale, tasso tecnico del 2,5 %
rh_acttbl_1950<-new("actuarialtable",x=rh_lifetable_1950@x,lx=rh_lifetable_1950@lx, interest=0.025,name="RH ActTbl")
##calcoliamo l'APV, actuarial present value, potizziamo una coorte di studio con stessa età di partenza per il contratto pensionistico (x=25), età di pensionamento (67).
pensionamento=c(68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101)
renditapensionistica=axn(actuarialtable = rh_acttbl_1950,x=pensionamento)
renditapensionistica
[1] 17.604182 17.197637 16.786008 16.372524 15.956155 15.530681 15.101439 14.667728
 [9] 14.225934 13.780988 13.331783 12.873051 12.412901 11.962722 11.510797 11.054927
[17] 10.590056 10.122982  9.653236  9.177053  8.702471  8.218166  7.733880  7.242784
[25]  6.757649  6.250330  5.730142  5.182871  4.608696  4.012319  3.356371  2.647126
[33]  1.869662  1.000000
##per calcolare i contributi che in età lavorativa gli impiegati devono pagare al fondo annualmente, per m=42 anni fino ad m=0,  partendo dall'eta x=25 fino ad arrivare ad x=67, con n=1 pagamento annuo.
##creiamo il vettore degli anni lavorativi
work=c(42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0) 
##creiamo il vettore delle età
workage=c(25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67)
pagamenti=axn(actuarialtable=rh_acttbl_1950, x=workage, m=work) / + axn(actuarialtable=rh_acttbl_1950, x=workage, n=1) 
pagamenti
[1]  5.610397  5.755236  5.903917  6.056208  6.212474  6.373020  6.537739  6.706844
 [9]  6.880465  7.058843  7.242132  7.430721  7.624345  7.823675  8.028792  8.239971
[17]  8.457402  8.681410  8.912441  9.150655  9.396528  9.650141  9.912424 10.183876
[25] 10.464750 10.755839 11.057450 11.369546 11.694214 12.030911 12.381023 12.745213
[33] 13.124314 13.519461 13.931584 14.361067 14.809863 15.280930 15.771582 16.287169
[41] 16.831458 17.402582 18.006033
##grafici andamenti prestazioni
plot(pagamenti)
plot(renditapensionistica)

 