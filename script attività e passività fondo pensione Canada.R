##calcolo attività e passività del fondo pensione utilizzando lifecontingencies.
##carico i pacchetti necessari
library(demography)
library(StMoMo)
library(lifecontingencies)
##carico dati canadesi
canada <- read.demogdata("CanadaMx_1x1.txt",
  "CanadaExposures_1x1.txt", type="mortality", label="Canada")
##come visto nell'analisi di sopravvivenza, il modello RH ha la migliore bontà di adattamento per la popolazione canadese, fitto i dati.
canada.=extract.years(canada,1950:2015)
canadast=StMoMoData(canada., series="total")
##inizializzo i dati
canadaini=central2initial(canadast)
##assumo età massima 100
ages.fit=0:100
wxt <- genWeightMat(ages = ages.fit, years = canadaini$years,clip = 3)
##fitting
LC <- lc(link = "logit")
RH <- rh(link = "logit", cohortAgeFun = "1")
LCfit <- fit(LC, data = canadaini, ages.fit = ages.fit, wxt = wxt)
RHfit <- fit(RH, data = canadaini, ages.fit = ages.fit, wxt = wxt,start.ax = LCfit$ax,start.bx = LCfit$bx, start.kt = LCfit$kt)
##proiettiamo i dati per un orizzonte temporale di 50 anni
RHfor <- forecast(RHfit, h = 50)
##per ottenere la lifetable delle proiezioni di mortalità data una coorte di nascita (nel nostro caso scegliamo l'anno 1950), unisco i tassi di mortalità storici con quelli proiettati.
rh_historical_rates <- extractCohort(fitted(RHfit, type = "rates"), cohort = 1950)
rh_forecasted_rates <- extractCohort(RHfor$rates, cohort = 1950)
rh_rates_1950 <- c(rh_historical_rates,rh_forecasted_rates)
##adesso converto i tassi di mortalità in probabilità di morte
rh_qx_1950<-mx2qx(rh_rates_1950)
##ottengo la lifetable
rh_lifetable_1950<-probs2lifetable(probs=rh_qx_1950,type = "qx",name = paste("RH","1950","lt",sep="_"))
##per verificare la correttezza delle proiezioni, estraggo l'aspettativa di vita di una persona 65enne (quanti anni gli restano da vivere).
exn(rh_lifetable_1950,x=65) ##33,49644, un dato accettabile.
##adesso creiamo la tavola attuariale necessario per il calcolo del bilancio del fondo pensione, ipotizzando un tasso tecnico del 2,5 %
rh_acttbl_1950<-new("actuarialtable",x=rh_lifetable_1950@x,lx=rh_lifetable_1950@lx, interest=0.025,name="RH ActTbl")
##prima di tutto calcoliamo l'APV, l'acturial preent value, cioè il valore attuariale dei flussi che il fondo si aspetta di dover erogare durante il pensionamento. Ipotizziamo una coorte di studio con stessa età di partenza per il contratto pensionistico (x=25), età di pensionamento (67).
##creiamo un vettore con le varie età fino a quando il fondo si aspetta di dover erogare le pensioni. 
pensionamento=c(68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101)
renditapensionistica=axn(actuarialtable = rh_acttbl_1950,x=pensionamento)
renditapensionistica
 [1] 22.250338 21.943537 21.613758 21.261469 20.887196 20.487008 20.063077 19.617011
 [9] 19.148208 18.657941 18.146975 17.614992 17.063411 16.493559 15.904374 15.295801
[17] 14.668017 14.021463 13.356317 12.672297 11.969613 11.248041 10.507368  9.747339
[25]  8.967753  8.168121  7.348126  6.507361  5.645398  4.761750  3.855915  2.927382
[33]  1.975598  1.000000
##per calcolare i contributi che in età lavorativa gli impiegati devono pagare al fondo annualmente, per m=42 anni fino ad m=0,  partendo dall'eta x=25 fino ad arrivare ad x=67, con n=1 pagamento annuo.
##creiamo il vettore degli anni lavorativi
work=c(42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0) 
##creiamo il vettore delle età
workage=c(25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67)
pagamenti=axn(actuarialtable=rh_acttbl_1950, x=workage, m=work) / + axn(actuarialtable=rh_acttbl_1950, x=workage, n=1) 
pagamenti
 [1]  6.973011  7.154266  7.340147  7.530723  7.726137  7.926757  8.132702  8.344025
 [9]  8.561191  8.784191  9.013296  9.249055  9.491320  9.740531  9.997172 10.261277
[17] 10.533191 10.813405 11.102248 11.400420 11.708343 12.026401 12.355421 12.695777
[25] 13.047955 13.413160 13.791227 14.183155 14.590298 15.012842 15.452518 15.910179
[33] 16.386041 16.883884 17.402506 17.943100 18.509342 19.104660 19.724871 20.376192
[41] 21.064796 21.790161 22.532261
##grafici andamenti prestazioni
plot(pagamenti)
plot(renditapensionistica)
