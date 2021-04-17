##dopo aver provato a calcolare attività e passività con lifecontingencies, applico la formula del premio unitario provando a calcolare in maniera diversa attività e passività applicando il contratto target benefit.
library(demography)
library(StMoMo)
##carico i dati canadesi
canada <- read.demogdata("CanadaMx_1x1.txt",
  "CanadaExposures_1x1.txt", type="mortality", label="Canada")
##estraggo i dati della popolazione totale
canaday=extract.years(canada, 1950:2015)
canadatot=StMoMoData(canaday, series="total")
##modello LC e RH
LC=lc(link="log")
ages.fit=25:100
wxt <- genWeightMat(ages = ages.fit, years = canadatot$years, clip = 3)
LCfit=fit(LC, data= canadatot, ages.fit= ages.fit, years.fit=canadatot$years)
LCfor=forecast(LCfit, h=100)
##calcolo i tassi di mortalità per il Canada
mxMale=LCfor$rates
##qx Canada
eMale=exp(-mxMale)
qxMale=1-eMale
##px Canada
pxMale=1-qxMale
##utilizzo la formula del premio unitario. 
##Calcolo il numeratore, dai 67 ai 100 anni il fondo eroga la pensione ai membri.
#calcolo il numeratore della formula del premio iniziando con il vettore delle px 
#dai 67 ai 100 anni si riceve la prestazione, dunque il fondo eroga ai pensionati
px67=pxMale[43,43]
px68=pxMale[44,44]
px69=pxMale[45,45]
px70=pxMale[46,46]
px71=pxMale[47,47]
px72=pxMale[48,48]
px73=pxMale[49,49]
px74=pxMale[50,50]
px75=pxMale[51,51]
px76=pxMale[52,52]
px77=pxMale[53,53]
px78=pxMale[54,54]
px79=pxMale[55,55]
px80=pxMale[56,56]
px81=pxMale[57,57]
px82=pxMale[58,58]
px83=pxMale[59,59]
px84=pxMale[60,60]
px85=pxMale[61,61]
px86=pxMale[62,62]
px87=pxMale[63,63]
px88=pxMale[64,64]
px89=pxMale[65,65]
px90=pxMale[66,66]
px91=pxMale[67,67]
px92=pxMale[68,68]
px93=pxMale[69,69]
px94=pxMale[70,70]
px95=pxMale[71,71]
px96=pxMale[72,72]
px97=pxMale[73,73]
px98=pxMale[74,74]
px99=pxMale[75,75]
px100=pxMale[76,76]

i=0.025 ##tasso tecnico garantito di rendimento
##estraggo il vettore delle px e lo moltiplico per i flussi attualizzati
pxn=c(px67,px68,px69,px70,px71,px72,px73,px74,px75,px76,px77,px78,px79,px80,px81,px82,px83,px84,px85,px86,px87,px88,px89,px90,px91,px92,px93,px94,px95,px96,px97,px98,px99,px100)
dr=1+i ##elevo a -t e ottengo il tasso di sconto
t=c(42:75)
flusatt=dr^-t
A=pxn*flusatt #è un prodotto vettoriale
numer=sum(A) ##numeratore premio unitario =8,225814
##per ottenere il denominatore calcolo le px da 25 a 66
px25=pxMale[1,1]
px26=pxMale[2,2]
px27=pxMale[3,3]
px28=pxMale[4,4]
px29=pxMale[5,5]
px30=pxMale[6,6]
px31=pxMale[7,7]
px32=pxMale[8,8]
px33=pxMale[9,9]
px34=pxMale[10,10]
px35=pxMale[11,11]
px36=pxMale[12,12]
px37=pxMale[13,13]
px38=pxMale[14,14]
px39=pxMale[15,15]
px40=pxMale[16,16]
px41=pxMale[17,17]
px42=pxMale[18,18]
px43=pxMale[19,19]
px44=pxMale[20,20]
px45=pxMale[21,21]
px46=pxMale[22,22]
px47=pxMale[23,23]
px48=pxMale[24,24]
px49=pxMale[25,25]
px50=pxMale[26,26]
px51=pxMale[27,27]
px52=pxMale[28,28]
px53=pxMale[29,29]
px54=pxMale[30,30]
px55=pxMale[31,31]
px56=pxMale[32,32]
px57=pxMale[33,33]
px58=pxMale[34,34]
px59=pxMale[35,35]
px60=pxMale[36,36]
px61=pxMale[37,37]
px62=pxMale[38,38]
px63=pxMale[39,39]
px64=pxMale[40,40]
px65=pxMale[41,41]
px66=pxMale[42,42]
## ottengo il vettore e attualizzo i flussi
pxd=c(px25,px26,px27,px28,px29,px30,px31,px32,px33,px34,px35,px36,px37,px38,px39,px40,px41,px42,px43,px44,px45,px46,px47,px48,px49,px50,px51,px52,px53,px54,px55,px56,px57,px58,px59,px60,px61,px62,px63,px64,px65,px66)

time=c(0:41)
flusattden=dr^-time
B=pxd*flusattden
denominatore=sum(B)
denominatore =24.43295

premio=numeratore/denominatore
premio #è pari a 0,3111954
 ##ipotizzando 1000 contratti, calcolo attività ipotizzando capitalizzazioni annuali con tasso fisso + i il flusso costante annuo
 c=1000
 N1=c*pxMale[2,2]
N2=c*pxMale[3,3]
N3=c*pxMale[4,4]
N4=c*pxMale[5,5]
N5=c*pxMale[6,6]
N6=c*pxMale[7,7]
N7=c*pxMale[8,8]
N8=c*pxMale[9,9]
N9=c*pxMale[10,10]
N10=c*pxMale[11,11]
N11=c*pxMale[12,12]
N12=c*pxMale[13,13]
N13=c*pxMale[14,14]
N14=c*pxMale[15,15]
N15=c*pxMale[16,16]
N16=c*pxMale[17,17]
N17=c*pxMale[18,18]
N18=c*pxMale[19,19]
N19=c*pxMale[20,20]
N20=c*pxMale[21,21]
N21=c*pxMale[22,22]
N22=c*pxMale[23,23]
N23=c*pxMale[24,24]
N24=c*pxMale[25,25]
N25=c*pxMale[26,26]
N26=c*pxMale[27,27]
N27=c*pxMale[28,28]
N28=c*pxMale[29,29]
N29=c*pxMale[30,30]
N30=c*pxMale[31,31]
N31=c*pxMale[32,32]
N32=c*pxMale[33,33]
N33=c*pxMale[34,34]
N34=c*pxMale[35,35]
N35=c*pxMale[36,36]
N36=c*pxMale[37,37]
N37=c*pxMale[38,38]
N38=c*pxMale[39,39]
N39=c*pxMale[40,40]
N40=c*pxMale[41,41]
N41=c*pxMale[42,42]
##ottengo vettore
N=c(N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,N38,N39,N40,N41)
#Ora calcolo l'evoluzione annua del fondo, cheviene investito ad un tasso maggiore da quello erogato ai futuri pensionati. Ipotizziamo un tasso annuo del 5,76, come calcolato nella tesina investendo in azioni e obbligazioni canadesi.
i=0.0576
a=1+i
F0=c*premio
F1=(F0*a^1)+(N1*premio)
F2=(F1*a^1)+(N2*premio)
F3=(F2*a^1)+(N3*premio)
F4=(F3*a^1)+(N4*premio)
F5=(F4*a^1)+(N5*premio)
F6=(F5*a^1)+(N6*premio)
F7=(F6*a^1)+(N7*premio)
F8=(F7*a^1)+(N8*premio)
F9=(F8*a^1)+(N9*premio)
F10=(F9*a^1)+(N10*premio)
F11=(F10*a^1)+(N11*premio)
F12=(F11*a^1)+(N12*premio)
F13=(F12*a^1)+(N13*premio)
F14=(F13*a^1)+(N14*premio)
F15=(F14*a^1)+(N15*premio)
F16=(F15*a^1)+(N16*premio)
F17=(F16*a^1)+(N17*premio)
F18=(F17*a^1)+(N18*premio)
F19=(F18*a^1)+(N19*premio)
F20=(F19*a^1)+(N20*premio)
F21=(F20*a^1)+(N21*premio)
F22=(F21*a^1)+(N22*premio)
F23=(F22*a^1)+(N23*premio)
F24=(F23*a^1)+(N24*premio)
F25=(F24*a^1)+(N25*premio)
F26=(F25*a^1)+(N26*premio)
F27=(F26*a^1)+(N27*premio)
F28=(F27*a^1)+(N28*premio)
F29=(F28*a^1)+(N29*premio)
F30=(F29*a^1)+(N30*premio)
F31=(F30*a^1)+(N31*premio)
F32=(F31*a^1)+(N32*premio)
F33=(F32*a^1)+(N33*premio)
F34=(F33*a^1)+(N34*premio)
F35=(F34*a^1)+(N35*premio)
F36=(F35*a^1)+(N36*premio)
F37=(F36*a^1)+(N37*premio)
F38=(F37*a^1)+(N38*premio)
F39=(F38*a^1)+(N39*premio)
F40=(F39*a^1)+(N40*premio)
F41=(F40*a^1)+(N41*premio)
fondo=c(F0,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18,F19,F20,F21,F22,F23,F24,F25,F26,F27,F28,F29,F30,F31,F32,F33,F34,F35,F36,F37,F38,F39,F40,F41)
plot(fondo)
#le passività del fondo sono con benefit costante ai sopravvissuti fino all'età omega (100 anni).
N42=c*pxMale[43,43]
N43=c*pxMale[44,44]
N44=c*pxMale[45,45]
N45=c*pxMale[46,46]
N46=c*pxMale[47,47]
N47=c*pxMale[48,48]
N48=c*pxMale[49,49]
N49=c*pxMale[50,50]
N50=c*pxMale[51,51]
N51=c*pxMale[52,52]
N52=c*pxMale[53,53]
N53=c*pxMale[54,54]
N54=c*pxMale[55,55]
N55=c*pxMale[56,56]
N56=c*pxMale[57,57]
N57=c*pxMale[58,58]
N58=c*pxMale[59,59]
N59=c*pxMale[60,60]
N60=c*pxMale[61,61]
N61=c*pxMale[62,62]
N62=c*pxMale[63,63]
N63=c*pxMale[64,64]
N64=c*pxMale[65,65]
N65=c*pxMale[66,66]
N66=c*pxMale[67,67]
N67=c*pxMale[68,68]
N68=c*pxMale[69,69]
N69=c*pxMale[70,70]
N70=c*pxMale[71,71]
N71=c*pxMale[72,72]
N72=c*pxMale[73,73]
N73=c*pxMale[74,74]
N74=c*pxMale[75,75]
N75=c*pxMale[76,76]
erogazioni=c(N42,N43,N44,N45,N46,N47,N48,N49,N50,N51,N52,N53,N54,N55,N56,N57,N58,N59,N60,N61,N62,N63,N64,N65,N66,N67,N68,N69,N70,N71,N72,N73,N74,N75)
plot(erogazioni)
##utilizziamo uguali modalità di calcolo per Alberta e Ontario. Sostituisco il dataset nelle prime funzioni.
alberta <- read.demogdata("AlbertaMx_1x1.txt", "AlbertaExposures_1x1.txt", type="mortality", label="Alberta")
albertay=extract.years(alberta, 1950:2015)
albertatot=StMoMoData(albertay, series="total")
ages.fit=25:100
wxt <- genWeightMat(ages = ages.fit, years = albertatot$years, clip = 3) 
 LCfit=fit(LC, data= albertatot, ages.fit= ages.fit, years.fit=albertatot$years)
 LCfor=forecast(LCfit, h=100)
 mxMale=LCfor$rates
 ##si prosegue ricopiando i comandi precedenti, alla fine si ottiene
fondo ##alberta
 [1]   336.6689   692.5179  1068.8785  1466.9250  1887.8671  2333.0585  2803.8998
 [8]  3301.8515  3828.4846  4385.4397  4974.4507  5597.4072  6256.2339  6952.9961
[15]  7689.8630  8469.1639  9293.3217 10164.9483 11086.7732 12061.6434 13092.6496
[22] 14182.9946 15336.0937 16555.6187 17845.3508 19209.3231 20651.7869 22177.3090
[29] 23790.6643 25496.8680 27301.2695 29209.6154 31227.7295 33362.0173 35619.2121
[36] 38006.2554 40530.6515 43200.4633 46023.9673 49009.9350 52167.8084 55507.3211
> erogazioni ##alberta
 [1] 993.5466 992.9004 992.3646 991.6202 990.5525 990.2914 989.4271 988.8261 987.9388
[10] 987.3192 985.4079 984.8091 981.6203 978.2842 975.1372 971.6729 969.4589 967.2234
[19] 958.4844 953.8773 948.6767 939.1815 926.3022 914.9863 896.6739 884.7057 865.3407
[28] 857.6520 832.2278 806.3202 798.8816 699.9552 715.8534 655.1101
##vale lo stesso per l'Ontario
ontario <- read.demogdata("OntarioMx_1x1.txt", "OntarioExposures_1x1.txt", type="mortality", label="Ontario")
ontarioy=extract.years(ontario, 1950:2015)
ontariotot=StMoMoData(ontarioy, series="total")
ages.fit=25:100
wxt <- genWeightMat(ages = ages.fit, years = ontariotot$years, clip = 3) 
 LCfit=fit(LC, data= ontariotot, ages.fit= ages.fit, years.fit=ontariotot$years)
 LCfor=forecast(LCfit, h=100)
