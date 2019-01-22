#DDD model Module based

#----------------------------------------------------------------------------------
#     Script name:  DDD_Modul_Eb_Oct2018_func.R
#     Author:     Thomas Skaugen
#     Revised: 15.10.2018
#              
#Purpose   : DDD (Distance Distribution Dynamics) hydrological model simulates
#               1)saturated soil water flow, states of saturated (2-D)) and unsaturates zones
#               2) snow accumulation, distribution and melt
#               3) river run off for a given catchment.
#               4)saves the states for each day where precip > 0, i. e where SCHADEX might insert simulated extreme precip
#
# Model input: 1) catchment physical/geographical description
#              2) weather data, temperature,precipitation and river discharge
#              3) model states variables (optional)
#
# Model output: simulated 1)saturated and unsaturated soil water
#                         2)snow water equivalent and snow distribution
#                         3)river discharge
#                         4)state files for each rainy day. Remember to include temperature for 5 days t-1, t=0, t+1, t+2, t+3 
#
# Running DDD:    It is possible to save model state variables and run the model
#                 starting from saved state variables
#----------------------------------------------------------------------------------


DDD_Modul_EB_Oct2018_func <- function(startsim, tprm,days,param.file,ptq.file,utfile,r2fil,modstate,savest,
                                        savestate,ebelements,UP,catchment,ebfile, OFDDD, kal)  
{

prm <- read.table(param.file)
ptqinn <- read.table(ptq.file)
navn <- list("yr","mnt","day","hr","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10","t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q") # i tilfelle av høydesoneverdier
names(ptqinn)[1:25] <-navn[1:25]                                                                                                                                # i tilfelle av høydesoneverdier

days <- length(ptqinn$yr)
savest <-  seq(544,1870,1)   #save state at this time (a value of i) , state file "state_dd_mm_yyyy_hh"
tell <- 1                    #counts number in savest vector.For each added state tell <- tell +1

#-------------------------------------preprocess start------------------------------------------
#-----------------------------------------------------------------------------------------------

#####################VARIABLES FOR 10 ELEVATION ZONES###########################################
isoil<- vector("numeric",10)# precipitation and snowmelt from the elevation zones
gwgt <- vector("numeric",10)# weights for glaciermelt pr each elevation zone
swgt <- vector("numeric",10)# weights for input to soils pr each elevation zone NB bogs are not part of this yet
gisoil <- vector("numeric",10)#glaciermelt from the elevation zones
spd<- vector("numeric",10)
swe_h <- vector("numeric",10)#SWE pr altitude level
wcd<- vector("numeric",10)
sca<- vector("numeric",10)
nsno<- vector("numeric",10)
alfa<- vector("numeric",10)
ny <- vector("numeric",10)
hprecip<- vector("numeric",10)
htemp<- vector("numeric",10)
hfelt<- vector("numeric",10)
snowfree <- vector("numeric",10)  #Indicator variable 1 is sca  less than gca 0 if sca > gca
tempvec <- vector("numeric",5) # temperature vector  to save in statefile
PR<- vector("numeric",10)
PS<- vector("numeric",10)
MW<- vector("numeric",10)
MWGLAC<- vector("numeric",10)
scaobx<- vector("numeric",10)
SWrad<- vector("numeric",10)
LA<- vector("numeric",10)
LT<- vector("numeric",10)
Pa<- vector("numeric",10)
#Addition for Enegy balance snowmelt

snro <- vector("numeric",10)            #Snow density
#sz <- vector("numeric",10)              #Snowdepth
melt <- vector("numeric",10)            #instead of degreeday melt (MW)
snowdepth <-vector("numeric",10)        #hmm..
sndens <-vector("numeric",10)           #hm ...
#initializing
snro[1:10]<-0.0 
#sz[1:10]<-0.0 
#melt[1:10]<-0.0 
snowdepth[1:10]<-0.0
sndens[1:10]<-0.0
# end # addition for EB snowmelt

#initialisering
isoil<- 0.0
gisoil[1:10] <-0.0
snowfree[1:10] <- 0.0
spd[1:10]<- 0.0
wcd[1:10]<- 0.0
sca[1:10]<- 0.0
nsno[1:10]<- 0.0
alfa[1:10]<- 0.0
ny[1:10] <- 0.0
hprecip[1:10]<- 0.0
htemp[1:10]<- 0.0
hfelt[1:10]<-0.0
PR[1:10]<-0.0
PS[1:10]<-0.0
MW[1:10]<-0.0
MWGLAC[1:10]<-0.0
scaobx[1:10]<-0.0
SWrad[1:10]<-0.0
LA[1:10]<-0.0
LT[1:10]<-0.0
Pa[1:10]<-0.0 

#Reading parameters
hfelt[1] <-prm$V2[3]+(prm$V2[4]-prm$V2[3])/2.0
hfelt[2] <-prm$V2[4]+(prm$V2[5]-prm$V2[4])/2.0
hfelt[3] <-prm$V2[5]+(prm$V2[6]-prm$V2[5])/2.0
hfelt[4] <-prm$V2[6]+(prm$V2[7]-prm$V2[6])/2.0
hfelt[5] <-prm$V2[7]+(prm$V2[8]-prm$V2[7])/2.0
hfelt[6] <-prm$V2[8]+(prm$V2[9]-prm$V2[8])/2.0
hfelt[7] <-prm$V2[9]+(prm$V2[10]-prm$V2[9])/2.0
hfelt[8] <-prm$V2[10]+(prm$V2[11]-prm$V2[10])/2.0
hfelt[9] <-prm$V2[11]+(prm$V2[12]-prm$V2[11])/2.0
hfelt[10] <-prm$V2[12]+(prm$V2[13]-prm$V2[12])/2.0
pro <-tprm[1]     #prm$V2[14]    #% liquid water content in snow
xcor <-prm$V2[15]    # x-coordinateUTM33
ycor <-prm$V2[16]    # y-coordinateUTM33
Aprim <-prm$V2[17]    # Initial Albedo 0.86 (Walter)
taux <-prm$V2[18]     # Age of snow  0.0 (UEB)
hst1 <-prm$V2[19]   #  used
u <-tprm[7]      # Wind speed [m/s]
hfeltmid <-prm$V2[21]  # mean elevation of catchment
tgrad <-prm$V2[22] # tgrad, pgrad   :Lapse rates pr 100 metres
pgrad <-prm$V2[23] # tgrad, pgrad   :Lapse rates pr 100 metres
pkorr <-tprm[2]   #prm$V2[24]# pskorr,prkorr   :Met corrections
skorr <-tprm[3]     #prm$V2[25]# pskorr,prkorr   :Met corrections
TX <-tprm[4] #prm$V2[26]   # TX              :Threshold temp for rain/snow
TS <-0.0 #prm$V2[27]   # TS              :Threshold temp. for melting/freezing
CX <-prm$V2[28]   # CX              :Degreedays factor for snowmelt
CFR <-prm$V2[29]  # CFR  	         :Correksjonsfactor of CX for refreezing
CGLAC <-tprm[5]   #prm$V2[30]# CGLAC          :Degreeday index for glaciermelt
a0 <-prm$V2[31]  #                 :Alfa null in snofall unit distribution
d <-prm$V2[32]   #                 :Decorrelation length (Measured in n) for precipitation
Timeresinsec <-prm$V2[33] # in seconds
MAD <-prm$V2[34]          # mean annual discharge
totarea <-prm$V2[35]      # total area in m2
maxLbog <- prm$V2[36]     #Maximum distance in wetlands
midLbog <-prm$V2[37]      # mean distance in wetlands
bogfrac <-prm$V2[38]      #areal fraction wetlands
zsoil <-prm$V2[39]        #frac zero dist soil
zbog <-prm$V2[40]         #frac zero dist bog
NoL <-prm$V2[41]          #Number of layers
cea <- -1000.0#tprm[6]            # prm$V2[42]         #Degreeday factor for evapotranspiration
R <-prm$V2[43]            #soil moisture content/100 for field capacity of soils 
Gshape <-prm$V2[44]       # shapeparametr celertity   Lower case lambda
Gscale <-prm$V2[45]       # scaleparameter celerity   Lower case lambda
GshInt <-prm$V2[46]       # shapeparameter saturation Capital lambda
GscInt <-prm$V2[47]       # scaleparameter saturation Capital lambda
cvHBV <- prm$V2[48]       # not used
dummy3 <- prm$V2[49]      # not used
rv <-tprm[6]              # prm$V2[50]           # celerity of water in rivernetwork
midFl <-prm$V2[51]        # mean of distance in river network
stdFl <-prm$V2[52]        # standard deviation of distances in river network
maxFl<-prm$V2[53]         # maximum distance in river network
maxDl <-prm$V2[54]        # maximum distance in hillslopes
CritFlux <- 250.0#tprm[8]       # prm$V2[55]    # Critical volume/time unit to create overland flow channels
midDL <- prm$V2[56]       # mean distance in hillslopes
glacfrac <- prm$V2[57]    # fraction of glaciers
midGl <- prm$V2[58]       # mean of distances of glaciers
stdGl <- prm$V2[59]       # standard deviation of distances of glaciers
maxGl <- prm$V2[60]       # maximum distance of glaciers
g1    <- prm$V2[61]       # areal fraction of glaciers in first elevation zone
g2    <- prm$V2[62]
g3    <- prm$V2[63]
g4    <- prm$V2[64]
g5    <- prm$V2[65]
g6    <- prm$V2[66]
g7    <- prm$V2[67]
g8    <- prm$V2[68]
g9    <- prm$V2[69]
g10   <- prm$V2[70]

#Merging assumed Glacier river network  with observed river network. We assume independent normal distributions
maxFl <- maxFl+ maxGl
midFl <-midFl + midGl  #mean of flowlength distribution
stdFl <-stdFl + stdGl

#tprm <- c(CX,TS,TX,pkorr, pro, tgrad, pgrad, cea) #parameters of your choosing to be written to the R2 file

#Constants
R <- 0.3            # Field capapcity from litterature (Skaugen and Onof, 2014)
hson <- 10          # Elevation zones
unitsnow <- 0.1     # mean of unit snow [mm]
n0 <- unitsnow*a0   # shape parameter of unit snow
gtcel <- 0.99       # threshold groundwater capacity -> Overland flow
CFR <- 2.5*(Timeresinsec/86400)*0.0833 #Fixed as 1/12 of estimate of CX for 24 hours 
len <-5*(86400/Timeresinsec)              #number og timestepes to estimate snowpack temperature. first factor is days
startsim <- startsim +len                #taking into account estimating snowpack temperature
STempvec <-vector("numeric",len)        #Temperature vector for estimating snowpack temperature
Pa[1:10] <-   Pa <-101.3*((293-0.0065*hfelt[1:10])/293)^5.26     # Air pressure as a function of height Zhu et al. 2012 and Stoy et al, 2018
MPa <- mean(Pa)
#vectors and matrixes
ddist <- vector("numeric", NoL)          # distribution of water in layers
ddistx <- vector("numeric", NoL)         # water present in layers
ddistxtemp <-vector("numeric", NoL)      # temporay of the above
aktMag <- vector("numeric", NoL)         # current water in Layers
qberegn <-vector("numeric",days)         # for calculation of skill scores etc. 
if(kal==0) simresult <- matrix(0.0, days,26)        # matrix which into the results are written
if(ebelements==1) energyel <- matrix(0.0, (days-startsim+1),194) #matrix to which the simulated energybalance elements are written
tell <-1                                 # counter for writning energyel properly
#grwpointresult <-matrix(0.0, days,29)    # matrix which into the groundwater point results are written

# Areas and weights
area <-(1-(bogfrac+glacfrac))*totarea    #area not covered by wetlands or glaciers. These three landscape types have different distancedistributions
area2 <- (1-(bogfrac))*totarea           #area in which we have hillslope process and rivernetwork processes
bogarea <- bogfrac*totarea               #area af wetlands
glacarea <-glacfrac *totarea             #area of glaciers
elevarea <- (totarea/hson)               #area pr elevationzone
gca <- c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10) #fraction of glacier pr elevation zone 
soilca <-1-gca
 # gwgt is the fraction of glaciated area pr elevation zone in relation to total glaciated area
ifelse(glacarea >0,gwgt[1:hson] <-gca[1:10]*elevarea/glacarea,gwgt[1:hson]<- 0.0) 
#Finds the fraction of soils (and bogs) in each elevation zone in relation to total soil (and bog) area
swgt[1:hson] <-soilca[1:10]*elevarea/(area + bogarea) 

#Subsurface celerities 
k <- Celerity_SubSurface(NoL, Gshape, Gscale, midDL, Timeresinsec) # Celerity of subsurface (and overland) flow
bogspeed <- k[1]*1.0                                                # Same celerity as overland flow

#Unit hydrgraphs for Wetlands
antBogsteps <- trunc(maxLbog/(bogspeed*Timeresinsec))+1 # Number of timsteps to drain wetlands
UHbog <- vector("numeric",antBogsteps)                  # Vector for Bog unit hydrograph
if (bogarea >0) UHbog <- SingleUH(bogspeed,Timeresinsec, midLbog, maxLbog, zbog)
ifelse(bogarea >0, UHbog, UHbog <-1)

#Unit hydrographs for hillslopes based on celerities and distance distributions,exponentially distributed
antHorlag <- trunc(maxDl/(k[NoL]*Timeresinsec))+1       # +1 includes the final day.  NB only in soils part
nodaysvector <- vector("numeric", NoL)                        # Number of days to drin each layer
layerUH <- matrix(0.0, ncol=antHorlag, nrow =NoL)
nodaysvector[1:NoL] <- trunc(maxDl/k[1:NoL]/Timeresinsec)+1   # integer time steps
for (i in 1: NoL)layerUH[i,1:nodaysvector[i]] <- SingleUH(k[i], Timeresinsec, midDL, maxDl, zsoil)

#UH for rivers, normal distributed 
#number of time units in river routing, rounded up. If set equal to one day, time is actually less
nodaysRiv <- trunc(maxFl/rv/Timeresinsec)+1
timeres <-seq(0,(nodaysRiv-1),1)                      # hourly resolution on riverrouting hydrogram
midFlscl <-midFl/rv/Timeresinsec
stdFlscl <-stdFl/rv/Timeresinsec
UHriver <- vector("numeric", nodaysRiv)
UHriver[1:nodaysRiv] <- 0.0
UHriver <- dnorm(timeres,midFlscl,stdFlscl)           # makes  the PDF
UHriver<-UHriver/sum(UHriver)                         # scale to give unity sum = 1
noDT <- length(UHriver)
if(nodaysRiv == 1) UHriver <-1.0                      # implies no river routing

#convolve, total UH for bogs and river, Fixed, does not change, is either on or off.
#BogRiverUH <- convolve(UHbog,UHriver, type="o")

vectorlengde <- antHorlag+noDT-1
qsimutx <-vector("numeric",noDT)  #vector for runoff i m3/s
qsimut <-vector("numeric",noDT)   #vector for runoff in mm
qBogut <-vector("numeric",noDT)    #vector for runoff from wetlands bogs
QRivx <- vector("numeric",noDT)   #all qsimutvectors are stored antHorlag timesteps back

#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#
Layers <- matrix(0.0, ncol=antHorlag,nrow=NoL)
BogLayers <- vector("numeric", antBogsteps)

FromLayerEst <- Layer_Estimation(GshInt,GscInt,Timeresinsec,maxDl,midDL, MAD, totarea, NoL, gtcel) 
     M <- FromLayerEst$M
Magkap <- FromLayerEst$Magkap

#Extracting areas and heights in pyramide-plot
FromAreas <- PyrAreas(NoL,totarea,maxDl,nodaysvector, layerUH, antHorlag) # A matrix of areas[in m2] for each time-step box, 
Areas <-FromAreas$Areas                                        # Important: same dimensions as Layers
delta_d <- FromAreas$delta_d                                   # The area drained pr time-step box for each layer


#States
if(modstate ==1) load(tilst.file)  # Run with state
if (modstate==0)                   # Run with initial values
{  
   sm <- 0.0                     # Soilmoisture
   smbog <-0.0                   # Saturation of Bogs
   
   # groundwater initial value
   for (i in 1:NoL) Layers[i,1:nodaysvector[i]] <- 0.0/nodaysvector[i]
   QRivx[1:noDT] <- qsimutx                             #discharge in m3/s
   totdef <-M
}
isoilx <- 0
Snowdepth <-0.0
snittT <- 0.0
fraconv <- ConvUTM_latlong(xcor,ycor) #Denne virker perfekt
phi <- fraconv$phi              #latitude radianer
thi <- fraconv$thi              #longitude radianer
############################################################################
#                                                                          #
#simulation starts here.............. days is no. time steps               #
############################################################################
for (i in startsim:days)
{

##  if ((i > 1) && (i< (days-3)) ) tempvec <- c(ptqinn$t1[i-1], ptqinn$t1[i], ptqinn$t1[i+1],ptqinn$t1[i+2],ptqinn$t1[i+3])# There have to be 10 of this one for EB calculations

dato <- paste(ptqinn$day[i],ptqinn$mnt[i], ptqinn$yr[i]) #date
if (Timeresinsec ==86400) ptqinn$hr[i] <- 12 
datohr <- paste(ptqinn$day[i],ptqinn$mnt[i], ptqinn$yr[i], ptqinn$hr[i]) #date 
DN <- strptime(dato,"%d%m%Y")$yday +1           #daynumber        
#print(paste("Timestep= ", i," max= ",days, " dato= ",datohr," DN= ", DN, sep=""))  
  
#Reads Precipitation and temperature for each elevation zone   
FromReadPrecipTemp <- ReadPrecipTemp(ptqinn[i,],hprecip,htemp)
htemp <-FromReadPrecipTemp$htemp
hprecip <- FromReadPrecipTemp$hprecip

meanprecip <-mean(hprecip[1:hson])
meantemp <- mean(htemp)
  
#Extracts SCA observations
#FromSCAobservation <- SCAobservation(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],scaob$dagsdato,hson,scaobx, scaob)
#scaobx <- FromSCAobservation$scaobx

for(idim in 1:hson)#elevation zones
{
  STempvec <- TemperatureVector(ptqinn[(i-len+1):i,15:24],idim, len, STempvec)
#browser()
#Calculates Rain vs Snow and Snow/Glaciermelt
  WaterIn <- nedb_eb(DN,hprecip[idim],htemp[idim],spd[idim],Timeresinsec,STempvec,TX,pkorr,skorr,ptqinn$hr[i],u,Pa[idim],CGLAC,taux,snittT,phi,thi)
#From WaterIn per elevation zone:         
  PR[idim] <- WaterIn$PR         # Precipitation as rain
  PS[idim] <- WaterIn$PS         # Precipitation as snow
  MW[idim] <- WaterIn$MW         # Snow melt potential
  MWGLAC[idim] <-WaterIn$MWGLA   # Glacier melt potential
  SWrad[idim] <-WaterIn$SWrad    # Short wave radiation (net)
  LA[idim] <- WaterIn$LA         # Long wave radiation atomspheric
  LT[idim] <-WaterIn$LT          # Long wave radiation terrestrial
  SH <-WaterIn$SH          # Turbulent heat
  LE <-WaterIn$LE          # Latent heat
  GH <-WaterIn$GH          # Ground heat
  PH <-WaterIn$PH          # Precipitation heat
  CC <-WaterIn$CC          # Cold content
  A <-WaterIn$A            # Albedo
  snittT <-WaterIn$snittT  # Snowpack temperature
  Tss <-WaterIn$Tss        # Snow surface temperature
  RH <- WaterIn$RH         # Relative humidity
  Cl <-WaterIn$Cl          # CloudCover
  taux <- WaterIn$taux     # Aging factor of snow for albedo calulations

# calculates Snow accumulation and Snowdistribution     
    FromSnow <-snogamma(PR[idim],PS[idim],MW[idim],sca[idim],scaobx[idim],spd[idim],wcd[idim],pro,nsno[idim],alfa[idim],ny[idim],a0,n0,a0,d,UP)
   gisoil[idim] <- MWGLAC[idim]
    isoil[idim] <- FromSnow$isoil  #precipitation and snowmelt
      spd[idim] <- FromSnow$spd 
      wcd[idim] <- FromSnow$wcd 
      sca[idim] <- FromSnow$sca 
     nsno[idim] <- FromSnow$nsno 
     alfa[idim] <- FromSnow$alfa 
       ny[idim] <- FromSnow$ny
     
    if (spd[idim] > 20000)
     {
      spd[idim] <-20000
      print("Glaciers are made. This speeds up calibration")
     }
    ifelse(sca[idim] < gca[idim],snowfree[idim] <- 1.0, snowfree[idim] <-0.0)#test for glaciermelt,  no melt if snowcovered, snowfree[idim] =0.0 

    #Snowdepths and snowdenities
    snowdepthx <- snowdepth[idim]
    snomag <-  (sca[idim]*spd[idim]) #snowreservir this timestep
    #Calling snowdepth routine
    FromSD <- newsnowsd(PR[idim],PS[idim],htemp[idim],TX,MW[idim],snomag,snomag,snowdepthx)
    #from newsnowdepth
    snowdepth[idim] <- FromSD$sndepth    #snowdepth
    sndens[idim] <- FromSD$dens          #snowdensity 
    
    #browser()
    
    if (ebelements==1)
    {
    if(idim==1)energyel[(i-startsim+1), 1:23] <-c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i],round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                  + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                  + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                  + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==2)energyel[(i-startsim+1), 24:42] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                   + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                   + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                   + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==3)energyel[(i-startsim+1), 43:61] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                   + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                   + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                   + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==4)energyel[(i-startsim+1), 62:80] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                   + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                   + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                   + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==5)energyel[(i-startsim+1), 81:99] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                   + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                   + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                   + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==6)energyel[(i-startsim+1), 100:118] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                     + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                     + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                     + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==7)energyel[(i-startsim+1), 119:137] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                     + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                     + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                     + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==8)energyel[(i-startsim+1), 138:156] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                     + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                     + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                     + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==9)energyel[(i-startsim+1), 157:175] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                     + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                     + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                     + round(Tss,3), round(RH,3), round(Cl,3))
    if(idim==10)energyel[(i-startsim+1), 176:194] <-c(round(sca[idim],3),round(sca[idim]*spd[idim],3),
                                                      + round(snowdepth[idim],3),round(sndens[idim],3),round(MW[idim],3), round(MWGLAC[idim],3),round(SWrad[idim],3),round(LA[idim],3),
                                                      + round(LT[idim],3),round(SH,3),round(LE,3),round(GH,3),round(PH,3),round(CC,3), round(A,3),round(snittT,3),
                                                      + round(Tss,3), round(RH,3), round(Cl,3))
    } # kal == 0   
    
} #for elevation zones
  MSWrad <- mean(SWrad)
  MLA <- mean(LA)
  MLT <- mean(LT)
  #Snowreservoir
  snomag <- mean(sca[1:hson]*spd[1:hson]) #mean catchment SWE ,must multiply with SCA to have arealvalues
  swe_h[1:hson] <- sca[1:hson]*spd[1:hson] #SWE pr. elevation zone
  middelsca <- sum(sca[1:hson])/hson
  snofritt <- 1-middelsca
  wcs <- mean(wcd[1:hson])

  #Estimate isoil from all elevation zones
  misoil <-sum(isoil*swgt)  #snowmelt and rain on soils. isoil is weighted by the fraction of soils  in relation to soil an bog area pr elevation band
 
  m_r_onglac <- sum(isoil*gwgt) #snowmelt and rain from glaciated area. isoil is weighted by the fraction of glaciers pr elevation band in relation to total glaciated area  
  #glacier melt (gisoil) in mm this timestep. m_r_onglac + outglac is total output from glacier
  outglac <-sum(gisoil[1:hson]*gwgt[1:hson]*snowfree[1:hson])#gwgt because it is going to be scaled by glacfrac later on

  #Estimating current capacity in Layers
  FromLayerCapUp <- Layer_capacity_update(Layers, nodaysvector, Magkap, NoL, ddistx, aktMag)
  ddistx <- FromLayerCapUp$ddistx
  
  # State of unsaturated zone, it can etiher be zero (complete saturation) or positive, the deficit. If negative we have overland flow
  totdef <- sum(ddistx[2:NoL])       # This may be negative if input outstrips deficits
  D <-max(totdef,0)
  
  toSoil <- misoil*(1-(glacfrac))+glacfrac*(outglac+m_r_onglac) # input from rain, snow and glaciers
  #if(toSoil < 0.0) browser() alltid positiv
  
  PotEvap <- Potential_evap_PT(meantemp, MSWrad, MLA, MLT, MPa)$ET
  #print(paste("PotEvap=",PotEvap, "eatemp =", eatemp))
  if(PotEvap < 0.0)PotEvap <- 0.0
  
  #calculating evapotranspiration directly from Layers and updating Layers
  FromLayerEvap <-  Layer_evap(Layers, nodaysvector, PotEvap, meantemp, D, sm,M, toSoil,layerUH,NoL)
  Layers <- FromLayerEvap$Layers
  ea <- FromLayerEvap$ea
  
  #Checking for capacity after evapotranspiration
  FromLayerCapUp <- Layer_capacity_update(Layers, nodaysvector, Magkap, NoL, ddistx, aktMag)
  ddistx <- FromLayerCapUp$ddistx
  
  # State of unsaturated zone, after Layers are updated (after evapotranspiration) 
  totdef <- sum(ddistx[2:NoL])       # This may be negative if input outstrips deficits
  D <-max(totdef,0)
  
# Call for the soilwater routine  calculates soil water (unsaturated)
  Fromsoil <- Unsaturated_ex_evap(toSoil, sm, M, D)
  outx <- Fromsoil$outx
    sm <- Fromsoil$sm

  #From Wetlands/Bogs
  Fromwetlands <- Wetlands_EB(misoil, meantemp,middelsca, smbog,M,PotEvap)
  outbog <- Fromwetlands$outbog
   smbog <- Fromwetlands$smbog
   eabog <- Fromwetlands$eabog
 
  #Estimating input to Layers according to capacity
  FromGrvInputDistribution <- Grv_input_distribution(outx,NoL,ddistx, ddist)
  ddist <- FromGrvInputDistribution$ddist
  
  #Calculating UH for overland flow Layer #1
  if (OFDDD ==1)
  {
  if(ddist[1]*outx > 0)# overland flow is confirmed for this event
    {
      FromOFLayer<- OverlandFlowDynamicDD(k,ddist,outx, layerUH, nodaysvector,NoL, midDL, CritFlux, Timeresinsec)
      layerUH[1,] <- FromOFLayer$OFlayer
    }
  }
  
  #if(i==3552)browser()
  GDT <- sum(Layers[,1])+ sum(ddist*outx*layerUH[,1])  #Groundwater to be discharged into the river network + this timesteps contribution
  GDT_Bog <- BogLayers[1] + outbog*UHbog[1]              #Bogwater to be discharged into the river network + this timesteps contribution
  
  #Updating the saturation Layers
  FromLayerUpdate <- Layer_update(ddist,outx, Layers, layerUH, nodaysvector, NoL)
  Layers <- FromLayerUpdate$Layers
  
  FromBogLayerUpdate <- BogLayer_update(outbog, BogLayers, UHbog, antBogsteps)#
  BogLayers <- FromBogLayerUpdate$BogLayers 
  
  lyrs <- sum(Layers)# gives the sum of layers after todays runoff has taken place
  boglyrs <- sum(BogLayers)
  
  #if(i==3566)browser()
  #GrWPoint <- GrW_Point(nodaysvector, i, NoL, Areas, Layers, totarea, mindist[RN,],delta_d)
  
  #Response from the soils part [m3/s] from todays input. This is a vector  giving todays runoff and antHorlag + nodaysRiv ahead
  qsimutx <- (((GDT/1000)*area2)/Timeresinsec)*UHriver           

  #Response from Wetlands
  ifelse(bogarea > 0, qBogut <- (((GDT_Bog/1000)*bogarea)/Timeresinsec)*UHriver, qBogut[1:noDT] <-0)  #m3/s
  
  #Total response
  qsimutx[1:noDT] <- qsimutx[1:noDT]+qBogut[1:noDT]#adding contribution from bogs

  #Updating the routing in river network
  FromRiverUpdate <- River_update(noDT,QRivx,qsimutx)
  QRivx <- FromRiverUpdate$QRivx
  QRD <- FromRiverUpdate$QRD
  qmm_state <- (sum(QRivx)-QRD)*(Timeresinsec*1000/totarea)       # This is also a reservoir [mm], water from todays event is stored for future runoff in the RN, relevant for WB.
    
  Qmm <- (QRD*Timeresinsec*1000/totarea)  #QRD in mm/day

  if (kal==0)
  {
  #Assigning outdata to vector, one vector for each timestep
   simresult[i, 1:26] <- c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i],round(meanprecip,4),round(meantemp,3),
                           +     ptqinn$q[i],round(QRD,6),round(middelsca,3),round(snomag,4),round((M-totdef),2),
                           round(totdef,2),round(sm,4),round(ea,4),round(outx,4),round(outbog,4),round(outglac,2),
                           round(m_r_onglac+outglac,2),round(lyrs,4), round(Qmm,4), round(smbog,4), round(eabog,4),round(boglyrs,4),round(qmm_state,4), round(wcs,4), round(PotEvap,4))
  }
  #grwpointresult[i, 1:29] <- c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i], round(GrWPoint[1:25],2))
  
  
                                
  qberegn[i] <-QRD      #Simulated runoff for current timestep. Used for skillscores etc.

  #---------------------------------------------------------------------------------------------------

  #Saving states, SWE, sm, Layers, etc in a RData file
  
  if(i == (savest[tell]) && savestate == 1) 
  {
    Monthlast<-ptqinn$mnt[savest[tell]]
    Daylast<-ptqinn$day[savest[tell]]
    Yearlast<-ptqinn$yr[savest[tell]]
    Hourlast<-ptqinn$hr[savest[tell]]
    tdate <- paste(Yearlast,".",Monthlast,".",Daylast,".",Hourlast,sep="") 
    tdate1 <- strptime(tdate, "%Y.%m.%d.%H")
    tilsdate <- format(tdate1,"%Y%m%d.%H")
    tilst.file2 <- paste(path2,"tilst\\","tilst_",i,"_",tilsdate,".rdata",sep="") 
    tilst <-c("nsno","alfa","ny","sm","spd","wcd","sca", "Layers", "QRivx","smbog", "totdef","tilsdate","k")
    save(list=tilst, file = tilst.file2)
     if(tell < 1870)tell <- tell +1
  }
  #---------------------------------------------------------------------------------------------

} # end of loop for number of timesteps in timeseries


 skillstart <- 365
 days2 <-days
 meansim <- mean(qberegn[skillstart:days2])
 meanobs <- mean(ptqinn$q[skillstart:days2])
# Computing skillscores NSE (R2), kge, bias
 KGE_crit <-KGE(qberegn[skillstart:days2],ptqinn$q[skillstart:days2],method="2012",out.type="full")
 kge1 <-round(KGE_crit$KGE.value,3)
 kge2 <-round(KGE_crit$KGE.elements[1],3)
 kge3 <-round(KGE_crit$KGE.elements[2],3)
 kge4 <-round(KGE_crit$KGE.elements[3],3)
 f2 <- sum((ptqinn$q[skillstart:days2]-qberegn[skillstart:days2])^2) #sum of square errors
 R2 <- round(1-((f2/length(ptqinn$q[skillstart:days2]))/var(ptqinn$q[skillstart:days2])),3)
#Recall, tprm is the parameter vector
 y <- c(R2,kge1,kge2, kge3,kge4, round(meansim/meanobs,3),tprm)
 cat(round(y,3),"\n",file=r2fil,append=TRUE,sep="\t")
 print(paste("R2=",ptq.file," =",round(R2,3),"KGE=",kge1,"msim/mobs=",round(meansim/meanobs,3),"days=",days2))

 if (kal==0)write.table(simresult, file=utfile,quote =FALSE, row.names=FALSE,na = "NA", col.names=FALSE,sep ="\t", eol="\n")
 #write.table(grwpointresult, file=grwutfile,quote =FALSE, row.names=FALSE,na = "NA", col.names=FALSE,sep ="\t", eol="\n")
 if (ebelements==1)write.table(energyel, file=ebfile,quote =FALSE, row.names=FALSE,na = "NA", col.names=FALSE,sep ="\t", eol="\n" )

if(kal==0)
{
 Monthlast<-ptqinn$mnt
 Daylast<-ptqinn$day
 Yearlast<-ptqinn$yr
 Hourlast<-ptqinn$hr
 tdate <- paste(Yearlast,".",Monthlast,".",Daylast,".",Hourlast,sep="") 
 tdate1 <- strptime(tdate, "%Y.%m.%d.%H")

 windows(19,12)
 par(mfrow=c(1,1))
 par(mar=c(1.5,1.5,1.5,1.5), oma=c(5,5,1,1))

##Plots Runoff
 start <-943
 stopp <- start + 210 #days2-16  #length(ptqinn$yr)
 plot(tdate1[start:stopp],simresult[start:stopp,7],type ="l",col=1,lwd=2.0, xaxt="n",xlab="",ylab="Q [m3/s]",ylim=c(0,1.0*max(simresult[,7]))) #Qsim
 axis.POSIXct(1, at=seq(tdate1[start],tdate1[stopp], by ="month"),format ="%Y%m%d.%H", labels=TRUE, cex.axis=0.9, las=2)
 lines(tdate1[start:stopp],simresult[start:stopp,8],type="l",lwd=2.0, col="red")
 #lines(tdate1[start:stopp],0.01*simresult[start:stopp,10],type="l",lwd=2.0, col="blue")# SWE
 #lines(tdate1[start:stopp],1.0*simresult[start:stopp,14],type="l",lwd=1.0, col="green") #ea
 #lines(tdate1[start:stopp],0.1*simresult[start:stopp,5],type="l",lwd=1.0, col="grey") #precip
 legend(tdate1[start],0.8*max(simresult[,8]), c("ObsQ", "SimQ"), col = c("black","red"),
       lty=1,lwd=2.5, text.col = "black", bg = "gray99")

}  
t11 <- Sys.time()
print(t11-t0)
t0
if(kal==0)
{
start2005 <- 365
psum <- sum(simresult[start2005:(start2005+365),5])#2005
print(psum)
easum <- sum(simresult[start2005:(start2005+365),14])#2005
print(easum)
print(paste ("ratio sum_ea/sum_p =",round(easum/psum,3)))
}

result3D <- NULL
result3D$qest <- qberegn
result3D$qobs <- ptqinn$q
result3D$NSE <- round(R2,3)
result3D$KGE1 <- kge1
result3D$KGE2 <- kge2
result3D$KGE3 <- kge3
result3D$KGE4 <- kge4
result3D$Mquant <- M
#result3D$Mmean <-mRes
#result3D$Mstd <- stdRes
result3D$bias <- round(meansim/meanobs,3)
result3D
} #end of func
