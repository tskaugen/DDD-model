# Runs DDD2018_Modul_DynOF_func in a loop for calibration
# 05.02.2018
#
#
rm(list=ls())    # removes old variables
t0 <- Sys.time() # start time
library(hydroPSO)
library(hydroGOF)
library(hydroTSM)

source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_pakke\\DDDGitHub\\DDD_Modul_EB_Oct2018_func.r") #generic DDD version as a function 
#runoff dynamics
##########################################################################################
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SnowVar.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SnowGamma_SCAupdate.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SCAbasedSWEupdate.R")
#source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Unsaturated_evap_EB.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Unsaturated_ex_evap.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Wetlands_EB.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Grv_input_distribution.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Celerity_SubSurface.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Layer_Estimation.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Layer_capacity_update.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Layer_update.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Layer_evap.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\River_update.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\BogLayer_update.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SingleUH.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\RainSnowSnowmelt.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SCAobservation.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\ReadPrecipTemp.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\PyramidAreas.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\GrWPoint.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\OverlandFlowDynamicDD.R")
#########################################################################################
#energybalance components
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\nedb_eb.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\densityage.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\newsnowsd_eb.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\newsnowdensity_eb.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Smelt_eb.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\ConvUTM_latlong.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\solrad_trans_albedo.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Long_wave_rad.r")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Sensible_lat_heat.r")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Ground_prec_CC.r")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\SnowpackTemp.r")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\AlbedoUEB.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\AlbedoWalter.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\CloudCover.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\TemperatureVector.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Tss_dewpoint.R")
source("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Modul\\Functions\\Potential_evap_PT.R")

path <-"F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_pakke\\DDDGitHub\\"
pathinn <-"F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_pakke\\DDDGitHub\\" #Input data

st_file <- paste(path,"\\vfst_list_test.txt",sep="")

st_list <- scan(st_file, character(0))
TR <- "24h"                                                       #temporal resolution

for(i in 1:1) # Loop of runof ststions. In this case there is only one station
{
catchment <- st_list[i]

utfile <- paste(path,"res_OF_",catchment,TR,"_.txt",sep="")       # output file of hydrological simulations
ebfile <- paste(path,"Energy_elements_",catchment,"_.txt",sep="") # file in which you store the estimated EB elements. BIG file

ptq.file <-paste(pathinn,"\\",catchment,"_",TR,"ptq_kal.txt",sep="") 
ptqinn <- read.table(ptq.file)
navn <- list("yr","mnt","day","hr","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10","t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q")
names(ptqinn)[1:25] <-navn[1:25]
days <- length(ptqinn$yr) #actually the number of time steps in the input file
qobs <- ptqinn$q[1:days]

param.file <-paste(path,"best_par_",catchment,"_",TR,".txt",sep="")# I første omgang har vi bare ptq filer for 24 timer med elevation
prm <-read.table(param.file)
tprm <- c(prm$V2[14], prm$V2[24:26],prm$V2[30],prm$V2[50], prm$V2[20]) # no 30 er glacial melt degreeday factor

modstate <- 0     #Run model with state: modstate = 1 , no state: modstate = 0
savestate <- 0    #Save  state: yes = 1, no = 0
UP <- 0           # if 1: update from sattelite derived SCA; if 0: do not
OFDDD <- 0        # if 1: we include dynamic distance distributions for overland flow. if 0: we do not
kal <- 0          # calibrate or not. 1:calibrate, 0: not
ebelements <-0    # choose if yoou want to write the energybalance elements to a file

r2fil <-  paste(path,"test_R2_OF_",catchment,"_",TR,"_",OFDDD,"_.txt",sep="") # Skill score file, NSE, KGE etc
startsim <- 1
savest <-  seq(544,1870,1)

     indata = list(
	    startsim=startsim,
      tprm=tprm,
      days=days,
      param.file=param.file,
      ptq.file=ptq.file,
      utfile=utfile,
	    r2fil=r2fil,
      modstate=modstate,
      savest=savest,
      savestate=savestate,
      ebelements=ebelements, 
      UP=UP,
      catchment=catchment,
	    ebfile=ebfile,
      OFDDD=OFDDD,
      kal=kal
     )

   calib_wrapper_model <- function(tprm, indata) { #startsim, tprm,days,param.file,ptq.file,utfile,r2fil,modstate,savest,savestate,ebelements,UP,catchment,param_eb.file,ebfile, OFDDD, kal)
   	startsim=indata$startsim
      days=indata$days
      param.file=indata$param.file
      ptq.file=indata$ptq.file
      utfile=indata$utfile
	    r2fil=indata$r2fil
      modstate=indata$modstate
      savest=indata$savest
      savestate=indata$savestate
      ebelements =indata$ebelements 
      UP=indata$UP
      catchment=indata$catchment
      ebfile=ebfile
      OFDDD=indata$OFDDD
      kal=indata$kal
	    res <- DDD_Modul_EB_Oct2018_func(startsim, tprm,days,param.file,ptq.file,utfile,r2fil,modstate,savest,savestate,ebelements,UP,catchment,ebfile,OFDDD,kal)
    return(res$KGE1)   
   }  
   
    calib_single_wsh <- function(indata) #startsim, tprm,days,param.file,ptq.file,utfile,r2fil,tilst.file,modstate,savest,savestate,sca.file,UP,catchment,mindistfile,param_eb.file,ebfile, OFDDD, kal)
   {	
     #lparam <- c(0.01,0.4,0.4,-2.0, 0.25,0.2,0.5)#3h
     #uparam <- c(0.1,2.0,2.0,2.0,2.5,3.0,4.5)#3h
     lparam <- c(0.01,0.4,0.4,-2.0, 2.5,0.2,0.5)#24h
	   uparam <- c(0.1,2.0,2.0,2.0,8.5,3.0,4.5)#24h
     sparam<- (lparam + uparam)/2
     control <- list(MinMax="max", write2disk = FALSE, verbose = FALSE, npart=8, maxit=100)
     res_calib <- hydroPSO::hydroPSO(indata$tprm, fn=calib_wrapper_model, indata, lower=lparam, upper=uparam, control=control)
     
     return(res_calib)
    }

    if(kal==0)calib_wrapper_model(tprm, indata) # a single run 
    if(kal==1)calib_single_wsh(indata)          # calibrate
}    