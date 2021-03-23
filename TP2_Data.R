rm(list=ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(openxlsx)
library(readxl)
library(seasonal)
library(readr)

# Data PBI y Gasto
Data_FM <- read_csv("Data_FM.csv")

# Consumer Price Index ####

# Historic CPI
ipc <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=178.1_NL_GENERAL_0_0_13&limit=5000&format=csv"))
ipc <- ts(ipc$nivel_general, start = c(1943, 01), frequency = 12)
ipc <- diff(log(ipc))
ipc <- window(ipc, end = c(2006, 12))

# CPI, Province of San Luis 
ipc.sl <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=197.1_NIVEL_GENERAL_2014_0_13&limit=5000&format=csv"))
ipc.sl <- ts(ipc.sl$nivel_general, start = c(2005, 10), frequency = 12)
ipc.sl <- diff(log(ipc.sl))
ipc.sl <- window(ipc.sl, start = c(2007, 01), end = c(2012, 07))

# CPI, City of Buenos Aires
ipc.ba <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=193.1_NIVEL_GENERAL_JULI_0_13&limit=5000&format=csv"))
ipc.ba <- ts(ipc.ba$nivel_general, start = c(2012, 07), frequency = 12)
ipc.ba <- diff(log(ipc.ba))
ipc.ba <- window(ipc.ba, end = c(2016, 04))

# CPI, Greater Buenos Aires (INDEC)
ipc.gba <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=101.1_I2NG_2016_M_22&limit=5000&format=csv"))
ipc.gba <- ts(ipc.gba$ipc_2016_nivel_general, start = c(2016, 04), frequency = 12)
ipc.gba <- diff(log(ipc.gba))
ipc.gba <- window(ipc.gba, end = c(2016, 12))

# CPI, National (INDEC)
ipc.nac <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=145.3_INGNACNAL_DICI_M_15&limit=5000&format=csv"))
ipc.nac <- ts(ipc.nac$ipc_ng_nacional, start = c(2016, 12), frequency = 12)
ipc.nac <- diff(log(ipc.nac))

pc <- c(ipc, ipc.sl, ipc.ba, ipc.gba, ipc.nac)
pc <- c(1, cumprod(exp(pc)))
pc <- ts(pc, start = c(1943, 01), frequency = 12)
pc <- window(pc, start = c(2004, 01),end = c(2019, 12))

#pc <- 100 * (pc / mean(tail(pc, 12)))

remove(ipc, ipc.sl, ipc.ba, ipc.gba, ipc.nac)
pc <- aggregate(pc, nfrequency = 4, mean) #Trimestralizo la serie


#Declaro TS
gc <- ts(Data_FM$CP_K / 4, start = c(2004, 01), frequency = 4)
yk <- ts(Data_FM$PIB_K / 4, start = c(2004, 01), frequency = 4)
gs <- ts(Data_FM$GS_K / 4, start = c(2004, 01), frequency = 4)
gk <- ts(Data_FM$GK_K / 4, start = c(2004, 01), frequency = 4)

#Ventana de estimacion
gc <- window(gc, end = c(2019, 04))
yk <- window(yk, end = c(2019, 04))
gs <- window(gs, end = c(2019, 04))
gk <- window(gk, end = c(2019, 04))


#Deflacto Series
a<- matrix(1, nrow = 64)
ruben<- gs %/% a


for (i in 1:64) {
  infla<- matrix(0, nrow = 64)
  infla<- (pc[i] - pc[i-1]) / pc[i-1]
  infla<- infla * 100
}


#Desestacionalizo series
seas.adj <- seas(gs)
seas.adj2 <- seas(gk)
gs <- seas.adj$series$s11
gk <- seas.adj2$series$s11
rm(seas.adj, seas.adj2)