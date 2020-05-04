## ==============================================
## Non-invasive monitoring of drug action (NiMDA)
## 
## Analyse the data
## Dose-response curves over time of drug exposure
## Author : Anna Fesser
## date : 22.05.2019
## last modified : 21.04.2020
## ==============================================

# This script is intended to be used to build dose-response models from parasite counts produced by time-lapse imaging (Fesser et al. 2020)

# It is intended to be used with time-lapse data from high-content imaging
# The different parameters should be adapted to your needs
# Copyright (C) <2020> <Anna F. Fesser, Olivier Braissant, Francisco Olmo, 
# John M. Kelly, Pascal Mäser, Marcel Kaiser>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


# Load packages
library(drc) # for dose-response curves
library(readxl) # I use for reading xlxs
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble

# tidy up
rm(list = ls())

# load the functions
source(file = "R-Scripts/DRM_3versions_functions.R")

# load the data
load("Data/NiMDA_00_raw_data.RData")

# drug incubation time
GC.site.fixed$hpd[GC.site.fixed$drug != "0"] <- GC.site.fixed$hpi[GC.site.fixed$drug != "0"] - 24
GC.site.fixed$hpd[is.na(GC.site.fixed$hpd)] <- 0

GC.site.live$hpd[GC.site.live$drug != "0"] <- GC.site.live$hpi[GC.site.live$drug != "0"] - 24
GC.site.live$hpd[is.na(GC.site.live$hpd)] <- 0

GC.site.live.summary$hpd[GC.site.live.summary$drug != "0"] <- GC.site.live.summary$hpi[GC.site.live.summary$drug != "0"] - 24
GC.site.live.summary$hpd[is.na(GC.site.live.summary$hpd)] <- 0


## modelling - DRM-loop

# loop indices
# pn <- 1
# rd <- 1
# drex <- 1

plates.used <- unique(GC.site.live$plate.name)

# data frames to store simulation results
DRC.fixed.Nr.par.simulated <- data.frame(plate.name = character(), read.date = character(), 
                                         drug.exp = numeric(), 
                                         pos.mean = numeric(),
                                         EC50 = numeric(),
                                         Hill = numeric(),
                                         Bottom = numeric(),
                                         Top = numeric(),
                                         Molar.sim = numeric(), Molar.sim.log10 = numeric(),
                                         sim.gro = numeric(), sim.par = numeric())

DRC.fixed.sim.key <- data.frame()

# data frames to store parameter results
results.fixed.combi.GF.par.Nr <- data.frame(read.date = character(), plate.name = character(), 
                                            drug.exp = numeric(), 
                                            EC50 = numeric(), EC50lo = numeric(), EC50hi = numeric(), 
                                            EC90 = numeric(), EC90lo = numeric(), EC90hi = numeric(), 
                                            Hill = numeric(), Bottom = numeric(), Top = numeric(),
                                            Hill.fixed = logical(), Bottom.fixed = logical(), Top.fixed = logical(),
                                            model.vers = character())

for (pn in 1:length(plates.used)) {
  
  for (rd in 1:length(read.dates)){
    
    drug.exp <- as.numeric(levels(as.factor(GC.site.live$hpd[GC.site.live$plate.name == plates.used[pn] & 
                                                               GC.site.live$read.date == read.dates[rd]])))
    drug.exp <- drug.exp[drug.exp > 0]
    
    for (drex in 1:length(as.factor(drug.exp))){
      
      hpi.current <- as.numeric(levels(as.factor(GC.site.live$hpi[GC.site.live$plate.name == plates.used[pn] & 
                                                                    GC.site.live$read.date == read.dates[rd] &
                                                                    GC.site.live$hpd == as.factor(drug.exp[drex])])))
      
      pos.mean <- mean(GC.site.live$GF.par.Nr[GC.site.live$plate.name == plates.used[pn] & 
                                                GC.site.live$read.date == read.dates[rd] &
                                                GC.site.live$hpi == as.factor(hpi.current) & 
                                                GC.site.live$inoc.trypo > 0 &
                                                GC.site.live$drug == "0"])
      
      rel.growth <- GC.site.live$GF.par.Nr[GC.site.live$plate.name == plates.used[pn] & 
                                             GC.site.live$read.date == read.dates[rd] &
                                             GC.site.live$hpd == as.factor(drug.exp[drex])] * 100 / pos.mean
      
      concentration <- GC.site.live$Molar[GC.site.live$plate.name == plates.used[pn] & 
                                                 GC.site.live$read.date == read.dates[rd] &
                                                 GC.site.live$hpd == as.factor(drug.exp[drex])]
      
      # Concentration range to simulate DRCs
      conc.min <- min(concentration)
      conc.max <- max(concentration)
      Molar.sim.log10 <- seq(log10(conc.min) - 1 , log10(conc.max) + 1, by = 0.01)
      Molar.sim <- 10^(Molar.sim.log10)
      
      # modelling
      results.DR.modelling <- DRM.LL4.3versions(rel.growth = rel.growth,
                                                concentration = concentration,
                                                LL.4.m1 = c(NA, 0, 100, NA),
                                                LL.4.m2 = c(NA, NA, 100, NA),
                                                LL.4.m3 = c(NA, NA, NA, NA),
                                                sim.conc.log10 = Molar.sim.log10)
      
      # save results in data frame
      results <- results.DR.modelling$results
      results$read.date <- read.dates[rd]
      results$plate.name <- plates.used[pn]
      results$drug.exp <- drug.exp[drex]
      
      # append results to full data frame
      results.fixed.combi.GF.par.Nr <- rbind(results.fixed.combi.GF.par.Nr, results[, colnames(results.fixed.combi.GF.par.Nr)])
      
      sim.gro <- results.DR.modelling$fit.drc
      sim.par <- sim.gro/100 * pos.mean
      
      # save fitting results in new column in data fra,e
      fit.results.temp <- data.frame(plate.name = plates.used[pn], read.date = read.dates[rd], 
                                     drug.exp = drug.exp[drex], 
                                     pos.mean, 
                                     EC50 = results$EC50,
                                     Hill = results$Hill,
                                     Bottom = results$Bottom,
                                     Top = results$Top,
                                     Molar.sim, Molar.sim.log10,
                                     sim.gro, sim.par)
      
      DRC.fixed.Nr.par.simulated <- rbind(DRC.fixed.Nr.par.simulated, fit.results.temp)
      
    }
    
  }  
  
}

# Determine which EC50 values have a reasonable confidence interval
results.fixed.combi.GF.par.Nr$CI.EC50.okay <- log10(results.fixed.combi.GF.par.Nr$EC50) - log10(results.fixed.combi.GF.par.Nr$EC50lo) <= 1 & 
  log10(results.fixed.combi.GF.par.Nr$EC50hi) - log10(results.fixed.combi.GF.par.Nr$EC50) <= 1

# no NAs
results.fixed.combi.GF.par.Nr$CI.EC50.okay <- !is.na(results.fixed.combi.GF.par.Nr$CI.EC50.okay) & results.fixed.combi.GF.par.Nr$CI.EC50.okay


# write necessary files into .RData file
save(DRC.fixed.Nr.par.simulated, 
     results.fixed.combi.GF.par.Nr,
     file = "Data/NiMDA_02_Analysis_DRM.RData")
