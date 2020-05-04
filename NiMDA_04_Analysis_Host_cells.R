## ==============================================
## Non-invasive monitoring of drug action (NiMDA)
## 
## Analyzing host cell numbers
##
## Author : Anna Fesser
## date : 30.03.2020
## last modified : 21.04.2020
## ==============================================

# This script is intended to be used to analyze the effects of host cell age and the drugs on host cell numbers (Fesser et al. 2020)

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
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble

# load the data
rm(list = ls())
load(file = "Data/NiMDA_00_raw_data.RData") 

# Statistical Analysis of host cell numbers (and parasites per host cells) over host cell age

# loop indices
# br <- "1"
# inf <- "8"

lm.hcage.par <- data.frame("bio.rep" = numeric(),
                           "dpi" = numeric(),
                           "read.date" = character(),
                           "hc.int" = numeric(),
                           "hc.int.p" = numeric(),
                           "hc.inc" = numeric(),
                           "hc.inc.p" = numeric(),
                           "hc.Rsq" = numeric(),
                           "ph.int" = numeric(),
                           "ph.int.p" = numeric(),
                           "ph.inc" = numeric(),
                           "ph.inc.p" = numeric(),
                           "ph.Rsq" = numeric())
sim.lm.hcage <- data.frame("bio.rep" = numeric(),
                           "dpi" = numeric(),
                           "plate.name" = character(),
                           "read.date" = character(),
                           "hc.age" = numeric(),
                           "hc.fit" = numeric(), "hc.lo"= numeric(), "hc.up"= numeric(),
                           "pc.fit"= numeric(), "pc.lo"= numeric(), "pc.up"= numeric())

for (br in levels(as.factor(GC.site.fixed.age$bio.rep))){
  for (inf in levels(as.factor(GC.site.fixed.age$dpi))){
    temp.df <- GC.site.fixed.age[GC.site.fixed.age$bio.rep == as.numeric(br) &
                                   GC.site.fixed.age$dpi == as.numeric(inf) &
                                   GC.site.fixed.age$drug == "0",
                                 c("plate.name", "read.date", "drug", "bio.rep", "dpi", "host.nuc", "Par.kin", "hc.age")]
    lm.hc.hcage <- lm(host.nuc ~ hc.age, temp.df)
    lm.pph.hcage <- lm((Par.kin + 0.1)/(host.nuc + 0.1) ~hc.age, temp.df)
    age.range <- seq(min(temp.df$hc.age[!is.na(temp.df$hc.age)]), max(temp.df$hc.age[!is.na(temp.df$hc.age)]), length.out = 21)
    sim.hc <- as.data.frame(predict(object = lm.hc.hcage, interval = "confidence"))
    colnames(sim.hc) <- c("hc.fit", "hc.lo", "hc.up")
    sim.ppc <- as.data.frame(predict(object = lm.pph.hcage, interval = "confidence"))
    colnames(sim.ppc) <- c("pc.fit", "pc.lo", "pc.up")
    temp.par <- data.frame("bio.rep" = as.numeric(br),
                           "dpi" = as.numeric(inf),
                           "read.date" = unique(temp.df$read.date[!is.na(temp.df$read.date)]),
                           "hc.int" = as.numeric(lm.hc.hcage$coefficients[1]),
                           "hc.int.p" = as.numeric(summary(lm.hc.hcage)$coefficients[1, 4]),
                           "hc.inc" = as.numeric(lm.hc.hcage$coefficients[2]),
                           "hc.inc.p" = as.numeric(summary(lm.hc.hcage)$coefficients[2, 4]),
                           "hc.Rsq" = as.numeric(summary(lm.hc.hcage)$r.squared),
                           "ph.int" = as.numeric(lm.pph.hcage$coefficients[1]),
                           "ph.int.p" = as.numeric(summary(lm.pph.hcage)$coefficients[1, 4]),
                           "ph.inc" = as.numeric(lm.pph.hcage$coefficients[2]),
                           "ph.inc.p" = as.numeric(summary(lm.pph.hcage)$coefficients[2, 4]),
                           "ph.Rsq" = as.numeric(summary(lm.pph.hcage)$r.squared))
    temp.sim <- data.frame("bio.rep" = temp.df$bio.rep[rownames(temp.df)%in% rownames(sim.ppc)],
                           "dpi" = temp.df$dpi[rownames(temp.df)%in% rownames(sim.ppc)],
                           "plate.name" = temp.df$plate.name[rownames(temp.df)%in% rownames(sim.ppc)],
                           "read.date" = temp.df$read.date[rownames(temp.df)%in% rownames(sim.ppc)],
                           "hc.age" = temp.df$hc.age[rownames(temp.df)%in% rownames(sim.ppc)],
                           "hc.fit" = sim.hc$hc.fit, "hc.lo"= sim.hc$hc.lo, "hc.up"= sim.hc$hc.up,
                           "pc.fit"= sim.ppc$pc.fit, "pc.lo"= sim.ppc$pc.lo, "pc.up"= sim.ppc$pc.up)
    lm.hcage.par <- rbind(lm.hcage.par, temp.par)
    sim.lm.hcage <- rbind(sim.lm.hcage, temp.sim)
    
  }
  
}


# save the data
save(GC.site.fixed.age, lm.hcage.par, sim.lm.hcage,
     file = "Data/NiMDA_04_Analysis_Host_cells.RData")

## dose-response modelling
# load the data
rm(list = ls())
load(file = "Data/NiMDA_04_Analysis_Host_cells.RData")

# load the functions
source(file = "R-Scripts/DRM_3versions_functions.R")

# drug incubation time
GC.site.fixed.age$hpd[GC.site.fixed.age$drug != "0"] <- GC.site.fixed.age$hpi[GC.site.fixed.age$drug != "0"] - 24
GC.site.fixed.age$hpd[is.na(GC.site.fixed.age$hpd)] <- 0
GC.site.fixed.age$drug.exp <- GC.site.fixed.age$hpd

## modelling - DRM-loop

# loop indices
# pn <- 1
# rd <- 1
# drex <- 1

plates.used <- unique(GC.site.fixed.age$plate.name)
read.dates <- unique(GC.site.fixed.age$read.date)

# data frames to store simulation results
DRC.fixed.Nr.hc.simulated <- data.frame(plate.name = character(), read.date = character(), 
                                         drug.exp = numeric(), 
                                         pos.mean = numeric(),
                                         EC50 = numeric(),
                                         Hill = numeric(),
                                         Bottom = numeric(),
                                         Top = numeric(),
                                         Molar.sim = numeric(), Molar.sim.log10 = numeric(),
                                         sim.gro = numeric(), sim.hc = numeric())

# data frames to store parameter results
results.fixed.combi.hc.Nr <- data.frame(read.date = character(), plate.name = character(), 
                                            drug.exp = numeric(), 
                                            EC50 = numeric(), EC50lo = numeric(), EC50hi = numeric(), 
                                            EC90 = numeric(), EC90lo = numeric(), EC90hi = numeric(), 
                                            Hill = numeric(), Bottom = numeric(), Top = numeric(),
                                            Hill.fixed = logical(), Bottom.fixed = logical(), Top.fixed = logical(),
                                            model.vers = character())

for (pn in 1:length(plates.used)) {
  
  for (rd in 1:length(read.dates)){
    
    drug.exp <- as.numeric(levels(as.factor(GC.site.fixed.age$hpd[GC.site.fixed.age$plate.name == plates.used[pn] & 
                                                               GC.site.fixed.age$read.date == read.dates[rd]])))
    drug.exp <- drug.exp[drug.exp > 0]
    
    for (drex in 1:length(as.factor(drug.exp))){
      
      hpi.current <- as.numeric(levels(as.factor(GC.site.fixed.age$hpi[GC.site.fixed.age$plate.name == plates.used[pn] & 
                                                                    GC.site.fixed.age$read.date == read.dates[rd] &
                                                                    GC.site.fixed.age$hpd == as.factor(drug.exp[drex])])))
      
      pos.mean <- mean(GC.site.fixed.age$host.nuc[GC.site.fixed.age$plate.name == plates.used[pn] & 
                                                GC.site.fixed.age$read.date == read.dates[rd] &
                                                GC.site.fixed.age$hpi == as.factor(hpi.current) & 
                                                GC.site.fixed.age$inoc.trypo > 0 &
                                                GC.site.fixed.age$drug == "0"])
      
      rel.growth <- GC.site.fixed.age$host.nuc[GC.site.fixed.age$plate.name == plates.used[pn] & 
                                             GC.site.fixed.age$read.date == read.dates[rd] &
                                             GC.site.fixed.age$hpd == as.factor(drug.exp[drex])] * 100 / pos.mean
      
      concentration <- GC.site.fixed.age$Molar[GC.site.fixed.age$plate.name == plates.used[pn] & 
                                            GC.site.fixed.age$read.date == read.dates[rd] &
                                            GC.site.fixed.age$hpd == as.factor(drug.exp[drex])]
      
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
      results.fixed.combi.hc.Nr <- rbind(results.fixed.combi.hc.Nr, results[, colnames(results.fixed.combi.hc.Nr)])
      
      sim.gro <- results.DR.modelling$fit.drc
      sim.hc <- sim.gro/100 * pos.mean
      
      # save fitting results in new column in data fra,e
      fit.results.temp <- data.frame(plate.name = plates.used[pn], read.date = read.dates[rd], 
                                     drug.exp = drug.exp[drex], 
                                     pos.mean, 
                                     EC50 = results$EC50,
                                     Hill = results$Hill,
                                     Bottom = results$Bottom,
                                     Top = results$Top,
                                     Molar.sim, Molar.sim.log10,
                                     sim.gro, sim.hc)
      
      DRC.fixed.Nr.hc.simulated <- rbind(DRC.fixed.Nr.hc.simulated, fit.results.temp)
      
    }
    
  }  
  
}

# # Determine which EC50 values have a reasonable confidence interval
results.fixed.combi.hc.Nr$CI.EC50.okay <- log10(results.fixed.combi.hc.Nr$EC50) - log10(results.fixed.combi.hc.Nr$EC50lo) <= 1 & 
  log10(results.fixed.combi.hc.Nr$EC50hi) - log10(results.fixed.combi.hc.Nr$EC50) <= 1

# no NAs
results.fixed.combi.hc.Nr$CI.EC50.okay <- !is.na(results.fixed.combi.hc.Nr$CI.EC50.okay ) & results.fixed.combi.hc.Nr$CI.EC50.okay 

# save the data
save(GC.site.fixed.age, lm.hcage.par, sim.lm.hcage, results.fixed.combi.hc.Nr, DRC.fixed.Nr.hc.simulated,
     file = "Data/NiMDA_04_Analysis_Host_cells.RData") 

rm(list = ls())

