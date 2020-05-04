## ==============================================
## Non-invasive monitoring of drug action (NiMDA)
## 
## Analyse the data
## Define DYING (tipping point of drug action)
## Author : Anna Fesser
## date : 25.05.2019
## last modified : 21.04.2020
## ==============================================

# This script is intended to be used to build define the tipping point of drug action from time-lapse parasite counts (Fesser et al. 2020)

# It is intended to be used with exponential model results from time-lapse data from high-content imaging
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
library(reshape2)
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble


# load the data
rm(list = ls())
load("Data/NiMDA_00_raw_data.RData") # this returns all necessary data for this analysis
load("Data/NiMDA_01_Analysis_fold-change.RData")

## modelling
# looping indices
#pn <- 1
#rd <- 1
#rw <- 2
#cl <- 2

# data frames to store results
expm.hpi.site.df$dpd <- expm.hpi.site.df$dpi - 1
results.slope.DYING <- data.frame(plate.name = character(0), read.date = character(0), column = numeric(0), row = character(0), 
                                  dpd = numeric(0), Molar = numeric(), slope.mean = numeric(), slope.up.CI = numeric())

expm.hpi.site.df$Molar[expm.hpi.site.df$dpd == 1] <- expm.hpi.site.df$Molar[expm.hpi.site.df$dpd == 2]


# loop through all parameter
for (pn in 1:length(levels(GC.site.live$plate.name))){
  
  for(rd in 1:length(read.dates)) {
    
    for (rw in 2:7){
      
      for (cl in 2:10){
        
        if (!(pn == 4 & rd == 4)){
          
          temp.slopes <- expm.hpi.site.df[expm.hpi.site.df$plate.name == levels(as.factor(expm.hpi.site.df$plate.name))[pn] &
                                            expm.hpi.site.df$read.date == read.dates[rd] &
                                            expm.hpi.site.df$row == levels(GC.site.live$row)[rw] &
                                            expm.hpi.site.df$column == cl, 
                                          "slope"]
          
          temp.n <- length(temp.slopes)
          temp.mean <- mean(temp.slopes)
          temp.se <- sd(temp.slopes)/sqrt(temp.n)
          temp.CI.upper <- temp.mean + qt(0.95, df = temp.n) * temp.se
          
          temp.dpd <- unique(expm.hpi.site.df[expm.hpi.site.df$plate.name == levels(as.factor(expm.hpi.site.df$plate.name))[pn] &
                                                expm.hpi.site.df$read.date == read.dates[rd] &
                                                expm.hpi.site.df$row == levels(GC.site.live$row)[rw] &
                                                expm.hpi.site.df$column == cl, 
                                              "dpd"])
          
          temp.Molar <- unique(expm.hpi.site.df[expm.hpi.site.df$plate.name == levels(as.factor(expm.hpi.site.df$plate.name))[pn] &
                                                       expm.hpi.site.df$read.date == read.dates[rd] &
                                                       expm.hpi.site.df$row == levels(GC.site.live$row)[rw] &
                                                       expm.hpi.site.df$column == cl, 
                                                     "Molar"])
          
          temp.df <- data.frame(plate.name = as.character(levels(as.factor(expm.hpi.site.df$plate.name))[pn]), 
                                read.date = as.character(read.dates[rd]), 
                                column = cl, row = as.character(levels(GC.site.live$row)[rw]), 
                                dpd = as.numeric(temp.dpd), Molar = as.numeric(temp.Molar), 
                                slope.mean = as.numeric(temp.mean), slope.up.CI = as.numeric(temp.CI.upper))
          
          if (temp.CI.upper < 0){
            
            results.slope.DYING <- rbind(results.slope.DYING, temp.df)
            
          }
          
        }
        
      }
      
    }
    
  }
  
}



## Find the first time to kill per concentration
# looping indices
#pn <- 1
#rd <- 1
#con <- 1

results.days.to.kill <- data.frame(plate.name = character(), read.date = character(),
                                  Molar = numeric(), dpd = numeric())

for (pn in 1:length(levels(GC.site.live$plate.name))){
  
  for(rd in 1:length(read.dates)){
    
    temp.conc <- unique(results.slope.DYING[results.slope.DYING$plate.name == levels(as.factor(expm.hpi.site.df$plate.name))[pn] &
                                              results.slope.DYING$read.date == read.dates[rd], 
                                            "Molar"])
    
    temp.conc <- temp.conc[!is.na(temp.conc)]
    
    if (length(temp.conc) > 0){
      
      for (con in 1:length(temp.conc)){
        
        temp.dpd <- min(results.slope.DYING[results.slope.DYING$plate.name == levels(as.factor(expm.hpi.site.df$plate.name))[pn] &
                                              results.slope.DYING$read.date == read.dates[rd] &
                                              results.slope.DYING$Molar == temp.conc[con], 
                                            "dpd"])
        
        temp.df <- data.frame(plate.name = levels(as.factor(expm.hpi.site.df$plate.name))[pn], read.date = read.dates[rd],
                              Molar = temp.conc[con], dpd = temp.dpd)
        
        results.days.to.kill <- rbind(results.days.to.kill, temp.df)
        
      }
      
    }
  
  }

}

#results.days.to.kill
results.days.to.kill <- merge(results.days.to.kill, expm.hpi.site.df.rep.sum, all.y = FALSE, all.x = TRUE)

# write necessary files into .RData file
save(expm.hpi.site.df, expm.hpi.site.df.rep.sum, 
     results.slope.DYING, results.days.to.kill,
     file = "Data/NiMDA_03_Analysis_SLOPES_to_death.RData")


