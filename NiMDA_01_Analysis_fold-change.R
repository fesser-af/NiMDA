## ==============================================
## Non-invasive monitoring of drug action (NiMDA)
## 
## Analyse the data - Exponential models
## of the daily change in parasite numbers (live imaging)
## Author : Anna Fesser
## date : 17.05.2019
## last modified : 21.04.2020
## ==============================================

# This script is intended to be used to build exponential models from parasite counts produced by time-lapse imaging (Fesser et al. 2020)

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
library(reshape2)
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble

# load the data
rm(list = ls())
load("Data/NiMDA_00_raw_data.RData")

## modelling
# looping indices
# rd <- 1
# pn <- 1
# cl <- 2
# rw <- 2
# st <- 1

# data frames to store parameter results
expm.hpi.site.df <- data.frame(plate.name = character(0), read.date = character(0),
                               column = integer(0), row = character(0), site = integer(0), 
                               intercept = numeric(0), slope = numeric(0), 
                               p.intercept = numeric(0), p.slope = numeric(0), Rsq = numeric(0))

# data frames to store simulation results
simulation <- data.frame(hpi = numeric(0), expm.hpi = numeric(0), lm.hpi = numeric(0),
                         plate.name = character(0), column = numeric(0), row = character(0), site = numeric(0))

# loop through all parameter
for (pn in 1:length(levels(GC.site.live$plate.name))){
  
  for(rd in 1:length(read.dates)) {
    
    for (cl in 2:11){
      
      for (rw in 2:7){
        
        for (st in 1:9){
          
          # reduce the data frame to work with
          
          temp.df.input <- GC.site.live[GC.site.live$plate.name == levels(GC.site.live$plate.name)[pn] &
                                          GC.site.live$read.date == read.dates[rd] &
                                          GC.site.live$column == cl &
                                          GC.site.live$row == levels(GC.site.live$row)[rw] &
                                          GC.site.live$site == st, 
                                        c("hpi", "GF.par.Nr")]
          
          # build the model
          temp.expm.hpi.site <- lm(log(GF.par.Nr + 0.1) ~ hpi, data = temp.df.input)
          
          # store the important parameters of the model temporarily
          temp.model.par <- data.frame(plate.name = levels(GC.site.live$plate.name)[pn], 
                                       read.date = read.dates[rd],
                                       column = cl, 
                                       row = levels(GC.site.live$row)[rw], 
                                       site = st, 
                                       intercept = as.numeric(temp.expm.hpi.site$coefficients[1]), 
                                       slope = as.numeric(temp.expm.hpi.site$coefficients[2]), 
                                       p.intercept = summary(temp.expm.hpi.site)$coefficient[1, 4], 
                                       p.slope = summary(temp.expm.hpi.site)$coefficient[2, 4], 
                                       Rsq = summary(temp.expm.hpi.site)$r.squared)
          
          # add them to the data frame
          expm.hpi.site.df <- rbind(expm.hpi.site.df, temp.model.par)

          # hpi for simulation
          hpi <- seq(min(temp.df.input$hpi), max(temp.df.input$hpi), length.out = 100)
          
          # fit 
          expm.hpi <- exp(as.numeric(temp.expm.hpi.site$coefficients[1]) + as.numeric(temp.expm.hpi.site$coefficients[2]) * hpi)
          
          # store the fitted model values
          temp.result.df <- data.frame(hpi, expm.hpi, 
                                       plate.name = rep(levels(GC.site.live$plate.name)[pn], 100), 
                                       read.date = rep(read.dates[rd], 100 ),
                                       column = rep(cl, 100), 
                                       row = rep(levels(GC.site.live$row)[rw], 100), 
                                       site = rep(st, 100))
          
          # put it into general dataframe
          simulation <- rbind(simulation, temp.result.df)
          
        }
        
      }
      
    }
    
  }
  
}

# add column about dpi
dpi.key <- data.frame(row = LETTERS[seq(7, 2, -1)], dpi = seq(1, 6, 1))
expm.hpi.site.df <- merge(expm.hpi.site.df, dpi.key, by = "row")

# add drug concentrations per well
concentration.key <- design.live[, c("plate.name", "read.date", "row", "column", "drug", "Molar", "drug.colours")]
expm.hpi.site.df <- merge(expm.hpi.site.df, concentration.key, by = c("plate.name", "read.date", "row", "column"))

# add days of drug exposure
expm.hpi.site.df$dpd <- expm.hpi.site.df$dpi - 1

# summarise expm
expm.hpi.site.df.sum <- expm.hpi.site.df %>% # over all replicates
  group_by(plate.name, drug, Molar, dpd, row, column) %>%
  summarise(n.sum = n(),
            slope.mean = mean(slope),
            slope.sd = sd(slope),
            slope.se = sd(slope)/sqrt(n())
  ) %>%
  as.data.frame()

expm.hpi.site.df.rep.sum <- expm.hpi.site.df %>% # of each replicate
  group_by(plate.name, read.date, drug, Molar, dpd, row, column) %>%
  summarise(n.sum = n(),
            slope.mean = mean(slope),
            slope.sd = sd(slope),
            slope.se = sd(slope)/sqrt(n())
  ) %>%
  as.data.frame()

# add an identifier per site column
simulation$ident <- paste(simulation$row, sprintf("%02d", simulation$column), simulation$site, sep = ".")
simulation <- merge(simulation, concentration.key, by = c("plate.name", "read.date", "row", "column"))

# write necessary files into .RData file
save(GC.site.live.summary, GC.site.live, design.live, 
     expm.hpi.site.df, 
     expm.hpi.site.df.sum, expm.hpi.site.df.rep.sum,
     simulation,
     file = "Data/NiMDA_01_Analysis_fold-change.RData")

# tidy up environment
rm(list = ls())

