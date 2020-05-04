## ==============================================
## Customized function
## 
## Dose-respnse curves
## iterative with different fixed values
## Author : Anna Fesser
## date : 22.05.2019
## last modified : 20.04.2020
## ==============================================

# This script is intended to be used to load a customized function to build dose-response models with three hierarchical parameter combinations m1-m2 (Fesser et al. 2020)

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



# load packages
library(drc) # for dose-response curves
library(reshape2)
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble


# those are the variables needed for the function

# rel.growth
# concentration

# LL.4.m1 <- c(NA, 0, 100, NA)
# LL.4.m2 <- c(NA, NA, 100, NA)
# LL.4.m3 <- c(-1, NA, 100, NA)

# sim.conc.log10


# and here comes the function
DRM.LL4.3versions<- function(rel.growth, concentration, 
                             LL.4.m1 = c(NA, 0, 100, NA),
                             LL.4.m2 = c(NA, NA, 100, NA),
                             LL.4.m3 = c(1, NA, 100, NA),
                             sim.conc.log10)
{
  ## Preparation: create storage objects
  # simulation results are stored in a vector
  fit.drc <- vector(mode = "numeric", length = length(sim.conc.log10))
  
  ## Modelling part
  
  # Coefficients and their naming in Prism:
  # b:Hill slope (in normal DRC this value is positive)
  # c:Bottom
  # d:Top
  # e:IC50
  
  # Fit LL4-model1 e.g. (fixed = c(NA, 0, 100, NA), names = c("b", "c", "d", "e")) [DRM.m1], IF POSSIBLE.
  DRM.m1 <- tryCatch( # tests first, whether DRM is calculatable
    {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m1, names = c("b", "c", "d", "e")))}, 
    error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
  
  # if this gives an error:
  if (DRM.m1 == "DRM not calculatable.")
  {
    
    # Fit LL4-model2 e.g. (fixed = c(NA, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m2], IF POSSIBLE.
    DRM.m2 <- tryCatch( # tests first, whether DRM is calculatable
      {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m2, names = c("b", "c", "d", "e")))}, 
      error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
    
    # if this gives an error:
    if (DRM.m2 == "DRM not calculatable.")
    {
      
      # Fit LL4-model3 e.g. (fixed = c(1, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m3], IF POSSIBLE.
      DRM.m3 <- tryCatch( # tests first, whether DRM is calculatable
        {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m3, names = c("b", "c", "d", "e")))}, 
        error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
      
      # if this gives an error:
      if (DRM.m3 == "DRM not calculatable.")
      {        
        
        # the otherwise determinable parameter are set to NA
        # take coefficients from DRM.m3
        
        # Are coefficients fixed?
        Hill.fixed <- !is.na(LL.4.m3[1])
        Bottom.fixed <- !is.na(LL.4.m3[2])
        Top.fixed <- !is.na(LL.4.m3[3])
        
        # e:EC50
        EC50 <- NA
        EC50lo <- NA
        EC50hi <- NA
        EC90 <- NA
        EC90lo <- NA
        EC90hi <- NA
        # b:Hill slope (absolute value)
        Hill <- NA
        # c:Bottom
        Bottom <- NA
        # d:Top
        Top <- NA
        
        # simulate the response level to the potential concentrations
        fit.drc <- rep(NA, length.out = length(sim.conc.log10))
        
        model.vers <- NA
        
      }else
      { # if DRM.32$b < 0:
        if (DRM.m3$coefficients[1] < 0)
        {        
          
          # the otherwise determinable parameter are set to NA
          # take coefficients from DRM.m3
          
          # Are coefficients fixed?
          Hill.fixed <- !is.na(LL.4.m3[1])
          Bottom.fixed <- !is.na(LL.4.m3[2])
          Top.fixed <- !is.na(LL.4.m3[3])
          
          # e:EC50
          EC50 <- NA
          EC50lo <- NA
          EC50hi <- NA
          EC90 <- NA
          EC90lo <- NA
          EC90hi <- NA
          # b:Hill slope (absolute value)
          Hill <- NA
          # c:Bottom
          Bottom <- NA
          # d:Top
          Top <- NA
          
          # simulate the response level to the potential concentrations
          fit.drc <- rep(NA, length.out = length(sim.conc.log10))
          
          model.vers <- NA
          
        }
        else
        {
          # take coefficients from DRM.m3
          
          # Are coefficients fixed?
          Hill.fixed <- !is.na(LL.4.m3[1])
          Bottom.fixed <- !is.na(LL.4.m3[2])
          Top.fixed <- !is.na(LL.4.m3[3])
          
          # Estimate effective doses
          ECs <- ED(DRM.m3,c(50,90), interval = "delta")
          # e:EC50
          EC50 <- ECs[1, 1]
          EC50lo <- ECs[1, 3]
          EC50hi <- ECs[1, 4]
          EC90 <- ECs[2, 1]
          EC90lo <- ECs[2, 3]
          EC90hi <- ECs[2, 4]
          # b:Hill slope (absolute value)
          if (Hill.fixed)
          {
            Hill <- LL.4.m3[1]
          }else
          {
            Hill <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "b:(Intercept)"])
          }
          # c:Bottom
          if (Bottom.fixed)
          {
            Bottom <- LL.4.m3[2]
          }else
          {
            Bottom <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "c:(Intercept)"])
          }
          # d:Top
          if (Top.fixed)
          {
            Top <- LL.4.m3[3]
          }else
          {
            Top <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "d:(Intercept)"])
          }
          
          # simulate the response level to the potential concentrations
          fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
          
          model.vers <- "m3"
        }
        
      }
      
    }else
    {
      
      # if DRM.m2$b < 0:
      if (DRM.m2$coefficients[1] < 0)
      {
        
        # Fit LL4-model3 e.g. (fixed = c(-1, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m3], IF POSSIBLE.
        DRM.m3 <- tryCatch( # tests first, whether DRM is calculatable
          {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m3, names = c("b", "c", "d", "e")))}, 
          error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
        
        # if this gives an error:
        if (DRM.m3 == "DRM not calculatable.")
        {        
          
          # the otherwise determinable parameter are set to NA
          # take coefficients from DRM.m3
          
          # Are coefficients fixed?
          Hill.fixed <- !is.na(LL.4.m3[1])
          Bottom.fixed <- !is.na(LL.4.m3[2])
          Top.fixed <- !is.na(LL.4.m3[3])
          
          # e:EC50
          EC50 <- NA
          EC50lo <- NA
          EC50hi <- NA
          EC90 <- NA
          EC90lo <- NA
          EC90hi <- NA
          # b:Hill slope (absolute value)
          Hill <- NA
          # c:Bottom
          Bottom <- NA
          # d:Top
          Top <- NA
          
          # simulate the response level to the potential concentrations
          fit.drc <- rep(NA, length.out = length(sim.conc.log10))
          
          model.vers <- NA
          
        }else
        { # if DRM.32$b < 0:
          if (DRM.m3$coefficients[1] < 0)
          {        
            
            # the otherwise determinable parameter are set to NA
            # take coefficients from DRM.m3
            
            # Are coefficients fixed?
            Hill.fixed <- !is.na(LL.4.m3[1])
            Bottom.fixed <- !is.na(LL.4.m3[2])
            Top.fixed <- !is.na(LL.4.m3[3])
            
            # e:EC50
            EC50 <- NA
            EC50lo <- NA
            EC50hi <- NA
            EC90 <- NA
            EC90lo <- NA
            EC90hi <- NA
            # b:Hill slope (absolute value)
            Hill <- NA
            # c:Bottom
            Bottom <- NA
            # d:Top
            Top <- NA
            
            # simulate the response level to the potential concentrations
            fit.drc <- rep(NA, length.out = length(sim.conc.log10))
            
            model.vers <- NA
            
          }
          else
          {
            # take coefficients from DRM.m3
            
            # Are coefficients fixed?
            Hill.fixed <- !is.na(LL.4.m3[1])
            Bottom.fixed <- !is.na(LL.4.m3[2])
            Top.fixed <- !is.na(LL.4.m3[3])
            
            # Estimate effective doses
            ECs <- ED(DRM.m3,c(50,90), interval = "delta")
            # e:EC50
            EC50 <- ECs[1, 1]
            EC50lo <- ECs[1, 3]
            EC50hi <- ECs[1, 4]
            EC90 <- ECs[2, 1]
            EC90lo <- ECs[2, 3]
            EC90hi <- ECs[2, 4]
            # b:Hill slope (absolute value)
            if (Hill.fixed)
            {
              Hill <- LL.4.m3[1]
            }else
            {
              Hill <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "b:(Intercept)"])
            }
            # c:Bottom
            if (Bottom.fixed)
            {
              Bottom <- LL.4.m3[2]
            }else
            {
              Bottom <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "c:(Intercept)"])
            }
            # d:Top
            if (Top.fixed)
            {
              Top <- LL.4.m3[3]
            }else
            {
              Top <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "d:(Intercept)"])
            }
            
            # simulate the response level to the potential concentrations
            fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
            
            model.vers <- "m3"
          }
          
        }
        
      }else
      {
        
        # DRM.m2$b >= 0:
        
        # take coefficients from DRM.m2 
        
        # Are coefficients fixed?
        Hill.fixed <- !is.na(LL.4.m2[1])
        Bottom.fixed <- !is.na(LL.4.m2[2])
        Top.fixed <- !is.na(LL.4.m2[3])
        
        # Estimate effective doses
        ECs <- ED(DRM.m2,c(50,90), interval = "delta")
        # e:EC50
        EC50 <- ECs[1, 1]
        EC50lo <- ECs[1, 3]
        EC50hi <- ECs[1, 4]
        EC90 <- ECs[2, 1]
        EC90lo <- ECs[2, 3]
        EC90hi <- ECs[2, 4]
        # b:Hill slope (absolute value)
        if (Hill.fixed)
        {
          Hill <- LL.4.m2[1]
        }else
        {
          Hill <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "b:(Intercept)"])
        }
        # c:Bottom
        if (Bottom.fixed)
        {
          Bottom <- LL.4.m2[2]
        }else
        {
          Bottom <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "c:(Intercept)"])
        }
        # d:Top
        if (Top.fixed)
        {
          Top <- LL.4.m2[3]
        }else
        {
          Top <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "d:(Intercept)"])
        }
        
        # simulate the response level to the potential concentrations
        fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
        
        model.vers <- "m2"
        
      }
      
    }
    
  }else
  {
    
    # if DRM.m1$b < 0:
    if (DRM.m1$coefficients[1] < 0)
    {
      
      # Fit LL4-model2 e.g. (fixed = c(NA, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m2], IF POSSIBLE.
      DRM.m2 <- tryCatch( # tests first, whether DRM is calculatable
        {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m2, names = c("b", "c", "d", "e")))}, 
        error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
      
      # if this gives an error:
      if (DRM.m2 == "DRM not calculatable.")
      {
        
        # Fit LL4-model3 e.g. (fixed = c(-1, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m3], IF POSSIBLE.
        DRM.m3 <- tryCatch( # tests first, whether DRM is calculatable
          {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m3, names = c("b", "c", "d", "e")))}, 
          error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
        
        # if this gives an error:
        if (DRM.m3 == "DRM not calculatable.")
        {        
          
          # the otherwise determinable parameter are set to NA
          # take coefficients from DRM.m3
          
          # Are coefficients fixed?
          Hill.fixed <- !is.na(LL.4.m3[1])
          Bottom.fixed <- !is.na(LL.4.m3[2])
          Top.fixed <- !is.na(LL.4.m3[3])
          
          # e:EC50
          EC50 <- NA
          EC50lo <- NA
          EC50hi <- NA
          EC90 <- NA
          EC90lo <- NA
          EC90hi <- NA
          # b:Hill slope (absolute value)
          Hill <- NA
          # c:Bottom
          Bottom <- NA
          # d:Top
          Top <- NA
          
          # simulate the response level to the potential concentrations
          fit.drc <- rep(NA, length.out = length(sim.conc.log10))
          
          model.vers <- NA
          
        }else
        { # if DRM.32$b < 0:
          if (DRM.m3$coefficients[1] < 0)
          {        
            
            # the otherwise determinable parameter are set to NA
            # take coefficients from DRM.m3
            
            # Are coefficients fixed?
            Hill.fixed <- !is.na(LL.4.m3[1])
            Bottom.fixed <- !is.na(LL.4.m3[2])
            Top.fixed <- !is.na(LL.4.m3[3])
            
            # e:EC50
            EC50 <- NA
            EC50lo <- NA
            EC50hi <- NA
            EC90 <- NA
            EC90lo <- NA
            EC90hi <- NA
            # b:Hill slope (absolute value)
            Hill <- NA
            # c:Bottom
            Bottom <- NA
            # d:Top
            Top <- NA
            
            # simulate the response level to the potential concentrations
            fit.drc <- rep(NA, length.out = length(sim.conc.log10))
            
            model.vers <- NA
            
          }
          else
          {
            # take coefficients from DRM.m3
            
            # Are coefficients fixed?
            Hill.fixed <- !is.na(LL.4.m3[1])
            Bottom.fixed <- !is.na(LL.4.m3[2])
            Top.fixed <- !is.na(LL.4.m3[3])
            
            # Estimate effective doses
            ECs <- ED(DRM.m3,c(50,90), interval = "delta")
            # e:EC50
            EC50 <- ECs[1, 1]
            EC50lo <- ECs[1, 3]
            EC50hi <- ECs[1, 4]
            EC90 <- ECs[2, 1]
            EC90lo <- ECs[2, 3]
            EC90hi <- ECs[2, 4]
            # b:Hill slope (absolute value)
            if (Hill.fixed)
            {
              Hill <- LL.4.m3[1]
            }else
            {
              Hill <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "b:(Intercept)"])
            }
            # c:Bottom
            if (Bottom.fixed)
            {
              Bottom <- LL.4.m3[2]
            }else
            {
              Bottom <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "c:(Intercept)"])
            }
            # d:Top
            if (Top.fixed)
            {
              Top <- LL.4.m3[3]
            }else
            {
              Top <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "d:(Intercept)"])
            }
            
            # simulate the response level to the potential concentrations
            fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
            
            model.vers <- "m3"
          }
          
        }
        
      }else
      {
        # if DRM.m2$b < 0:
        if (DRM.m2$coefficients[1] < 0)
        {
          
          # Fit LL4-model3 e.g. (fixed = c(-1, NA, 100, NA), names = c("b", "c", "d", "e")) [DRM.m3], IF POSSIBLE.
          DRM.m3 <- tryCatch( # tests first, whether DRM is calculatable
            {drm(rel.growth ~ concentration, fct = LL.4(fixed = LL.4.m3, names = c("b", "c", "d", "e")))}, 
            error=function(var.plus){return("DRM not calculatable.")}) # if yes, model is saved as such, if not the error message is saved
          
          # if this gives an error:
          if (DRM.m3 == "DRM not calculatable.")
          {        
            
            # the otherwise determinable parameter are set to NA
            # take coefficients from DRM.m3
            
            # Are coefficients fixed?
            Hill.fixed <- !is.na(LL.4.m3[1])
            Bottom.fixed <- !is.na(LL.4.m3[2])
            Top.fixed <- !is.na(LL.4.m3[3])
            
            # e:EC50
            EC50 <- NA
            EC50lo <- NA
            EC50hi <- NA
            EC90 <- NA
            EC90lo <- NA
            EC90hi <- NA
            # b:Hill slope (absolute value)
            Hill <- NA
            # c:Bottom
            Bottom <- NA
            # d:Top
            Top <- NA
            
            # simulate the response level to the potential concentrations
            fit.drc <- rep(NA, length.out = length(sim.conc.log10))
            
            model.vers <- NA
            
          }else
          { # if DRM.32$b < 0:
            if (DRM.m3$coefficients[1] < 0)
            {        
              
              # the otherwise determinable parameter are set to NA
              # take coefficients from DRM.m3
              
              # Are coefficients fixed?
              Hill.fixed <- !is.na(LL.4.m3[1])
              Bottom.fixed <- !is.na(LL.4.m3[2])
              Top.fixed <- !is.na(LL.4.m3[3])
              
              # e:EC50
              EC50 <- NA
              EC50lo <- NA
              EC50hi <- NA
              EC90 <- NA
              EC90lo <- NA
              EC90hi <- NA
              # b:Hill slope (absolute value)
              Hill <- NA
              # c:Bottom
              Bottom <- NA
              # d:Top
              Top <- NA
              
              # simulate the response level to the potential concentrations
              fit.drc <- rep(NA, length.out = length(sim.conc.log10))
              
              model.vers <- NA
              
            }
            else
            {
              # take coefficients from DRM.m3
              
              # Are coefficients fixed?
              Hill.fixed <- !is.na(LL.4.m3[1])
              Bottom.fixed <- !is.na(LL.4.m3[2])
              Top.fixed <- !is.na(LL.4.m3[3])
              
              # Estimate effective doses
              ECs <- ED(DRM.m3,c(50,90), interval = "delta")
              # e:EC50
              EC50 <- ECs[1, 1]
              EC50lo <- ECs[1, 3]
              EC50hi <- ECs[1, 4]
              EC90 <- ECs[2, 1]
              EC90lo <- ECs[2, 3]
              EC90hi <- ECs[2, 4]
              # b:Hill slope (absolute value)
              if (Hill.fixed)
              {
                Hill <- LL.4.m3[1]
              }else
              {
                Hill <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "b:(Intercept)"])
              }
              # c:Bottom
              if (Bottom.fixed)
              {
                Bottom <- LL.4.m3[2]
              }else
              {
                Bottom <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "c:(Intercept)"])
              }
              # d:Top
              if (Top.fixed)
              {
                Top <- LL.4.m3[3]
              }else
              {
                Top <- as.numeric(DRM.m3$coefficients[names(DRM.m3$coefficients) == "d:(Intercept)"])
              }
              
              # simulate the response level to the potential concentrations
              fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
              
              model.vers <- "m3"
            }
            
          }
          
        }else
        {
          
          # DRM.m2$b >= 0:
          
          # take coefficients from DRM.m2 
          
          # Are coefficients fixed?
          Hill.fixed <- !is.na(LL.4.m2[1])
          Bottom.fixed <- !is.na(LL.4.m2[2])
          Top.fixed <- !is.na(LL.4.m2[3])
          
          # Estimate effective doses
          ECs <- ED(DRM.m2,c(50,90), interval = "delta")
          # e:EC50
          EC50 <- ECs[1, 1]
          EC50lo <- ECs[1, 3]
          EC50hi <- ECs[1, 4]
          EC90 <- ECs[2, 1]
          EC90lo <- ECs[2, 3]
          EC90hi <- ECs[2, 4]
          # b:Hill slope (absolute value)
          if (Hill.fixed)
          {
            Hill <- LL.4.m2[1]
          }else
          {
            Hill <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "b:(Intercept)"])
          }
          # c:Bottom
          if (Bottom.fixed)
          {
            Bottom <- LL.4.m2[2]
          }else
          {
            Bottom <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "c:(Intercept)"])
          }
          # d:Top
          if (Top.fixed)
          {
            Top <- LL.4.m2[3]
          }else
          {
            Top <- as.numeric(DRM.m2$coefficients[names(DRM.m2$coefficients) == "d:(Intercept)"])
          }
          
          # simulate the response level to the potential concentrations
          fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
          
          model.vers <- "m2"
          
        }
        
      }  
      
    }else
    {
      # DRM.m1$b >= 0:
      # take coefficients from DRM.m1
      
      # Are coefficients fixed?
      Hill.fixed <- !is.na(LL.4.m1[1])
      Bottom.fixed <- !is.na(LL.4.m1[2])
      Top.fixed <- !is.na(LL.4.m1[3])
      
      # Estimate effective doses
      ECs <- ED(DRM.m1,c(50,90), interval = "delta")
      # e:EC50
      EC50 <- ECs[1, 1]
      EC50lo <- ECs[1, 3]
      EC50hi <- ECs[1, 4]
      EC90 <- ECs[2, 1]
      EC90lo <- ECs[2, 3]
      EC90hi <- ECs[2, 4]
      # b:Hill slope (absolute value)
      if (Hill.fixed)
      {
        Hill <- LL.4.m1[1]
      }else
      {
        Hill <- as.numeric(DRM.m1$coefficients[names(DRM.m1$coefficients) == "b:(Intercept)"])
      }
      # c:Bottom
      if (Bottom.fixed)
      {
        Bottom <- LL.4.m1[2]
      }else
      {
        Bottom <- as.numeric(DRM.m1$coefficients[names(DRM.m1$coefficients) == "c:(Intercept)"])
      }
      # d:Top
      if (Top.fixed)
      {
        Top <- LL.4.m1[3]
      }else
      {
        Top <- as.numeric(DRM.m1$coefficients[names(DRM.m1$coefficients) == "d:(Intercept)"])
      }
      
      # simulate the response level to the potential concentrations
      fit.drc <- Bottom + (Top - Bottom) / ( 1 + exp(Hill * (log(10^sim.conc.log10) - log(EC50))))
      
      model.vers <- "m1"
      
    }
    
  }
  
  # data frame to write parameters & coefficients into
  results <- data.frame(EC50, EC50lo, EC50hi, EC90, EC90lo, EC90hi,
                        Hill, Bottom, Top,
                        Hill.fixed, Bottom.fixed, Top.fixed, model.vers)
  
  toReturn <- list(results = data.frame(EC50, EC50lo, EC50hi, EC90, EC90lo, EC90hi,
                                        Hill, Bottom, Top,
                                        Hill.fixed, Bottom.fixed, Top.fixed, model.vers), 
                   fit.drc = fit.drc)
  
  # this will be returned
  return(toReturn)
}
