## ==============================================
## Non-invasive monitoring of drug action (NiMDA)
## 
## figures
## 
## Author : Anna Fesser
## date : 18.05.2019
## last modified : 21.04.2020
## ==============================================

# This script is intended to be used to reproduce the figures from Fesser et al. 2020

# It is intended to be used with the data from all the following scripts
# NiMDA_01_Analysis_fold-change.R
# NiMDA_02_Analysis_DRM.R
# NiMDA_03_Analysis_SLOPES_to_death.R
# NiMDA_04_Analysis_Host_cells.R 
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
library(tidyverse) # contains ggplot2, dplyr, tidyr, readr, purrr, tibble
library(viridis) # colorblind-friendly palette
library(gtable)

# load the data
rm(list = ls())
load("Data/NiMDA_00_raw_data.RData") 
load("Data/NiMDA_01_Analysis_fold-change.RData")
load("Data/NiMDA_02_Analysis_DRM.RData")
load("Data/NiMDA_03_Analysis_SLOPES_to_death.RData")
load("Data/NiMDA_04_Analysis_Host_cells.RData") 

# grid labelling
# labels for plate names
plate_names_drugs <- list(
  "A" = "Benznidazole",
  "E" = "Posaconazole"
)
# labels for replicates
replicate_labels_read_dates <- list(
  "18-08-15" = "rep1",
  "18-09-22" = "rep2",
  "18-09-27" = "rep3"
)
# labels for dpi
dpi_labels <- list(
  "1" = "24 hpi",
  "2" = "48 hpi",
  "3" = "72 hpi",
  "4" = "96 hpi",
  "5" = "120 hpi",
  "6" = "144 hpi"
)

# labels for bio.rep
bio_rep_labels <- list(
  "1" = "rep1",
  "2" = "rep2",
  "3" = "rep3"
)

# labelling functions for plate names
labeller_plate_names <- function(variable,value){
  return(plate_names[value])
}
labeller_rd_drugs <- function(variable,value){
  if (variable=="read.date"){
    return(replicate_labels_read_dates[value])
  }else{
    return(plate_names_drugs[value])
  }
}

labeller_rep_dpi <- function(variable,value){
  if (variable=="dpi"){
    return(dpi_labels[value])
  }else{
    return(bio_rep_labels[value])
  }
}

## FIGURES
# Fig 2 B
Fig2B <- ggplot(data = subset(design.live, 
                     subset = read.date == "18-08-15" & 
                       plate.name == "E"),
       aes(column, row, colour = strain, fill = log10(Molar))) +
  geom_point(shape = 21, size = 12, stroke = 1.5) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
  scale_fill_viridis(name = "drug", limits = c(-11, -3), 
                     breaks = c(-9, -6, -3), labels = c("1 nM", "1 uM", "1 mM"), 
                     na.value = "white", direction = -1) +
  scale_colour_manual(name = "par", labels = c("no", "yes"), values = c("black", "green3")) +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        legend.position = "left", legend.key.height =  unit(1, "cm"),
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig2B_plate_design_legend.png",
       plot = Fig2B,
       width = 220, height = 120,
       unit = "mm",
       dpi = 320)

#Fig 3 A
# all parasites - live (GFP) with legend
Fig3A <- ggplot(data = subset(GC.site.live, subset = drug == "0"),
       aes(hpi, log10(GF.par.Nr), shape = read.date)) +
  geom_jitter(width = 1.2, alpha = 0.2, colour = "green4", size = 2) +
  theme_bw() +
  scale_y_continuous(name = "GF parasites", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  ggtitle(label = "Untreated parasites detected in live imaging") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file ="figure/pub/Fig3A_GFP_par_all_image.png", 
       plot = Fig3A,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)


# Fig 3 B
# all parasites - kDNA - with legend
Fig3B <- ggplot(data = subset(GC.site.fixed, subset = drug == "0"),
       aes(hpi, log10(Par.kin), shape = read.date)) +
  geom_jitter(width = 6, alpha = 0.3, size = 2, colour = "black") +
  theme_bw() +
  scale_y_continuous(name = "Parasites detected as kDNA", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  ggtitle(label = "Untreated parasites detected in fixed imaging") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig3B_kDNA_par.png",
       plot = Fig3B,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)

## mean parasite numbers with drugs
# Fig 3 C
# Benznidazole
Fig3C <- ggplot(data = subset(GC.site.live.summary, subset = plate.name %in% c("A")), 
       aes(hpi, log10(mean.GF.Par), colour = log10(Molar), shape = read.date)) + 
  theme_bw()+
  geom_point(alpha = 1, size = 2) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  scale_y_continuous(name = "GF parasites", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  scale_colour_viridis(name = "drug \nconc", limits = c(-11, -3), 
                       breaks = c(-9, -6, -3), labels = c("1 nM", "1 uM", "1 mM"), 
                       na.value = "black", direction = -1) +
  geom_vline(xintercept = 24, 
             colour = "#39568CFF") +
  geom_text(x = 24, y = 3.2, label = "drug \nadded", hjust = -0.1, 
            colour = "black", size = 7.5) +
  ggtitle(label = "Benznidazole-treated parasites") +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig3C_GFP_par_benz_image.png",
       plot = Fig3C,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)

# FIg 3 D
# Posaconazole
Fig3D <- ggplot(data = subset(GC.site.live.summary, subset = plate.name %in% c("E")), 
       aes(hpi, log10(mean.GF.Par), colour = log10(Molar), shape = read.date)) +
  theme_bw()+
  geom_point(alpha = 1, size = 2) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  scale_y_continuous(name = "GF parasites", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  scale_colour_viridis(name = "drug \nconc", limits = c(-11, -3), 
                       breaks = c(-9, -6, -3), labels = c("1 nM", "1 uM", "1 mM"), 
                       na.value = "black", direction = -1) +
  geom_vline(xintercept = 24, 
             colour = "#39568CFF") +
  geom_text(x = 24, y = 3.2, label = "drug \nadded", hjust = -0.1, 
            colour = "black", size = 7.5) +
  ggtitle(label = "Posaconazole-treated parasites") +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig3D_GFP_par_posa_image.png",
       plot = Fig3D,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)

# Fig 4 A
Fig4A <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                     subset = plate.name == "A"),
       aes(as.numeric(Molar.sim.log10), sim.gro, colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_line(lwd = 1) +
  scale_y_continuous(name = "relative growth [%]", breaks = c(0, 50, 100, 150), limits = c(-5, 155)) +
  scale_x_continuous(name = "drug concentration",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig4A_simulation_benz_gro.png", 
       plot = Fig4A,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)

# Fig 4 B
Fig4B <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                     subset = plate.name == "E"),
       aes(as.numeric(Molar.sim.log10), sim.gro, colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_line(lwd = 1) +
  scale_y_continuous(name = "relative growth [%]", breaks = c(0, 50, 100, 150), limits = c(-5, 155)) +
  scale_x_continuous(name = "drug concentration", limits = c(-11.5, -6.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig4B_simulation_posa_gro.png", 
       plot = Fig4B,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)

# prepare to get a good legend for Fig 4 A and B
Fig4Al <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                               subset = plate.name == "A"),
                 aes(as.numeric(Molar.sim.log10), sim.gro, colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_line(lwd = 0.1) +
  scale_y_continuous(name = "relative growth [%]", breaks = c(0, 50, 100, 150), limits = c(-5, 155)) +
  scale_x_continuous(name = "drug concentration",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of drug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        legend.position = "bottom",legend.key.width = unit(2, "cm"),
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
# transform into a gtable
Fig4Al.gtable <- ggplot_gtable(ggplot_build(Fig4Al)) 
p.Fig4Al.legend <- gtable(unit(8, "cm"), unit(1, "cm"))
# take out the legend element from the gtable
p.Fig4Al.legend <- gtable::gtable_add_grob(p.Fig4Al.legend, Fig4Al.gtable$grobs[[which(Fig4Al.gtable$layout$name == "guide-box")]],
                                                    1, 1, 1, 1)
ggsave(file = "figure/pub/Fig4A_legend.png",
       plot = p.Fig4Al.legend,
       width = 160, height = 20,
       unit = "mm",
       dpi = 320)


# EC50 values
# Fig 4 C
Fig4C <- ggplot(data = subset(results.fixed.combi.GF.par.Nr,
                     subset = plate.name == "A"),
       aes(as.numeric(drug.exp))) +
  geom_errorbar(aes(ymin = log10(EC50lo), 
                    ymax = log10(EC50hi)), width=.02, colour = "darkgrey") +
  geom_point(aes(as.numeric(drug.exp), log10(EC50), colour = CI.EC50.okay, shape = read.date), size = 2) +
  scale_y_continuous(name = "EC50 (log10)",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_x_continuous(name = "time of drug exposure (h)", breaks = seq(0, 144, by = 24)) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  scale_colour_manual(name = "95% CI meets inclusion criteria", 
                      breaks = c(TRUE, FALSE), labels = c("yes", "no"), values = c("red", "black"), guide = FALSE) +
  ggtitle(label = "Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig4C_EC50s_benz.png",
       plot = Fig4C,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)

Fig4Cl <- ggplot(data = subset(results.fixed.combi.GF.par.Nr,
                               subset = plate.name == "A"),
                 aes(as.numeric(drug.exp))) +
  geom_errorbar(aes(ymin = log10(EC50lo), 
                    ymax = log10(EC50hi)), width=.02, colour = "darkgrey") +
  geom_point(aes(as.numeric(drug.exp), log10(EC50), colour = CI.EC50.okay, shape = read.date), size = 5) +
  scale_y_continuous(name = "EC50 (log10)",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_x_continuous(name = "time of drug exposure (h)", breaks = seq(0, 144, by = 24)) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  scale_colour_manual(name = "95% CI meets inclusion criteria", 
                      breaks = c(TRUE, FALSE), labels = c("yes", "no"), values = c("red", "black")) +
  ggtitle(label = "Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        legend.position = "bottom", legend.key.width = unit(2, "cm"),
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
# transform into a gtable
Fig4Cl.gtable <- ggplot_gtable(ggplot_build(Fig4Cl)) 
p.Fig4Cl.legend <- gtable(unit(8, "cm"), unit(1, "cm"))
# take out the legend element from the gtable
p.Fig4Cl.legend <- gtable::gtable_add_grob(p.Fig4Cl.legend, Fig4Cl.gtable$grobs[[which(Fig4Cl.gtable$layout$name == "guide-box")]],
                                           1, 1, 1, 1)
ggsave(file = "figure/pub/Fig4C_legend.png",
       plot = p.Fig4Cl.legend,
       width = 200, height = 20,
       unit = "mm",
       dpi = 320)

# Fig 4 D
Fig4D <- ggplot(data = subset(results.fixed.combi.GF.par.Nr,
                     subset = plate.name == "E"),
       aes(as.numeric(drug.exp))) +
  geom_errorbar(aes(ymin = log10(EC50lo), 
                    ymax = log10(EC50hi)), width=.02, colour = "darkgrey") +
  geom_point(aes(as.numeric(drug.exp), log10(EC50), colour = CI.EC50.okay, shape = read.date), size = 2) +
  scale_y_continuous(name = "EC50 (log10)",  limits = c(-11.5, -6.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_x_continuous(name = "time of drug exposure (h)", breaks = seq(0, 144, by = 24)) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  scale_colour_manual(name = "95% CI meets inclusion criteria", 
                      breaks = c(TRUE, FALSE), labels = c("yes", "no"), values = c("red", "black"), guide = FALSE) +
  ggtitle(label = "Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig4D_EC50s_posa.png",  
       plot = Fig4D,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)


# Fig 5 A
# in one well on 
Fig5A <- ggplot(data = subset(GC.site.live, subset = plate.name == "A" & 
                       read.date == "18-09-27"),
       aes(hpi, log10(GF.par.Nr))) +
  geom_point(data = subset(GC.site.live, subset = plate.name == "A" & 
                             read.date == "18-09-27"),
             aes(hpi, log10(GF.par.Nr)), alpha = 0.3, colour = "grey", shape = 17, size = 2) +
  geom_smooth(data = subset(simulation, subset = plate.name == "A" & 
                              read.date == "18-09-27" & 
                              column == 3 & 
                              row == "E"),
              aes(hpi, log10(expm.hpi), group = ident), method = "lm", colour = "black", alpha = 0.5, lwd = 1) +
  geom_point(data = subset(GC.site.live, subset = plate.name == "A" & 
                             read.date == "18-09-27" & 
                             column == 3 & 
                             row == "E"),
             aes(hpi, log10(GF.par.Nr)), alpha = 1, colour = "green4", shape = 17, size = 2) +
  theme_bw() +
  scale_y_continuous(name = "GF parasites", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  ggtitle(label = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig5A_GFP_par_1well_highlighted_point.png", 
       plot = Fig5A,
       width = 240, height = 160,
       unit = "mm",
       dpi = 320)

# Fig 5 B
Fig5B <- ggplot(data = subset(expm.hpi.site.df.rep.sum, 
                     subset = Molar == 0),
       aes(dpd + 1, slope.mean, shape = read.date)) +
  geom_jitter(alpha = 0.5, width = 0.1, colour = "green4", size = 2) +
  scale_x_continuous(name = "time after infection [d]", breaks = 1:6) +
  scale_y_continuous(name = "fold \nchange/h", limits = c(-0.15, 0.15),
                     breaks = c(-0.1, 0, 0.1), 
                     labels = c( round(exp(-0.1), digits = 1), 1, round(exp(0.1), digits = 1))
  ) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  geom_hline(yintercept = 0, colour = "black", lty = "dashed") +
  ggtitle("without drug") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig5B_slopes_0drug.png", 
       plot = Fig5B,
       width = 200, height = 160,
       unit = "mm",
       dpi = 320)

# Fig 6 A
expm.hpi.site.df.sum$Molar[expm.hpi.site.df.sum$plate.name == "A" & expm.hpi.site.df.sum$row == "F"] <- expm.hpi.site.df.sum$Molar[expm.hpi.site.df.sum$plate.name == "A" & expm.hpi.site.df.sum$row == "E"]
Fig6A <- ggplot(data = subset(expm.hpi.site.df.sum,
                     subset = plate.name %in% c("A") &
                       dpd > 0),
       aes(log10(Molar), dpd, fill = slope.mean)) +
  geom_tile(data = subset(expm.hpi.site.df.sum,
                          subset = plate.name %in% c("A") &
                            dpd > 0),
            aes(log10(Molar), dpd, fill = slope.mean)) +
  scale_fill_gradient2(name = "fold \nchange/h", limits = c(-0.15, 0.15),
                       breaks = c(-0.1, 0, 0.1), 
                       labels = c( round(exp(-0.1), digits = 1), 1, round(exp(0.1), digits = 1)),
                       low = "red", high = "blue", mid = "white", midpoint = 0, guide = FALSE) +
  geom_point(data = subset(results.days.to.kill, 
                           subset = plate.name %in% c("A") &
                             Molar > 0),
             aes(log10(Molar), dpd, shape = read.date), alpha = 1, size = 5, position = position_dodge(width = 0.3)) +
  scale_y_continuous(name = "time of \ndrug exposure [d]") +
  scale_x_continuous(name = "drug concentration", limits = c(-7, -3), 
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig6A_heat_benz_sum_t2k.png", 
       plot = Fig6A,
       width = 210, height = 160,
       unit = "mm",
       dpi = 320)


# Fig 6 B
Fig6B <- ggplot(data = subset(expm.hpi.site.df.sum,
                     subset = plate.name %in% c("E") &
                       dpd > 0),
       aes(log10(Molar), dpd, fill = slope.mean)) +
  geom_tile() +
  geom_point(data = subset(results.days.to.kill, 
                           subset = plate.name %in% c("E") & 
                             Molar > 0),
             aes(log10(Molar), dpd, shape = read.date), alpha = 1, size = 5, position = position_dodge(width = 0.3)) +
  scale_y_continuous(name = "time of \ndrug exposure [d]") +
  scale_x_continuous(name = "drug concentration", limits = c(-11, -7.5), 
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_fill_gradient2(name = "fold \nchange/h", limits = c(-0.15, 0.15),
                       breaks = c(-0.1, 0, 0.1), 
                       labels = c( round(exp(-0.1), digits = 1), 1, round(exp(0.1), digits = 1)),
                       low = "red", high = "blue", mid = "white", midpoint = 0, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  ggtitle("Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/Fig6B_heat_posa_sum_t2k.png", 
       plot = Fig6B,
       width = 210, height = 160,
       unit = "mm",
       dpi = 320)

Fig6Bl <- ggplot(data = subset(expm.hpi.site.df.sum,
                              subset = plate.name %in% c("E") &
                                dpd > 0),
                aes(log10(Molar), dpd, fill = slope.mean)) +
  geom_tile() +
  geom_point(data = subset(results.days.to.kill, 
                           subset = plate.name %in% c("E") & 
                             Molar > 0),
             aes(log10(Molar), dpd, shape = read.date), alpha = 1, size = 5, position = position_dodge(width = 0.3)) +
  scale_y_continuous(name = "time of \ndrug exposure [d]") +
  scale_x_continuous(name = "drug concentration", limits = c(-11, -7.5), 
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_fill_gradient2(name = "fold \nchange/h", limits = c(-0.15, 0.15),
                       breaks = c(-0.1, 0, 0.1), 
                       labels = c( round(exp(-0.1), digits = 1), 1, round(exp(0.1), digits = 1)),
                       low = "red", high = "blue", mid = "white", midpoint = 0) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  ggtitle("Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
# transform into a gtable
Fig6Bl.gtable <- ggplot_gtable(ggplot_build(Fig6Bl)) 
Fig6B.legend <- gtable(unit(8, "cm"), unit(1, "cm"))
# take out the legend element from the gtable
Fig6B.legend <- gtable::gtable_add_grob(Fig6B.legend, Fig6Bl.gtable$grobs[[which(Fig6Bl.gtable$layout$name == "guide-box")]],
                                        1, 1, 1, 1)
ggsave(file = "figure/pub/Fig6B_legend.png", 
       plot = Fig6B.legend,
       width = 50, height = 160,
       unit = "mm",
       dpi = 320)


# Supp Fig 1
SuppFig1 <- ggplot(data = subset(GC.site.live, subset = column %in% c(2, 3, 4)),
       aes(hpi)) +
  geom_jitter(data = subset(GC.site.live, subset = column %in% c(2, 3, 4) &
                              read.date != "18-08-10" 
                            & plate.name %in% c("A", "E")),
              aes(hpi, log10(GF.par.Nr), shape = read.date),
              width = 1.2, alpha = 0.3, colour = "green4", size = 2) +
  geom_jitter(data = subset(GC.site.fixed, subset = column %in% c(2, 3, 4)), 
              aes(hpi, log10(Par.kin), shape = read.date),
              width = 1.2, alpha = 0.3, colour = "black", size = 2) +
  scale_y_continuous(name = "Parasites", breaks = c(0, 1, 2, 3, 4), limits = c(0, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time since infection [h]", breaks = seq(0, 264, by = 24), limits = c(-1, 148))+
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(plate.name ~ read.date, labeller = labeller_rd_drugs) +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig1_GFP_and_DAPI.png", 
       plot = SuppFig1,
       width = 500, height = 240,
       unit = "mm",
       dpi = 320)

# Supp Fig 2
SuppFig2A <- ggplot(data = subset(GC.site.fixed.age, subset = drug == "0" &
                       bio.rep > 0 &
                       dpi > 0),
       aes(hc.age, log10(host.nuc))) +
  geom_jitter(alpha = 0.25, width = 0.2, size = 2) +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(hc.fit)), colour = "red") +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(hc.lo)), colour = "red", lty = "dashed") +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(hc.up)), colour = "red", lty = "dashed") +
  scale_y_continuous(name = "Host cells detected", breaks = c(0, 1, 2, 3, 4), limits = c(1, 3.5),
                     labels = c("   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time from plating the host cells to fixation [d]")+
  ggtitle(label = "Untreated host cells detected in fixed imaging") +
  facet_grid(bio.rep ~ dpi, labeller = labeller_rep_dpi) +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig2A_host-cell_per_age.png", 
       plot = SuppFig2A,
       width = 480, height = 280,
       unit = "mm",
       dpi = 320)

SuppFig2B <- ggplot(data = subset(GC.site.fixed.age, subset = drug == "0" &
                       bio.rep > 0 &
                       dpi > 0),
       aes(hc.age, log10((Par.kin +0.1)/(host.nuc+0.1)))) +
  geom_jitter(alpha = 0.25, width = 0.2, size = 2) +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(pc.fit)), colour = "red") +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(pc.lo)), colour = "red", lty = "dashed") +
  geom_line(data = subset(sim.lm.hcage, subset = bio.rep >0),
            aes(hc.age, log10(pc.up)), colour = "red", lty = "dashed") +
  scale_y_continuous(name = "Parasites per host cell", breaks = c(-3, -2, -1, 0, 1, 2, 3, 4), #limits = c(1, 3.5),
                     labels = c("0.001", "0.01", "0.1", "   1", "  10", " 100", "1000", "10000")) +
  scale_x_continuous(name = "time from plating the host cells to fixation [d]")+
  ggtitle(label = "Untreated parasites per host cells detected in fixed imaging") +
  facet_grid(bio.rep ~ dpi, labeller = labeller_rep_dpi) +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig2B_par-cell_per_age.png", 
       plot = SuppFig2B,
       width = 480, height = 280,
       unit = "mm",
       dpi = 320)

# Supp Fig 3
SuppFig3A <- ggplot(data = subset(DRC.fixed.Nr.hc.simulated,
                     subset = plate.name == "A" & read.date != "18-08-10"),
       aes(as.numeric(Molar.sim.log10), sim.hc, colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_jitter(data = subset(GC.site.fixed.age,
                            subset = plate.name == "A" & read.date != "18-08-10"), 
              aes(log10(Molar), host.nuc, colour = drug.exp), alpha = 0.3, width = 0.15, size = 2)+
  geom_line(lwd = 1) +
  scale_y_continuous(name = "host cell numbers", breaks = c(0, 200, 400, 600, 800, 1000), limits = c(-5, 1000)) +#labels = c("0.1", "1", "10", "100", "1000")) + #, limits = c(-1, 3.5)
  scale_x_continuous(name = "drug concentration",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3),
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig3A_simulation_data_hc_benz.png", 
       plot = SuppFig3A,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)

SuppFig3B <- ggplot(data = subset(DRC.fixed.Nr.hc.simulated,
                     subset = plate.name == "E" & read.date != "18-08-10"),
       aes(as.numeric(Molar.sim.log10), sim.hc, colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_jitter(data = subset(GC.site.fixed.age,
                            subset = plate.name == "E" & read.date != "18-08-10"), 
              aes(log10(Molar), host.nuc, colour = drug.exp), alpha = 0.3, width = 0.15, size = 2)+
  geom_line(lwd = 1) +
  scale_y_continuous(name = "host cell numbers", breaks = c(0, 200, 400, 600, 800, 1000), limits = c(-5, 1000)) +#labels = c("0.1", "1", "10", "100", "1000")) + #, limits = c(-1, 3.5)
  scale_x_continuous(name = "drug concentration", limits = c(-11.5, -6.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig3B_simulation_data_hc_posa.png", 
       plot = SuppFig3B,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)

# Supp Fig 4A
SuppFig4A <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                     subset = plate.name == "A"),
       aes(as.numeric(Molar.sim.log10), log10(sim.par), colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_jitter(data = subset(GC.site.live,
                            subset = plate.name == "A"), 
              aes(log10(Molar), log10(GF.par.Nr), colour = drug.exp), alpha = 0.3, width = 0.15, size = 2)+
  geom_line(lwd = 1) +
  scale_y_continuous(name = "parasite numbers", breaks = c(-1, 0, 1, 2, 3), labels = c("0.1", "1", "10", "100", "1000"), limits = c(-0.5, 3.5)) + 
  scale_x_continuous(name = "drug concentration",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig4A_simulation_data_log_benz.png", 
       plot = SuppFig4A,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)

SuppFig4Al <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                                   subset = plate.name == "A"),
                     aes(as.numeric(Molar.sim.log10), log10(sim.par), colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_jitter(data = subset(GC.site.live,
                            subset = plate.name == "A"), 
              aes(log10(Molar), log10(GF.par.Nr), colour = drug.exp), alpha = 0.3, width = 0.15, size = 5)+
  geom_line(lwd = 0.1) +
  scale_y_continuous(name = "parasite numbers", breaks = c(-1, 0, 1, 2, 3), labels = c("0.1", "1", "10", "100", "1000"), limits = c(-0.5, 3.5)) + 
  scale_x_continuous(name = "drug concentration",  limits = c(-7.5, -2.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of drug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3")) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Benznidazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        legend.position = "bottom",legend.key.width = unit(2, "cm"),
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
# transform into a gtable
SuppFig4A.gtable <- ggplot_gtable(ggplot_build(SuppFig4Al)) 
SuppFig4A.legend <- gtable(unit(8, "cm"), unit(1, "cm"))
# take out the legend element from the gtable
SuppFig4A.legend <- gtable::gtable_add_grob(SuppFig4A.legend, SuppFig4A.gtable$grobs[[which(SuppFig4A.gtable$layout$name == "guide-box")]],
                                        1, 1, 1, 1)
ggsave(file = "figure/pub/SuppFig4A_legend.png", 
       plot = SuppFig4A.legend,
       width = 300, height = 20,
       unit = "mm",
       dpi = 320)

SuppFig4B <- ggplot(data = subset(DRC.fixed.Nr.par.simulated,
                     subset = plate.name == "E"),
       aes(as.numeric(Molar.sim.log10), log10(sim.par), colour = drug.exp, group = drug.exp, shape = read.date)) +
  geom_jitter(data = subset(GC.site.live,
                            subset = plate.name == "E"), 
              aes(log10(Molar), log10(GF.par.Nr), colour = drug.exp), alpha = 0.3, width = 0.15, size = 2)+
  geom_line(lwd = 1) +
  scale_y_continuous(name = "parasite numbers", breaks = c(-1, 0, 1, 2, 3), 
                     labels = c("0.1", "1", "10", "100", "1000"), limits = c(-1, 3.5)) + 
  scale_x_continuous(name = "drug concentration", limits = c(-11.5, -6.5),
                     breaks = c(-12, -11, -10, -9, -8, -7, -6,  -5, -4, -3), 
                     labels = c("1 pM", "10 pM", "100 pM", "1 nM", "10 nM", "100 nM", "1 uM", "10 uM", "100 uM", "1 mM")) +
  scale_colour_viridis(name = "time of \ndrug \nexposure (h)", breaks = seq(0, 144, by = 24), option = "magma", direction = -1, guide = FALSE) +
  scale_shape_manual(name = "biological \nreplicate", values = c(15, 16, 17), labels = c("1", "2", "3"), guide = FALSE) +
  facet_grid(read.date ~ ., labeller = labeller_rd_drugs) +
  ggtitle("Posaconazole") +
  theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.title = element_text(size = 24, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size = 20, vjust = -1), axis.text.y = element_text(size = 20, hjust = 1), 
        panel.border = element_rect(size = 1),
        legend.text = element_text(size = 20, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
        legend.position = "bottom",legend.key.width = unit(2, "cm"),
        strip.text.x = element_text(size = 24, angle = 0), strip.text.y = element_text(size = 24, angle = 270))
ggsave(file = "figure/pub/SuppFig4B_simulation_data_log_posa.png", 
       plot = SuppFig4B,
       width = 240, height = 400,
       unit = "mm",
       dpi = 320)
