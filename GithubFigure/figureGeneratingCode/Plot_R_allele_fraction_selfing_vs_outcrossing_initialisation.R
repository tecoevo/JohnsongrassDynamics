library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(cowplot)
library(scales)
library(latex2exp)
library(ggh4x)

# Read files with simulation results
dynamics_outcrossing <- read.table("../data/Table_Deterministic_Initialisation_control3_selfing0_sexualReproduction1_seedbank1_cost300_kCost5_kHerb5.txt",
                                   header = TRUE, sep = ",", dec = ".")
dynamics_selfing <- read.table("../data/Table_Deterministic_Initialisation_control3_selfing1_sexualReproduction1_seedbank1_cost300_kCost5_kHerb5.txt",
                               header = TRUE, sep = ",", dec = ".")


# Dataframe contaning R allele frequencies for the different control regimes
R_allele_fraction <- data.frame(Season = 0:30,
                                R = c((dynamics_outcrossing$RSplants + 2*dynamics_outcrossing$RRplants) / (2*(dynamics_outcrossing$SSplants + dynamics_outcrossing$RSplants + dynamics_outcrossing$RRplants)),
                                      (dynamics_selfing$RSplants + 2*dynamics_selfing$RRplants) / (2*(dynamics_selfing$SSplants + dynamics_selfing$RSplants + dynamics_selfing$RRplants))),
                                Type = rep(c("outcrossing", "selfing"),
                                           each = 31))


# Properties for plotting
# Theme for the plots
th <- theme_ipsum(plot_margin = margin(3, 3, 3, 3), 
                  base_family = "Arial",
                  base_size = 10.5, axis_title_size = 11,
                  grid = FALSE,
                  ticks = TRUE,
                  axis_col = "black") +
  theme(legend.position="none",
        legend.text=element_text(size=rel(1)),
        axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(colour = "black",
                                   margin=unit(c(7,7,3.5,7), "pt")), 
        axis.text.y = element_text(colour = "black",
                                   margin=unit(c(7,7,7,3.5), "pt")),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        axis.ticks = element_line(size=2),
        axis.ticks.length = unit(-3.5, "pt"))

# Size of lines
s <- 1


# Dynamics of R allele fraction under different control regimes
# Plot R allele fraction
Plot_R_allele_fraction <- ggplot(data = R_allele_fraction, 
                                 aes(x=Season, y = R, group = Type, 
                                     color = Type)) +
  geom_line(size=0.5*s) +
  geom_point(size=s) +
  scale_x_continuous(limits = c(0,14), breaks = seq(from = 0, to = 14, by = 2), 
                     minor_breaks = seq(0, 30, by = 1),
                     guide = "axis_minor") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  scale_color_manual(name = "Proportion of self-plollination: ",
                     values = c("#3C4C59",  "#99ABBA"),
                     breaks = c("outcrossing", "selfing"),
                     labels = c("0 %","95 %" ), 
                     guide = "legend") +
  th +
  theme(legend.position = "bottom") +
  ylab(TeX("R allele fraction")) +
  xlab("Time (Years)")



ggsave("R_allele_fraction_selving_vs_outcrossing_initialisation.jpeg",
       plot = Plot_R_allele_fraction, scale = 0.75)


