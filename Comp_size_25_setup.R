# Set up file 


# Project organization------
library(renv)
# renv::init() # Initatlize the local library options for the R project. Do not run
# renv::install("ggplot2")  # Command to localy install packages

# R packages-----

library(conflicted) # To check for conflicts between packaged
library(ggplot2) # Data vidualization
library(plyr) # Data manipulation
library(dplyr) # Data manipulation
library(tidyr) # For data tiding
library(parallel) # parallize process
library(multcompView) # letters above boxplots
library(vegan) # Community analysis
library(stringr) # String manipulation in R
library(MASS) # Supportive statistics package
library(tibble)  # alternative to data.frames
library(VennDiagram) # Venn-diagram
library(gridExtra) # Arrange multiple figures
library(indicspecies) # Indicator analysis
library(ALDEx2) # Differential abundance analysis
library(RColorBrewer) # Color pallette
library(readxl) # Read excle files
library(ggpubr) # ggplots with common legends
library(corrplot) # Correlation analysis & heat maps
library(patchwork) # putting together several graphics
library(lme4) # linear mixed effect models
library(emmeans) # estimating marginal means, post-hoc analysis
library(lmerTest) # linear mixed effect models
library(car) # regression analysis
library(multcomp) # Pairwise comparisons
library(colorBlindness) # Colorblind friendly color palettes
library(cowplot) # Formating of graphics, multiple figures
library(colorspace) # Color generator
library(UpSetR) # Visualizing intersections


# Conflicts

conflicted::conflicts_prefer(dplyr::arrange)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)


bg_theme <- theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.3, vjust = 0, size = 14, color ="black"),
    axis.text.y = element_text(size = 14, color ="black"),
    text = element_text(size = 14, color ="black"),
    axis.title = element_text(size=14),
    legend.text = element_text(size=14)
  )
color.compost <- c("grey80", "grey40", "grey10")



# Maybe needed
#library(RColorBrewer)

#library(factoextra)
#library(rstatix)
#library(Hmisc)
#library(reshape2)
#library(ggpubr) # Graphic
#library(car) 
#library(gtools)
#library(emmeans)
#library(factoextra)

#library(multcomp)