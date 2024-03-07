#This is the code to create 2D microbial communities and reproduce figures published in Schmidt et al paper. See comments in .R files for details.
# Please not that calculating the %HSA covered can tale a long time for dense communities. 
# Comment line 76 and 109 if this info is not needed. 

library(spatstat) 
library(dplyr)
library(pals) # This is to have nice colour scales
library(ggplot2) 
library(ggforce) # Needed to plot interaction distance discs in subplots
library(patchwork)


source("R/populate_species.R")
source("R/make_community.R")

#Parameters for simulations
density = c(649, 577, 1313, 493, 216) # cell density in images
seed= c(11,2,1,1,2) # random seeds to obtain the same images as in the Schmidt et al

dimxy = 250 # size of main window (in µm)

props = c(Acidobacteria = 0.21,Actinobacteria = 0.04, Proteobacteria = 0.34, Firmicutes=0.01, Nitrospirae = 0.03, Other = 0.09,  Bacteroidetes = 0.1, Chloroflexi = 0.02,  Thaumarcheota = 0.02, Verrucomicrobia = 0.14) # Proportion of microbial groups in simulations. Must add to 1. 

cell_diameter = 0.75 # microbial cell diameter

prop_types = c(0.55, 0.4, 0.05) # proportions of cells in microcolonies, individual cells, filamentous  

# Microcolonies
microcolonies_size_range = c(3,8) # cell number range in microcolonies
microcolonies_diameter = 10 # microcolony diameter
microcolonies_distribution = "NeymanScott" # Spatial distribution for microcolonies
filamentous_size_range = c(3,8) # cell number range in filamentous colonies

for (s in seq_along(density)) { 
  
  print(paste("Image", s))
  set.seed(seed[s])

  # Create parameter list for make_community function
  parms = list(dimxy = dimxy,
               props = props,
               prop_types = prop_types,
               cell_density = density[s]/(100*100),
               cell_diameter = cell_diameter,
               microcolonies_size_range = microcolonies_size_range,
               microcolonies_diameter = microcolonies_diameter,
               microcolonies_distribution = microcolonies_distribution,
               filamentous_size_range = filamentous_size_range
  )
  
  # Create community
  cells = populate_species(make_community(parms), parms$props)
  
  # Make plot
  ## plot annotations
  length = 50
  posx = parms$dimxy -20 -length
  posy = 10
  zones = bind_rows(c(x=130, y=130, size = 100), c(x = 20, y = 20, size = 100))
  
  ## Plot colors
  colors = pals::kelly(length(parms$props)+1)[-1] #here we define the colours for each groups (the [-1] is because the first colour of this palette is white...)
  legend_labels = cbind(as.data.frame(table(marks(cells)$phylum)/cells$n), colour = colors[1:length(unique(marks(cells)$phylum))]) %>% arrange(desc(Freq)) %>% mutate(name = paste0(Var1, " (",round(100*Freq),"%)"))
  
  ## Main plot
  global = ggplot(as.data.frame(cells),aes(x = x, y=y, colour= phylum)) + 
    geom_point( size  = parms$cell_diameter) + 
    annotate("segment", x = posx, y = posy, xend = posx+length, yend = posy, linewidth = 2) + 
    annotate("text", x = posx + length/2, posy -5, label = paste(length, "µm")) + 
    annotate("rect", xmin = zones$x, ymin = zones$y, xmax = zones$x+zones$size, ymax = zones$y+zones$size, fill=NA, colour = "black")+
    scale_x_continuous(expand=c(0,0), limits= c(0,parms$dimxy)) +  
    scale_y_continuous(expand=c(0,0), limits= c(0,parms$dimxy)) +
    scale_color_manual(values = with(legend_labels, setNames(colour, Var1)), breaks = legend_labels$Var1, name = NULL, labels = with(legend_labels, setNames(name, Var1))) + 
    scale_alpha_discrete(range = c(1,1), guide = "none") + 
    guides(color = guide_legend(override.aes = list(size = 2), nrow = 3)) +
    coord_fixed() + 
    theme_void() + 
    theme(panel.background = element_rect(), legend.position = "bottom") +
    labs(tag = "A", title = paste(cells$n, "cells in 250x250 µm window"), subtitle = paste0(paste0(format(area.owin(discs(cells, cell_diameter/2))*100/area.owin(as.owin(cells)), digits = 2, nsmall=2),"% of HSA covered\n"), paste0(format(mean(rowSums(pairdist(cells) < 20)-1),digits = 1,nsmall =1), " neighbouring cells"))) +
    NULL

  ## Subplots
    subplots = list()
  for (z in 1:dim(zones)[1]) {
    X = cells[owin(c(zones[z,]$x, zones[z,]$x+ zones[z,]$size),c(zones[z,]$y, zones[z,]$y+ zones[z,]$size))]
    radius = 20
    cell_ref = X[X$x> min(X$x)+1.5*radius & X$x< max(X$x)-1.5*radius & X$y > min(X$y)+ 1.5*radius & X$y< max(X$y)-1.5*radius ][1]
    # Calculate the circle
    range = disc(centre = c(cell_ref$x,cell_ref$y), radius = radius)
    marks(X)$inrange = F
    marks(X[range])$inrange = T
    
    length = 20
    posx = zones[z,]$x+zones$size - 10 -length
    posy = zones[z,]$y + 10
    
    legend_labels = left_join( legend_labels[,c('Var1', 'colour')], as.data.frame(table(marks(X)$phylum)/X$n)) %>% mutate(name = paste0(Var1, " (",round(100*Freq),"%)"))
    
    subplots[[z]] = ggplot(as.data.frame(X),aes(x = x, y=y, colour= phylum)) + 
      geom_point(aes(alpha = inrange), size  = parms$cell_diameter) + 
      geom_circle(data = as.data.frame(cell_ref), aes(x0= x, y0 = y, r = radius), colour = "black", linetype = "dashed", size = 0.1)+
      annotate("segment", x = posx, y = posy, xend = posx+length, yend = posy, linewidth = 1) + 
      annotate("text", x = posx + length/2, posy -3, label = paste(length, "µm")) + 
      scale_x_continuous(expand=c(0,0)) +  
      scale_y_continuous(expand=c(0,0)) +
      scale_color_manual(values = with(legend_labels, setNames(colour, Var1)), breaks = legend_labels$Var1, name = NULL, labels = with(legend_labels, setNames(name, Var1))) + 
      scale_alpha_discrete(range = c(0.5,1), guide = "none") + 
      guides(color = guide_legend(override.aes = list(size = 2))) +
      coord_fixed() + 
      theme_void() + 
      theme(panel.background = element_rect(), legend.position = "right" ) +#
      labs(tag = c("B","C")[z], title = paste(X$n, "cells in 100x100 µm window"), subtitle = paste0(paste0(format(area.owin(discs(X, parms$cell_diameter/2))*100/area.owin(as.owin(X)), digits = 2, nsmall=2),"% of HSA covered\n"), paste0(format(mean(rowSums(pairdist(X) < 20)-1),digits = 2,nsmall =1)), " neighbouring cells")) +
      NULL
  }
  
  ## Make Figure
  left_panel = (subplots[[1]]+ theme(plot.margin = unit(c(10,0,20,0), "pt")))/ (subplots[[2]] + theme(plot.margin = unit(c(20,0,0,0), "pt")))
  (global + theme(plot.margin = unit(c(0,0,0,0), "pt"), plot.subtitle = element_text(margin = unit(c(0,0,-15,0), "pt")), plot.title = element_text(margin = unit(c(0,0,-15,0), "pt")))) + (left_panel + theme(plot.margin = unit(c(20,25,20,25), "pt")))+ plot_layout(widths = c(250, 100))
  ggsave(paste0(density[s],".pdf"), width = 11, height=8.27)
}

