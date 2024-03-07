populate_species = function(X, props) {
  # This function adds a mark to each point: the group (phylum) to which it belong .
  # Depending on the total number of cells and the proportion in each group, we may not obtain round numbers. The if is here to deal with this problem. If  cell_number*props do not give integers, we SAMPLE the groups with given proportion. THIS MEANS THAT THE PROPORTION IN THIS SAMPLE CAN DEVIATES FROM YOUR PROPORTIONS (especially when cell density is low)
  # cells within a colonie are assigned to the same phyla
  
  # We create a vector of names with the correct proportions  
  phylum = rep(names(props), spatstat.geom::npoints(X)*props) #This create a vector of names according to given proportions
  if (spatstat.geom::npoints(X) - length(phylum) >0) {
    phylum = factor(c(phylum, sample(names(props), size = spatstat.geom::npoints(X) - length(phylum), p=props, replace=T)), levels = names(props)) # If the proportions*ncells does not sum up to the total number of cells we resample the proportions to have a complete set.
  }
  
  # clusters
  ncolonies = c(spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "microcolony",]$phylum, spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "filamentous",]$phylum)
  
  retry = T
  while(retry) {
    retry = F
    microcol_phylum = factor(rep(sample(names(props), p=props, size = length(unique(ncolonies)), replace = T), table(ncolonies)), levels = names(props))
    if (any(table(phylum)-table(microcol_phylum)<0) ) {
      retry = T
    }
  }
  
  spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "microcolony",]$phylum = as.character(microcol_phylum[1:spatstat.geom::npoints(X[spatstat.geom::marks(X)$type == "microcolony",])])
  spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "filamentous",]$phylum = as.character(microcol_phylum[spatstat.geom::npoints(X[spatstat.geom::marks(X)$type == "microcolony",])+1:spatstat.geom::npoints(X[spatstat.geom::marks(X)$type == "filamentous",])])
  
  #  phylum = unlist(lapply(union(phylum, microcol_phylum), function(i) 
  #    rep(i, length(phylum[phylum == i]) - length(microcol_phylum[microcol_phylum == i]))))
  
  phylum = unlist(lapply(names(props), function(p) rep(p, length(phylum[phylum==p])-length(microcol_phylum[microcol_phylum==p]))))
  
  #individual cells
  spatstat.geom::marks(X[spatstat.geom::marks(X)$type == "individual"])$phylum = phylum[1:spatstat.geom::npoints(X[spatstat.geom::marks(X)$type == "individual"])]
  
  return(X)
}