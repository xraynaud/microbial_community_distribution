#' Add phylum marks to microbila community made with `make_community()`
#' @description This function adds a mark to each point: the group (phylum) to which it belong.
#' @details
#' Depending on the total number of cells and the proportion in each group, we may not obtain round numbers. If the total number of cells does not allow to have integer groups, the function samples the groups with given proportion. This means that realized proportions can deviate from the given proportions. This is especially the case when cell density is low.
#' In the case of colonies, cells within a colonie are assigned to the same phylum
#' @param X a simulated community produced by `make_community()`
#' @param props a named vector of proportions for the different phyla. Must add to 1.
#' @examples
#' parms = list(dimxy = 100,
#'              prop_types = c(40,50,10),
#'              cell_density = 0.0002,
#'              cell_diameter = 1,
#'              microcolonies_size_range = c(8,10),
#'              microcolonies_diameter = 10,
#'              microcolonies_distribution = "Matern",
#'              filamentous_size_range = c(8,10)
#'             )
#' populate_species(make_community(parms) , c(phylum1 = 0.20,
#'                                            phylum2 = 0.40, 
#'                                            phylum3 = 0.10,
#'                                            phylum4 = 0.25, 
#'                                            phylum5 = 0.05)
#'                                            )
#' @export
populate_species = function(X, props) {
  
  # We create a vector of names with the correct proportions  
  phylum = rep(names(props), spatstat.geom::npoints(X)*props) #This create a vector of names according to given proportions
  if (spatstat.geom::npoints(X) - length(phylum) >0) {
    phylum = factor(c(phylum, sample(names(props), size = spatstat.geom::npoints(X) - length(phylum), prob=props, replace=T)), levels = names(props)) # If the proportions*ncells does not sum up to the total number of cells we resample the proportions to have a complete set.
  }
  
  # clusters
  ncolonies = c(spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "microcolony",]$phylum, spatstat.geom::marks(X)[spatstat.geom::marks(X)$type == "filamentous",]$phylum)
  
  retry = T
  while(retry) {
    retry = F
    microcol_phylum = factor(rep(sample(names(props), prob=props, size = length(unique(ncolonies)), replace = T), table(ncolonies)), levels = names(props))
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