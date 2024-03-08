#' Simulate microbial community on a 2D space
#' @description The function `make_community` distributes points on a 2D space to simulate the distribution of microbial cells on pore surfaces. 
#' @param parameters contains the list of parameters to build the community. See `Details`. 
#' @details
#' The `parameters` list should contain:
#' \describe{
#' \item{`dimxy`}{the dimension of the square window} 
#' \item{`prop_types`}{the proportions of microcolonies, filamentous and isolated cells}
#' \item{`cell_diameter`}{the diameters of cells (identical for all groups)}
#' \item{`microcolonies_size_range`}{the range for cell numbers in microcolonies} 
#' \item{`microcolonies_diameter`}{the range for microcolony diameter}
#' \item{`microcolonies_distribution`}{the type of spatial distribution for microcolonies. Choose between *NeymanScott*, *Matern* or *Thomas*}
#' \item{`filamentous_size_range`}{the the range for cell numbers in filamentous colonies} 
#' }
#' @return a marked point process with type (isolated, microcolony, filamentous) and parent cluster id (column phylum)
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
#' make_community(parms)
#' @export
make_community = function(parameters) {
  
  cells = with(parameters,{
    cell_number = cell_density*dimxy^2
    microcolonies_number = cell_number*prop_types[1]/mean(microcolonies_size_range)
    
    microcolonies = switch(microcolonies_distribution,
           NeymanScott = {
            nclust <- function(x0, y0, radius, nmin, nmax) {
              off = spatstat.random::runifdisc(floor(runif(1, nmin, nmax+1)), radius, centre=c(x0, y0))
            }
            spatstat.random::rNeymanScott(microcolonies_number/dimxy^2,0.2, nclust, radius = microcolonies_diameter/2,nmin = microcolonies_size_range[1], nmax = microcolonies_size_range[2], win=spatstat.geom::owin(c(0,dimxy), c(0,dimxy)), saveparents = T)}, 
           Matern = spatstat.random::rMatClust(microcolonies_number/dimxy^2, microcolonies_diameter/2,mean(microcolonies_size_range), win=spatstat.geom::owin(c(0,dimxy), c(0,dimxy)), saveparents = T),
           Thomas = spatstat.random::rThomas(microcolonies_number/dimxy^2, microcolonies_diameter/4, mean(microcolonies_size_range), win=spatstat.geom::owin(c(0,dimxy), c(0,dimxy)), saveparents = T)
    )
    
    spatstat.geom::marks(microcolonies) = data.frame(type = "microcolony", phylum = attr(microcolonies,"parentid"))
    
    # Filamentous
    filamentous_number = round(cell_number*prop_types[3]/mean(filamentous_size_range))
    
    x_mothers = runif(filamentous_number, min =0, max = dimxy)
    y_mothers = runif(filamentous_number, min =0, max = dimxy)
    filamentous_alpha = runif(filamentous_number, min = 0, max = 2*pi)
    filamentous_curve = rnorm(filamentous_number, sd = 0.05)
    filamentous_size = sample(filamentous_size_range[1]:filamentous_size_range[2],filamentous_number, replace = T)
    filamentous_length = filamentous_size*cell_diameter
    filamentous_parent_id = rep(1:filamentous_number, filamentous_size)
    
    filamentous = dplyr::bind_rows(
      lapply(
        1:filamentous_number,
        function(x) lapply(
          0:filamentous_size[x]-1, 
          function(c) data.frame(
            x = x_mothers[x]+ c*cell_diameter*cos(filamentous_alpha[x]+c*2.5*filamentous_curve[x])+ rnorm(1, sd=0.1),
            y = y_mothers[x]+ c*2.5*cell_diameter*sin(filamentous_alpha[x]+c*filamentous_curve[x])+ rnorm(1, sd = 0.1), 
            parent_id = x))))
    
    filamentous = filamentous[filamentous$x >= 0 &filamentous$y >= 0 & filamentous$x <= dimxy &filamentous$y <= dimxy, ]
    
    filamentous = spatstat.geom::ppp(filamentous$x, filamentous$y, marks = data.frame(type = "filamentous", phylum = filamentous$parent_id), window = spatstat.geom::owin(c(0,dimxy), c(0,dimxy)), check = T)

    # isolated cells
    random = spatstat.random::rpoispp(cell_number*prop_types[2]/dimxy^2, win=spatstat.geom::owin(c(0,dimxy), c(0,dimxy)))
    spatstat.geom::marks(random) = data.frame(type = "individual", phylum =  NA)
    
    cells = spatstat.geom::superimpose(microcolonies, filamentous, random)
  }
  )
  spatstat.geom::marks(cells[spatstat.geom::marks(cells)$type == 'filamentous'])$phylum = spatstat.geom::marks(cells[spatstat.geom::marks(cells)$type == 'filamentous'])$phylum + max(spatstat.geom::marks(cells[spatstat.geom::marks(cells)$type == 'microcolony'])$phylum)
  return(cells)
}

