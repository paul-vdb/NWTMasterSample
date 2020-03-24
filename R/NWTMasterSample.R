#===========================================
# Code specific to NWT Master Sample:
# Currently this accepts lakes but it is 
# minimal effort to expand it to rivers.
#===========================================

#' @import leaflet
#' @import sp
#' @import data.table
NULL

#' Function to create Master Sample from points.
#' 
#' @param pts SpatialPointsDataFrame or SpatialPoints for pts that are to be sampled.
#' @param B Number of HIP Boxes to create for partitioning
#' @param P1 Permutation for base 2 (0,1)
#' @param P2 Permutation for base 3 (0,1,2)
#' @param n Sample size. Keep NULL if you want to return a shapefile with all points given the Master Sample Ordering.
#'
#'
#' @export
MasterSample <- function(pts, B = 12, P1 = 0:1, P2 = 0:2, n = NULL)
{
	# Set up the Master Sample components
	# What is the NWT MS Bounding Box:
	bb <- data.frame(min = c(-1026000, 8121500), max = c(558000, 9375000), row.names = c("x", "y"))
	# Random Seed
	u <- c(5857, 8700)
	# Projection
	gnwt_proj <- CRS("+init=epsg:3580")

	# Shift things to unit box.
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]

	# Currently assuming 100 m boxes.
	J <- round(log((bb[,2] - bb[,1])/c(100,100))/log(c(2,3)),0)	# Roughly 100 m grid
	base <- c(2,3)
	Bx <- base[1]^J[1]
	By <- base[2]^J[2]
	
	old.proj <- proj4string(pts)
	
	# Ensure the points are in the correct projection.
	pts <- spTransform(pts, gnwt_proj)
	XY <- data.frame(coordinates(pts))
	names(XY) <- c("X", "Y")
	# These permutations and seeds are fixed for the master sample: DO NOT CHANGE INPUT TO THIS FUNCTION!!!
	# Unless of course you are writing a master sample for a  new location such as BC :)
	MSI <- getIndividualBoxIndices(XY, J = J, bb = bb, base = c(2,3), seed = u, s1 = 1:0, s2 = c(0,2,1))
	XY$index <- MSI

	hip <- getHipSample(X = XY$X, Y = XY$Y, index = XY$index, N = B, bb = bb, seed = u, hipS1 = P1, hipS2 = P2)
	Bi <- max(hip$HIPIndex) + 1 
	hip[ , SampleIndex := (as.integer(factor(index)) - 1)*Bi + HIPIndex, by = "HIPIndex"]
	hip <- hip[!duplicated(index),]
	hip <- hip[,.(X,Y, MasterSampleIndex = index, HIPOrder = HIPIndex, SampleIndex)]
	smp.pts <- SpatialPointsDataFrame(SpatialPoints(cbind(hip$X, hip$Y), proj4string = gnwt_proj), data = hip[,.(MasterSampleIndex, HIPOrder, SampleIndex)])	
	pts.wgs <- spTransform(smp.pts, CRS(old.proj))

	if(is.null(n)) {
		return(smp.pts)
	}else{
		return(smp.pts[smp.pts$SampleIndex < n, ])
	}
}

#' Plots a hip dataframe for sample Index into a leaflet map for easy visualization.
#'
#' @param pts Spatial Points Data Fram exported from the master sample function.
#' @param n Number of points you want to plot.
#'
#' @export
plotMS <- function(pts, n)
{
	pts.wgs <- spTransform(pts, CRS("+proj=longlat"))
	pts.wgs <- pts.wgs[pts.wgs$SampleIndex < n, ]
	leaflet(pts.wgs) %>% addProviderTiles(providers$Esri.WorldImagery) %>%  
			addMarkers(popup = ~paste0("SampleIndex: ", SampleIndex, "<br> HIPOrder: ", HIPOrder, 
						"<br> Master Sample Index: ", MasterSampleIndex))
}
