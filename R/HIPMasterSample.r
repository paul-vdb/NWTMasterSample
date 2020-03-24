#' @useDynLib NWTMasterSample
#' @importFrom Rcpp sourceCpp
NULL


#' This is an internal Master Sample function to assign
#' indices to points based on the discrete Halton Box overlay.
#'
#' @param input An sp or sf spatial points. Accepts either. Or even a data frame with X, Y names.
#' @param J Integer for number of Halton Boxes to make. If not set it defaults to 100m roughly.
#' @param bb Master Sample bounding box.
#' @param base Generally 2,3. If you change it read the literature.
#' @param seed Master Sample Seed for two bases.
#' @param s1 Permutation of x orderings (0, 1)
#' @param s2 Permutation of y orderings (0, 1, 2)
#'
#' @export
getIndividualBoxIndices <- function(input, J = NULL, bb, base = c(2,3), seed = c(587, 8750), s1 = 0:1, s2 = 0:2)
{
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	
	if(is.null(J)) J <- round(log((bb[,2] - bb[,1])/100)/log(base), 0)

	B <- prod(base^J)
	
	if(any(class(input) == "sf"))
	{
		dat <- data.table(st_coordinates(input))
	} else if(ncol(coordinates(input[1,])) == 2){
		# Assume sp input:
		dat <- data.table(coordinates(input))
		setnames(dat, c("X", "Y"))
	}else{
		dat <- data.table(input$X, input$Y)
		setnames(dat, c("X", "Y"))
	}
	Bx <- base[1]^J[1]
	By <- base[2]^J[2]
	# Shift points to the lower box coordinate.
	dat <- dat[,c("x", "y") := list((X - shift.bas[1])/scale.bas[1], (Y - shift.bas[2])/scale.bas[2])]
	# Put it in 0-1 terms, with an adjustment to deal with "zeros" and machine precision.
	dat <- dat[,c("Ax", "Ay") := list(floor((x + 2*.Machine$double.eps)*Bx), 
		floor((y + 2*.Machine$double.eps)*By))]
		
	# This code seems to be doing something wrong. Switched to a more explicit merge above instead...
	dat[, c("Ax.new", "Ay.new") := list(permutation_B(Ji = J[1], bi = 2, si = s1)[Ax+1], permutation_B(Ji = J[2], bi = 3, si = s2)[Ay+1])] # Add one since indices include 0

	haltonIndex <- SolveCongruence(cbind(dat$Ax.new, dat$Ay.new), base = c(2,3), J = J)

	# Adjust everything for the Master Sample random seed.
	a1 <- seed[1:2] %% base^J
	boxInit <- SolveCongruence(matrix(a1, ncol = 2), base = c(2,3), J = J)
	
	# Adjusted index:
	haltonIndex <- ifelse(haltonIndex < boxInit, B + (haltonIndex - boxInit), haltonIndex - boxInit)
	# Return the Halton Index for all "pts" in dat that are passed to this function.
	return(haltonIndex)
}

#' Find the cut points for HIP loop.
#'
#' @param n an integer for number of points
#' @param xi is the coordinates
#' @param index Unique point ID
#' @param slice Which base to cut up.
#'
#' @export
subset_dist <- function(n, xi, index, slice = 2)
{
	K <- n %% slice

	# If they don't have an index then randomly remove points...
	if(is.null(index))
	{
		index <- sample.int(n, n, replace = FALSE)
	}
	
	out <- vector("integer", length = n)

	if(K == 0){
		out <- (1:n > floor(n/slice))*1 + (1:n > floor(2*n/slice))*(slice == 3)
		return(out)
	}else if(K == 1)
	{
		max.i <- which.max(index)
		out[-max.i] <- (1:(n-1) > floor((n - 1)/slice))*1 + (1:(n - 1) > floor(2*(n - 1)/slice))*(slice == 3)
		out[max.i] <- -1
		return(out)
	}else if(K == 2){
		maxi.1 <- max(index)
		maxi.2 <- max(index[index != maxi.1])
		rm.max <- which(index %in% c(maxi.1, maxi.2))
		out[-(rm.max)] <- (1:(n - 2) > floor((n - 2)/slice))*1 + (1:(n - 2) > floor(2*(n - 2)/slice))*(slice == 3)
		out[rm.max] <- -1
	}
	return(out)
}

	
#' Permute every level of Js but with the same random 0,1,2.
#' @param Ji integer for power of the base
#' @param bi the co-prime base
#' @param si random permutation of 0,1,2
#'
#' @export
permutation_B <- function(Ji = 3, bi = 2, si = c(1, 0))
{
	if(Ji == 1) return(si)
	#Permute Base bi
	I <- si	
	for(k in 1:(Ji - 1))
	{
		v <- c()
		for(j in 1:bi^k)
		{
			v <- c(v, rep(I[j], bi) + bi^k*si)
		}
	I <- v
	}
	return(I)
}


#' This is the Main Master Sample function that processes the data and then calls HIP and orders based the MSI and per stratum orderings.
#' Only  mess with this if you are writing your own master sample. Otherwise, we will have it in a NWT or NZ wrapper function for set
#' Master Sample orderings etc. It works on centroids of any discrete Halton space (rivers and lakes or just boxes from continuous space).
#'
#' @param X Centroid coordinate for points be ordered x-axis.
#' @param Y Centroid coordinate for points be ordered y-axis.
#' @param index Master Sample Index for each point. If not supplied it will be calculated but slow things down and currently attemps to create 100 m boxes.
#' @param bb Master Sample bounding box.
#' @param base Generally 2,3. If you change it read the literature.
#' @param seed Master Sample Seed for two bases used to define the Master Sample Index.
#' @param quite if you don't want to show each iteration of HIP keep this TRUE
#' @param Ps1 Permutation for base 2 for calculating Master Sample Index.
#' @param Ps2 Permutation for base 3 for calculating Master Sample Index.
#' @param hipS1 Permutation of x orderings (0, 1) for this HIP (per stratum or individual sample)
#' @param hipS2 Permutation of y orderings (0, 1, 2) for this HIP (per stratum or individual sample)
#'
#' @export
getHipSample <- function(X, Y, index = NULL, N = NULL, bb,  base = c(2,3), seed = c(587, 8750), quiet = TRUE, Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = 0:2)
{
	if(is.null(N))
	{
		N <- 50
		print("Making 50 sites.")
	}
	
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	 
	dat <- data.table(X, Y)
	if(!is.null(index)) dat$index <- index

	# Shift points to the lower box coordinate.
	dat <- dat[,c("x", "y") := list((X - shift.bas[1])/scale.bas[1], (Y - shift.bas[2])/scale.bas[2])]
	
	dat[, ID := 1:.N]
	
	# Choosing boxes slightly as close to 100 m^2 as possible but maybe a
	# bit on the bigger side.
	J.index <- ceiling(log((bb[,2] - bb[,1])/100)/log(base))

	#Produces individual indices
	if(is.null(index)){
	# Put it in 0-1 terms, with an adjustment to deal with "zeros" and machine precision.
	Bx <- base[1]^J.index[1]
	By <- base[2]^J.index[2]

	dat <- dat[,c("Ix", "Iy") := list(floor((x + 2*.Machine$double.eps)*Bx), 
		floor((y + 2*.Machine$double.eps)*By))]

		# Permute the boxes for the new method of adding "extra" randomness.
		dat[, c("Ix.new", "Iy.new") := list(permutation_B(Ji = J[1], bi = 2, si = Ps1)[Ix+1], permutation_B(Ji = J[2], bi = 3, si = Ps2)[Iy+1])] # Add one since indices include 0

		haltonIndex <- SolveCongruence(cbind(dat$Ix.new, dat$Iy.new), base = c(2,3), J = J.index)

		dat[, index := haltonIndex]

		# Adjust everything for the Master Sample random seed.
		a1 <- seed %% base^J.index
		boxInit <- SolveCongruence(matrix(a1, ncol = 2), base = c(2,3), J = J.index)
		dat[, index := ifelse(index < boxInit, prod(base^J.index) + (index - boxInit), index - boxInit)]
	}
	
	dat[, c("Ax", "Ay") := list(0, 0)]
	
	# Before we do HIP I want to subset for all boxes that might have more than one observation...
	# The reason to do this is that we are considering the box itself as the sample unit to "partition" and not
	# the actual point. We can add it back later...
	dat.rm <- dat[duplicated(index)]
	dat <- dat[!duplicated(index)]
	
	# Partition options for 
	J.opts <- data.table(J1 = c(1,2,3,2,3,4,3,5,4,5,6,5,7,6,5,7,8,7,8) , J2 = c(1,1,1,2,2,2,3,2,3,3,3,4,3,4,5,4,4,5,5) )
	J.opts[, "B" := 2^J1*3^J2]
	Jmax <- J.opts[B >= N, ][which.min(B),]
	
	dat <- doHIP(dat, n = N, quiet = quiet, base = base, Jmax = Jmax)	# This is the iterative while loop now in a pretty function :)

	J <- attributes(dat)$J
	B <- prod(base^J)
	
	# HIP specific permutatations that are random each time
	# s1 <- sample(0:1, 2, replace = FALSE)
	# s2 <- sample(0:2, 3, replace = FALSE)
	s1 <- hipS1
	s2 <- hipS2
	
	# Add permutation from HIP Paper to add more "randomness" to the seed selection.
	dat[, c("Ax.new", "Ay.new") := list(permutation_B(Ji = J[1], bi = 2, si = s1)[Ax+1], permutation_B(Ji = J[2], bi = 3, si = s2)[Ay+1])] # Add one since indices include 0
	
	# Use the new permutation to solve the congruence.
	dat[, "HIPIndex" := SolveCongruence(cbind(Ax.new, Ay.new), base = c(2,3), J = J)]
	
	#Find initial box.
	# boxSeed <- seed %% base^J
	# boxSeed <- SolveCongruence(matrix(boxSeed, ncol = 2, nrow = 1), base = c(2,3), J = J)
	boxSeed <- dat[which.min(index),]$HIPIndex
	dat[, "HIPIndex" := ifelse(HIPIndex < boxSeed, B + (HIPIndex - boxSeed), HIPIndex - boxSeed)]
			
	tmp <- dat[, c("index", colnames(dat)[!names(dat) %in% names(dat.rm)]), with = FALSE]
	
	dat.tmp <- merge(dat.rm, tmp, by = "index", all.x = TRUE, all.y = FALSE)
	dat <- rbind(dat, dat.tmp)
	setorder(dat, ID)
	setattr(dat, "J", J)
	setattr(dat, "B", B)
	setattr(dat, "Permutation1", s1)
	setattr(dat, "Permutation2", s2)
	
	# Return the Halton Index for all "pts" in dat that are passed to this function.
	return(dat)
}

#' Run the loop for Halton Iterative Partitioning
#' @param dat2 internal data table passed to the function.
#' @param n Sample size required
#' @param quiet whether or not to be print on each iteration.
#' @param base Co-prime bases 2,3
#' @param Maximum number of cuts by x,y J values.
#'
#' @export
doHIP <- function(dat2, n = NULL, quiet = TRUE, base = base, Jmax)
{
	if(is.null(n)) n = nrow(dat2)*2

	grp.size <- dat2[,.N, by = c("Ax", "Ay")]

	i <- 0
	j <- 0
	while(all(grp.size$N > 1) & prod(base^c(i,j)) < n)
	{
		if(all(grp.size$N > 1) & prod(base^c(i,j)) < n & (i + 1 <= Jmax$J1)){
		setorder(dat2, x)
		
		dat2[, "Ax2" := subset_dist(.N, x, index, 2), by = c("Ax","Ay")]
		dat2 <- dat2[Ax2 != -1]		
		dat2[, "Ax" := Ax*2 + Ax2]

		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		i <- i + 1
		}
		if(all(grp.size$N > 2) & ((base[1]^i > base[2]^j) | (i >= Jmax$J1)) & (prod(base^c(i,j)) < n) & (j + 1 <= Jmax$J2)){	# Don't divide by 3 every time so we can achieve reasonably square divisions. Need to play with this.
			setorder(dat2, y)
			dat2[, "Ay2" := subset_dist(.N, y, index, 3), by = c("Ax","Ay")]
			dat2 <- dat2[Ay2 != -1]
			dat2[, "Ay" := Ay*3 + Ay2]			
			j <- j + 1
		}
		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		
		if(quiet == FALSE){
			print(paste0("J1 = ", i, " and J2 = ", j))
			print(min(grp.size$N))
		}
	}

	J <- c(i,j)
	setattr(dat2, "J", J)
	return(dat2)
}

#' Select a sample given a HIP exported data.table.
#'
#' @param dat HIP Exported data.table. Requires certain columns.
#' @param n sample size
#' @param lks If lakes then change how spatial points are processed.
#' @param proj.string Projection if needing to turn points into spatial points data frame for export.
#'
#' @export
getSamples <- function(dat, n, lks = FALSE, proj.string = NULL)
{

	B <- max(dat$HIPIndex) + 1 
	dat[ , samplei := (as.integer(factor(index)) - 1)*B + HIPIndex, by = "HIPIndex"]
	dat <- dat[!duplicated(index),]
	smp <- dat[samplei < n, ]
	if(lks == FALSE) return(smp)
	smp.pts <- SpatialPointsDataFrame(SpatialPoints(cbind(smp$X, smp$Y), proj4string = proj.string), data = smp)
	return(smp.pts)
}

#' I honestly can't remember why I wrote this. I think it's an older version of getSamples. HB must stand for something afterall...
#'
#' @param dat HIP Exported data.table. Requires certain columns.
#' @param n sample size
#' @param proj.string Projection if needing to turn points into spatial points data frame for export.
#'
#' @export
getHBSamples <- function(dat, n, proj.string)
{
	dat <- dat[!duplicated(index),]
	dat[ , samplei := (as.integer(factor(index)) -1)]
	smp <- dat[samplei < n, ]
	smp.pts <- SpatialPointsDataFrame(SpatialPoints(cbind(smp$X, smp$Y), proj4string = proj.string), data = smp)
	return(smp.pts)
}


#' Written for Ollie to get over samples within a HIP box.
#'
#' @param dat HIP Exported data.table. Requires certain columns.
#' @param n sample size
#' @param ni Number of oversamples per HIP Box
#'
#' @export
getOverSamples <- function(dat, n, ni)
{

	B <- max(dat$HIPIndex) + 1 
	dat[ , sampleRank := (rank(index) - 1), by = "HIPIndex"]
	return(dat[HIPIndex < n & sampleRank < ni, ])
}

#' Do HIP but track the split points for making boxes.
#'
#' @param subdat HIP Exported data.table. Requires certain columns.
#' @param boxes passing the existing boxes to partition.
#' @param slice Which co-prime base are we partitioning meow.
#'
#' @export
subset_boxes <- function(subdat, boxes, slice = 2)
{
	bb.new <- list()
	for(i in sort(unique(subdat$Ax)))
	{
		for(j in sort(unique(subdat$Ay)))
		{
			bb.i <- boxes[[paste(c(i, j), collapse = " ")]]
			tmp <- subdat[Ax == i & Ay == j]
			tmp <- setorder(tmp, -index)			
			if(slice == 2)
			{
				miss <- nrow(tmp) %% 2

				if(miss >= 1 ) tmp <- tmp[-(1:miss)]
				md <- median(tmp$x)
				
				bb1 <- bb.i
				bb2 <- bb.i
				bb1[1,2] <- md
				bb2[1,1] <- md
				bb.new[[paste(i*2, j, collapse = " ")]] <- bb1
				bb.new[[paste(i*2 + 1, j, collapse = " ")]] <- bb2
			}
			if(slice == 3)
			{
				miss <- nrow(tmp) %% 3
				if(miss >= 1 ) tmp <- tmp[-(1:miss)]
				md <- quantile(tmp$y, c(1/3, 2/3))
				
				bb1 <- bb.i
				bb2 <- bb.i
				bb3 <- bb.i
				
				bb1[2,2] <- md[1]
				bb2[2,1] <- md[1]
				bb2[2,2] <- md[2]
				bb3[2,1] <- md[2]

				bb.new[[paste(i, j*3, collapse = " ")]] <- bb1
				bb.new[[paste(i, j*3 + 1, collapse = " ")]] <- bb2
				bb.new[[paste(i, j*3 + 2, collapse = " ")]] <- bb3
				
			}
		}
	}
	return(bb.new)
}

#' This function returns a bunch of polygons that represent the HIP splits of the population.
#' It really just does HIP again and could be combined with doHIP in the future with just the "boxes" 
#' added and create Polygon = TRUE or something. For now it's separate since it's plenty fast.
#'
#' @param hip HIP Exported data.table. Requires certain columns.
#' @param bb bounding box for this region.
#' @param n Number of HIP boxes to make.
#'
#' @export
getHIPBoxes <- function(hip, bb, n)
{
	# Partition options for 
	J.opts <- data.table(J1 = c(1,2,3,2,3,4,3,5,4,5,6,5,7,6,5,7,8,7,8) , J2 = c(1,1,1,2,2,2,3,2,3,3,3,4,3,4,5,4,4,5,5) )
	J.opts[, "B" := 2^J1*3^J2]
	Jmax <- J.opts[B >= n, ][which.min(B),]

	boxes <- list("0 0" = bb)

	dat2 <- hip[,.(x = X, y = Y, index, Ax = 0, Ay = 0)]
	
	grp.size <- dat2[,.N, by = c("Ax", "Ay")]

	i <- 0
	j <- 0
	while(all(grp.size$N > 1) & length(boxes) < n)
	{
		if(all(grp.size$N > 1) & prod(base^c(i,j)) < n & (i+1 <= Jmax$J1)){
		setorder(dat2, x)
		boxes <- subset_boxes(dat2, boxes, slice = 2)
		
		dat2[, "Ax2" := subset_dist(.N, x, index, 2), by = c("Ax","Ay")]
		dat2 <- dat2[Ax2 != -1]	
		dat2[, "Ax" := Ax*2 + Ax2]

		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		i <- i + 1
		}
		if(all(grp.size$N > 2) & ((base[1]^i > base[2]^j) | (i >= Jmax$J1)) & length(boxes) < n & (j+1 <= Jmax$J2)){	# Don't divide by 3 every time so we can achieve reasonably square divisions. Need to play with this.		
			setorder(dat2, y)
			boxes <- subset_boxes(dat2, boxes, slice = 3)

			dat2[, "Ay2" := subset_dist(.N, y, index, 3), by = c("Ax","Ay")]
			dat2 <- dat2[Ay2 != -1]
			dat2[, "Ay" := Ay*3 + Ay2]			
			j <- j + 1
		}
		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		
	}

	polys <- lapply(1:length(boxes), FUN = function(i){as(extent(as.matrix(boxes[[i]])), 'SpatialPolygons')})
	for(i in 1:length(polys))
	{
		polys[[i]] <- polys[[i]]@polygons[[1]]
		slot(polys[[i]], "ID") <- names(boxes)[i]
	}
	
	polys <- SpatialPolygonsDataFrame(SpatialPolygons(polys), data.frame(ID = names(boxes), row.names = names(boxes)))
	return(polys)
}

#' Create the Sample Raster:
#' This is for linear features. You pass the lines as a lines feature as well as the spatial sample smp. Then match by 
#' what was clipped in the Halton Boxes. Then draw a Halton sample along the lines. For continuous sampling inside a Halton box on a line.
#'
#' @param shp Linear features for rivers or some sort of lines.
#' @param smp The chosen samples as data.table.
#' @param bb Master Sample bounding box.
#' @param J Master Sample cuts to figure out what the raster is.
#' @param base co-prime base. Almost always 2,3. Might really bugger things up if it weren't.
#'
#' @export
getPolyRast <- function(shp, smp, bb, J = c(5,3), base = c(2,3))
{
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	box.size <- scale.bas/base^J

	pts <- list()
	for(i in 1:nrow(smp))
	{
		# Find the box, clip the shape and get the dang sample!
		x <- c(smp[i]$X - box.size[1]/2, smp[i]$X + box.size[1]/2)
		y <- c(smp[i]$Y - box.size[2]/2, smp[i]$Y + box.size[2]/2)		
		b_poly <- as(extent(c(x,y)), "SpatialPolygons")
		if(any(class(shp) %in% "sf"))
		{
			rv <- st_intersection(shp, st_set_crs(st_as_sf(b_poly), st_crs(shp)))	
			rv <- as_Spatial(rv)
		}else{
			proj4string(b_poly) <- proj4string(shp)
			rv <- gIntersection(shp, b_poly, byid = T)
		}
		pts[[i]] <- sampRiv(n = 1, x = rv, seed = smp[i]$Ax.new)
	}
	pts <- do.call("rbind", pts)
	pts <- SpatialPointsDataFrame(pts, data = smp)
	return(pts)
}

#' Sample along a line using BAS.
#'
#' @param n Sample size
#' @param x Linear shape file clipped by a raster.
#' @param seed Random starting seed in one dimension for BAS linear features.
#'
#' @export
sampRiv <- function(n = 10, x, seed = 0)
{
	cc <- coordinates(x)
	cc <- do.call("rbind", cc)
	cc.mat <- as.matrix(do.call("rbind", cc))
	lengths = LineLength(cc.mat, longlat = FALSE, sum = FALSE)
	csl = c(0, cumsum(lengths))
	maxl = csl[length(csl)]

	pts = HaltonSeq(seed, 2, n)* maxl
	int = findInterval(pts, csl, all.inside = TRUE)
	where = (pts - csl[int])/diff(csl)[int]
	xy = cc.mat[int, , drop = FALSE] + where * (cc.mat[int + 1, , drop = FALSE] - cc.mat[int, , drop = FALSE])
	SpatialPoints(xy, proj4string = CRS(proj4string(x)))
}

#' Sample along a line using BAS.
#'
#' @param hip HIP exported data.table
#' @param n Number of samples
#' @param smp Spatial Points Sample, assumed sp.
#'
#' @export
getBalance <- function(hip, n, smp = NULL)
{
	if(is.null(smp)){
		smp <- getSamples(hip, n)
		smp <- SpatialPoints(cbind(smp$X, smp$Y))
	}	
	pts <- SpatialPoints(cbind(hip$X, hip$Y)) 
	vp <- getVoronoi(smp, bbox(pts))
	proj4string(pts) <- proj4string(vp)
	v <- var(colSums(gContains(vp, pts, byid = TRUE)))
	return(v)
}

# Bunch of stuff specific to NZ Freshwater Master Sample.
# Ignore it for now.	
# getFWSeed <- function(island = "South")
# {
	# if(island == "South"){
		# seed <-  c(7815, 699)
		# s1 <- c(0,1)
		# s2 <- c(2,0,1)		
	# }
	# if(island == "North"){
		# seed <- c(601, 5024)
		# s1 <- c(0,1)
		# s2 <- c(2,1,0)
	# }
	# return(list(seed = seed, s1 = s1, s2 = s2))
# }

# getStratumPermutation <- function(island = "South", stratum = "1-2")
# {
	# if(island == "South"){
		# if(stratum == "1-2") return(list(s1 = c(1,0), s2 = c(2,1,0)))
		# if(stratum == "3-4") return(list(s1 = c(1,0), s2 = c(0,2,1)))
		# if(stratum == "5-8") return(list(s1 = c(1,0), s2 = c(2,0,1)))
	# }
	# if(island == "North"){
		# if(stratum == "1-2") return(list(s1 = c(0,1), s2 = c(1,0,2)))
		# if(stratum == "3-4") return(list(s1 = c(0,1), s2 = c(0,2,1)))
		# if(stratum == "5-8") return(list(s1 = c(0,1), s2 = c(0,1,2)))	
	# }
	# return("ERROR")
# }

#' This is a modified SDraw voronoi polygons code for the specific purposes of testing Spatial Balance with
#' objects typically built in this package.
#'
#' @param x is an sp spatial points
#' @param bbox the bounding box for the whole region.
#'
#' @export
getVoronoi <- function (x, bbox) 
{
    if (!inherits(x, "SpatialPoints")) {
        stop("Must pass a SpatialPoints* object to voronoi.polygons.")
    }
    crds = coordinates(x)
    z = deldir::deldir(crds[, 1], crds[, 2], rw = bbox)
    w = deldir::tile.list(z)
    polys = vector(mode = "list", length = length(w))
    for (i in seq(along = polys)) {
        pcrds = cbind(w[[i]]$x, w[[i]]$y)
        pcrds = rbind(pcrds, pcrds[1, ])
        polys[[i]] = Polygons(list(Polygon(pcrds)), ID = as.character(i))
    }
    SP = SpatialPolygons(polys, proj4string = CRS(proj4string(x)))
    voronoi = SpatialPolygonsDataFrame(SP, data = data.frame(x = crds[, 
        1], y = crds[, 2], area = sapply(slot(SP, "polygons"), 
        slot, "area"), row.names = sapply(slot(SP, "polygons"), 
        slot, "ID")))
    return(voronoi)
}

