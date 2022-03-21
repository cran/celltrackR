## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(dpi=72)

## ----pack, warning = FALSE, message = FALSE---------------------------------------------------------------------------
library( celltrackR )
library( ggplot2 )

## ---------------------------------------------------------------------------------------------------------------------
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
load( system.file("extdata", "BCellsRaw.rda", package="celltrackR" ) )
load( system.file("extdata", "NeutrophilsRaw.rda", package="celltrackR" ) )

## ---------------------------------------------------------------------------------------------------------------------
# nrow on a track gives # coordinates; number of steps is this minus one
minStepsT <- min( sapply( TCellsRaw, nrow ) - 1 )
minStepsB <- min( sapply( BCellsRaw, nrow ) - 1 )
minStepsN <- min( sapply( NeutrophilsRaw, nrow ) - 1 )
c( "T cells" = minStepsT, "B cells" = minStepsB, "Neutrophils" = minStepsN )

## ---------------------------------------------------------------------------------------------------------------------
veryShort <- sum( sapply( NeutrophilsRaw, nrow ) < 4 )
100 * veryShort / length( NeutrophilsRaw )


## ---------------------------------------------------------------------------------------------------------------------
Neutrophils <- filterTracks( function(t) nrow(t) >= 4, NeutrophilsRaw )
TCells <- TCellsRaw
BCells <- BCellsRaw

## ----fig.width = 8, fig.height = 4------------------------------------------------------------------------------------
hotellingsTest( TCells, step.spacing = 10 )
hotellingsTest( BCells, step.spacing = 10 )
hotellingsTest( Neutrophils, step.spacing = 10 )

## ---- fig.width = 5, fig.height = 4-----------------------------------------------------------------------------------
par( mfrow=c(3,1), mar = c(0,0,0,0) + 0.1 )
plot( TCells, dims = c("x","z"), xaxt='n', yaxt = 'n', ann=FALSE )
plot( BCells, dims = c("x","z"), xaxt='n', yaxt = 'n', ann=FALSE )
plot( Neutrophils, dims = c("x","z"), xaxt='n', yaxt = 'n', ann=FALSE )

## ---------------------------------------------------------------------------------------------------------------------
# zoom in on border cells
plot( TCells, xlim = c(400, 420), ylim = c(250,350))

## ---------------------------------------------------------------------------------------------------------------------
# Checks angle of a cell's steps to the borders
# (bb is the bounding box of all tracks, used to define those borders)
# returns the fraction of a cell's steps that are aligned with 
# one of the borders
angleCheck <- function( steps, bb, thresholdAngle = 0.1 ){
	# only consider x and y borders since filtering on the z-border would
  # remove too many cells (we'll later project on the xy plane instead):
	minx <- bb["min","x"]
	miny <- bb["min","y"]
	
	angles <- matrix( 0, nrow = length(steps), ncol = 2 )
	angles[,1] <- sapply( steps, angleToPlane, 
	                      p1 = c(minx,0,1), p2=c(minx,1,0), p3 = c(minx,0,0) )
	angles[,2] <- sapply( steps, angleToPlane, 
	                      p1 = c(0,miny,1), p2=c(1,miny,0), p3 = c(0,miny,0) )

	minAng <- apply( angles, 1, min, na.rm = TRUE )
	maxAng <- apply( angles, 1, max, na.rm = TRUE )

	# Steps are suspect if they are at angle ~0 or ~180 to the border plane.
	return( sum( minAng < thresholdAngle | maxAng > (180-thresholdAngle) )/length(steps) )
}

# Checks distance of a cell's steps to the borders; returns the fraction of steps
# that are closer than a certain threshold to one of the borders.
distanceCheck <- function( steps, bb, threshold = 1 ){
	total <- numeric( length(steps) )
	for( d in c("x","y") ){
		# distance to the lower border
		minDist <- sapply( steps, function(x) min( abs( x[,d] - bb["min",d] ) ) )
		
		# distance to the higher border
		maxDist <- sapply( steps, function(x) min( abs( x[,d] - bb["max",d] ) ) )
		
		# suspect if one of these distances is below threshold
		nearBorder <- ( minDist < threshold ) | ( maxDist < threshold )		
		total[nearBorder] <- 1
	}
	return( sum(total)/length(total) )	
}

# Remove tracks that have more than maxFrac steps that are aligned with the border AND 
# within a certain distance to the border:
notAtBorder <- function( tracks, angleThreshold = 0.1, distanceThreshold = 1, maxFrac = 0.2 ){
	bb <- boundingBox( tracks )
	stepsByCell <- lapply(tracks, function(x){ subtracks(x,1) })
	atBorderAngle <- sapply( stepsByCell, angleCheck, bb, threshold = angleThreshold ) > maxFrac
	atBorderDistance <- sapply( stepsByCell, distanceCheck, bb, threshold = distanceThreshold ) > maxFrac
	atBorder <- atBorderAngle & atBorderDistance

	return( tracks[!atBorder] )
}


## ---- fig.width = 7, fig.height = 2.5---------------------------------------------------------------------------------
par( mfrow=c(1,3) )
old <- TCells
TCells <- notAtBorder( TCells )
TRemoved <- old[ !is.element( names(old), names(TCells) ) ]
plot( TRemoved, col = "red" )

old <- BCells
BCells <- notAtBorder( BCells )
BRemoved <- old[ !is.element( names(old), names(BCells) ) ]
plot( BRemoved, col = "red" )

old <- Neutrophils
Neutrophils <- notAtBorder( Neutrophils )
NRemoved <- old[ !is.element( names(old), names(Neutrophils) ) ]
plot( NRemoved, col = "red" )

# show how many removed:
c( paste0( "T cells : ", length( TRemoved), " of ", 
           length( TRemoved ) + length( TCells ), " tracks removed"),
   paste0( "B cells : ", length( BRemoved), " of ", 
           length( BRemoved ) + length( BCells ), " tracks removed"),
   paste0( "Neutrophils : ", length( NRemoved), " of ", 
           length( NRemoved ) + length( Neutrophils ), " tracks removed")
   )


## ---------------------------------------------------------------------------------------------------------------------
TCells <- projectDimensions( TCells, c("x","y") )
BCells <- projectDimensions( BCells, c("x","y") )
Neutrophils <- projectDimensions( Neutrophils, c("x","y") )

## ---------------------------------------------------------------------------------------------------------------------
bicNonMotile <- function( track, sigma ){
  
  # we'll use only x and y coordinates since we saw earlier that the z-dimension was
  # not so reliable
  allPoints <- track[,c("x","y")]
  
  # Compute the log likelihood under a multivariate gaussian.
  # For each point, we get the density under the Gaussian distribution 
  # (using dmvnorm from the mvtnorm package).
  # The product of these densities is then the likelihood; but since we need the 
  # log likelihood, we can also first log-transform and then sum:
  Lpoints <- mvtnorm::dmvnorm( allPoints, 
                      mean = colMeans(allPoints), # for a Gaussian around the mean position
                      sigma = sigma*diag(2), # sd of the Gaussian (which we should choose)
                      log = TRUE )
  logL <- sum( Lpoints )
  
  # BIC = k log n - 2 logL; here k = 3 ( mean x, mean y, sigma )
  return( 3*log(nrow(allPoints)) - 2*logL )
}
# the BIC for a given cutoff m
bicAtCutoff <- function( track, m, sigma ){
  
  # we'll use only x and y coordinates since we saw earlier that the z-dimension was
  # not so reliable
  allPoints <- track[,c("x","y")]
  
  # Split into two coordinate sets based on the cutoff m:
  firstCoords <- allPoints[1:m, , drop = FALSE]
  lastCoords <- allPoints[(m+1):nrow(allPoints), , drop = FALSE ]
  
  # Compute log likelihood under two separate Gaussians:
  Lpoints1 <- mvtnorm::dmvnorm( firstCoords, 
    mean = colMeans(firstCoords), 
    sigma = sigma*diag(2), 
    log = TRUE )
  Lpoints2 <- mvtnorm::dmvnorm( lastCoords, 
    mean = colMeans(lastCoords), 
    sigma = sigma*diag(2), 
    log = TRUE )
  logL <- sum( Lpoints1 ) + sum( Lpoints2 )
  
  # BIC = k log n - 2 logL; here k = 6 ( 2*mean x, 2*mean y, sigma, and m )
  return( 6*log(nrow(allPoints)) - 2*logL )
}

# We'll try all possible cutoffs m, and choose best model (minimal BIC)
# to compare to our non-motile "null hypothesis":
bicMotile <- function( track, sigma ){
  
  # cutoff anywhere from after the first two coordinates to 
  # before the last two (we want at least 2 points in each Gaussian,
  # to prevent fitting of a single point)
  cutoffOptions <- 2:(nrow(track)-2)
  
    min( sapply( cutoffOptions, function(m) bicAtCutoff(track,m,sigma) ) )
}

# Delta BIC between the two models
deltaBIC <- function( x, sigma ){
  b1 <- bicNonMotile( x, sigma )
  b2 <- bicMotile( x, sigma )
  d <- b1 - b2
  d
}


## ---- fig.width = 7, fig.height = 10----------------------------------------------------------------------------------
TCellsBIC <- sapply( TCells, deltaBIC, 7 )
BCellsBIC <- sapply( BCells, deltaBIC, 7 )
NeutrophilsBIC <- sapply( Neutrophils, deltaBIC, 7 )

# Keep only the motile cells; BIC > 6 means reasonable evidence for motility
TNonMotile <- TCells[ TCellsBIC < 6 ]
BNonMotile <- BCells[ BCellsBIC < 6 ]
NNonMotile <- Neutrophils[ NeutrophilsBIC < 6 ]

TCells <- TCells[ TCellsBIC >= 6 ]
BCells <- BCells[ BCellsBIC >= 6 ]
Neutrophils <- Neutrophils[ NeutrophilsBIC >= 6 ]

# Check how many removed:
c( paste0( "T cells : ", length( TNonMotile), " of ", 
           length( TNonMotile ) + length( TCells ), " tracks removed"),
   paste0( "B cells : ", length( BNonMotile), " of ", 
           length( BNonMotile ) + length( BCells ), " tracks removed"),
   paste0( "Neutrophils : ", length( NNonMotile), " of ", 
           length( NNonMotile ) + length( Neutrophils ), " tracks removed")
   )

# Plot for comparison:
par(mfrow=c(3,2))
plot( TCells, main = "T cells (motile)" )
plot( TNonMotile, main = "T cells (non-motile)" )
plot( BCells, main = "B cells (motile)" )
plot( BNonMotile, main = "B cells (non-motile)" )
plot( Neutrophils, main = "Neutrophils (motile)" )
plot( NNonMotile, main = "Neutrophils (non-motile)" )

## ----eval = FALSE-----------------------------------------------------------------------------------------------------
#  checkPotentialDoubles <- function( tracks, distanceThreshold = 10, angleThreshold = 10 ){
#    # na.omit because when cells do not share time points, their distance is NA.
#    pairs <- na.omit( analyzeCellPairs( tracks ) )
#    check <- pairs[ pairs$dist <= distanceThreshold & pairs$angle <= angleThreshold, ]
#  
#    # return if no pairs to check
#    if( nrow(check) == 0 ){
#      message("No suspicious pairs found!")
#      return(NULL)
#    }
#  
#    # Plot suspicious pairs; let user navigate with keystrokes:
#    oldpar <- par()
#    par( mfrow=c(2,2), mar=c(0, 0, 4, 0))
#    for( i in 1:nrow( check ) ) {
#      c1 <- pairs$cell1[i]; c2 <- pairs$cell2[i]
#      plot( tracks[c(c1,c2)], main = paste0( c1,"-",c2),axes=FALSE,
#          frame.plot=TRUE, xlab=NA, ylab=NA )
#      if( i %% 4 == 0 ) invisible(readline(prompt="Press [enter] to continue"))
#    }
#    par( oldpar )
#  
#    return(check)
#  }
#  
#  checkPotentialDoubles( TCells )
#  checkPotentialDoubles( BCells )
#  checkPotentialDoubles( Neutrophils )

## ---------------------------------------------------------------------------------------------------------------------
# Check median dt for all datasets:
all.data <- list( TCells = TCells, BCells = BCells, Neutrophils = Neutrophils )
lapply( all.data, timeStep )

## ---------------------------------------------------------------------------------------------------------------------
# Find durations of all steps in each dataset
step.dt <- lapply( all.data, function(x) {
  steps <- subtracks(x,1)
  sapply( steps, duration)
})

# Check the range of these durations:
range.dt <- lapply( step.dt, range )
range.dt

## ---------------------------------------------------------------------------------------------------------------------
percentage.missing <- lapply( step.dt, function(x) 100*sum( x != 24 ) / length(x) )
percentage.missing

## ---------------------------------------------------------------------------------------------------------------------
# Split tracks when there is a gap; after splitting, keep only tracks of at least length 4.
TCells <- repairGaps( TCells, how = "split", split.min.length = 4 )
BCells <- repairGaps( BCells, how = "split", split.min.length = 4 )
Neutrophils <- repairGaps( Neutrophils, how = "split", split.min.length = 4 )

# check that it has worked:
corrected.data <- list( TCells = TCells, BCells = BCells, Neutrophils = Neutrophils )
step.dt <- lapply( corrected.data, function(x) {
  steps <- subtracks(x,1)
  sapply( steps, duration)
})
lapply( step.dt, range )

