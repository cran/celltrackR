## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(dpi=72)

## ----pack, warning = FALSE, message = FALSE---------------------------------------------------------------------------
library( celltrackR )
library( ggplot2 )

## ---- echo = FALSE----------------------------------------------------------------------------------------------------
# Save current par() settings
oldpar <- par( no.readonly =TRUE )

## ---------------------------------------------------------------------------------------------------------------------
str( TCells, list.len = 1 )

## ---------------------------------------------------------------------------------------------------------------------
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
load( system.file("extdata", "BCellsRaw.rda", package="celltrackR" ) )

## ----Tdata------------------------------------------------------------------------------------------------------------
str( TCellsRaw, list.len = 1 )

## ---------------------------------------------------------------------------------------------------------------------
head( TCellsRaw[[1]] )

## ---------------------------------------------------------------------------------------------------------------------
# for reproducibility
set.seed(2021)

# sample names, sort them in numeric order, and use them to subset the tracks
Bsample <- sample( names(BCellsRaw),50)
Bsample <- Bsample[ order(as.integer(Bsample)) ]
BCellsRaw <- BCellsRaw[ Bsample ]

# same for the T cells, but now we include a few specific IDs we'll need for the
# analysis below.
Tsample <- sample( names(TCellsRaw),44)
Tsample <- unique( c( Tsample, "5","9658","83","6080","7832","8352" ) )
Tsample <- Tsample[ order(as.integer(Tsample)) ]
TCellsRaw <- TCellsRaw[ Tsample ]

## ----tracklength-check------------------------------------------------------------------------------------------------
# Each track has a coordinate matrix with one row per coordinate;
# The number of steps is the number of rows minus one.
track.lengths <- sapply( TCellsRaw, nrow ) - 1
hist( track.lengths, xlab = "Track length (#steps)", breaks = seq(0,40 ) )
summary( track.lengths )

## ----max-tracklength--------------------------------------------------------------------------------------------------
# This is the number of coordinates, so the number of steps is one less.
maxTrackLength( TCellsRaw )

## ----filter-short-----------------------------------------------------------------------------------------------------
# nrow() of a track is always the number of steps plus one.
# For steps >= min.steps, we can substitute nrow > min.steps:
filter.criterion <- function( x, min.steps ){
  nrow(x) > min.steps
}

TCells.filtered <- filterTracks( filter.criterion, TCellsRaw, min.steps = 11 )
# Or shorthand: filterTracks( function(x) nrow(x) > min.steps, TCellsRaw )

## ----check-filtered---------------------------------------------------------------------------------------------------
# Find lengths of the filtered dataset
track.lengths.filtered <- sapply( TCells.filtered, nrow ) - 1

# Histograms of track lengths before and after
c( "min length before filter" = min( track.lengths ),
   "min length after filter" = min( track.lengths.filtered ) )

# Check how many tracks are left:
length( TCells.filtered )

## ---- fig.width=7, fig.height = 3.5-----------------------------------------------------------------------------------
# Drift is 0.05 micron/sec in each dimension
drift.speed <- c( 0.05, 0.05, 0.05 )
add.drift <- function( x, drift.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  drift.matrix <- matrix( rep( drift.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  drift.matrix <- drift.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords + drift.matrix
  return(x)
}

# Create data with drift
TCells.drift <- as.tracks( lapply( TCellsRaw, add.drift, drift.speed ) )

# Plot tracks and 'rose plot' with overlayed starting point
par(mfrow=c(1,2) )
plot(TCells.drift, main = "drift", cex = 0 )
plot( normalizeTracks(TCells.drift), main = "drift (rose plot)", cex = 0, col = "gray" )
abline( h = 0 )
abline( v = 0 )


## ---- fig.width = 4, fig.height = 3-----------------------------------------------------------------------------------
hotellingsTest( TCellsRaw, col = "gray" )

## ---- fig.width = 6, fig.height = 6-----------------------------------------------------------------------------------
hotellingsTest( BCellsRaw, col = "gray" )

## ---- fig.width = 4, fig.height = 3.5---------------------------------------------------------------------------------
# Compute autocovariance, normalize so the first point lies at 1
Tacov <- aggregate( TCellsRaw, overallDot )
Tacov$value <- Tacov$value / Tacov$value[1]

# The same for B cells:
Bacov <- aggregate( BCellsRaw, overallDot )
Bacov$value <- Bacov$value / Bacov$value[1]

# Compare autocovariances:
plot( Tacov )
points( Bacov, col = "red" )

## ---- fig.width = 6, fig.height = 6-----------------------------------------------------------------------------------
hotellingsTest( BCellsRaw, col = "gray", step.spacing = 10 )

## ---- fig.width = 4, fig.height = 4-----------------------------------------------------------------------------------
hotellingsTest( TCells.drift, plot = TRUE, col = "gray", step.spacing = 10 )

## ---- fig.width = 6, fig.height = 6-----------------------------------------------------------------------------------
hotellingsTest( TCellsRaw, dim = c("x","y","z"), step.spacing = 10 )
hotellingsTest( TCells.drift, dim = c("x","y","z"), step.spacing = 10 )

## ---- echo = FALSE, fig.width=7---------------------------------------------------------------------------------------
# sp <- seq(0,10)
# 
# htest.nodrift <- lapply( sp, function(x) hotellingsTest( TCells, 
#                                                  dim = c("x","y","z"),
#                                                  step.spacing = x ))
# htest.drift <- lapply( sp, function(x) hotellingsTest( TCells.drift, 
#                                                  dim = c("x","y","z"),
#                                                  step.spacing = x ))
# 
# pval.nodrift <- sapply( htest.nodrift, function(x) x$p.value )
# pval.drift <- sapply( htest.drift, function(x) x$p.value )
# stat.nodrift <- sapply( htest.nodrift, function(x) unname(x$statistic) )
# stat.drift <- sapply( htest.drift, function(x) unname(x$statistic) )
# 
# d.drift <- data.frame( step.spacing = sp,
#                        exp = "drift",
#                        pval = pval.drift,
#                        Tstatistic = stat.drift )
# d.nodrift <- data.frame( step.spacing = sp,
#                        exp = "no drift",
#                        pval = pval.nodrift,
#                        Tstatistic = stat.nodrift )
# d <- rbind( d.drift, d.nodrift )
# 
# p1 <- ggplot( d, aes( x = step.spacing, y = pval, color = exp ) ) +
#   geom_line() +
#   geom_hline( yintercept = 0.05, color = "red" ) +
#   scale_color_manual( values = c( "drift" = "black", "no drift" = "gray")  ) +
#   scale_y_log10() +
#   theme_classic()
# 
# p2 <- ggplot( d, aes( x = step.spacing, y = Tstatistic, color = exp ) ) +
#   geom_line() +
#   scale_y_continuous( limits = c(0,NA) ) +
#   scale_color_manual( values = c( "drift" = "black", "no drift" = "gray")  ) +
#   theme_classic()
# 
# gridExtra::grid.arrange(p1,p2,ncol=2)

## ----cellPairs, warning = FALSE, message = FALSE, fig.width=6, fig.heigth = 2.5---------------------------------------
# compute for both original as drift data
df.drift <- analyzeCellPairs( TCells.drift )
df.norm <- analyzeCellPairs( TCellsRaw )

# Plot
p.norm <- ggplot( df.norm, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray70", size = 0.5 ) +
  stat_smooth( span = 1, fill = "blue" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs",
        title = "original data") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

p.drift <- ggplot( df.drift, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray70", size = 0.5 ) +
  stat_smooth( span = 1, fill = "blue" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs",
        title = "data with drift") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

gridExtra::grid.arrange( p.norm, p.drift, ncol = 2 )

## ----stepPairs, warning = FALSE, message = FALSE, fig.width=6, fig.heigth = 2.5---------------------------------------
# compute for both original as drift data.
df.drift <- analyzeStepPairs( TCells.drift, filter.steps = function(x) displacement(x)>2  )
df.norm <- analyzeStepPairs( TCellsRaw, filter.steps = function(x) displacement(x)>2  )

# Plot
p.norm <- ggplot( df.norm, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray70", size = 0.5 ) +
  stat_smooth( span = 0.75, fill="blue" )+
  labs( x = "distance between step pairs",
        y = "angle between step pairs",
        title = "original data") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

p.drift <- ggplot( df.drift, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray70", size = 0.5 ) +
  stat_smooth( span = 0.75,  fill="blue")+
  labs( x = "distance between step pairs",
        y = "angle between step pairs",
        title = "data with drift") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

gridExtra::grid.arrange( p.norm, p.drift, ncol = 2 )

## ---------------------------------------------------------------------------------------------------------------------
# Get steps and find their displacement vectors
steps.drift <- subtracks( TCells.drift, 1 )
step.disp <- t( sapply( steps.drift, displacementVector ) )

# Get the mean
mean.displacement <- colMeans( step.disp )

# Divide this by the mean timestep to get a drift speed
drift.speed <- mean.displacement/timeStep( TCells.drift )
drift.speed

## ---- fig.width = 6, fig.height = 3.5---------------------------------------------------------------------------------
correct.drift <- function( x, drift.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  drift.matrix <- matrix( rep( drift.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  drift.matrix <- drift.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords - drift.matrix
  return(x)
}

# Create data with drift
TCells.corrected <- as.tracks( lapply( TCells.drift, correct.drift, drift.speed ) )

# Compare (zoom in on part of the field to see better)
par( mfrow = c(1,2) )
plot( TCellsRaw, col = "blue", main = "uncorrected", cex = 0, pch.start = NA, xlim = c(200,400), ylim = c(0,200) )
plot( TCells.drift, col = "red", add = TRUE, cex = 0, pch.start = NA )
plot( TCellsRaw, col = "blue", main = "corrected", cex = 0, pch.start = NA, xlim = c(200,400), ylim = c(0,200)  )
plot( TCells.corrected, col = "red", add = TRUE, cex = 0, pch.start = NA  )


## ---------------------------------------------------------------------------------------------------------------------
# Take the track with id "5"
dup.track <- TCellsRaw[["5"]]

# Add some noise to coordinates
dup.track[,"x"] <- dup.track[,"x"] + rnorm( nrow(dup.track), sd = 0.5 )
dup.track[,"y"] <- dup.track[,"y"] + rnorm( nrow(dup.track), sd = 0.5 )
dup.track[,"z"] <- dup.track[,"z"] + rnorm( nrow(dup.track), sd = 0.5 )

# Wrap the track in a tracks object and add it to the TCell data with
# a unique id number
dup.track <- wrapTrack( dup.track )
new.id <- max( as.numeric( names( TCellsRaw ) ) ) + 1
names(dup.track) <- as.character( new.id )
TCells.dup <- c( TCellsRaw, dup.track )

## ---- warning = FALSE, fig.width = 3, fig.height = 2.5----------------------------------------------------------------
df <- analyzeCellPairs( TCells.dup )

# label cellpairs that have both angle and distance below threshold
angle.thresh <- 5 # in degrees
dist.thresh <- 10 # this should be the expected cell radius
df$id <- paste0( df$cell1,"-",df$cell2 )
df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- "" 

# Plot; zoom in on the region with small angles and distances
ggplot( df, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  geom_text( aes( label = id ), color = "red", hjust = -0.1 ) +
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  coord_cartesian( xlim=c(0,30), ylim=c(0,20) ) +
  geom_hline( yintercept = angle.thresh, col = "blue",lty=2 ) +
  geom_vline( xintercept = dist.thresh, col = "blue", lty=2) +
  theme_classic()


## ---- fig.width = 7, fig.height = 2.5---------------------------------------------------------------------------------
par( mfrow=c(1,3))
plot( TCells.dup[c("5","9659") ] )
plot( TCells.dup[c("83","6080") ] )
plot( TCells.dup[c("7832","8352") ] )

## ---------------------------------------------------------------------------------------------------------------------
corrected.dup <- TCells.dup[ names( TCells.dup ) != "9659" ]

## ---------------------------------------------------------------------------------------------------------------------
tracks <- TCellsRaw
bb <- boundingBox( tracks )
bb

## ---- fig.width = 3, fig.height = 2-----------------------------------------------------------------------------------
# Define points:
lower1 <- c( bb["min","x"], bb["min","y"], bb["min","z"] )
lower2 <- c( bb["max","x"], bb["min","y"], bb["min","z"] )
lower3 <- c( bb["max","x"], bb["max","y"], bb["min","z"] )
zsize <- bb["max","z"] - bb["min","z"]

# Compute angles and distances of steps to this plane.
single.steps <- subtracks( tracks, 1 )
angles <- sapply( single.steps, angleToPlane, p1 = lower1, 
                  p2 = lower2, p3 = lower3 )
distances <- sapply( single.steps, distanceToPlane, p1 = lower1,
                     p2 = lower2, p3 = lower3 )
df <- data.frame( angles = angles,
                  distances = distances )

# Plot
ggplot( df, aes( x = distances,
                 y = angles ) ) +
  geom_point( color = "gray70", size = 0.5 ) +
  stat_smooth( method = "loess", span = 1, fill="blue" ) +
  geom_hline( yintercept = 32.7, color = "red" ) +
  scale_x_continuous( limits=c(0,zsize ), expand = c(0,0) ) +
  theme_classic()

## ---- fig.width = 4, fig.height = 3-----------------------------------------------------------------------------------
par( mar = c(5, 4, 0.5, 2) + 0.1 )
plot( TCellsRaw, dims = c("x","z" ) )

## ---------------------------------------------------------------------------------------------------------------------
boundingBox( TCellsRaw )

## ----get-steps--------------------------------------------------------------------------------------------------------
# Extract all subtracks of length 1 (that is, all "steps")
single.steps <- subtracks( TCellsRaw, 1 )

# The output is a new tracks object with a unique track for each step
# in the data (no longer grouped by the original cell they came from):
str( single.steps, list.len = 3 )

## ----check-avdt-------------------------------------------------------------------------------------------------------
median.dt <- timeStep( TCellsRaw )
median.dt

## ----check-dt---------------------------------------------------------------------------------------------------------
step.dt <- sapply( single.steps, duration )
str(step.dt)

## ----check-dt-hist----------------------------------------------------------------------------------------------------
dt.diff.perc <- (step.dt - median.dt) * 100 / median.dt
hist( dt.diff.perc, xlab = "dt (percentage difference from median)" )

## ---------------------------------------------------------------------------------------------------------------------
range( step.dt )
unique( step.dt )

## ---------------------------------------------------------------------------------------------------------------------
sum( step.dt == 48 )

## ---------------------------------------------------------------------------------------------------------------------
# This function randomly removes coordinates from a track dataset with probability "prob"
remove.points <- function( track, prob=0.1 ){
  
    tlength <- nrow( track )
    remove.rows <- sample( c(TRUE,FALSE), tlength, replace=TRUE,
                           prob = c(prob, (1-prob) ) )
    track <- track[!remove.rows,]
  return(track)
}

# Apply function to dataset to randomly remove coordinates in the data
TCells.gap <- as.tracks( lapply( TCellsRaw, remove.points ) )

## ---------------------------------------------------------------------------------------------------------------------
# duration of the individual steps
steps.gap <- subtracks( TCells.gap, 1 )

T1.step.disp <- sapply( single.steps, displacement )
T1.gap.disp <- sapply( steps.gap, displacement )

lapply( list( original = T1.step.disp, gaps = T1.gap.disp ), summary )

## ---------------------------------------------------------------------------------------------------------------------
T1.norm.disp <- sapply( single.steps, normalizeToDuration( displacement ) )
T1.norm.gap.disp <- sapply( steps.gap, normalizeToDuration( displacement ) )
lapply( list( original = T1.norm.disp, gaps = T1.norm.gap.disp ), summary )

## ---- fig.width=6-----------------------------------------------------------------------------------------------------

# Repair gaps by splitting or interpolation, the number of tracks is 
# different after each fix
split.gap <- repairGaps( TCells.gap, how = "split" )
interpolate.gap <- repairGaps( TCells.gap, how = "interpolate" )

c( "after splitting" = length( split.gap),
   "after interpolation" = length( interpolate.gap ) )

## ---------------------------------------------------------------------------------------------------------------------
T2 <- subsample( TCellsRaw, k = 2 )

## ---------------------------------------------------------------------------------------------------------------------
# displacement
T1.steps <- subtracks( TCellsRaw, 1 )
T1.disp <- sapply( T1.steps, displacement )

T2.steps <- subtracks( T2, 1 )
T2.disp <- sapply( T2.steps, displacement )

lapply( list( T1 = T1.disp, T2 = T2.disp ), summary )


## ---------------------------------------------------------------------------------------------------------------------
# interpolate both datasets at the time resolution of the neutrophils
dt <- timeStep( TCellsRaw )
interpolate.dt <- function( x, dt, how = "spline" ){
  trange <- range( timePoints( wrapTrack( x ) ) )
  tvec <- seq( trange[1], trange[2], by = dt )
  x <- interpolateTrack( x, tvec, how = how )
  return(x)
}

T1.corrected <- as.tracks( lapply( TCellsRaw, interpolate.dt, dt = dt ) )
T2.corrected <- as.tracks( lapply( T2, interpolate.dt, dt = dt ) )

# Check the effect on the displacement statistics:
T1.corr.steps <- subtracks( T1.corrected, 1 )
T1.corr.disp <- sapply( T1.corr.steps, displacement )

T2.corr.steps <- subtracks( T2.corrected, 1 )
T2.corr.disp <- sapply( T2.corr.steps, displacement )

lapply( list( T1 = T1.disp, T2 = T2.disp, 
              T1.corr = T1.corr.disp, T2.corr = T2.corr.disp ), 
        summary )

## ---- echo = FALSE----------------------------------------------------------------------------------------------------
# Reset par() settings
par(oldpar)

## ---- dpi = 100-------------------------------------------------------------------------------------------------------
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )

par( mar = c(2, 2, 1, 1) + 0.1 )
plot( TCellsRaw )

## ---------------------------------------------------------------------------------------------------------------------
cellSpeeds <- sapply( TCellsRaw, speed )
hist( cellSpeeds, breaks = 30 )

## ---------------------------------------------------------------------------------------------------------------------
cutoff <- 0.1 # 0.1 micron/sec = 6 micron/min
loSpeed <- filterTracks( function(t) speed(t) < cutoff, TCellsRaw )
hiSpeed <- filterTracks( function(t) speed(t) >= cutoff, TCellsRaw )

## ---- fig.width = 6, fig.height = 3.5---------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot( loSpeed, main = "Below cutoff" )
plot( hiSpeed, main = "Above cutoff" )

## ----bic-nonmotile----------------------------------------------------------------------------------------------------
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

## ---------------------------------------------------------------------------------------------------------------------
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

## ---- fig.width = 6, fig.height = 3.5---------------------------------------------------------------------------------
deltaBIC <- function( x, sigma ){
  b1 <- bicNonMotile( x, sigma )
  b2 <- bicMotile( x, sigma )
  d <- b1 - b2
  d
}

dBIC <- sapply( TCellsRaw, deltaBIC, 7 )
motile <- TCellsRaw[ dBIC >= 6 ]
nonMotile <- TCellsRaw[ dBIC < 6 ]


# Plot again for comparison:
par(mfrow=c(1,2))
plot( motile, main = "Motile" )
plot( nonMotile, main = "Non-Motile" )

## ---- fig.width = 4, fig.height = 3-----------------------------------------------------------------------------------
plot( dBIC, cellSpeeds, xlim = c(0,20) )
abline( v = 6, col = "red" ) # delta BIC cutoff
abline( h = max(cellSpeeds), col = "blue", lty = 2 ) # max speed of all cells

