## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi=72)

## ----pack, warning = FALSE, message = FALSE-----------------------------------
library( celltrackR )
library( ggplot2 )

## ---- echo = FALSE------------------------------------------------------------
# Save current par() settings
oldpar <- par( no.readonly =TRUE )

## ----Tdata--------------------------------------------------------------------
str( TCells, list.len = 1 )

## -----------------------------------------------------------------------------
head( TCells[[1]] )

## ----bdata--------------------------------------------------------------------
str( BCells, list.len = 1 )
str( Neutrophils, list.len = 1 )

## -----------------------------------------------------------------------------
# for reproducibility
set.seed(2021)

# sample names, sort them in numeric order, and use them to subset the tracks
Bsample <- sample( names(BCells),50)
BCells2 <- BCells[ Bsample ]

# same for the T cells
Tsample <- sample( names(TCells),50)
TCells2 <- TCells[ Tsample ]

# and the neutrophils
Nsample <- sample( names(Neutrophils),50)
Neutrophils2 <- Neutrophils[ Nsample ]


## ---- fig.width = 6, fig.height=2---------------------------------------------
par( mfrow = c(1,3), mar = c(2,2,3,1)+0.1 )
plot( TCells2, main = "T cells" )
plot( BCells2, main = "B cells" )
plot( Neutrophils2, main = "Neutrophils" )

## ---- fig.width = 8, fig.height=2.5-------------------------------------------
# load original, 3D tracks as an example, since the processed tracks are 2D
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
load( system.file("extdata", "BCellsRaw.rda", package="celltrackR" ) )
load( system.file("extdata", "NeutrophilsRaw.rda", package="celltrackR" ) )

# omit axes and labels to save space here;
par( mfrow = c(1,3), mar = c(0,0,2,0) + 0.1 )
plot3d( TCellsRaw, main = "T cells", tick.marks = FALSE )
plot3d( BCellsRaw, main = "B cells", tick.marks = FALSE )
plot3d( NeutrophilsRaw, main = "Neutrophils", tick.marks = FALSE )

## ---- fig.width = 6, fig.height=2---------------------------------------------
par( mfrow = c(1,3), mar = c(2,2,3,1)+0.1 )
plot( normalizeTracks(TCells2), main = "T cells" )
plot( normalizeTracks(BCells2), main = "B cells" )
plot( normalizeTracks(Neutrophils2), main = "Neutrophils" )

## -----------------------------------------------------------------------------
# Obtain mean speeds for each track using sapply
Tcell.speeds <- sapply( TCells2, speed )
str( Tcell.speeds )
summary( Tcell.speeds )

## ----eval = FALSE-------------------------------------------------------------
#  hist( Tcell.speeds )

## -----------------------------------------------------------------------------
Bcell.speeds <- sapply( BCells2, speed )
Nphil.speeds <- sapply( Neutrophils2, speed )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", speed = Tcell.speeds )
dB <- data.frame( cells = "B cells", speed = Bcell.speeds )
dN <- data.frame( cells = "Neutrophils", speed = Nphil.speeds )
d <- rbind( dT, dB, dN )

# Compare:
ggplot( d, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom() +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

## -----------------------------------------------------------------------------
Tcell.steps <- subtracks( TCells2, i = 1 )

## -----------------------------------------------------------------------------
# Obtain instantaneous speeds for each step using sapply
Tcell.step.speeds <- sapply( Tcell.steps, speed )
str( Tcell.step.speeds )
summary( Tcell.step.speeds )

## -----------------------------------------------------------------------------
d.step <- data.frame( method = "step-based", speed = Tcell.step.speeds )
d.cell <- data.frame( method = "cell-based", speed = Tcell.speeds )
d.method <- rbind( d.step, d.cell )

ggplot( d.method, aes( x = method, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## -----------------------------------------------------------------------------
Bcell.step.speeds <- sapply( subtracks(BCells2,1), speed )
Nphil.step.speeds <- sapply( subtracks(Neutrophils2,1), speed )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", speed = Tcell.step.speeds )
dB <- data.frame( cells = "B cells", speed = Bcell.step.speeds )
dN <- data.frame( cells = "Neutrophils", speed = Nphil.step.speeds )
dstep <- rbind( dT, dB, dN )

# Compare (violin plot now because there are many individual steps)
ggplot( dstep, aes( x = cells, y = speed ) ) +
  geom_violin() +
  #ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## -----------------------------------------------------------------------------
Tcell.step.angle <- sapply( subtracks(TCells2,2), overallAngle )
Bcell.step.angle <- sapply( subtracks(BCells2,2), overallAngle )
Nphil.step.angle <- sapply( subtracks(Neutrophils2,2), overallAngle )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", turn.angle = Tcell.step.angle )
dB <- data.frame( cells = "B cells", turn.angle = Bcell.step.angle )
dN <- data.frame( cells = "Neutrophils", turn.angle = Nphil.step.angle )
dangle <- rbind( dT, dB, dN )

# convert radians to degrees
dangle$turn.angle <- pracma::rad2deg( dangle$turn.angle )

# Compare:
ggplot( dangle, aes( x = cells, y = turn.angle ) ) +
  #ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  geom_violin() +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## -----------------------------------------------------------------------------
x <- TCells2[[4]]
x

## ---------------------------------------------------------------------------------------------------------------------
options(width = 120)
stag <- applyStaggered( x, speed, matrix = TRUE )
stag

## ---------------------------------------------------------------------------------------------------------------------
image( stag )

## ---- fig.width=4-----------------------------------------------------------------------------------------------------
filled.contour( stag )

## ---- fig.width=4-----------------------------------------------------------------------------------------------------
filled.contour( applyStaggered( TCells2[[1]], speed, matrix = TRUE ) )

## ---------------------------------------------------------------------------------------------------------------------
# These are the same:
applyStaggered( TCells2[[1]], speed )
mean( stag, na.rm = TRUE )

## ---------------------------------------------------------------------------------------------------------------------
# Compare the three different methods
Tcell.stag.speed <- sapply( TCells2, staggered( speed ) )
d.staggered <- data.frame( method = "staggered", speed = Tcell.speeds )
d.method <- rbind( d.step, d.cell, d.staggered )

ggplot( d.method, aes( x = method, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

# cell-based and staggered cell-based are slightly different:
summary( Tcell.speeds - Tcell.stag.speed )

## ---------------------------------------------------------------------------------------------------------------------
aggregate( TCells2, squareDisplacement, subtrack.length = 1 )

## ---------------------------------------------------------------------------------------------------------------------
aggregate( TCells2, squareDisplacement, subtrack.length = 1, FUN = "mean.se" )

## ---------------------------------------------------------------------------------------------------------------------
Tcell.msd <- aggregate( TCells2, squareDisplacement, FUN = "mean.se" )
str( Tcell.msd )

## ---------------------------------------------------------------------------------------------------------------------
ggplot( Tcell.msd, aes( x = i, y = mean ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (steps)") ),
        y = "mean square displacement") +
  theme_classic()

## ---------------------------------------------------------------------------------------------------------------------
num10 <- length( subtracks( TCells2, 10 ) )
num35 <- length( subtracks( TCells2, 35 ) )
c( "10" = num10, "35" = num35 )

## ---- fig.width = 7---------------------------------------------------------------------------------------------------
# Combine into a single dataframe with one column indicating the celltype
# To truly compare them, report subtrack length not in number of steps but
# in their duration (which may differ between different datasets)
Tcell.msd$cells <- "T cells"
Tcell.msd$dt <- Tcell.msd$i * timeStep( TCells2 )
Bcell.msd <- aggregate( BCells2, squareDisplacement, FUN = "mean.se" )
Bcell.msd$cells <- "B cells"
Bcell.msd$dt <- Bcell.msd$i * timeStep( BCells2 )
Nphil.msd <- aggregate( Neutrophils2, squareDisplacement, FUN = "mean.se" )
Nphil.msd$cells <- "Neutrophils"
Nphil.msd$dt <- Nphil.msd$i * timeStep( Neutrophils2 )
msddata <- rbind( Tcell.msd, Bcell.msd, Nphil.msd )
head(msddata)

# Plot
p1 <- ggplot( msddata, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "mean square displacement") +
  theme_classic() + theme(
    legend.position = "top"
  )
# Also make a zoomed version to look only at first part of the plot
pzoom <- p1 + coord_cartesian( xlim = c(0,500), ylim = c(0,3000) ) 
gridExtra::grid.arrange( p1, pzoom, ncol = 2 ) 

## ---- fig.width=6-----------------------------------------------------------------------------------------------------
Tcell.acor <- aggregate( TCells2, overallNormDot, FUN = "mean.se"  )
Tcell.acor$dt <- Tcell.acor$i * timeStep(TCells2)
Tcell.acor$cells <- "T cells"

ggplot( Tcell.acor, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = expression(cos(alpha))) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width=6-----------------------------------------------------------------------------------------------------
# Returns TRUE if first and last step of a track are of the minimum length of 1 micron
check <- function(x){ 
  all( sapply( list(head(x,2),tail(x,2)), trackLength ) >= 1.0 )
}
# repeat analysis
Tcell.acor2 <- aggregate( TCells2, overallNormDot, 
                     FUN = "mean.se", filter.subtracks = check )
Tcell.acor2$dt <- Tcell.acor2$i * timeStep(TCells2)
Tcell.acor2$cells <- "T cells (filtered)"

d <- rbind( Tcell.acor, Tcell.acor2 )

# compare
ggplot( d, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = expression(cos(alpha))) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width=6-----------------------------------------------------------------------------------------------------
# autocovariance
Tcell.acov <- aggregate( TCells2, overallDot, 
                     FUN = "mean.se", filter.subtracks = check )
Tcell.acov$dt <- Tcell.acov$i * timeStep(TCells2)
Tcell.acov$cells <- "T cells"


# compare
ggplot( Tcell.acov, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width = 7---------------------------------------------------------------------------------------------------
# Normalized autocovariance
Tcell.acov[,2:4] <- Tcell.acov[,2:4] / Tcell.acov$mean[1]

# Other cells
Bcell.acov <- aggregate( BCells2, overallDot, FUN = "mean.se" )
Bcell.acov$cells <- "B cells"
Bcell.acov$dt <- Bcell.acov$i * timeStep( BCells2 )
Bcell.acov[,2:4] <- Bcell.acov[,2:4] / Bcell.acov$mean[1]

Nphil.acov <- aggregate( Neutrophils2, overallDot, FUN = "mean.se" )
Nphil.acov$cells <- "Neutrophils"
Nphil.acov$dt <- Nphil.acov$i * timeStep( Neutrophils2 )
Nphil.acov[,2:4] <- Nphil.acov[,2:4] / Nphil.acov$mean[1]

acovdata <- rbind( Tcell.acov, Bcell.acov, Nphil.acov )

# compare;
p1 <- ggplot( acovdata, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank(), 
         legend.position = "top" )

pzoom <- p1 + coord_cartesian( xlim = c(0,500), ylim  = c(-0.1,1) )
gridExtra::grid.arrange( p1, pzoom, ncol = 2 )

## ---- fig.width = 6, fig.height = 3-----------------------------------------------------------------------------------
par( mfrow=c(1,2) )
plot( Neutrophils2 )
hotellingsTest( Neutrophils2, plot = TRUE, col = "gray", step.spacing = 10 )

## ---- fig.width = 7---------------------------------------------------------------------------------------------------
# Directional movement is 0.05 micron/sec in each dimension
directional.speed <- c( 0.05, 0.05 )
add.dir <- function( x, speed.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  speed.matrix <- matrix( rep( speed.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  speed.matrix <- speed.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords + speed.matrix
  return(x)
}

# Create data with directionality
TCells2.dir <- as.tracks( lapply( TCells2, add.dir, directional.speed ) )

# Plot both for comparison
# par(mfrow=c(1,2) )
# plot(TCells2, main = "original data" )
# plot(TCells2.dir, main = "with directional bias" )

## ---- fig.width = 6---------------------------------------------------------------------------------------------------
step.angles <- sapply( subtracks( TCells2.dir, 1), angleToDir, dvec = c(1,1) )
step.angles.original <- sapply( subtracks( TCells2, 1 ), angleToDir, dvec = c(1,1) )
par(mfrow=c(1,2) )
hist( step.angles.original, main = "original data" )
hist( step.angles, main = "with directional bias" )

## ---------------------------------------------------------------------------------------------------------------------
# Normalize the endpoint by subtracting its coordinates from all other coordinates:
normalize.endpoint <- function( x ){
  coords <- x[,-1]
  coords.norm <- t( apply( coords, 1, function(a) a - coords[nrow(coords),] ) )
  x[,-1] <- coords.norm
  return(x)
}
TCells2.point <- as.tracks( lapply( TCells2, normalize.endpoint ) )

# Plot for comparison:
plot( TCells2.point, main = "with point attraction" )

## ---------------------------------------------------------------------------------------------------------------------
hotellingsTest( TCells2.point, step.spacing = 10, col = "gray" )

## ---- fig.width = 7---------------------------------------------------------------------------------------------------
step.angles <- sapply( subtracks( TCells2.point, 1), angleToPoint, p = c(0,0) )
step.angles.original <- sapply( subtracks( TCells2, 1 ), angleToPoint, p = c(0,0) )
par(mfrow=c(1,2) )
hist( step.angles.original, main = "original data" )
hist( step.angles, main = "with point attraction" )

## ---------------------------------------------------------------------------------------------------------------------
step.distances <- sapply( subtracks( TCells2.point, 1), distanceToPoint, p = c(0,0) )

df <- data.frame( dist = step.distances,
                  angles = step.angles )
ggplot( df, aes( x = dist, y = angles ) ) +
  geom_point( size = 0.2, color = "gray") +
  geom_hline( yintercept = 90, color = "black" ) +
  stat_smooth( color = "red", method = "loess" ) +
  labs( x = "distance to reference point",
        y = "angle to reference point" ) +
  theme_classic()

## ---------------------------------------------------------------------------------------------------------------------
df.cells <- analyzeCellPairs( TCells2 )

# Plot
ggplot( df.cells, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray" ) +
  stat_smooth( span = 1, color = "red" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = 90, color = "black" ) +
  theme_classic()

## ---- warning = FALSE, message = FALSE--------------------------------------------------------------------------------
df.steps <- analyzeStepPairs( TCells2, filter.steps = function(x) displacement(x)>2  )

# Plot
ggplot( df.steps, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray", size = 0.5 ) +
  stat_smooth( color = "red" )+
  labs( x = "distance between step pairs",
        y = "angle between step pairs") +
  geom_hline( yintercept = 90, color = "black" ) +
  theme_classic()


## ---- echo = FALSE----------------------------------------------------------------------------------------------------
# Reset par() settings
par(oldpar)

