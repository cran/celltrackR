## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(dpi=72)

## ----pack, warning = FALSE, message = FALSE---------------------------------------------------------------------------
library( celltrackR )
library( ggplot2 )

## ----savepar, echo = FALSE--------------------------------------------------------------------------------------------
# Save current par() settings
oldpar <- par( no.readonly =TRUE )

## ----Tdata------------------------------------------------------------------------------------------------------------
str( TCells, list.len = 3 )

## ----showTdata--------------------------------------------------------------------------------------------------------
head( TCells[[1]] )

## ----brownian1--------------------------------------------------------------------------------------------------------
brownian <- brownianTrack( nsteps = 20, dim = 3 )
str(brownian)

## ----plotSingleBrownian-----------------------------------------------------------------------------------------------
plot( wrapTrack(brownian) )

## ----plotBrownian, fig.width = 7--------------------------------------------------------------------------------------
brownian.tracks <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3 ) )
par(mfrow=c(1,2))
plot( brownian.tracks, main = "simulated random walk" )
plot( normalizeTracks( TCells ), main = "real T cell data" )

## ----matchDisp, fig.width = 7-----------------------------------------------------------------------------------------
# get displacement vectors
step.displacements <- t( sapply( subtracks(TCells,1), displacementVector ) )

# get mean and sd of displacement in each dimension
step.means <- apply( step.displacements, 2, mean )
step.sd <- apply( step.displacements, 2, sd )

# simulate brownian motion with the same statistics
brownian.tracks.matched <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = step.means,
                                                              sd = step.sd ) )

# compare displacement distributions
data.displacement <- sapply( subtracks( TCells,1), displacement )
matched.displacement <- sapply( subtracks( brownian.tracks.matched, 1 ), displacement )

df <- data.frame( disp = data.displacement,
                  data = "TCells" )
df2 <- data.frame( disp = matched.displacement,
                   data = "model" )
df <- rbind( df, df2 )
ggplot( df, aes( x = data, y = disp ) ) +
  geom_boxplot() +
  theme_classic()

# Plot new simulation versus real data
par(mfrow=c(1,2))
plot( brownian.tracks.matched, main = "simulated random walk" )
plot( normalizeTracks( TCells ), main = "real T cell data" )


## ----biasedWalk-------------------------------------------------------------------------------------------------------

# simulate brownian motion with bias
brownian.tracks.bias <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = c(1,1,1) ) )

plot( brownian.tracks.bias, main = "biased random walk" )


## ----beauchemin-------------------------------------------------------------------------------------------------------
beauchemin.tracks <- simulateTracks( 10, beaucheminTrack(sim.time=20) )
plot( beauchemin.tracks )

## ----bootstrap--------------------------------------------------------------------------------------------------------
bootstrap.tracks <- simulateTracks( 10, bootstrapTrack( nsteps = 20, TCells ) )
plot( bootstrap.tracks )

## ----distributionComp, fig.width = 7, warning = FALSE, message = FALSE------------------------------------------------
# Simulate more tracks to reduce noice
bootstrap.tracks <- simulateTracks( 50, bootstrapTrack( nsteps = 20, TCells ) )

# Compare step speeds in real data to those in bootstrap data
real.speeds <- sapply( subtracks( TCells,1 ), speed )
bootstrap.speeds <- sapply( subtracks( bootstrap.tracks,1), speed )
dspeed <- data.frame( tracks = c( rep( "data", length( real.speeds ) ),
                                  rep( "bootstrap", length( bootstrap.speeds ) ) ),
                      speed = c( real.speeds, bootstrap.speeds ) )

# Same for turning angles
real.angles <- sapply( subtracks( TCells,2 ), overallAngle, degrees = TRUE )
bootstrap.angles <- sapply( subtracks( bootstrap.tracks,2), overallAngle, degrees = TRUE )
dangle <- data.frame( tracks = c( rep( "data", length( real.angles ) ),
                                  rep( "bootstrap", length( bootstrap.angles ) ) ),
                      angle = c( real.angles, bootstrap.angles ) )

# plot
pspeed <- ggplot( dspeed, aes( x = tracks, y = speed ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

pangle <- ggplot( dangle, aes( x = tracks, y = angle ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

gridExtra::grid.arrange( pspeed, pangle, ncol = 2 )


## ----msdComp, fig.width=6---------------------------------------------------------------------------------------------
# Simulate more tracks
brownian.tracks <- simulateTracks( 50, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = step.means,
                                                              sd = step.sd ) )
bootstrap.tracks <- simulateTracks( 50, bootstrapTrack( nsteps = 20, TCells ) )

msd.data <- aggregate( TCells, squareDisplacement, FUN = "mean.se" )
msd.data$data <- "data"
msd.brownian <- aggregate( brownian.tracks, squareDisplacement, FUN = "mean.se" )
msd.brownian$data <- "brownian"
msd.bootstrap <- aggregate( bootstrap.tracks, squareDisplacement, FUN = "mean.se" )
msd.bootstrap$data <-"bootstrap"

msd <- rbind( msd.data, msd.brownian, msd.bootstrap )
ggplot( msd, aes( x = i, y = mean, ymin = lower, ymax = upper, color = data, fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "t (steps)",
        y = "square displacement" ) +
  scale_x_log10(limits= c(NA,10) ) +
  scale_y_log10() +
  theme_bw()

## ----acovComp, fig.width=6--------------------------------------------------------------------------------------------
# compute autocorrelation
acor.data <- aggregate( TCells, overallDot, FUN = "mean.se" )
acor.data$data <- "data"
acor.brownian <- aggregate( brownian.tracks, overallDot, FUN = "mean.se" )
acor.brownian$data <- "brownian"
acor.bootstrap <- aggregate( bootstrap.tracks, overallDot, FUN = "mean.se" )
acor.bootstrap$data <-"bootstrap"

acor <- rbind( acor.data, acor.brownian, acor.bootstrap )
ggplot( acor, aes( x = i, y = mean, ymin = lower, ymax = upper, color = data, fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "dt (steps)",
        y = "autocovariance" ) +
  scale_x_continuous(limits= c(0,10) ) +
  theme_bw()

## ----bb---------------------------------------------------------------------------------------------------------------
boundingBox( TCells )

## ----project----------------------------------------------------------------------------------------------------------
Txy <- projectDimensions( TCells, c("x","y") )

## ----msdXY, fig.width = 6---------------------------------------------------------------------------------------------
# Compute MSD
msdData <- aggregate( Txy, squareDisplacement, FUN = "mean.se" )

# Scale time axis: by default, this is in number of steps; convert to minutes.
tau <- timeStep( Txy ) / 60 # in minutes
msdData$dt <- msdData$i * tau

# plot
ggplot( msdData, aes( x = dt, y = mean ) ) +
  geom_ribbon( alpha = 0.3, aes( ymin = lower, ymax = upper ) ) +
  geom_line() +
  labs( x = expression( Delta*"t (min)" ), 
        y = expression( "displacement"^2*"("*mu*"m"^2*")") ) +
  coord_cartesian( xlim = c(0,NA), ylim = c(0,NA), expand = FALSE ) +
  theme_bw()

## ----fitThreshold-----------------------------------------------------------------------------------------------------
fitThreshold <- 5 # minutes

## ----furthFit---------------------------------------------------------------------------------------------------------
fuerthMSD <- function( dt, D, P, dim ){
  return( 2*dim*D*( dt - P*(1-exp(-dt/P) ) ) )
}

# Fit this function using nls. We fit only on the data where 
# dt < fitThreshold (see above), and need to provide reasonable starting
# values or the fitting algorithm will not work properly. 
model <- nls( mean ~ fuerthMSD( dt, D, P, dim = 2 ), 
              data = msdData[ msdData$dt < fitThreshold, ], 
              start = list( D = 10, P = 0.5 ), 
              lower = list( D = 0.001, P = 0.001 ), 
              algorithm = "port" 
)
D <- coefficients(model)[["D"]] # this is now in units of um^2/min
P <- coefficients(model)[["P"]] # persistence time in minutes
D
P

## ----fittedBrownian, fig.width = 5, fig.height = 5--------------------------------------------------------------------
# simulate tracks using the estimated D. We simulate with time unit seconds
# to match the data, so divide D by 60 to go from um^2/min to um^2/sec.
# The dimension-wise variance of brownian motion is 2*D*tau (with tau the step
# duration), so use this to parametrise the model:
tau <- timeStep( Txy )
brownianScaled <- function( nsteps, D, tau, dim = 2 ){
  tr <- brownianTrack( nsteps = nsteps, dim = 2, sd = sqrt( 2*(D/60)*tau ) )
  tr[,"t"] <- tr[,"t"]*tau # scale time
  return(tr)
}
brownian.tracks <- simulateTracks( 1000,  
                                   brownianScaled( nsteps = maxTrackLength(Txy ),
                                                   D = D, tau = tau ) )
plot( brownian.tracks[1:10], main = "brownian motion")


## ----msdFittedBrownian, fig.width = 6---------------------------------------------------------------------------------
# get msd of the simulated tracks and scale time to minutes
msdBrownian <- aggregate( brownian.tracks, squareDisplacement, FUN = "mean.se" )
tau <- timeStep( brownian.tracks ) / 60 # in minutes
msdBrownian$dt <- msdBrownian$i * tau

# Get Furth MSD using fitted P and D for comparison
msdFuerth <- data.frame( dt = msdBrownian$dt, 
                         mean = fuerthMSD( msdBrownian$dt, D, P, dim = 2 ))

# plot
ggplot( msdBrownian, aes( x = dt, y = mean) ) +
  geom_point( data = msdData, shape = 21, aes( fill = dt < fitThreshold ), color = "blue", show.legend = FALSE ) +
  geom_line( data = msdFuerth, color = "red", lty = 2 ) +
  geom_line( )+
  scale_fill_manual( values = c( "TRUE" = "blue", "FALSE" = "transparent" ) ) +
  scale_x_log10( expand = c(0,0) ) +
  scale_y_log10( expand = c(0,0) ) +
  labs( color = NULL, x = expression( Delta*"t (min)"), 
        y = expression( "displacement"^2*"("*mu*"m"^2*")") , fill = NULL ) +
  theme_bw()

## ----fittedBeauchemin-------------------------------------------------------------------------------------------------
# analytical formula
beauchemin.msd <- function( x, M=60, t.free=2, t.pause=0.5, dim=2 ){
  v.free <- sqrt( 6 * M * (t.free+t.pause) ) / t.free
  M <- (v.free^2 * t.free^2) / 6 / (t.free + t.pause)
  multiplier <- rep(1/3, length(x))
  xg <- x[x<t.free]
  multiplier[x<t.free] <- 1/3 * (xg/t.free)^3 - (xg/t.free)^2 + (xg/t.free)
  2 * M * x * dim - 2 * M * dim * t.free * multiplier
}
# Again, fit on the first fitThreshold min, before confinement becomes an issue.
beaucheminFit <- nls( mean ~ beauchemin.msd( dt, M, t.free ), 
                      msdData[ msdData$dt < fitThreshold, ], 
                      start=list(M=30, t.free=1))

# extract parameters M and t.free
M <- coef(beaucheminFit)[["M"]] # in um^2/min
t.free <- coef(beaucheminFit)[["t.free"]] # in min

## ----compareFittedParms-----------------------------------------------------------------------------------------------
c( D = D, M = M )

## ----beaucheminMSD, fig.width = 6-------------------------------------------------------------------------------------
# Get the fitted MSD using these values
msdBeauchemin <- data.frame(
  dt = msdBrownian$dt,
  mean = beauchemin.msd( msdBrownian$dt, M, t.free )
)

# plot
ggplot( msdBeauchemin, aes( x = dt, y = mean) ) +
  geom_point( data = msdData, shape = 21, aes( fill = dt < fitThreshold ), color = "blue", show.legend = FALSE ) +
  geom_line( )+
  scale_fill_manual( values = c( "TRUE" = "blue", "FALSE" = "transparent" ) ) +
  scale_x_log10( expand = c(0,0) ) +
  scale_y_log10( expand = c(0,0) ) +
  labs( color = NULL, x = expression( Delta*"t (min)"), 
        y = expression( "displacement"^2*"("*mu*"m"^2*")") , fill = NULL ) +
  theme_bw()

## ----beaucheminParameters, fig.width=5, fig.height = 5----------------------------------------------------------------
# t.free has been fitted; set t.pause to its fixed value and 
# compute v.free from the fitted M:
t.pause <- 0.5 # fixed
v.free <- sqrt( 6 * M * (t.free+t.pause) ) / t.free

# simulate 3 tracks using the fitted parameters.
# we don't use the arguments p.persist, p.bias, bias.dir, and taxis.mode,
# since these are not parameters of the original Beauchemin model.
tau <- timeStep( Txy )
beauchemin.tracks <- simulateTracks( 3,  
                                   beaucheminTrack( sim.time = 10*max( timePoints( Txy) ),
                                                    delta.t = tau,
                                                    t.free = t.free,
                                                    v.free = v.free,
                                                    t.pause = t.pause ) )

plot( beauchemin.tracks, main = "Beauchemin tracks")

## ----resetpar, echo = FALSE-------------------------------------------------------------------------------------------
# Reset par() settings
par(oldpar)

