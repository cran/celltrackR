## ----setup, include = FALSE-------------------------------------------------------------------------------------------
matrix_check <- ( packageVersion( "Matrix" ) > "1.6.1" )
irlba_check <- ( packageVersion( "irlba" ) <= "2.3.5.1" )

matrix_irlba_clash <- TRUE #( matrix_check & irlba_check )
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4.5,
  fig.height = 3
)

## ----savepar, echo = FALSE--------------------------------------------------------------------------------------------
# Save current par() settings
oldpar <- par( no.readonly =TRUE )

## ----pack, warning = FALSE, message = FALSE---------------------------------------------------------------------------
library( celltrackR )

## ----Tdata------------------------------------------------------------------------------------------------------------
str( TCells, list.len = 2 )

## ----showTdata--------------------------------------------------------------------------------------------------------
head( TCells[[1]] )

## ----bdata------------------------------------------------------------------------------------------------------------
str( BCells, list.len = 2 )
str( Neutrophils, list.len = 2 )

## ----sampleData-------------------------------------------------------------------------------------------------------
# Take a sample
set.seed(1234)
TCells <- TCells[ sample( names(TCells), 30 ) ]
BCells <- BCells[ sample( names(BCells), 30 ) ]
Neutrophils <- Neutrophils[ sample( names(Neutrophils), 30 ) ]

## ----combineData------------------------------------------------------------------------------------------------------
T2 <- TCells
names(T2) <- paste0( "T", names(T2) )
tlab <- rep( "T", length(T2) )

B2 <- BCells
names(B2) <- paste0( "B", names(B2) )
blab <- rep( "B", length(B2) )

N2 <- Neutrophils
names(N2) <- paste0( "N", names(Neutrophils) )
nlab <- rep( "N", length( N2) )

all.tracks <- c( T2, B2, N2 )
real.celltype <- c( tlab, blab, nlab )


## ----featurematrix----------------------------------------------------------------------------------------------------
m <- getFeatureMatrix( all.tracks, 
                       c(speed, meanTurningAngle, 
                         outreachRatio, squareDisplacement) )

# We get a matrix with a row per track and one column for each metric:
head(m)

## ----plotFM, fig.width = 4, fig.height = 3----------------------------------------------------------------------------
plot( m, xlab = "speed", ylab = "mean turning angle" )

## ----pca--------------------------------------------------------------------------------------------------------------
pca <- trackFeatureMap( all.tracks, 
               c(speed,meanTurningAngle,squareDisplacement,
                 maxDisplacement,outreachRatio ), method = "PCA", 
               labels = real.celltype, return.mapping = TRUE )

## ----inspectPC--------------------------------------------------------------------------------------------------------
pc1 <- pca[,1]
pc2 <- pca[,2]
track.speed <- sapply( all.tracks, speed )
cor.test( pc1, track.speed )

## ----mds--------------------------------------------------------------------------------------------------------------
trackFeatureMap( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "MDS",
               labels = real.celltype )

## ----rspec, echo = FALSE, eval = matrix_irlba_clash-------------------------------------------------------------------
library( RSpectra )

## ----umap-------------------------------------------------------------------------------------------------------------
trackFeatureMap( all.tracks,
        c(speed,meanTurningAngle,squareDisplacement,
          maxDisplacement,outreachRatio ), method = "UMAP",
          labels = real.celltype )

## ----hclus, fig.width = 7.5, fig.height = 3---------------------------------------------------------------------------
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "hclust", labels = real.celltype )

## ----kmean, fig.height = 8, fig.width = 7-----------------------------------------------------------------------------
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "kmeans", 
               labels = real.celltype, centers = 3 )

## ----resetpar, echo = FALSE-------------------------------------------------------------------------------------------
# Reset par() settings
par(oldpar)

