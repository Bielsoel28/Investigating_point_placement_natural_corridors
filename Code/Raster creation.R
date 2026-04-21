##########################################################################
# This script is Raster creation for Investigating the influence         # 
# of point placement strategies in the stability of modeled movement     #   
# corridors: The case of the middle neolithic in the North-East          #                                                   
# Biel Soriano-Elias, Anna Bach G?mez, Miquel Molist Montany?            #                                                   
#                                                                        #                                                   
# Author: Biel Soriano Elias                                             #
# Affiliation : Autonomous University of Barcelona                       #
# Creation date : 8/1/2026                                               #
# E-mail: biel.soriano@uab.cat                                           #
##########################################################################

# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Set Working directory
setwd("XXXXX/XXXXX/XXXXX")

# Clean up workspace
rm(list=ls())

# 0.2 Install packages =========================================================

# Needed packages
packages <- c("raster","leastcostpath","terra", "sf", "ggplot2","dplyr", "spatstat", "stars","corrplot","tidyr")

#Optional, Run this if the pacakges are not already installed
for (package in packages) {
  install.packages(package, character.only = TRUE)
}

#Load the packages
for (package in packages) {
  library(package, character.only = TRUE)
}

# 0.3 Show session infos =======================================================

sessionInfo()

# 1 Setting up needed functions ################################################

#Function to generate fractal DEM using Fourier Filtering
Fractal_dem <- function(n, hurst, seed = 1234) {
  set.seed(seed)
  
  # Define parameters
  beta <- 2 * hurst + 1
  size <- n
  
  # Initialize arrays for Fourier coefficients
  A <- matrix(0, nrow = size, ncol = size)
  B <- matrix(0, nrow = size, ncol = size)
  
  # Generate Fourier coefficients
  for (i in 1:(size / 2)) {
    for (j in 1:(size / 2)) {
      radius <- (sqrt(i^2 + j^2))^(-beta / 2) * rnorm(1)
      phase <- 2 * pi * runif(1)
      A[i, j] <- radius * cos(phase)
      B[i, j] <- radius * sin(phase)
    }
  }
  
  # Create complex matrix for FFT
  complex_matrix <- A + 1i * B
  
  # Perform 2D inverse FFT
  terrain <- Re(fft(complex_matrix, inverse = TRUE))
  
  # Normalize values to range 0-255 for DEM
  terrain <- (terrain - min(terrain)) / (max(terrain) - min(terrain)) * 10
  
  # Convert to terra raster
  dem <- rast(terrain)
  ext(dem) <- c(0, n, 0, n) # Set arbitrary extent for spatial dimensions
  return(dem)
}

# 2 Raster creation ############################################################
# H defines roughness
hurst_exponent <- 0.8  

# Raster A
Raster_A <- generate_fractal_dem_terra(1000, hurst_exponent, 300)
writeRaster(dem, filename = "Raster_A.tif") # Saving the results

# Raster B
Raster_B <- generate_fractal_dem_terra(1000, hurst_exponent, 600)
writeRaster(dem, filename = "Raster_B.tif") # Saving the results

# Raster C
Raster_C <- generate_fractal_dem_terra(1000, hurst_exponent, 900)
writeRaster(dem, filename = "Raster_C.tif") # Saving the results

# Raster un sis km
Raster_un_sis <- generate_fractal_dem_terra(1600, hurst_exponent, 300)
writeRaster(dem, filename = "Raster_un_sis_km.tif") # Saving the results

# Raster dos
Raster_dos <- generate_fractal_dem_terra(2000, hurst_exponent, 300)
writeRaster(dem, filename = "Raster_dos_km.tif") # Saving the results  


