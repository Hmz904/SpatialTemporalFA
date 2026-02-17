rm(list = ls())
# Load necessary libraries
library(FITSio)  # For reading FITS files
library(ggplot2)  # For advanced plotting
library(gridExtra)  # For combining multiple ggplot plots
library(reshape2)  # For data reshaping
library(dplyr)     # For data manipulation
library(stringr)   # For string manipulation
library(urca)      # For cointegration tests
library(fields)    # For image plotting

# Set random seed for reproducibility
set.seed(123)

# Read FITS data - update the path as needed
datah <- (readFITS("D:/Google Download/datacube_128x128x296.fits"))$imDat

# Extract data dimensions
dim_data <- dim(datah)
N <- dim_data[1]  # Assuming 128x128 image
T <- dim_data[3]  # Number of time points
q <- 10  # Number of principal components

# Restructure 3D array into matrix format
# Each row is a pixel, each column is a time point
Y <- matrix(c(datah), nrow = N^2, ncol = T, byrow = FALSE)

# Calculate spatial covariance matrix
S1 <- t(Y) %*% Y

# Eigendecomposition
eigen_result <- eigen(S1)
lambda <- eigen_result$values
L.hat <- eigen_result$vectors[, 1:q]
Lam <- diag(lambda[1:q])

# Calculate loading matrix
L.hat <- Y %*% L.hat %*% sqrt(solve(Lam))

# Calculate factor scores
F.hat <- t(Y) %*% L.hat / N

# Reshape loading matrix into 3D array
L.hat3 <- array(L.hat, dim = c(N, N, q))

#------------------------------------------------------------
# Find extreme pixels in components 2-5 (top 0.5% by absolute value)
#------------------------------------------------------------

# Create a list to store extreme pixel coordinates for each component
extreme_pixels <- list()

# Find extreme pixels for components 2-5
for (comp in 2:5) {
  # Get component loadings
  component_loadings <- L.hat3[,,comp]
  
  # Calculate threshold for top 0.5% of absolute values
  threshold <- quantile(abs(component_loadings), 0.995)
  
  # Find pixels that exceed the threshold
  extreme_indices <- which(abs(component_loadings) > threshold, arr.ind = TRUE)
  
  # Store the coordinates and values
  extreme_pixels[[comp-1]] <- data.frame(
    row = extreme_indices[,1],
    col = extreme_indices[,2],
    value = component_loadings[extreme_indices],
    component = comp
  )
  
  # Print summary
  cat(sprintf("Component %d: Found %d extreme pixels (threshold = %.4f)\n", 
              comp, nrow(extreme_indices), threshold))
}

# Combine all extreme pixels
all_extreme_pixels <- do.call(rbind, extreme_pixels)
cat("Total extreme pixels across components 2-5:", nrow(all_extreme_pixels), "\n")

#------------------------------------------------------------
# MODIFIED: Filter extreme pixels to only include those with at least one neighbor
# This is now done BEFORE the cointegration tests
#------------------------------------------------------------

# Count neighbors for extreme pixels
count_neighbors <- function(pixels_df) {
  # Create a map of pixel positions
  pixel_map <- matrix(0, nrow = N, ncol = N)
  for (i in 1:nrow(pixels_df)) {
    pixel_map[pixels_df$row[i], pixels_df$col[i]] <- 1
  }
  
  # Count neighbors for each pixel
  neighbor_counts <- numeric(nrow(pixels_df))
  for (i in 1:nrow(pixels_df)) {
    row <- pixels_df$row[i]
    col <- pixels_df$col[i]
    
    # Define 3x3 window boundaries
    row_start <- max(1, row - 1)
    row_end <- min(N, row + 1)
    col_start <- max(1, col - 1)
    col_end <- min(N, col + 1)
    
    # Count neighbors (excluding self)
    neighbors <- sum(pixel_map[row_start:row_end, col_start:col_end]) - 1
    neighbor_counts[i] <- neighbors
  }
  
  return(neighbor_counts)
}

# Filter extreme pixels
filter_extreme_pixels <- function(pixels_df, min_neighbors = 1) {
  # Count neighbors
  pixels_df$neighbors <- count_neighbors(pixels_df)
  
  # Filter pixels with at least min_neighbors
  filtered_pixels <- pixels_df[pixels_df$neighbors >= min_neighbors, ]
  
  cat("Original extreme pixels count:", nrow(pixels_df), "\n")
  cat("Filtered extreme pixels count:", nrow(filtered_pixels), "\n")
  cat("Removed", nrow(pixels_df) - nrow(filtered_pixels), "isolated points\n")
  
  return(filtered_pixels)
}

# Filter extreme pixels to only include those with at least one neighbor
# This is now done BEFORE cointegration tests
filtered_extreme_pixels <- filter_extreme_pixels(all_extreme_pixels, min_neighbors = 1)



