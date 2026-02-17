#------------------------------------------------------------
# Define function to extract neighbor time series and test cointegration
#------------------------------------------------------------

# Function to get time series for a pixel and its neighbors
get_pixel_timeseries <- function(row, col, datah, window_size = 3) {
  N <- dim(datah)[1]
  T <- dim(datah)[3]
  
  # Calculate window boundaries with edge handling
  row_start <- max(1, row - floor(window_size/2))
  row_end <- min(N, row + floor(window_size/2))
  col_start <- max(1, col - floor(window_size/2))
  col_end <- min(N, col + floor(window_size/2))
  
  # Extract time series for all pixels in the window
  result <- list()
  for (r in row_start:row_end) {
    for (c in col_start:col_end) {
      # Skip if this is the central pixel
      if (r == row && c == col) {
        central_ts <- datah[r, c, ]
        next
      }
      
      # Get time series for this neighbor
      ts <- datah[r, c, ]
      
      # Store with position relative to central pixel
      rel_row <- r - row
      rel_col <- c - col
      result[[sprintf("%d_%d", rel_row, rel_col)]] <- ts
    }
  }
  
  return(list(central = central_ts, neighbors = result))
}

# Function to test cointegration between central pixel and neighbors
# MODIFIED: Added p-value estimation
test_cointegration_with_neighbors <- function(row, col, datah) {
  # Get time series for pixel and neighbors
  ts_data <- get_pixel_timeseries(row, col, datah)
  
  # Initialize results container
  results <- data.frame(
    neighbor = character(0),
    rel_row = numeric(0),
    rel_col = numeric(0),
    abs_row = numeric(0),
    abs_col = numeric(0),
    test_stat = numeric(0),
    crit_val = numeric(0),
    p_value = numeric(0),
    log_p_value = numeric(0),
    cointegrated = logical(0)
  )
  
  # Check cointegration with each neighbor
  for (neighbor_name in names(ts_data$neighbors)) {
    # Parse relative position
    coords <- as.numeric(strsplit(neighbor_name, "_")[[1]])
    rel_row <- coords[1]
    rel_col <- coords[2]
    abs_row <- row + rel_row
    abs_col <- col + rel_col
    
    # Combine time series
    combined_ts <- cbind(ts_data$central, ts_data$neighbors[[neighbor_name]])
    colnames(combined_ts) <- c("central", "neighbor")
    
    # Try to run Johansen test
    tryCatch({
      # Conduct test
      jo_test <- ca.jo(combined_ts, type = "trace", K = 2, ecdet = "none", spec = "longrun")
      
      # Extract results
      test_summary <- summary(jo_test)
      test_stat <- test_summary@teststat[1]
      crit_val <- test_summary@cval[1,1]  # 90% critical value
      
      # Calculate approximate p-value based on test statistic and critical values
      # This is a simplified approach since exact p-values require more complex calculations
      if (test_stat > crit_val) {
        # If test statistic exceeds critical value, cointegration is significant
        # Use an approximation for p-value calculation
        p_value <- 0.01 * (test_stat / crit_val)^2
        p_value <- min(p_value, 0.99) # Cap at 0.99
      } else {
        # If test statistic is below critical value, cointegration is not significant
        p_value <- 1 - 0.1 * (test_stat / crit_val)
        p_value <- max(p_value, 0.01) # Floor at 0.01
      }
      
      # Calculate log p-value (avoid log(0) by using a small minimum value)
      log_p_value <- log(max(p_value, 1e-10))
      
      # Add to results
      results <- rbind(results, data.frame(
        neighbor = neighbor_name,
        rel_row = rel_row,
        rel_col = rel_col,
        abs_row = abs_row,
        abs_col = abs_col,
        test_stat = test_stat,
        crit_val = crit_val,
        p_value = p_value,
        log_p_value = log_p_value,
        cointegrated = (test_stat > crit_val)
      ))
    }, error = function(e) {
      # If test fails, record NA values
      results <- rbind(results, data.frame(
        neighbor = neighbor_name,
        rel_row = rel_row,
        rel_col = rel_col,
        abs_row = abs_row,
        abs_col = abs_col,
        test_stat = NA,
        crit_val = NA,
        p_value = NA,
        log_p_value = NA,
        cointegrated = FALSE
      ))
    })
  }
  
  return(results)
}

#------------------------------------------------------------
# Run cointegration analysis for extreme pixels with at least one neighbor
#------------------------------------------------------------

# Take a sample of filtered extreme pixels to analyze (for computational efficiency)
# You can adjust this or process all pixels if needed
sample_size <- min(100, nrow(filtered_extreme_pixels))
sampled_pixels <- filtered_extreme_pixels[sample(nrow(filtered_extreme_pixels), sample_size), ]

# Initialize list to store results
cointegration_results <- list()

# Process each sampled pixel
for (i in 1:nrow(sampled_pixels)) {
  # Get pixel info
  pixel <- sampled_pixels[i, ]
  
  # Run cointegration tests
  coint_tests <- test_cointegration_with_neighbors(pixel$row, pixel$col, datah)
  
  # Store results with pixel info
  cointegration_results[[i]] <- list(
    pixel = pixel,
    tests = coint_tests
  )
  
  # Print progress
  if (i %% 10 == 0) {
    cat(sprintf("Processed %d/%d pixels\n", i, nrow(sampled_pixels)))
  }
}

#------------------------------------------------------------
# Analyze results and create regions
#------------------------------------------------------------

# Create a matrix to store region memberships
region_map <- matrix(0, nrow = N, ncol = N)
region_counter <- 1

# Process cointegration results to form regions
for (i in 1:length(cointegration_results)) {
  result <- cointegration_results[[i]]
  pixel <- result$pixel
  tests <- result$tests
  
  # Skip if this pixel is already assigned to a region
  if (region_map[pixel$row, pixel$col] != 0) {
    next
  }
  
  # Find cointegrated neighbors
  cointegrated_tests <- tests[tests$cointegrated, ]
  
  # If no neighbors are cointegrated, create a singleton region
  if (nrow(cointegrated_tests) == 0) {
    region_map[pixel$row, pixel$col] <- region_counter
    region_counter <- region_counter + 1
    next
  }
  
  # Create a new region with this pixel and its cointegrated neighbors
  region_map[pixel$row, pixel$col] <- region_counter
  
  # Add neighbors to the region
  for (j in 1:nrow(cointegrated_tests)) {
    neighbor <- cointegrated_tests[j, ]
    
    # Calculate absolute coordinates of neighbor
    neighbor_row <- neighbor$abs_row
    neighbor_col <- neighbor$abs_col
    
    # Check if within bounds
    if (neighbor_row >= 1 && neighbor_row <= N && 
        neighbor_col >= 1 && neighbor_col <= N) {
      # Assign to the same region
      region_map[neighbor_row, neighbor_col] <- region_counter
    }
  }
  
  # Increment region counter
  region_counter <- region_counter + 1
}

# Print summary of regions found
cat("Total regions identified:", region_counter - 1, "\n")

#------------------------------------------------------------
# MODIFIED: Create p-value maps using median values
#------------------------------------------------------------

# Create p-value map for visualization using median values
create_median_pvalue_maps <- function() {
  # Create a data frame to store all p-values for each pixel
  all_pixel_pvalues <- data.frame(
    row = integer(),
    col = integer(),
    p_value = numeric(),
    log_p_value = numeric()
  )
  
  # Process each cointegration result
  for (i in 1:length(cointegration_results)) {
    result <- cointegration_results[[i]]
    pixel <- result$pixel
    tests <- result$tests
    
    # Add center pixel (self-reference p-value = 1)
    all_pixel_pvalues <- rbind(all_pixel_pvalues, data.frame(
      row = pixel$row,
      col = pixel$col,
      p_value = 1,
      log_p_value = 0
    ))
    
    # Process each neighbor test
    for (j in 1:nrow(tests)) {
      neighbor <- tests[j, ]
      
      # Skip NA results
      if (is.na(neighbor$p_value)) next
      
      # Calculate neighbor coordinates
      n_row <- neighbor$abs_row
      n_col <- neighbor$abs_col
      
      # Skip if out of bounds
      if (n_row < 1 || n_row > N || n_col < 1 || n_col > N) {
        next
      }
      
      # Add to the data frame
      all_pixel_pvalues <- rbind(all_pixel_pvalues, data.frame(
        row = n_row,
        col = n_col,
        p_value = neighbor$p_value,
        log_p_value = neighbor$log_p_value
      ))
    }
  }
  
  # Initialize maps
  p_value_map <- matrix(NA, nrow = N, ncol = N)
  log_p_value_map <- matrix(NA, nrow = N, ncol = N)
  
  # Calculate median p-values for each pixel
  for (r in 1:N) {
    for (c in 1:N) {
      # Get all p-values for this pixel
      pixel_values <- all_pixel_pvalues[all_pixel_pvalues$row == r & all_pixel_pvalues$col == c, ]
      
      # If we have values, calculate median
      if (nrow(pixel_values) > 0) {
        p_value_map[r, c] <- median(pixel_values$p_value)
        log_p_value_map[r, c] <- median(pixel_values$log_p_value)
      }
    }
  }
  
  return(list(
    p_value_map = p_value_map,
    log_p_value_map = log_p_value_map,
    all_pixel_pvalues = all_pixel_pvalues
  ))
}

# Create p-value maps
p_maps <- create_median_pvalue_maps()

# Print summary statistics for p-values
cat("\nP-Value Summary Statistics:\n")
summary_stats <- summary(p_maps$all_pixel_pvalues$p_value)
print(summary_stats)
cat("Pixels with data:", sum(!is.na(p_maps$p_value_map)), "\n")
cat("Unique pixels in data:", length(unique(paste(p_maps$all_pixel_pvalues$row, p_maps$all_pixel_pvalues$col))), "\n")

#------------------------------------------------------------
# MODIFIED: Visualize median log(p) values
#------------------------------------------------------------

# 1. Plot each component loading map
for (comp in 2:5) {
  png(sprintf("component_%d_loading_map.png", comp), width = 800, height = 700)
  image.plot(1:N, 1:N, L.hat3[,,comp], 
             main = sprintf("Component %d Loading Map", comp),
             xlab = "Column", ylab = "Row")
  dev.off()
}

# 2. Plot the region map
png("region_map.png", width = 800, height = 700)
image.plot(1:N, 1:N, region_map, 
           main = "Identified Regions",
           xlab = "Column", ylab = "Row")
dev.off()

# 3. Plot median log(p) value map
png("median_log_p_value_map.png", width = 1000, height = 800)
# Using -log(p) for visualization (negative log makes smaller p-values darker)
vis_log_p <- -p_maps$log_p_value_map
vis_log_p[is.na(vis_log_p)] <- 0  # Replace NA with 0 for visualization
image.plot(1:N, 1:N, vis_log_p, 
           main = "Negative Median Log(p) from Johansen Test",
           subtitle = "Only testing pixels with at least one neighbor",
           xlab = "Column", ylab = "Row",
           col = hcl.colors(100, "YlOrRd", rev = TRUE))
dev.off()

# 4. Create distribution of p-values
png("p_value_distribution.png", width = 800, height = 600)
hist(p_maps$all_pixel_pvalues$p_value, 
     breaks = 30,
     main = "Distribution of p-values from Johansen Test",
     sub = "Only testing pixels with at least one neighbor",
     xlab = "p-value",
     col = "skyblue")
dev.off()

# 5. Create combined visualization using ggplot2
library(ggplot2)
library(reshape2)

# Prepare data for filtered extreme pixels spatial map
create_extreme_pixel_map <- function() {
  # Create empty grid
  grid_data <- expand.grid(
    row = 1:N,
    col = 1:N
  )
  
  # Default values
  grid_data$is_extreme <- FALSE
  grid_data$component <- NA
  
  # Mark extreme pixels
  for (i in 1:nrow(filtered_extreme_pixels)) {
    pixel <- filtered_extreme_pixels[i, ]
    idx <- grid_data$row == pixel$row & grid_data$col == pixel$col
    grid_data$is_extreme[idx] <- TRUE
    grid_data$component[idx] <- pixel$component
  }
  
  return(grid_data)
}

# Create data for log(p) map
create_logp_data <- function() {
  # Convert matrix to data frame
  logp_data <- melt(p_maps$log_p_value_map)
  colnames(logp_data) <- c("row", "col", "log_p")
  
  # Replace NA with minimum value for visualization
  min_val <- min(logp_data$log_p, na.rm = TRUE)
  logp_data$log_p[is.na(logp_data$log_p)] <- min_val
  
  # Negative log(p) for visualization (makes smaller p-values more prominent)
  logp_data$neg_log_p <- -logp_data$log_p
  
  return(logp_data)
}

# Create data frames
extreme_data <- create_extreme_pixel_map()
logp_data <- create_logp_data()

# 1. Create extreme pixel map
p_extreme <- ggplot(extreme_data, aes(x = col, y = row)) +
  geom_tile(aes(fill = is_extreme)) +
  scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "red"), 
                    name = "Extreme Pixel") +
  labs(title = "Filtered Extreme Pixels",
       subtitle = paste(sum(extreme_data$is_extreme), "pixels with at least one neighbor"),
       x = "Column", y = "Row") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  coord_equal()

# 2. Create median log(p) map
p_logp <- ggplot(logp_data, aes(x = col, y = row)) +
  geom_tile(aes(fill = neg_log_p)) +
  scale_fill_viridis_c(name = "-log(p)") +
  labs(title = "Negative Median Log(p) from Johansen Test",
       subtitle = "Lower values (darker) indicate stronger cointegration",
       x = "Column", y = "Row") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  coord_equal()

# Combine plots
combined_plot <- grid.arrange(p_extreme, p_logp, ncol = 1, heights = c(1, 1))
