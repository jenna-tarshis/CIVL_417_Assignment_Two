# read the free surface elevation data
data <- read.table("data/FreeSurfaceElevationTS.txt", header=TRUE)

# plot for visual ref...
plot(data$Time, data$Elevation, type="l", col="blue",
     xlab="Time (s)", ylab="Elevation (m)",
     main="Free Surface Elevation Time Series")
abline(h=mean_elevation, col="red", lty=2)  # Add mean elevation line

sink("results/Q3_results.txt")


# -------------------------------- PART A --------------------------------
mean_elevation <- mean(data$Elevation)
variance_elevation <- var(data$Elevation)
cat("Mean Free Surface Elevation:", mean_elevation, "\n")
cat("Variance of Free Surface Elevation:", variance_elevation, "\n")


# -------------------------------- PART B --------------------------------
signs <- sign(data$Elevation - mean_elevation)
next_sign <- c(signs[-1], NA)
up_crossings <- which(signs == -1 & next_sign == 1)
cat("Number of up-crossings:", length(up_crossings), "\n")


# -------------------------------- PART C --------------------------------
# ... see document


# -------------------------------- PART D --------------------------------
wave_heights <- numeric(length(up_crossings) - 1)

for (i in 1:(length(up_crossings) - 1)) {
  start_idx <- up_crossings[i]
  end_idx <- up_crossings[i + 1]
  
  # Find max (crest) and min (trough) in this wave
  wave_heights[i] <- max(data$Elevation[start_idx:end_idx]) - min(data$Elevation[start_idx:end_idx])
}

# Sort wave heights in descending order
sorted_wave_heights <- sort(wave_heights, decreasing = TRUE)

# Compute significant wave height (Hs) - mean of highest 1/3 of waves
num_significant_waves <- ceiling(length(sorted_wave_heights) / 3)
Hs <- mean(sorted_wave_heights[1:num_significant_waves])

# Compute wave periods - time difference between zero up-crossings
wave_periods <- diff(data$Time[up_crossings])

# Select corresponding periods for significant waves
sorted_wave_periods <- sort(wave_periods, decreasing = TRUE)
Ts <- mean(sorted_wave_periods[1:num_significant_waves])

# Print results
print(paste("Significant Wave Height (Hs):", round(Hs, 2), "m"))
print(paste("Significant Wave Period (Ts):", round(Ts, 2), "s"))


# -------------------------------- PART E --------------------------------

# Find wave heights (crest - trough between two up-crossings)
for (i in 1:(length(up_crossings) - 1)) {
  start_idx <- up_crossings[i]
  end_idx <- up_crossings[i + 1]
  
  # Find max (crest) and min (trough) in this wave
  wave_heights[i] <- max(data$Elevation[start_idx:end_idx]) - min(data$Elevation[start_idx:end_idx])
}

# Compute mean wave height (Hm)
Hm <- mean(wave_heights)

# Compute wave periods - time difference between zero up-crossings
wave_periods <- diff(data$Time[up_crossings])

# Compute mean wave period (Tm)
Tm <- mean(wave_periods)

# Print results
print(paste("Mean Wave Height (Hm):", round(Hm, 2), "m"))
print(paste("Mean Wave Period (Tm):", round(Tm, 2), "s"))


# -------------------------------- PART F --------------------------------
# Compute root-mean-square wave height (Hrms)
Hrms <- sqrt(mean(wave_heights^2))

# Print result
print(paste("Root-Mean-Square Wave Height (Hrms):", round(Hrms, 2), "m"))


# -------------------------------- PART G --------------------------------


sink()

