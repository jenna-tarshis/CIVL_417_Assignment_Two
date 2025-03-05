library(ggplot2)

#read wave spectrum data
data <- read.table("data/WaveSpectrumFromFFT.txt", header=TRUE)
# Accesses the first column dynamically cause the column names have special characters like []
frequency <- data[[names(data)[1]]] 
spectrum <- data[[names(data)[2]]]


# -------------------------------- PART A --------------------------------

peak_index <- which.max(spectrum)
fp <- frequency[peak_index]  # Peak frequency
Tp <- 1 / fp                 # Peak wave period

cat("Peak Frequency (fp):", fp, "Hz\n")
cat("Peak Wave Period (Tp):", Tp, "seconds\n")

ggplot(data, aes(x = frequency, y = spectrum)) +
  geom_line(color = "blue", size = 1) +
  geom_point(aes(x = fp, y = max(spectrum)), color = "red", size = 3) +
  ggtitle("Wave Frequency Spectrum") +
  xlab("Frequency (Hz)") +
  ylab("Spectral Density (m^2/Hz)") +
  theme_minimal()


# -------------------------------- PART B --------------------------------

f_diffs <- diff(frequency)
print(f_diffs)
df <- f_diffs[1]
print(df)

m0 <- sum(spectrum) * df
cat("Zeroth Moment (m0):", m0, "\n")


# -------------------------------- PART C --------------------------------
Hs_spectral <- 4 * sqrt(m0)
cat("Spectral significant wave height:", Hs_spectral, "\n")


# -------------------------------- PART C --------------------------------
PM_spectrum <- (5 * Hs_spectral^2) / (16 * fp) * (frequency / fp)^(-5) * 
  exp(-5/4 * (frequency / fp)^(-4))


# -------------------------------- PART D --------------------------------
ggplot() +
  geom_line(data = data, aes(x = frequency, y = spectrum), color = "blue", size = 1, linetype = "solid") +
  geom_line(aes(x = frequency, y = PM_spectrum), color = "red", size = 1, linetype = "dashed") +
  ggtitle("Wave Spectrum with Pierson-Moskowitz Fit") +
  xlab("Frequency (Hz)") +
  ylab("Spectral Density (m^2/Hz)") +
  theme_minimal() +
  scale_color_manual(name = "Legend", values = c("Measured" = "blue", "P-M Model" = "red"))

# -------------------------------- PART E --------------------------------
#... see doc

