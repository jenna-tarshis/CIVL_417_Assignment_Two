library(ggplot2)

# -------------------------------- PART A --------------------------------
T <- 10  # Wave period (s)
g <- 9.81  # Gravity (m/s^2)
d <- 25  # Water depth (m)
H <- 2  # Wave height (m)
a <- H / 2  # Amplitude (m)

omega <- 2 * pi / T  #angular frequency
lambda_0 <- (g * T^2) / (2 * pi)  # deep-water wavelength approx

# Function to solve dispersion relation iteratively
dispersion_relation <- function(k, d, omega, g) {
  return(omega^2 - g * k * tanh(k * d))
}

# Newton-Raphson method to find k
k_guess <- 2 * pi / lambda_0  #initial guess
tolerance <- 1e-6
max_iterations <- 100

for (i in 1:max_iterations) {
  f_k <- dispersion_relation(k_guess, d, omega, g)
  df_k <- -g * (tanh(k_guess * d) + k_guess * d * (1 / cosh(k_guess * d))^2)
  k_new <- k_guess - f_k / df_k
  
  if (abs(k_new - k_guess) < tolerance) {
    break
  }
  
  k_guess <- k_new
}

k <- k_new  # final computed wavenumber
lambda_exact <- 2 * pi / k  # Exact wavelength

# Define x values over one wavelength
x <- seq(0, lambda_exact, length.out = 200)  
phi <- k * x  # Phase variable

# Compute eta1 and eta2
eta1 <- a * cos(phi)  # First-order component
eta2 <- (0.5 * a^2 * k * cosh(2 * k * d) * cos(2 * phi))  # Second-order component
eta <- eta1 + eta2  # Total free-surface elevation

# Create a data frame for plotting
wave_data <- data.frame(
  x = rep(x, 3),
  elevation = c(eta1, eta2, eta),
  component = rep(c("First-Order (η1)", "Second-Order (η2)", "Total Elevation (η)"), each = length(x))
)

# Plot using ggplot2
ggplot(wave_data, aes(x = x, y = elevation, color = component)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("blue", "red", "black")) +
  labs(title = "Stokes Second-Order Wave Profile",
       x = "Position (m)",
       y = "Wave Elevation (m)",
       color = "Component") +
  theme_minimal() +
  theme(legend.position = "bottom")


# -------------------------------- PART B --------------------------------
# ... see doc


# -------------------------------- PART C --------------------------------
# ... see part A, same script but change initial water depth to 25m


# -------------------------------- PART D --------------------------------
H_s <- 2  # in meters
H_max <- 1.86 * H_s
cat("The conservative estimate of H_max is:", round(H_max, 2), "meters\n")


# -------------------------------- PART E --------------------------------
H_max <- 3.72  # Maximum wave height (m)
T <- 10  # Wave period (s)
g <- 9.81  # Gravity (m/s^2)
d <- 25  # Water depth (m)
D <- 1.1  # Pipe diameter (m)
rho <- 1025  # Seawater density (kg/m^3)
C_D <- 1  # Drag coefficient
C_m <- 2  # Inertia coefficient

omega <- 2 * pi / T 
lambda_0 <- (g * T^2) / (2 * pi)
k_guess <- 2 * pi / lambda_0  

dispersion_relation <- function(k, d, omega, g) {
  return(omega^2 - g * k * tanh(k * d))
}

tolerance <- 1e-6
max_iterations <- 100

for (i in 1:max_iterations) {
  f_k <- dispersion_relation(k_guess, d, omega, g)
  df_k <- -g * (tanh(k_guess * d) + k_guess * d * (1 / cosh(k_guess * d))^2)
  k_new <- k_guess - f_k / df_k
  
  if (abs(k_new - k_guess) < tolerance) {
    break
  }
  
  k_guess <- k_new
}

k <- k_new  
lambda_exact <- 2 * pi / k 

# Compute maximum orbital velocity u_m at pipe centerline (z = -d + D/2)
z <- -d + D/2
u_m <- (omega * H_max / 2) * (cosh(k * z) / sinh(k * d))

# Time vector for one wave period (0 ≤ ωt ≤ 2π)
t <- seq(0, T, length.out = 200)
u <- u_m * cos(omega * t)  # Orbital velocity (Table A.2)
u_dot <- -u_m * omega * sin(omega * t)  # Orbital acceleration (Table A.2)

#drag force F_D and inertia force F_I
F_D <- (0.5 * rho * D * C_D) * u * abs(u)
F_I <- (rho * (pi * D^2 / 4) * C_m) * u_dot
F_total <- F_D + F_I # total in-line force

# Create data frame for plotting
force_data <- data.frame(
  Time = rep(t, 3),
  Force = c(F_D, F_I, F_total),
  Component = rep(c("Drag Force", "Inertia Force", "Total Force"), each = length(t))
)

# Plot forces over one wave period
ggplot(force_data, aes(x = Time, y = Force, color = Component)) +
  geom_line(size = 1) +
  labs(title = "Time-Varying Wave Forces on Submarine Pipeline",
       x = "Time (s)",
       y = "Force per Unit Length (N/m)",
       color = "Force Component") +
  theme_minimal() +
  theme(legend.position = "bottom")



# -------------------------------- PART F --------------------------------
# Compute maximum orbital velocity u_m at pipe centerline (z = -d + D/2)
z <- -d + D/2
u_m <- (omega * H_max / 2) * (cosh(k * z) / sinh(k * d))

# Compute second-order correction factor
second_order_factor <- (3 * omega * pi * H_max^2) / (8 * lambda_exact * sinh(k * d)^4)

# Time vector for one wave period (0 ≤ ωt ≤ 2π)
t <- seq(0, T, length.out = 200)

# Compute linear velocity and acceleration
u_linear <- u_m * cos(omega * t)  
u_dot_linear <- -u_m * omega * sin(omega * t)  

# Compute second-order velocity and acceleration corrections
u_second_order <- second_order_factor * cosh(2 * k * (z + d)) * cos(2 * omega * t)
u_dot_second_order <- -(3 * omega^2 * pi * H_max^2) / (4 * lambda_exact * sinh(k * d)^4) * cosh(2 * k * (z + d)) * sin(2 * omega * t)

# Compute total velocity and acceleration (Stokes 2nd-order)
u_stokes <- u_linear + u_second_order
u_dot_stokes <- u_dot_linear + u_dot_second_order

# Compute drag force F_D using Stokes' velocity
F_D_stokes <- (0.5 * rho * D * C_D) * u_stokes * abs(u_stokes)

# Compute inertia force F_I using Stokes' acceleration
F_I_stokes <- (rho * (pi * D^2 / 4) * C_m) * u_dot_stokes

# Compute total in-line force
F_total_stokes <- F_D_stokes + F_I_stokes

# Create data frame for plotting
force_data_stokes <- data.frame(
  Time = rep(t, 3),
  Force = c(F_D_stokes, F_I_stokes, F_total_stokes),
  Component = rep(c("Drag Force", "Inertia Force", "Total Force"), each = length(t))
)

# Plot forces over one wave period
ggplot(force_data_stokes, aes(x = Time, y = Force, color = Component)) +
  geom_line(size = 1) +
  labs(title = "Time-Varying Wave Forces (Stokes 2nd-Order)",
       x = "Time (s)",
       y = "Force per Unit Length (N/m)",
       color = "Force Component") +
  theme_minimal() +
  theme(legend.position = "bottom")



