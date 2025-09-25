devtools::document()
devtools::load_all()
?predict_RECKMON

set.seed(123)
n <- 100
p <- 10
x_data <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(x_data) <- paste0("V", 1:p)

# Create a response variable y
y_data <- 2 * x_data[, 1] - 1.5 * x_data[, 2] +      # Linear effects
          x_data[, 3] * x_data[, 4] +               # Interaction effect
          sin(pi * x_data[, 5]) +                   # Non-linear effect
          rnorm(n, 0, 0.5)                          # Noise (CORRECTED: using 'n')

# Run the model with parameters suitable for a quick example
model_fit <- main_RECKMON(                          # CORRECTED: Function name
  x = x_data,
  y = y_data,
  poly_features = 5,
  step_gaussian = 1,
  cv = 3, # Use 3 folds for a quick example
  mc.cores = 1
)


# View the selected features
cat("Selected interaction features:\n")
print(model_fit$poly_feature)

cat("Selected Gaussian kernel features:\n")
print(model_fit$gaussian_feature)

# Predict new observation
predict_RECKMON(newx = matrix( rnorm(p*2), nrow = 2), model_fit)