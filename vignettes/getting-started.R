## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval = FALSE----------------------------------------------------
# # Install from GitHub
# # remotes::install_github("YourUsername/DynamicAlphaHull")
# 
# # Load the package
# library(DynamicAlphaHull)

## ----load, echo = FALSE-------------------------------------------------------
library(DynamicAlphaHull)

## ----quickstart---------------------------------------------------------------
# Create a small occurrence dataset
points <- data.frame(
  Longitude = c(116.0, 116.5, 116.4, 116.1, 115.8, 116.6),
  Latitude = c(39.0, 39.0, 39.4, 39.3, 39.1, 39.2)
)

# Build a range polygon
result <- getDynamicAlphaHull(
  points,
  buff = 10000,        # 10km buffer
  clipToCoast = "no"   # don't clip to coastline
)

# Plot the result
plot(result$hull, col = "lightblue", border = "blue")
points(points$Longitude, points$Latitude, pch = 16, cex = 1.5)

## ----quickstart_inspect-------------------------------------------------------
cat("Selected alpha:", result$alpha, "\n")
cat("Range is a", class(result$hull)[1], "\n")

## ----coastline_demo, eval = FALSE---------------------------------------------
# # Terrestrial species - clip to land
# terrestrial_range <- getDynamicAlphaHull(
#   points,
#   clipToCoast = "terrestrial"
# )

## ----example_data-------------------------------------------------------------
data(rosales_example)
head(rosales_example)

