library(spatialbp)
library(ggplot2)

# Load and parse NIfTI
message("Downloading MNI structural template from CRAN...")
db <- available.packages(repos = "https://cloud.r-project.org")
ver <- db["oro.nifti", "Version"]
package_url <- paste0("https://cloud.r-project.org/src/contrib/oro.nifti_", ver, ".tar.gz")
temp_file <- tempfile(fileext = ".tar.gz")
download.file(package_url, temp_file, mode = "wb", quiet = TRUE)
untar(temp_file, files = "oro.nifti/inst/nifti/mniRL.nii.gz", exdir = tempdir())
nii_path <- file.path(tempdir(), "oro.nifti/inst/nifti/mniRL.nii.gz")

con <- gzfile(nii_path, "rb")
header <- readBin(con, raw(), 352)
dims <- readBin(header[41:56], integer(), size = 2, n = 8, endian = "swap")
nx <- dims[2]; ny <- dims[3]; nz <- dims[4]
datatype <- readBin(header[71:72], integer(), size = 2, n = 1, endian = "swap")
voxels_flat <- readBin(con, integer(), n = nx * ny * nz, size = 2, signed = FALSE, endian = "swap")
close(con)
voxels_3d <- array(voxels_flat, dim = c(nx, ny, nz))

# Extract ROI
x_range <- 35:55
y_range <- 45:65
z_range <- 43:47
grid_3d <- expand.grid(x = x_range, y = y_range, z = z_range)
grid_3d$T1_MRI <- apply(grid_3d, 1, function(row) {
  voxels_3d[row["x"], row["y"], row["z"]]
})

grid_3d$x <- grid_3d$x - min(grid_3d$x)
grid_3d$y <- grid_3d$y - min(grid_3d$y)

# 1. Define an elongated, organic ventricular core skeleton
centers <- matrix(c(
  10, 5,
  11, 8,
  12, 11,
  11, 14,
  10, 17
), ncol = 2, byrow = TRUE)

# Calculate minimum distance to this elongated ventricular core
dists <- apply(grid_3d, 1, function(row) {
  min(sqrt((row["x"] - centers[, 1])^2 + (row["y"] - centers[, 2])^2))
})

# Construct a spatially continuous organic latent field
grid_3d$Latent_WM <- 10 - dists

# 2. Simulate corresponding DTI FA channel matching the real T1 structures
set.seed(42)
grid_3d$DTI_FA <- 0.15 + 0.55 * (grid_3d$T1_MRI / max(grid_3d$T1_MRI)) + rnorm(nrow(grid_3d), 0, 0.05)

# Fit the boundary
fit <- fit_spatial_boundary(
  df = grid_3d,
  y1 = "T1_MRI",
  y2 = "DTI_FA",
  latent_col = "Latent_WM",
  iterations = 3000
)

# Plot central slice
p <- plot(fit, z_slice = 3)
ggsave("C:/Users/bindoffa/.gemini/antigravity/brain/93f67760-4ac4-41dc-91ef-e9758c4812b2/artifacts/mri_organic_boundary.png", p, width = 6, height = 6, dpi = 150)
message("Saved mri_organic_boundary.png successfully!")
