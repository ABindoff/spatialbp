# spatialbp

**2D Continuous Spatial Piecewise Regression**

`spatialbp` is an R package with Rust bindings that extends 1D piecewise regression into two dimensions. It fits 2D spatial surfaces stitched together at 1D boundary curves to detect spatial regimes, affinity, and repulsion in continuous spatial fields.

> [!WARNING]
> **Work In Progress**
> This project is currently under active development. The API is not yet stable and features may change without notice. Please use with caution in production environments.

## Installation

Because this package contains Rust extensions, you will need Rust installed on your system (via [rustup](https://rustup.rs/)).

Once Rust is installed, you can install the development version from GitHub using `remotes` or `devtools`:

```r
# install.packages("remotes")
remotes::install_github("ABindoff/spatialbp")
```

## Overview

The `spatialbp` package provides tools to discover spatial relationships and boundaries between interacting variables or species. By fitting piecewise surfaces across a 2D spatial domain, you can estimate abrupt transition zones (ecotones) and map out regions of competitive co-dependence or biological affinity.
