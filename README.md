# FEISTY R package

The FEISTY model (FishErIes Size and functional TYpe model) is documented as an R package. The package includes four FEISTY model setups: two published FEISTY model setups ([Petrik et al., 2019](https://doi.org/10.1016/j.pocean.2019.102124); [van Denderen et al., 2021](https://doi.org/10.1111/geb.13348)) and their modified versions. Also, it allows researchers to customize and simulate new FEISTY model setups, for model development and other marine science research.  You can test FEISTY in the [online simulator](http://oceanlife.dtuaqua.dk/FEISTY). 

<img width="806" height="527" alt="image" src="https://github.com/user-attachments/assets/946884ec-cd9a-43a9-a7b4-54279c30ef99" />

---

# FEISTY Installation Guide

FEISTY is an R package with compiled Fortran code. Most users only need to install and load the package. Developers who modify the source code need additional tools for documentation, testing, and package checks.

---

## For Users

Use this route if you want to run FEISTY simulations but do not plan to edit the package source code.

### Requirements

- R, preferably R 4.3.x
- R package `remotes`
- A working package-compilation toolchain:
  - Windows: Rtools matching your R version
  - Linux: `gfortran` and `make`
  - macOS: Xcode command line tools and `gfortran`

Install `remotes` in R:

```r
install.packages("remotes")
```

Install the latest development version:

```r
remotes::install_github("Kenhasteandersen/FEISTY", upgrade = "never")
library(FEISTY)
```

Install a stable release instead:

```r
remotes::install_url(
  "https://github.com/Kenhasteandersen/FEISTY/archive/refs/tags/v1.1.1.tar.gz",
  upgrade = "never"
)
library(FEISTY)
```

Test the installation:

```r
webFEISTY()
```

---

## For Developers

Use this route if you want to modify FEISTY source code, rebuild documentation, run tests, or contribute changes.

### Additional Requirements

- Git
- R packages `devtools`, `roxygen2`, `testthat`, `knitr`, and `rmarkdown`

Install the recommended developer packages:

```r
install.packages(c(
  "devtools",
  "roxygen2",
  "testthat",
  "knitr",
  "rmarkdown"
))
```

Clone the repository:

```bash
git clone https://github.com/Kenhasteandersen/FEISTY.git
cd FEISTY
```

Install the local source package from R:

```r
remotes::install_local(".", upgrade = "never")
library(FEISTY)
```

If you edit exported functions or roxygen comments, rebuild the documentation:

```r
devtools::document()
```

Run tests:

```r
devtools::test()
```

Run a package check:

```r
devtools::check()
```

---

## Platform Notes

### Windows

Install Rtools matching your R version. For example, use Rtools43 with R 4.3.x. Rtools provides the `make` and Fortran tools needed to install FEISTY from source.

### Linux

On Ubuntu/Debian:

```bash
sudo apt update
sudo apt install gfortran make
```

### macOS

Install Xcode command line tools:

```bash
xcode-select --install
```

Install `gfortran`, for example with Homebrew:

```bash
brew install gcc
```

If R cannot find `gfortran`, check its location:

```bash
which gfortran
```

and configure `~/.R/Makevars` if needed:

```bash
mkdir -p ~/.R
nano ~/.R/Makevars
```

Add the compiler path returned by `which gfortran`. For example:

```make
FC=<gfortran_location>/gfortran
F77=<gfortran_location>/gfortran
FLIBS=-L<gfortran_location>/lib
```

Replace `<gfortran_location>` with the actual directory on your system. This step is only needed if R cannot locate the compiler automatically.

---
