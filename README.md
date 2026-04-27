# FEISTY R package

The FEISTY model (FishErIes Size and functional TYpe model) is documented as an R package. The package includes four FEISTY model setups: two published FEISTY model setups ([Petrik et al., 2019](https://doi.org/10.1016/j.pocean.2019.102124); [van Denderen et al., 2021](https://doi.org/10.1111/geb.13348)) and their modified versions. Also, it allows researchers to customize and simulate new FEISTY model setups, for model development and other marine science research.  You can test FEISTY in the [online simulator](http://oceanlife.dtuaqua.dk/FEISTY). 

<img width="806" height="527" alt="image" src="https://github.com/user-attachments/assets/946884ec-cd9a-43a9-a7b4-54279c30ef99" />

---

# FEISTY Installation Guide

The following provides a comprehensive guide to installing and setting up FEISTY, a Fortran-based marine ecosystem model, on Linux, Windows, and macOS systems.

Before installing FEISTY, ensure that the following tools are installed on your system:

- **gfortran**: The GNU Fortran compiler.
- **make**: A build automation tool.
- **git**: A version control system.

If these tools are not already installed, follow the instructions below for your operating system to install them.

---

## Installation by Operating System

### 1. Linux

1. **Install Dependencies**  
   Open a terminal and run the following commands to install the necessary tools:
   ```bash
   sudo apt update
   sudo apt install gfortran make git
   ```

---

### 2. Windows

1. **Install Dependencies**  
   - Download and install [MinGW-w64](https://sourceforge.net/projects/mingw/), which includes `gfortran` and `make`.
   - Download and install [Git for Windows](https://git-scm.com/).

   During the MinGW-w64 installation, ensure `gfortran` and `make` are selected.

2. **Add to path**  
Add the MinGW-w64 and make installation directory to your system's PATH environment variable. 
To do so go to: Environenment variables > user variables > path > edit; then paste `C:\MinGW\bin`


---

### 3. macOS

1. **Install Dependencies**  
   Use Homebrew to install the required tools. If you do not have Homebrew installed, visit [brew.sh](https://brew.sh/) to set it up. Then, run:
   ```bash
   brew install gcc make git
   ```

2. **Error compilation failed**  
   If you get an error like:
   ``make: /opt/gfortran/bin/gfortran: no such file or directory [...] compilation failed for package 'FEISTY' `` 
   this is due to make trying to search for the fortran compiler in the wrong path. To solve the problem you need to 
   add a Makevars file with the right directoy for the fortran compiler and the libraries, this directory depends on where 
   Homebrew installed gfortran.
    ```bash
   cd ~/
   mkdir .R
   cd .R/

   ```
   Find gfortran location:
   ```
   which gfortran
   <gfortran_location>/gfortran
   ```
   Then add the path to the Makevars file
   ```
   nano Makevars
   ```
   copy in Makevars the following by changing <gfortran_location> by the path you obtained from `which gfortran`
   ```
   FC=<gfortran_location>/gfortran
   F77=<gfortran_location>/gfortran
   FLIBS=-L<gfortran_location>/lib
   ```
---

## Download and build FEISTY

Clone FEISTY from the Rstudio console and load the FEISTY library.

The latest development version:
```bash
   remotes::install_github("Kenhasteandersen/FEISTY")
   library(FEISTY)
```

The latest stable version:
```bash
   remotes::install_url("https://github.com/Kenhasteandersen/FEISTY/archive/refs/tags/v1.0.0.tar.gz")
   library(FEISTY)
```
You should get something like: 
```
   ==> Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source FEISTY
   
   * installing to library 'C:/Users/rdenechere/AppData/Local/R/win-library/4.3'
   * installing *source* package 'FEISTY' ...
   ...
   * DONE (FEISTY)
```
Now you can try to run the FEISTY web app by typing in the console:
```
webFEISTY()
```

---
