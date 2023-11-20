# Repository for NEURON Model of Striatal Projection Neurons

This repository contains a NEURON + Python model of striatal projection neurons (or SPNs) designed to simulate the interaction between GABAergic and glutamatergic synaptic inputs. 

It also provides all the R code used to produce the graph output (in particular as `svg`) from the resultant NEURON + Python output. 

The NEURON + Python model is built on top of the 'striatal_SPN_lib' repository created by Lindroos and Kotaleski, 2020:

[Lindroos R, Kotaleski JH. Predicting complex spikes in striatal projection neurons of the direct pathway following neuromodulation by acetylcholine and dopamine. Eur J Neurosci. 2020](https://doi.org/10.1111/ejn.14891)

The original model can be found here [modelDB](https://senselab.med.yale.edu/ModelDB/ShowModel?model=266775&file=/lib/params_dMSN.json#tabs-2) or [GitHub](https://github.com/ModelDBRepository/266775)

## Table of Contents
- [Initial Set Up](#Initial-Set-Up)
- [Running the Models](#running-the-models)
  - [Getting Started](#getting-started)
  - [Simulations](#running-simulations-in-jupyter-notebook)
- [Data Analysis](#data-analysis)
  - [Setting up](#setting-up)
  - [Using R to analyse a simulation](#using-r-to-analyse-a-simulation)
- [Anaconda vs Miniconda](#anaconda-vs-miniconda)
- [Virtual Environments](#virtual-environments)
- [GitHub](#using-github)
- [References](#references)
- [Contact](#contact)

## Initial Set Up

### Prerequisites
- Conda [Conda offical website](https://docs.conda.io/projects/conda/en/stable)
- NEURON (tested on versions 8.1 and 8.2) [NEURON official website](https://www.neuron.yale.edu)

  [Hines ML, Carnevale NT. The NEURON Simulation Environment. Neural Comput. 1997](https://doi.org/10.1162/neco.1997.9.6.1179)
- Python (tested using version 3.9.16)

### Steps
1. **Install [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)** (Python package manager)
   
   Conda should include a version of Python. 

   Conda is simply a package manager and environment management system that is used to install, run and update packages and their dependencies.

   Whether to install [Anaconda or Miniconda](#anaconda-vs-miniconda) largely depends on personal preferences.

3. **Install [Jupyter Notebook](https://jupyter.org)**

   The simplest method to install Jupyter Notebook is via Conda using the command in Terminal:

   ```bash
   conda install -c conda-forge notebook
   ```

4. **Install [NEURON](https://www.neuron.yale.edu)**

  Follow the guide at [NEURON Installation](https://www.neuron.yale.edu/neuron/download)

* [Quickstart](https://www.neuron.yale.edu/ftp/neuron/2019umn/neuron-quickstart.pdf)

* [(Old) Mac instructions](https://www.neuron.yale.edu/ftp/neuron/nrn_mac_install.pdf)
  

## Running the Models

The following sections explain the initial set up required and instructions to create simulations subsequently used to generate figures.

### Getting Started

1. **Install NEURON with Python support** (see setup instructions)

2. **Open Terminal**:
   - On MacOS: Press `cmd + space` to open spotlight search and type `terminal`.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for `command prompt` or `PowerShell` in the start menu.

3. **Compile mechanisms**:

   Navigate to directory containing NEURON mechanisms.

   For instance if `SPNfinal` is in documents folder in MacOS, mechanisms are located in `cd/documents/SPNfinal/mechanisms/single`.

   Mechanisms are then compiled by entering `nrnivmodl` or `mknrndll`.

   On MacOS / Linux:
   ```bash
   cd documents/Repositories/SPNfinal/mechanisms/single
   nrnivmodl
   ```

   On Windows:
   ```bash
   cd C:\Users\YourUsername\documents\Repositories\SPNfinal\mechanisms\single
   nrnivmodl
   ```

4. **Create a conda environment**

   There is a `YAML` file in the main directory called `environment.yml` for MacOS/Linux. This can be used to create a conda environment called `neuron`. For further information see [Virtual Environments](#virtual-environments).

   Ensure make sure to navigate back to the main directory after step 3 above.

   Check installed correctly using 'conda list'.

   On MacOS / Linux:
   ```bash
   cd ../.. 
   conda env create -f environment.yml
   conda list
   ```

   On Windows:
   ```bash
   cd ..\.. 
   conda env create -f environment_pc.yml
   conda list
   ```
  
   Creating the environment on Windows is slightly different. `NEURON` cannot be installed via the terminal using `pip install neuron`. 
  
   Instead, `NEURON` must be installed via a downloaded setup.exe. 
  
   The `NEURON` version in   `environment_pc.yml` must match the installed version (for instance, the Windows laptop used for testing had `NEURON 8.2.0` installed). The `YAML` file may need editing in a text editor.
  
   As a result, a separate `environment_pc.yml` is provided for Windows. Limited testing on a Windows laptop showed the simulations working but suggested that, despite reasonable specs, `Intel(R) Core(TM) i5-8350U CPU @ 1.7 GHz 1.9 GHz 32GB`, the code ran extremely slowly in the Windows environment (~ 9-fold slower) when compared to a `MacBook M2 pro 32GB`. In fact, it was slower (~ 4-fold) than a 2015 `MacBook Pro 2.7 GHz Dual-Core Intel core i5`.     
     
   
5. **Quit Terminal**

   ```bash
   exit
   ```
   
### Running simulations in Jupyter Notebook

  The following steps 1-4 must be performed every time a new Jupyter Notebook session is started.

1. **Open Terminal**:
   - On MacOS: Press `cmd + space` to open spotlight search and type `terminal`.
   - On Linux: Search for `terminal` in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for `command prompt`, `PowerShell` or the appropriate `Miniconda/Anaconda prompt` in the start menu.

2. **Activate conda environment 'neuron'**

   Navigate back to the main directory

   On MacOS / Linux:
   ```bash
   cd documents/Repositories/SPNfinal
   conda activate neuron
   ```

   On Windows:
   ```bash
   cd C:\Users\YourUsername\Repositories\documents\SPNfinal
   conda activate neuron
   ```
3. **Run Jupyter notebook**

   Add `neuron` environment then open Jupyter Notebook
   ```bash
   python -m ipykernel install --user --name neuron --display-name "Python (neuron)"
   jupyter notebook
   ```

4. **Run a simulation**

   Jupyter Notebook should now be open in the default browser.

   Choose a Notebook to open (by clicking on any notebook - `*.ipynb`).
   
   Ensure kernel is set to Python 3 (ipykernel).

   From kernel dropdown menu choose `Restart Run All` (if running again then it's good practice to run `Restart and Clear Output` first).

   Code should run and generate raw data used to generate figures.

   If option `save = True` in the Notebook then the raw figures and pickled data is stored in a subdirectory within the main one.


## Data Analysis

The final analysis and figures presented in the manuscript were generated using R. 

The analyses were conducted in the R graphical user interface (GUI):
  - R version 4.3.1 – "Beagle Scouts"
  - [R Statistical Software](https://www.R-project.org/)

  Refer to the `R analysis` directory for the code.

  ### Setting up
  
  Only the R console was used for analysis. 
  
  If you prefer to work with `RStudio`, it can be downloaded [here](https://posit.co/products/open-source/rstudio/). The provided code should work although this has not been tested.
  
  In order for the R code to work, it is necessary to load various packages within the R environment.
  
  The following code should be executed in R prior to running any of this code. It checks if the required packages are present and, if they are not, it will install them.
  
  In addition, the first time the code is executed, it will install a Miniconda environment using the `reticulate` package in R: 

  The following steps 1-3 must be performed once.
  
  1. **Open R gui**
  2. **Run this code once**
  ```R
  rm( list=ls(all=TRUE ) )
  # Load and install necessary packages
  install_required_packages <- function(packages){
      new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
      if(length(new.packages)) install.packages(new.packages)    
  }
  
  required.packages <- c('reticulate', 'stringr', 'svglite', 'quantmod', 'xts', 'zoo')
  install_required_packages(required.packages)
  
  # Install miniconda if necessary
  if (!reticulate::py_available(initialize = TRUE)) {
      # Install Miniconda through Reticulate (if Miniconda is not already installed)
      reticulate::install_miniconda()
  
      # Load the reticulate library
      library(reticulate)
  
      # Create and activate a new Conda environment
      conda_create(envname = "myenv")
      use_condaenv("myenv", required = TRUE)
  
      # Install pandas in the new environment
      py_install("pandas", envname = "myenv")
  }  
  ```
  3. **Exit R session**

  **Once this code is run, it will do the necessary installations and load the necessary packages for the analysis.**
  
  A useful guide can be found [here](https://rstudio.github.io/reticulate/).
  
  ### Using R to analyse a simulation
  
  To run the analysis code, simply execute all the code for that particular R analysis file after having run the relevant `*.ipynb` file.
  
  Each simulation has a unique identifier in its `*.ipynb`; for instance, `Fig5_EF.ipynb` is `sim4`. 
  
  Once the Jupyter Notebook is executed with `save = True`, the outputs are stored automatically in a subfolder using the identifier. 
  
  In this case, raw trace data is stored as pickled files in the subdirectory `dspn/model1/physiological/simulations/sim4`. 
  
  Any images generated are found in `dspn/model1/physiological/images/sim4`. 
  
  The R code to analyse the output from `Fig5_EF.ipynb` is found in `Fig5_EF.R` in the `R analysis` directory. 
  
  1. **Open R gui**
  2. **Open `Fig5_EF.R` and run the code**
  
  ```R
  rm( list=ls(all=TRUE ) )
  # Load and install necessary packages
  load_required_packages <- function(packages){
      new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
      if(length(new.packages)) install.packages(new.packages)
      invisible(lapply(packages, library, character.only = TRUE))
  }
  
  required.packages <- c('reticulate', 'stringr', 'svglite', 'quantmod', 'xts', 'zoo')
  load_required_packages(required.packages)
  
  use_condaenv("myenv", required = TRUE)
  # Import pandas
  pd <- reticulate::import("pandas")
  
  # Metadata settings
  plotsave <- TRUE
  spn <- 'dspn'
  sim <- 4
  model <- 1
  path <- paste0('/Documents/Repositories/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')
  
  # Set working directory
  cd <- sub("/Documents.*", "", getwd())
  wd <- paste0(cd, path)
  setwd(wd)
  
  # Extract folder data
  folders1 <- dir()
  folders <- folders1[!grepl('GABA0|GLUT0', folders1)]
  
  # function to read data
  read_data_from_folders <- function(folders1, wd) {
  	folders <- folders1[!grepl('GABA0|GLUT0', folders1)]
    	Vtree.path <- list()
    	Vsoma <- list()
    	n.folders <- length(folders)
    
    	for (iii in 1:n.folders) {
      	XX <- folders[iii]
      	setwd(XX)
      
      	V_soma <- pd$read_pickle(paste0(XX, '_v_data.pickle'))
      	V_dend <- pd$read_pickle(paste0(XX, '_vdend_data.pickle'))
      	V_tree <- pd$read_pickle(paste0(XX, '_v_tree.pickle'))
      	time <- pd$read_pickle(paste0(XX, '_t_data.pickle'))
      	m <- length(time)
      	col.names <- colnames(V_tree[[1]])
      
      	if (iii == 1) {
        		t <- time[[1]]
        		gt <- rep(0, n.folders)
        		n <- length(time[[1]])
        		mat <- matrix(rep(0, n.folders*n), ncol=n.folders)
        		Vsoma <- Vdend <- rep(list(mat), m)
      	}
      
      	gt[iii] <- pd$read_pickle(paste0(XX, '_gt.pickle'))[[1]][1]
      	Vtree.path[[iii]] <- V_tree 
      	for (ii in 1:m) {
        		Vsoma[[ii]][,iii] <- V_soma[[ii]]
        		Vdend[[ii]][,iii] <- V_dend[[ii]]
      	}
      	setwd(wd)
    	}
    
    	if (sum(grepl("GABA0", folders1)) == 1) { # contains 0 GABA folder
      	XX <- folders[grepl("GABA0", folders1)]
      	setwd(XX)
      	Vsoma_GABA0 <- pd$read_pickle(paste0(XX, '_v_data.pickle'))
      	Vdend_GABA0 <- pd$read_pickle(paste0(XX, '_vdend_data.pickle'))
      	i_mechs_GABA0 <- pd$read_pickle(paste0(XX, '_i_mechanisms.pickle'))
      	setwd(wd)
      	# Add this to return data for GABA0 as well
      	return(list(Vtree.path = Vtree.path, Vsoma = Vsoma, Vdend = Vdend, gt = gt, t = t, m = m,  n.folders = n.folders, Vsoma_GABA0 =      Vsoma_GABA0, Vdend_GABA0 = Vdend_GABA0, i_mechs_GABA0 = i_mechs_GABA0))
    	}
    	# Return data
    	return(list(Vtree.path = Vtree.path, Vsoma = Vsoma, Vdend = Vdend, gt = gt, t = t, m = m, n.folders = n.folders))
  }

# read data
data <- read_data_from_folders(folders1, wd)
list2env(data, envir = .GlobalEnv)

# Clean environment
to_keep <- c('sim', 'n.folders', 't', 'm', 'gt','Vdend','Vsoma','Vsoma_GABA0','Vdend_GABA0', 'Vtree.path', 'plotsave', 'wd')
to_remove <- setdiff(ls(), to_keep)
rm(list = to_remove)

# Generate plots

# Sort and reorder based on 'gt' values
gt <- unlist(gt)
ind <- sort(gt, index.return=TRUE)$ix
gt <- gt[ind]
Vsoma <- lapply(1:m, function(ii) Vsoma[[ii]][, ind])
Vdend <- lapply(1:m, function(ii) Vdend[[ii]][, ind])

# Color palette function
fun_color_range <- colorRampPalette(c('slateblue', 'indianred', 'aquamarine4', 'darkorchid', 'sienna1'))

# Plotting function for graphical representation
grph1a <- function(lst, t, x.lim, y.lim, main = '', tcl = -0.3, lwd = 1, xaxis = TRUE, y.int = 10) {
  	n.folders <- dim(lst)[2]
  	color_range <- fun_color_range(n.folders-1)
  
  	ind <- sapply(1:length(x.lim), function(ii) which.min(abs(t - x.lim[ii])))
  
  	if (xaxis) x.lab <- 'time (ms)' else x.lab <- ''
  
  	for (ii in 1:n.folders) {
    	if (ii == 1) {
      		plot(t[ind[1]:ind[2]], lst[ind[1]:ind[2], ii], xlab = x.lab, ylim = y.lim, xlim = x.lim, ylab = 'mV', main = main, type = 'l', bty = 'n', col = 'gray', axes = FALSE, lwd = 2 * lwd, lty = 1, frame = FALSE)
      		axis(2, at = seq(y.lim[1], y.lim[2], y.int), labels = seq(y.lim[1], y.lim[2], y.int), 
           	las = 1, tcl = tcl, lwd = lwd)
    	} else {
      		lines(t[ind[1]:ind[2]], lst[ind[1]:ind[2], ii], col = color_range[ii-1], lwd = lwd, lty = 1)
    	}
 	}
  	if (xaxis) axis(1, at = seq(x.lim[1], x.lim[2], 50), labels = seq(x.lim[1], x.lim[2], 50) - x.lim[1], las = 1, tcl = tcl, lwd = lwd)
}

# Main plotting script
x.lim <- c(125, 375)
y.lim1 <- c(-90, -60)
y.lim2 <- c(-100, -20)
lwd <- 1

Vdend.plot <- cbind(Vdend_GABA0[[1]], Vdend[[1]])
Vsoma.plot <- cbind(Vsoma_GABA0[[1]], Vsoma[[1]])

# Create plots
dev.new(width = 6, height = 10 / 3, noRStudioGD = TRUE)
par(mar = c(1, 1, 1, 1), mfrow = c(1, 2), oma = c(2, 2, 2, 0), ps = 10, cex = 0.9, cex.main = 0.9)
grph1a(Vdend.plot, t, x.lim, y.lim2, main = 'dendritic voltage', lwd = lwd, xaxis = TRUE, y.int = 20)
grph1a(Vsoma.plot, t, x.lim, y.lim1, main = 'somatic voltage', lwd = lwd, xaxis = TRUE, y.int = 10)

if (plotsave) {	
    # Determine the save directory
    save_dir <- str_replace(wd, 'simulations', 'images')
    setwd(save_dir)
    
    # Create a timestamp string for file names
    timestamp <- gsub(':', '-', Sys.time())
    
    # Save the combined plots as a PDF for records
    quartz.save(paste0('sim', sim, ' ', timestamp, '.pdf'), type='pdf')
    
    # Save dendritic voltage plot as an SVG
    svglite(paste0('sim', sim, 'a', ' ', timestamp, '.svg'), width=3.5, height=3.75, pointsize=10)
    grph1a(Vdend.plot, t, x.lim, y.lim2, main='dendritic voltage', lwd=lwd, xaxis=TRUE, y.int=20)
    # Add reference lines to the plot
    abline(h=-20, lty=3, lwd=lwd)
    abline(h=-85, lty=3, lwd=lwd)
    dev.off()
    
    # Save somatic voltage plot as an SVG
    svglite(paste0('sim', sim, 'b', ' ',  timestamp, '.svg'), width=3.5, height=3.75, pointsize=10)
    grph1a(Vsoma.plot, t, x.lim, y.lim1, main='somatic voltage', lwd=lwd, xaxis=TRUE, y.int=10)
    # Add reference lines to the plot
    abline(h=-60, lty=3, lwd=lwd)
    abline(h=-85, lty=3, lwd=lwd)
    dev.off()
    
    # Return to original working directory
    setwd(wd)
}
```
  
  In order to locate the raw data to recreate the final images for the ms it may be necessary to alter the line (depending on where the original directory was created):

  On MacOS/Linux: 
  
  ```R
  path <- paste0('/Documents/Repositories/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')
  ```
 
  On Windows: 
  
  ```R
    path <- paste0('C:\\Users\\YourUsername\\Documents\\Repositories\\SPNfinal\\', spn, '\\model', model, '\\physiological\\simulations\\sim', sim, '\\')
  ```

  If running code on MacOS/Linux, this line should be the only one that it is necessary to change in order to execute the R code.
  On Windows, other lines would need to be altered:
  
  ```R
  # replace quartz.save(paste0('sim', sim, ' ', timestamp, '.pdf'), type='pdf')
  pdf(paste0('sim', sim, ' ', timestamp, '.pdf')) # replace quartz.save(paste0('sim', sim, ' ', timestamp, '.pdf'), type='pdf') 
  # ...plotting code...
  dev.off()
  ```
  ### Please note:
  
  The code provided in `Fig7_BD2.R` gives the graphical output for `Fig7_BD2.ipynb`. However, because I forgot to run gbar_gaba = 0 in `Fig7_BD2.ipynb` which is necessary for the R analysis,  the R code will retrieve this from  `Fig7_BD.ipynb`. This means both `Fig7_BD.ipynb` and `Fig7_BD2.ipynb` need to be run prior to analysis with `Fig7_BD2.R`. A similar arrangement applies to `Fig7_E.ipynb`, `Fig7_E2.ipynb` and `Fig7_E2.R`.
 
  ### Key Points for Windows:
  - The path should be changed to reflect the typical Windows file structure. Make sure to replace `YourUsername` with your actual username.
  - File path separators are changed to `\\`.
  - The script assumes the R working directory aligns with the structure of the path. Adjust as needed.
  - Replace any MacOS-specific functions (like `quartz.save`) with their Windows-compatible equivalents.

## Anaconda vs Miniconda

Anaconda and Miniconda are both popular distributions for Python and R programming in data science. They include the Conda package manager and aim to simplify package management and deployment.

**Anaconda** is a full-featured distribution that includes:

- Python and R language
- Conda package manager
- Over 1,500 pre-installed scientific packages
- Tools like Jupyter Notebook, RStudio, etc.

Anaconda provides an out-of-the-box setup for data science and scientific computing.

**Miniconda** offers a minimalistic approach:

- Python and R language
- Conda package manager
- No pre-installed packages

Miniconda provides a lightweight base to start with and but packages must be installed if needed.

**Advantages of Anaconda**:

- Quick, easy setup with a comprehensive suite of scientific packages and tools.
- Wide array of data science tools readily available within a single application.

**Advantages of Miniconda**:

- Lightweight, minimal base installation.
- Control over which packages are installed.
- Requires limited disk space or bandwidth.
- Clean environment that only includes packages required.

**Installation**

- For Anaconda, download the installer from the [Anaconda download page](https://www.anaconda.com/products/individual#Downloads).
- For Miniconda, visit the [Miniconda download page](https://docs.conda.io/en/latest/miniconda.html).

Follow the installation instructions provided on the respective download pages.

**Additional Resources**

- [Anaconda Documentation](https://docs.anaconda.com/)
- [Miniconda Documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Conda Package Management](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html)


## Virtual Environments

The `environment.yml` file is configured for NEURON 8.2.2 and Python 3.9.16. Use this file to create a consistent environment for running the models.

In brief:     

* **`YAML` Environment File**: The file `environment.yml` is a `YAML` file commonly used in Conda

  This file specifies the dependencies and settings for a particular virtual environment.

* **Environment Name** - name Key: In the `environment.yml` file, there is a key called name.

  The value associated with this key is the name of the Conda environment to be created.

  In this case, this name is `neuron`. This is the name that will subsequently be used

  to refer to the environment when activating it or installing additional packages into it.

* **Creating the Environment**: When the command `conda env create -f environment.yml` is executed

  Conda reads the `environment.yml` file and creates a new environment based on the specifications in that file.

  The environment will have the name given by the name key in the `YAML` file i.e. `neuron`.
  

## GitHub

For beginners, the [GitHub Desktop GUI](https://desktop.github.com/) is recommended. 

Instructions for cloning a repository using GitHub Desktop can be found [here](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop).

## References

[Du K, Wu Y-W, Lindroos R, Liu Y, Rózsa B, Katona G, et al. Cell-type–specific inhibition of the dendritic plateau potential in striatal spiny projection neurons. Proceedings of the National Academy of Sciences. 2017;114: E7612 E7621](https://doi.org/10.1073/pnas.1704893114)

[Hines ML, Carnevale NT. The NEURON Simulation Environment. Neural Comput. 1997;9: 1179–1209](https://doi.org/10.1162/neco.1997.9.6.1179)

[Lindroos R, Dorst MC, Du K, Filipović M, Keller D, Ketzef M, et al. Basal Ganglia Neuromodulation Over Multiple Temporal and Structural Scales-Simulations of Direct Pathway MSNs Investigate the Fast Onset of Dopaminergic Effects and Predict the Role of Kv4.2. Frontiers in neural circuits. 2018;12: 3](https://doi.org/10.3389/fncir.2018.00003) 

[Lindroos R, Kotaleski JH. Predicting complex spikes in striatal projection neurons of the direct pathway following neuromodulation by acetylcholine and dopamine. Eur J Neurosci. 2020](https://doi.org/10.1111/ejn.14891)



## Contact

The model was adapted from the publicly available code by [Clay Surmeier](mailto:clay-surmeier@northwestern.edu) and [Vernon Clarke](mailto:vernon.clarke@northwestern.edu).

The provided code was executed on a `MacBook M2 pro 32GB`. We have tried to ensure that the code works on other operating systems but it's inevitable that some errors and bugs exist. 

In order to make this code accessible for publication, it is necessary to create a permanent public repository with a citable DOI using, for example using `Zenodo` to archive a version of this `GitHub` package.

If any bug fixes are necessary (most likely related to providing help on other operating systems), it will be provided as an update on the parent [`GitHub` page](https://github.com/vernonclarke/SPNfinal).

For queries related to this repository, please [open an issue](https://github.com/vernonclarke/SPNfinal/issues) or [email](mailto:vernon.clarke@northwestern.edu) directly 

---


   

