# Repository for NEURON Model of Striatal Projection Neurons

This repository contains a NEURON model of striatal projection neurons (or SPNs) designed to simulate the interaction between GABAergic and glutamatergic synaptic inputs. 

It is built on top of the 'striatal_SPN_lib' repository created by Lindroos and Kotaleski, 2020:

[Lindroos R, Kotaleski JH. Predicting complex spikes in striatal projection neurons of the direct pathway following neuromodulation by acetylcholine and dopamine. Eur J Neurosci. 2020](https://doi.org/10.1111/ejn.14891)


The original model can be found here [modelDB](https://senselab.med.yale.edu/ModelDB/ShowModel?model=266775&file=/lib/params_dMSN.json#tabs-2) or [GitHub](https://github.com/ModelDBRepository/266775)

## Table of Contents
- [Initial Set Up](#Initial-Set-Up)
- [Running the Models](#running-the-models)
  - [Getting Started](#getting-started)
  - [Simulations](#running-simulations-in-jupyter-notebook)
- [Data Analysis](#data-analysis)
- [Anaconda vs Miniconda](#anaconda-vs-miniconda)
- [Virtual Environments](#virtual-environments)
- [GitHub](#using-github)
- [Contact](#contact)

## Initial Set Up

### Prerequisites
- Conda [Conda offical website](https://docs.conda.io/projects/conda/en/stable)
- NEURON (tested on versions 8.1 and 8.2) [NEURON official website](https://www.neuron.yale.edu)

  [Hines ML, Carnevale NT. The NEURON Simulation Environment. Neural Comput. 1997](https://doi:10.1162/neco.1997.9.6.1179)
- Python (tested using version 3.9.16)

### Steps
1. **Install [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)** (python package manager)
   
   Installing Conda should include a version of Python. 

   Conda is simply a package manager and environment management system that is used to install, run, and update packages and their dependencies.

   Whether to install [Anaconda or Miniconda](#anaconda-vs-miniconda) largely depends on personal preferences.

3. **Install [Jupyter Notebook](https://jupyter.org)**

   The simplest method is to install Jupyter Notebook is via Conda using the command in Terminal:

   ```bash
   conda install -c conda-forge notebook
   ```

   [For further info](https://jupyter.org)

5. **Install [NEURON](https://www.neuron.yale.edu)**

   doi:10.1162/neco.1997.9.6.1179

  Follow the guide at [NEURON Installation](https://www.neuron.yale.edu/neuron/download)

* [Quickstart](https://www.neuron.yale.edu/ftp/neuron/2019umn/neuron-quickstart.pdf)

* [(Old) Mac instructions](https://www.neuron.yale.edu/ftp/neuron/nrn_mac_install.pdf)
  

## Running the Models

The following sections explain the inital set up required and instructions to create simulations subsequently used to generate figures.

### Getting Started

1. **Install NEURON with Python support** (see setup instructions)

2. **Open Terminal**:
   - On MacOS: Press `cmd + space` to open spotlight search and type 'terminal'.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for command prompt or PowerShell in the start menu.

3. **Compile mechanisms**:

   Navigate to directory containing NEURON mechanisms.

   For instance if `SPNfinal` is in documents folder in OS, mechanisms are located in `cd/documents/SPNfinal/mechanisms/single`

   Mechanisms are then compiled by entering `nrnivmodl` or `mknrndll`

   On MacOS / Linux:
   ```bash
   cd documents/SPNfinal/mechanisms/single
   nrnivmodl
   ```

   On Windows:
   ```bash
   cd C:\Users\YourUsername\documents\SPNfinal\mechanisms\single
   nrnivmodl
   ```

5. **Create a conda environment**

   There is a yaml file in the main directory called `environment.yml`. This can be used to create a conda environment called `neuron`. For further information see [Virtual Environments](#virtual-environments).

   Ensure make sure to navigate back to the main directory after step 3 above

   Check installed correctly using 'conda list'

   On MacOS / Linux:
   ```bash
   cd ../.. 
   conda env create -f environment.yml
   conda list
   ```
   On Windows:
   ```bash
   cd ..\.. 
   conda env create -f environment.yml
   conda list
   ```
   
6. **Quit Terminal**

   ```bash
   exit
   ```
   
### Running simulations in Jupyter Notebook

  The following steps 1-4 must be every time a new Jupyter Notebook session is started

1. **Open Terminal**:
   - On MacOS: Press `cmd + space` to open spotlight search and type 'terminal'.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for command prompt or PowerShell in the start menu.

2. **Activate conda environment 'neuron'**

   Navigate back to the main directory

   On MacOS / Linux:
   ```bash
   cd documents/SPNfinal
   conda activate neuron
   ```

   On Windows:
   ```bash
   cd C:\Users\YourUsername\documents\SPNfinal
   conda activate neuron
   ```
4. **Run Jupyter notebook**

   Add 'neuron' environment then open Jupyter Notebook

   ```bash
   python -m ipykernel install --user --name neuron --display-name "Python (neuron)"
   jupyter notebook
   ```

5. **Run a simulation**

   Jupyter Notebook should be open in the default browser

   Choose a Notebook to open (by clicking on any notebook - *.ipynb)

   Ensure kernel is set to Python 3 (ipykernel)

   From kernel dropdown menu choose `Restart Run All` (if running again then it's good practice to run `Restart and Clear Output` first)

   Code should run and generate raw data used to generate figures

   If option `save = True` in the Notebook then the raw figures and pickled data is stored in a subdirectory within the main one


## Data Analysis

  The final analysis and figures used in the manuscript were made using R:
  - R version 4.3.1 – "Beagle Scouts"
  - [R Statistical Software](https://www.R-project.org/)

  Refer to the `R analysis` directory for the code.

  Each simulation has a unique identifier; for instance, `Fig5_EF.ipynb` is `sim4`. 
  
  Once the Jupyter Notebook is executed with `save = True`, the outputs are stored automatically. 
  
  In this case, raw trace data is stored as pickled files in the subdirectory `dspn/model1/physiological/simulations/sim4`. 
  
  Any images generated are found in `dspn/model1/physiological/images/sim4`. 
  
  The R code to analyse the output from `Fig5_EF.ipynb` is found in `Fig5_EF.R` in the `R analysis` directory. 
  
  In order to locate the raw data to recreate the final images for the ms it may be necessary to alter the line:

  On MacOS/Linux: 
    
    path <- paste0('/Documents/GitHub/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')

  On Windows: 
    
    path <- paste0('C:\\Users\\YourUsername\\Documents\\GitHub\\SPNfinal\\', spn, '\\model', model, '\\physiological\\simulations\\sim', sim, '\\')

  This line should be the only one that it is necessary to change in order to execute the R code.






## Anaconda vs Miniconda

Anaconda and Miniconda are both popular distributions for Python and R programming in data science. They include the Conda package manager and aim to simplify package management and deployment.

**Anaconda** is a full-featured distribution that includes:

- Python and R language
- Conda package manager
- Over 1,500 pre-installed scientific packages
- Tools like Jupyter Notebook, Spyder, etc.

It's ideal for those who prefer an out-of-the-box setup for data science and scientific computing.

**Miniconda** offers a minimalistic approach:

- Python or R language
- Conda package manager
- No pre-installed packages

It's suitable for users who want to have a lightweight base to start with and prefer to install packages as needed.

**Advantages of Anaconda**:

- Quick, easy setup with a comprehensive suite of scientific packages and tools.
- Wide array of data science tools readily available within a single application.

**Advantages of Miniconda**::

- Lightweight, minimal base installation.
- Control over which packages are installed.
- Requires limited disk space or bandwidth.
- Clean environment that only includes packages required

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

* **YAML Environment File**: The file `environment.yml` is a YAML file commonly used in Conda

  This file specifies the dependencies and settings for a particular virtual environment.

* **Environment Name** - name Key: In the `environment.yml` file, there is a key called name.

  The value associated with this key is the name of the Conda environment to be created.

  In this case, this name is `neuron`. This is the name that will subsequently be used

  to refer to the environment when activating it or installing additional packages into it.

* **Creating the Environment**: When the command `conda env create -f environment.yml` is executed

  Conda reads the `environment.yml` file and creates a new environment based on the specifications in that file.

  The environment will have the name given by the name key in the YAML file i.e. `neuron`.

  

## GitHub

For beginners, the [GitHub Desktop GUI](https://desktop.github.com/) is recommended. 

Instructions for cloning a repository using GitHub Desktop can be found [here](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop).

## Contact

The model was adapted from the publicly available code by Clay Surmeier and Vernon Clarke

vernon.clarke@northwestern.edu

clay-surmeier@northwestern.edu


For queries related to this repository, please [open an issue](https://github.com/your-repo-link/issues) or [email](vernon.clarke@northwestern.edu) directly 

---


   

