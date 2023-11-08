# SPNfinal

This repository contains a library of striatal projection neurons.


vernon.clarke@northwestern.edu

## General
The code is adapted from the work and published code found in the publication:

    Lindroos and Hellgren Kotaleski 2020. 
    Predicting complex spikes in striatal projection neurons of the direct pathway 
    following neuromodulation by acetylcholine and dopamine. EJN



Translated versions of the models in: 
    
    Hjorth et al., 2020. 
    The microcircuits of striatum in silico. PNAS
    
are also included. See the included example file and the specific branch for these models.



## Model software

NEURON+python: https://www.neuron.yale.edu/neuron/download
(tested on version 8.1 and 8.2)










## Instructions for setting up and running Jupyter Notebook

**I recommend the Anaconda/Miniconda method for running Jupyter Notebook**. 

On my machine (iOS Sonoma, M2 pro), there is a massive speed advantage in running this code in the conda environment.
The package versions seem the same so there is presumably a big difference in the underlying optimization. 
Since conda was designed with data science libraries in mind,the fact that conda packages are often more optimized than pip for scientific computing, reduced execution is expected. 
Conda has MKL (Math Kernel Library) optimized packages for many scientific libraries like NumPy, SciPy, etc., 
which can significantly improve performance for certain operations. If you wish to follow this advice then skip straight to B2

## A. Check if Python is Installed

To determine if Python is installed on your system and to check its version, follow these steps:

1. **Open your Terminal**:
   - On Windows: Search for command prompt or PowerShell in the start menu.
   - On MacOS: Press `cmd + cpace` to open spotlight search and type 'terminal'.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.

2. **Check Python version**:
   ```bash
   python --version
   ```

If Python is installed, the version number will be displayed

## Install Python (if not already installed)

Follow these instructions based on your operating system:

### On Windows

- **Download Python**: Navigate to the [official Python website](https://python.org) and download the latest version for Windows.
- **Install Python**: Open the downloaded installer. Ensure you select the "Add Python to PATH" option during installation for easier command-line access.

### On MacOS

- **Download Python**: Visit the [official Python website](https://python.org) to download Python for MacOS.
- Alternatively, use **Homebrew** to install Python by opening Terminal and running:
    ```bash
    brew install python3
    ```

### On Linux

Python is usually pre-installed on Linux. If you need to install or update it, use the package manager for your distribution:

- **On Ubuntu** (and Debian-based systems):
    ```bash
    sudo apt update
    sudo apt install python3
    ```
- **On Fedora** (and RHEL-based systems):
    ```bash
    sudo dnf install python3
    ```
- **On Linux**:
    ```bash
    sudo pacman -S python3
    ```

## B1. Setting up Jupyter Notebook without Anaconda/Miniconda

#### In this example a directory called Jupyter (this can be named whatever you like) has been created within the documents folder

### Create 'fitting' environment <span style="font-size: 80%;">(do once to set up)</span>:

1. **Open Terminal**

2. **Navigate to the Jupyter Directory**
    ```bash
    cd documents
    cd Jupyter
    ```

3. **Create a new virtual environment named 'fitting'**
    ```bash
    python3 -m venv fitting
    ```

4. **Activate the 'fitting' environment**
   - On MacOS and Linux:
     ```bash
     source fitting/bin/activate
     ```
   - On Windows:
     ```bash
     .\fitting\Scripts\activate
     ```

5. **Install the required packages**
    ```bash
    pip install jupyter numpy pandas matplotlib openpyxl plotly scipy tqdm ipywidgets
    ```

6. **Quit terminal**
    ```bash
    deactivate
    exit
    ```
    
7. **Close the terminal window/tab**

### To Run  <span style="font-size: 80%;">(do this every time you want to run the code)</span>:

1. **Open terminal**

2. **Navigate back to the Jupyter directory and activate the 'fitting' environment**
    ```bash
    cd documents
    cd Jupyter
    source fitting/bin/activate  # On Windows use: .\fitting\Scripts\activate
    ```

3. **Launch Jupyter Notebook**
    ```bash
    jupyter notebook
    ```

4. **Exit Jupyter Notebook when finished**
    ```bash
    ctrl+C
    ```

5. **Deactivate the virtual environment**
    ```bash
    deactivate
    ```

6. **Close the terminal window/tab**
   

## B2. Setting up Jupyter Notebook with Anaconda/Miniconda

**An alternative to the above is to use Anaconda or Miniconda**

Anaconda and Miniconda are Python distributions that include the Python interpreter along with a suite of other tools and libraries.

Even if Python is already installed, you may prefer to use Anaconda/Miniconda.

- **Anaconda**: Download the installer from the [Anaconda website](https://www.anaconda.com/products/individual).
- **Miniconda**: Download the installer from the [Miniconda website](https://docs.conda.io/en/latest/miniconda.html).

Run the downloaded installer and follow the instructions to set up your Python environment.


### create 'fitting' environment <span style="font-size: 80%;">(do once to set up)</span>:

1. **Open Terminal**

2. **Navigate to Jupyter directory**
    ```bash
    cd documents
    cd Jupyter
    ```

3. **Create a new conda environment named 'fitting'**
    ```bash
    conda create -n fitting
    ```

4. **Activate the 'fitting' environment**
    ```bash
    conda activate fitting
    ```

5. **Install the required packages**
    ```bash
    conda install -n fitting sqlite jupyter numpy pandas matplotlib openpyxl plotly scipy tqdm ipywidgets
    ```

6. **Quit Terminal**
    ```bash
    exit
    ```

7. **Close the terminal window/tab**

### To run:

1. **Open Terminal**

2. **Navigate back to the Jupyter directory and activate the 'fitting' environment**
    ```bash
    cd documents
    cd Jupyter
    conda activate fitting
    ```

3. **Launch Jupyter Notebook**
    ```bash
    jupyter notebook
    ```

4. **Exit Jupyter Notebook when finished**
    ```bash
    ctrl+C
    ```

5. **Close the terminal window/tab**





How to run the models (would need neuron+python installed, see below)
------------------------------------------------------------------------------

1) compile the mechanisms
2) run example.py


e.g. from a terminal:

    cd mechanisms/single/
    nrnivmodl
    cd ../../
    python3 example.py


Model software
------------------------------------------------------------------------------

NEURON+python: https://www.neuron.yale.edu/neuron/download
(tested on version 8.1 and 8.2)


General
------------------------------------------------------------------------------

Some of the models are used in the publication:

    Lindroos and Hellgren Kotaleski 2020. 
    Predicting complex spikes in striatal projection neurons of the direct pathway 
    following neuromodulation by acetylcholine and dopamine. EJN

Code for simulating, analysing and plotting of ispn is also included (in the ispn branch).
Code for dspn is located in another repo:
https://bitbucket.org/rlindroos/neuron/src/cholinergic_modulation/Complex_spike/


Simulation/analysis/plotting is included in the ispn branch under "Simulations".


Translated versions of the models in: 
    
    Hjorth et al., 2020. 
    The microcircuits of striatum in silico. PNAS
    
are also included. See the included example file and the specific branch for these models.
