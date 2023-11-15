# SPNfinal
vernon.clarke@northwestern.edu

clay-surmeier@northwestern.edu

This repository contains a Neuron model of striatal projection neurons built on top of 'striatal_SPN_lib' created by Lindroos and Hellgren Kotaleski 2020
    
    Lindroos and Hellgren Kotaleski 2020. 
    Predicting complex spikes in striatal projection neurons of the direct pathway 
    following neuromodulation by acetylcholine and dopamine. EJN

### NEURON

* The models are built in NEURON+python: https://www.neuron.yale.edu

  (tested on version 8.1 and 8.2)

* Install Neuron: https://www.neuron.yale.edu/neuron/download

* Quickstart: https://www.neuron.yale.edu/ftp/neuron/2019umn/neuron-quickstart.pdf
* (Old) Mac instructions: https://www.neuron.yale.edu/ftp/neuron/nrn_mac_install.pdf

* Install Python
https://python-docs.readthedocs.io/en/latest/starting/install3/osx.html




### Virtual Environments
There is a yaml environment file that is set to work with NEURON 8.2.1 and python 3.9.16 named 'environment'. 

Be sure that your versions of NEURON and python are compatible if using a different distribution.

* For setting up Conda (python package manager): https://conda.io/projects/conda/en/latest/user-guide/getting-started.html
* For setting up Jupyter Notebook (interactive code notebooks): https://jupyter.org/install

### Github
If you are unfamiliar with Github, the desktop app is a useful interface to use.

* For installing Github Desktop GUI: https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/installing-and-authenticating-to-github-desktop/installing-github-desktop
* Cloning a repository using Github Desktop: https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop

How to run the models
---------------------

The following steps 1-6 must be run once prior to running the code in Jupyter Notebook

1. **Install NEURON with python support (see setup instructions)**

2. **Open Terminal**:
   - On MacOS: Press `cmd + cpace` to open spotlight search and type 'terminal'.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for command prompt or PowerShell in the start menu.

3. **Compile mechanisms**:

   Navigate to directory containing NEURON mechanisms. For instance if 'SPNfinal' is in documents folder in OS, mechanisms are located in cd/documents/SPNfinal/mechanisms/single

   Mechanisms are then compiled by entering 'nrnivmodl' or 'mknrndll'

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
   There is a yaml file in the main directory called environment.yml. This can be used to create a conda environment called 'neuron'

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

   
7. **Quit Terminal**
   ```bash
   exit
   ```

The following steps 1-6 must be every time a new Jupyter Notebook session is started

1. **Open Terminal**:
   - On MacOS: Press `cmd + cpace` to open spotlight search and type 'terminal'.
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
3. **Run Jupyter notebook'**

   Adds environment then open Jupyter Notebook

   ```bash
   python -m ipykernel install --user --name neuron --display-name "Python (neuron)"
   jupyter notebook
   ```

4. **Run a simulation**

   Jupyter Notebook should be open in default browser

   
   
using command 'conda env create -f environment.yml' (make sure to 'cd ../..' back to main directory) then check this installed correctly with 'conda list'

4. Activate conda environment 'conda activate neuron'
5. Add environment to Jupyter notebook ipykernel using command 'python -m ipykernel install --user --name neuron --display-name "Python (neuron)"'
6. Run Jupyter notebook using command 'jupyter notebook'
7. Open notebook (i.e. click on 'example.ipynb')
8. Choose "Python (neuron)" kernel under Kernel > Change Kernel > Python (neuron)
9. Run code!







## Model software

NEURON+python: https://www.neuron.yale.edu/neuron/download
(tested on version 8.1 and 8.2)

Install Neuron
https://www.neuron.yale.edu/neuron/download

Install Python
https://python-docs.readthedocs.io/en/latest/starting/install3/osx.html










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

