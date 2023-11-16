# SPNfinal

clay-surmeier@northwestern.edu

vernon.clarke@northwestern.edu

This repository contains a NEURON model of striatal projection neurons built 

on top of 'striatal_SPN_lib' created by Lindroos and Hellgren Kotaleski 2020:

Lindroos R, Kotaleski JH. Predicting complex spikes in striatal projection neurons of the 

direct pathway following neuromodulation by acetylcholine and dopamine. Eur J Neurosci. 2020. 

[doi:10.1111/ejn.14891](https://doi.org/10.1111/ejn.14891)


### NEURON

* The models are built in NEURON+python: https://www.neuron.yale.edu

  (tested on version 8.1 and 8.2)

* Install Neuron: https://www.neuron.yale.edu/neuron/download

* Quickstart: https://www.neuron.yale.edu/ftp/neuron/2019umn/neuron-quickstart.pdf
* (Old) Mac instructions: https://www.neuron.yale.edu/ftp/neuron/nrn_mac_install.pdf

* Install Python
https://python-docs.readthedocs.io/en/latest/starting/install3/osx.html


### Virtual Environments
There is a yaml environment file ('environment.yml') that is set to work with NEURON 8.2.1 and python 3.9.16 (name 'neuron'). 

In brief:     

* **YAML Environment File**: The file environment.yml is a YAML file commonly used in Conda

  This file specifies the dependencies and settings for a particular virtual environment.

* **Environment Name** - name Key: In the environment.yml file, there is a key called name.

  The value associated with this key is the name of the Conda environment to be created.

  In this case, this name is 'neuron'. This is the name that will subsequently be used

  to refer to the environment when activating it or installing additional packages into it.

* **Creating the Environment**: When the command 'conda env create -f environment.yml' is executed (see later)

  Conda reads the environment.yml file and creates a new environment based on the specifications in that file.

  The environment will have the name given by the name key in the YAML file i.e. 'neuron'.

* For setting up Conda (python package manager): https://conda.io/projects/conda/en/latest/user-guide/getting-started.html
* For setting up Jupyter Notebook (interactive code notebooks): https://jupyter.org/install

### Github
If you are unfamiliar with Github, the desktop app is a useful interface to use.

* For installing Github Desktop GUI: https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/installing-and-authenticating-to-github-desktop/installing-github-desktop
* Cloning a repository using Github Desktop: https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop

### Running the models

The following steps 1-6 must be run once prior to running the code in Jupyter Notebook

1. **Install NEURON with python support** (see setup instructions)

2. **Open Terminal**:
   - On MacOS: Press `cmd + space` to open spotlight search and type 'terminal'.
   - On Linux: Search for terminal in your applications menu or press `ctrl + alt + T`.
   - On Windows: Search for command prompt or PowerShell in the start menu.

3. **Compile mechanisms**:

   Navigate to directory containing NEURON mechanisms.

   For instance if 'SPNfinal' is in documents folder in OS, mechanisms are located in cd/documents/SPNfinal/mechanisms/single

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

   There is a yaml file in the main directory called environment.yml.

   This can be used to create a conda environment called 'neuron'

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

   From kernel dropdown menu choose 'Restart Run All' (if running again then it's good practice to run 'Restart and Clear Output' first)

   Code should run and generate raw data used to generate figures

   If option save = True in the Notebook then the raw figures and pickled data is stored in a subdirectory within the main one



### Graphs and analysis

The final analysis and figures used in the ms were made using R (R version 4.3.1 (2023-06-16) -- "Beagle Scouts")

    R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

This code is provided in the 'R analysis' directory

Each simulation has a unique identifier; for instance, Fig5_EF.ipynb is sim4. 

Once the Jupyter Notebook is executed with save = True, the outputs are stored automatically. 

In this case, raw trace data is stored as pickled files in the subdirectory dspn/model1/physiological/simulations/sim4; any images generated are found in dspn/model1/physiological/images/sim4

The R code to analyse the output from Fig5_EF.ipynb is found in Fig5_EF.R in the 'R analysis' directory.

In order to locate the raw data to recreate the final images for the ms it may be necessary to alter the line

On MacOS/Linux: **path <- paste0('/Documents/GitHub/SPNfinal/', spn, '/model', model, '/physiological/simulations/sim', sim,'/')**

On Windows: **path <- paste0('C:\\Users\\YourUsername\\Documents\\GitHub\\SPNfinal\\', spn, '\\model', model, '\\physiological\\simulations\\sim', sim, '\\')**

This line should be the only one that it is necessary to change in order to execute the R code.
   

