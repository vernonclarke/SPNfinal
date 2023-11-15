# SPNfinal
vernon.clarke@northwestern.edu

clay-surmeier@northwestern.edu

This repository contains a NEURON model of striatal projection neurons built on top of 'striatal_SPN_lib' created by Lindroos and Hellgren Kotaleski 2020
    
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

### Running the models

The following steps 1-6 must be run once prior to running the code in Jupyter Notebook

1. **Install NEURON with python support** (see setup instructions)

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

The following steps 1-4 must be every time a new Jupyter Notebook session is started

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
3. **Run Jupyter notebook**

   Adds environment then open Jupyter Notebook

   ```bash
   python -m ipykernel install --user --name neuron --display-name "Python (neuron)"
   jupyter notebook
   ```

4. **Run a simulation**

   Jupyter Notebook should be open in the default browser

   Choose a Notebook to open (by clicking on any notebook - *.ipynb)

   Ensure kernel is set to Python 3 (ipykernel)

   From kernel dropdown menu choose 'Restart Run All' (if running again then it's good practice to run 'Restart and Clear Output' first)

   Code should run and generate raw data used to generate figures

   If option save = True in the Notebook then the raw figures and pickled data is stored in a subdirectory within the main one

The final analysis and figures used in the ms were made using R

R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

This code is provided in the R analysis directory


   

