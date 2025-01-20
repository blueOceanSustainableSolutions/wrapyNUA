# Introduction
**wrapyNUA** is a pre- and postprocessing script for the Python wrapper around the Fortran numerical uncertainty tools available [here](https://www.marin.nl/en/research/free-resources/verification-and-validation/verification-tools).

It was designed to provide a similar experience to the original tools but with enhanced plotting capabilities. Compared to the original tools, the separate input configuration file *input_nua.ini* is not needed any more, the settings are located at the beginning of the wrapyNUA script. The names of the setting variables are mostly similar, the parameters are documented directly in the script. The required structure of the data file is also similar. (The main difference is that the number of solutions to use for the fit is **not** indicated in the data file. Instead, the wrapyNUA script provides an option to specify a list with solutions to be ignored.)

The wrapyNUA repository comes with a tutorial for a steady and an unsteady case. The data files of both cases are commented with explanations of the required structure. The corresponding settings to configure the script are documented within the script. The settings for the MARIN examples supplied with the original tools are also documented. It is a good exercise to try to replicate those examples and check if you can achieve the same results. (The estimated uncertainties **must** be identical, the Fortran procedure working under the hood is the same.)

**If you are new to numerical uncertainty analysis in general, it is recommended that you first get familiar with the theory through the papers by Eca et al. referenced on the Verification Tools homepage.** The `doc` folder in the downloadable archive contains a selection of publications in PDF format.

# System requirements
- Linux operating system
- Python 3.6 or newer with **Matplotlib**, **NumPy** and **pandas** installed (Matplotlib $\geq$ 3.9 is recommended, which fixes a bug trying to avoid that the legend is plotted over text)
- Libraries shipped with the **Verification Tools**: *libifcoremt.so.5*, *libifport.so.5*, *libimf<area>.so*, *libintlc.so.5*, *libsvml<area>.so*

# Installation
- Create an installation directory and place the *wrapyNUA<area>.py* script inside it
- Download the Linux version of the Verification Tools and extract the archive
- Copy following three files from the `/bin` folder of the Verification Tools into the installation directory: *numerical_uncertainty*, *numerical_uncertainty_wrapper*, *wrap<area>.py*
- Copy the libraries into a directory of your choice, it may be the same as the installation directory

# Usage
 - Open a terminal in the installation directory
 - Export the location of the libraries (adapt the path as needed):
 `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/libraries`
 - If you installed the required Python packages in your system's Python you can skip this step. Otherwise load your Python environment. In the case you use Anaconda:
 `conda activate name-of-your-environment`
 - Prepare your uncertainty data file and adapt the settings in the wrapyNUA script to your needs. If you are new to wrapyNUA, check the tutorial examples supplied in the repository.
 - Launch the script:
 `python wrapyNUA.py`
