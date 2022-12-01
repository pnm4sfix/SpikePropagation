# SpikePropagation
Code and data repository for paper - "Dorsal root ganglia control nociceptive input to the central nervous system"


Code and instructions for the Spike sorting and analysis are provided in GitHub repository DOI:10.5281/zenodo.7350612

Fig.2A-E. These were created using the matplotlib and seaborn libraries and can be reproduced by running the scripts or Jupyter notebook “Analysis (060720).ipynb” in the GitHub repository DOI:10.5281/zenodo.7350612. Panels 2A and 2C are from raw data and spiketrains that are saved in Neuroscience Information Exchange (NIX) format as “Fig2.blk” in folder NIX_data (DOI:10.5281/zenodo.7350612). These can be accessed using python and by following the jupyter notebook in the repository (DOI:10.5281/zenodo.7350612). Panel 2C was created by using the command  spike_analysis() in sa_core.py on Github. Panles B, D, E were created in python from data stored in “Fig2.blk”  with scripts in the Github repository.

Fig.S5 A-B. Panel A was produced by running the spike_analysis() function in sa_core.py during the spike matching process. Panel B was produced by running the simulation notebook (1.5cm inter electrode sim.ipynb) in the github repository.

Fig. S6 A-D. Panel A was produced by the python function spike_analysis() in sa_core.py (in Github repository DOI:10.5281/zenodo.7350612) from FigS6.blk in NIX_data. Panels B-D were created in python from data stored in “FigS6.blk” with scripts in the Github repository (DOI:10.5281/zenodo.7350612). Data containing the waveforms in Fig. S6B, C, D  was also saved in excel sheets now in FigS6BCD.xlsx.

Fig. S7. Panel A was produced in Jupyter notebook “Analysis (060720).ipynb by extracting the waveforms from the NIX file from each recording. Panel B was created in python using matplotlib from the metadata in FigS7B.xlsx. Panel C was created in python from the instantaneous firing rates calculated from the NIX recording file and saved in FigS7C.xlsx.