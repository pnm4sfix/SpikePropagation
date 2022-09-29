# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:25:48 2020

@author: pierc
"""

import numpy as np
import pandas as pd

#import h5py
import os
import seaborn as sns
import matplotlib as mpl

from tkinter import filedialog
from tkinter import Tk
import spikeanalysis as sa
#import scipy.io as sio
    
mpl.rcParams['lines.linewidth'] = 0.3
mpl.rcParams['font.size'] = 18


def multiple_simulation(new_method, new_new_method, distance, poisson):
    """Multiple repeats. Problem exists in plotting alpha etc where it is delta etc that is changing"""
    final_df = pd.DataFrame()
    final_failure_df = pd.DataFrame()
    failures2 = []
    for trial in range(5):
        if poisson ==True:
            core_df, core_failure_df, failures = poisson_simulation(new_method, new_new_method, distance)
        else:
            core_df, core_failure_df, failures = real_simulation(new_method, new_new_method, distance)
            
        final_df = pd.concat([final_df, core_df], axis=0)
        final_failure_df = pd.concat([final_failure_df, core_failure_df], axis=0)
        failures2.append(failures)
        
    final_df.deleted = (final_df.deleted*100).astype("int64")
    
    if poisson == True:
        alpha=sns.relplot(data = final_df, x="freq", y="accuracy", hue = "deleted", 
                      col = "alpha", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        
        beta=sns.relplot(data = final_df, x="freq", y="accuracy", hue = "deleted", 
                      col = "beta", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        
        delta=sns.relplot(data = final_df, x="freq", y="accuracy", hue = "deleted", 
                      col = "delta", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        
        c=sns.relplot(data = final_df, x="freq", y="accuracy", hue = "deleted", 
                      col = "c", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
    else:
    
        alpha = sns.relplot(data = final_df, x="deleted", y="accuracy", 
                      col = "alpha", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        beta = sns.relplot(data = final_df, x="deleted", y="accuracy", 
                      col = "beta", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        delta = sns.relplot(data = final_df, x="deleted", y="accuracy", 
                      col = "delta", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
        c =sns.relplot(data = final_df, x="deleted", y="accuracy", 
                      col = "c", col_wrap = 1, kind="line", palette=sns.color_palette("colorblind", n_colors=5), linewidth=1)
    
    
    alpha.fig.tight_layout()
    beta.fig.tight_layout()
    delta.fig.tight_layout()
    c.fig.tight_layout()
    
    alpha.savefig("./"+str(new_method)+str(new_new_method)+str(distance)+str(poisson)+"_alpha.svg")
    beta.savefig("./"+str(new_method)+str(new_new_method)+str(distance)+str(poisson)+"_beta.svg")
    delta.savefig("./"+str(new_method)+str(new_new_method)+str(distance)+str(poisson)+"_delta.svg")
    c.savefig("./"+str(new_method)+str(new_new_method)+str(distance)+str(poisson)+"_c.svg")
    
    with pd.ExcelWriter("./"+str(new_method)+str(new_new_method)+str(distance)+str(poisson)+".xlsx") as writer:
        final_df.to_excel(writer, sheet_name = "final_df")
        final_failure_df.to_excel(writer, sheet_name = "final_failure_df")
    
    return final_df, final_failure_df, failures, alpha, beta, delta, c 



def real_simulation(new_method, new_new_method, distance):
    """This function varies velocities using real st1. 
    Maybe run this a few times because of the inherent probability variability.
    Seems to fail everynow and then"""
    core_df = pd.DataFrame()
    core_failure_df = pd.DataFrame()
    base_velocities = np.array([100., 50., 20., 1.])
    new_velocities =[]
    failures = []
    for vel in base_velocities:
        index = np.where(base_velocities == vel)[0]
        very_fast = vel * 4.
        fast = vel*2. 
        medium = vel
        slow = vel/2.
        very_slow = vel/5.
        new_vels = [very_fast, fast, medium, slow, very_slow]
        for new_vel in new_vels:
            new_velocity = np.copy(base_velocities)
            new_velocity[index] = new_vel
            new_velocities.append(new_velocity)
    freq = None        
    nixfile = filedialog.askopenfilename(initialdir = os.getcwd(), title = "Please select nix file")
    bl, st1, st2 = sa.read_block(nixfile)
    
    for velocities in new_velocities:
        try:
            success_df, failure_df = sa.create_simulation(st1, freq, velocities[0], velocities[1], velocities[2], velocities[3],
                                                    new_method, new_new_method, distance, False)
            core_df = pd.concat([core_df, success_df], axis=0)
            core_failure_df = pd.concat([core_failure_df, failure_df], axis=0)
        except:
            failures.append([freq, velocities])
            print("failed")
    return core_df, core_failure_df, failures



multiple_simulation(False, False, 2e-3, False)