# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:15:52 2020

@author: pierc
"""


import pandas as pd
import neo
import elephant as el
import numpy as np
from matplotlib import pyplot as plt
import quantities as pq
import sa_core as sa
import os
import seaborn as sns
#%matplotlib notebook

## TO DO
#Simulate 80% reduction
#Test new correlation
#Fix inst firing frequency excel save
#Check method of threshold with Quiroga
#Check normalisation of new correlation as its only doing so for the window snippet
sns.set()
sns.set_context("notebook", font_scale=2.5, rc={"lines.linewidth": 0.3})
sns.set(style="white")

def subset(condition, time, DR, SN):
    h, m, s = time.split(':')
    """Time in seconds from mins"""
    t = (int(h) * 3600 + int(m) * 60 + int(s)) * pq.s

    if condition == "drug":
        subtract = -35.001 * pq.s
        add = 70.001 * pq.s

        con_st = t + (-35 * pq.s)
        con_end = t + (-5 * pq.s)
        drug_st = t
        drug_end = t + (30 * pq.s)
        wash_st = t + (35 * pq.s)
        wash_end = t + (65.00 * pq.s)
        time_slices = {"control": (con_st, con_end),
                     "drug": (drug_st, drug_end),
                      "wash":(wash_st, wash_end)}

        t0, t1 = (t + subtract), (t + add)

    elif condition == "long":
        subtract = -60.001 * pq.s
        add = 160.001 * pq.s
        con_st = t + (-60 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (60 * pq.s)
        wash_st = t + (90 * pq.s)
        wash_end = t + (150.00 * pq.s)

        time_slices = {"control": (con_st, con_end),
                     "drug": (drug_st, drug_end),
                      "wash":(wash_st, wash_end)}

        t0, t1 = (t + subtract), (t + add)

    elif condition == "baseline":
        add = 60.01 * pq.s
        st = t
        end = t+(60*pq.s)

        time_slices = {"control": (st, end)}
        t0, t1 = (t+(0.01*pq.s)), (t+add)

    DR_subset = DR.time_slice(t0, t1)
    SN_subset = SN.time_slice(t0, t1)

    return DR_subset, SN_subset, time_slices

def han_io(time, condition):
    """Conditions to choose from include drug and baseline.
    Input time as hh:mm:ss"""
    filepath = sa.file_paths()

    time_split =  ",".join(time.split(":"))
    new_path = filepath[:-4]+"_"+time_split
    if not os.path.exists(new_path):
        new_folder = os.makedirs(filepath[:-4]+"_"+time_split)
        print(new_path)
    else:
        new_path = new_path

    df=pd.read_csv(filepath, header = 2, sep ="\t")
    df = df.iloc[:,2:4]
    df.columns =["DR", "SN"]

    DR = df.DR.dropna()
    SN = df.SN.dropna()

    """Create neo AnalogSignals"""

    DR_asig = neo.AnalogSignal(np.array(DR)*pq.uV, t_start = 0*pq.s,

                              sampling_rate=(5000)*pq.Hz,
                              name = "DR sig")

    SN_asig = neo.AnalogSignal(np.array(SN)*pq.uV, t_start = 0*pq.s,

                              sampling_rate=(5000)*pq.Hz,
                              name = "SN sig")

    DR_sig, SN_sig, time_slices = subset(condition, time, DR_asig, SN_asig)
    return DR_sig, SN_sig, time_slices, filepath


def cap_io(subdir):
    if subdir ==None:
        control = sa.file_paths()
        new_path = os.path.dirname(control)
    else:
        new_path = subdir
    files = os.listdir(new_path)
    files = [file.lower() for file in files]
    files = [os.path.join(new_path, file) for file in files]
    control = [file for file in files if "control" in file][0]
    cap = [file for file in files if "cap" in file][0]
    gaba = [file for file in files if "gaba" in file][0]

    filepaths = {"control":control, "cap":cap, "gaba":gaba}
    samp_freq = 20e3
    dt = 1/samp_freq
    core_df = pd.DataFrame()
    #time_slices = {}
    for key in filepaths.keys():
        path = filepaths[key]
        df = pd.read_csv(path, header = 2, sep ="\t")
        df = df.iloc[:,2:4]
        df.columns =["DR", "SN"]
        df["Condition"] = key


        core_df = pd.concat([core_df, df], axis=0)
    core_df["Time"] = np.arange(0, core_df.shape[0]*dt, dt)
    DR = core_df.DR.dropna()
    SN = core_df.SN.dropna()

    DR_asig = neo.AnalogSignal(np.array(DR)*pq.uV, t_start = 0*pq.s,

                              sampling_rate=(20000)*pq.Hz,
                              name = "DR sig")

    SN_asig = neo.AnalogSignal(np.array(SN)*pq.uV, t_start = 0*pq.s,

                              sampling_rate=(20000)*pq.Hz,
                              name = "SN sig")

    time_slices = {}
    for condition in core_df.Condition.unique():
            subset = core_df[core_df.Condition == condition]
            t0 = subset.iloc[0] ["Time"]*pq.s
            t1 = subset.iloc[-1]["Time"]*pq.s
            time_slices[condition] = (t0, t1)

    time = str(time_slices[next(iter(time_slices))][1])
    return DR_asig, SN_asig, time_slices, new_path, time


def get_spikes(condition, time, min_std=3, max_std=20, io=1, subdir=None):
    if io == 1:
        DR_sig, SN_sig, time_slices, new_path = han_io(time, condition)
    elif io == 2:
        DR_sig, SN_sig, time_slices, new_path, time = cap_io(subdir)



    DR_peaks, DR_amps, SN_peaks, SN_amps, DR_sig, SN_sig = sa.spike_detection(DR_sig, SN_sig, 1/DR_sig.sampling_rate,
                                                               DR_sig, SN_sig, time_slices,
                                                               DR_sig, SN_sig, min_std, max_std)
    DR_df = sa.extract_spikes(DR_sig, DR_peaks)
    SN_df = sa.extract_spikes(SN_sig, SN_peaks)

    try:
            DR_clusters, DR_pca_df, DR_sse = sa.clusters(DR_df.transpose(), True, 4)
            SN_clusters, SN_pca_df, SN_sse = sa.clusters(SN_df.transpose(), True, 4)
    except:
        print("Couldn't cluster-probs not enough peaks")
        DR_clusters, DR_pca_df, DR_sse = [None, None, None]
        SN_clusters, SN_pca_df, SN_sse = [None, None, None]

    """Spectra analysis"""
    DR_spectra = sa.spectra(DR_sig, time_slices)
    SN_spectra = sa.spectra(SN_sig, time_slices)

    DR_freq = sa.calc_freq(DR_peaks, time_slices)
    SN_freq = sa.calc_freq(SN_peaks, time_slices)

    sa.main_plot(DR_sig, DR_sig, DR_sse, DR_freq, DR_clusters, DR_pca_df, DR_spectra, DR_peaks, new_path,  "DR", "0", "noxious", time_slices)
    sa.main_plot(SN_sig, SN_sig, SN_sse, SN_freq, SN_clusters, SN_pca_df, SN_spectra, SN_peaks, new_path,  "SN", "0", "noxious", time_slices)

    sa.save(DR_freq, DR_clusters, DR_pca_df, DR_spectra, "0", "dorsal root", new_path)
    sa.save(SN_freq, SN_clusters, SN_pca_df, SN_spectra, "0", "spinal nerve", new_path)

    """Create bl and export"""
    bl =sa.create_neo_block(DR_sig, DR_sig, SN_sig, SN_sig, DR_peaks, SN_peaks, new_path, time_slices,
                            DR_df, SN_df)
    sa.spike_train_plot(bl)
    sa.wave_clus_export(DR_df, SN_df, DR_peaks, SN_peaks, DR_sig, SN_sig, new_path)
    try:
        sa.save_block(bl, "0", new_path, True)
    except:
        try:
            sa.save_block(bl, "0", new_path, True)
        except:
            print("failed to save blk")
def spike_matching():
    sa.spike_analysis(None, None, 5e-4)

def teased_fibre_io():
    """Get positive peaks which is stim, get negative peaks which should be spike.
    Df read slightly different"""
    filepath = sa.file_paths()
    df=pd.read_csv(filepath, header = 1, sep ="\t")
    plt.figure()
    df.plot()
    pass

#get_spikes("long", "00:01:40", 5, 30)

#spike_matching()
#teased_fibre()
