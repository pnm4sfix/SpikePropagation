# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:30:44 2019

@author: fbspmul
"""
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, find_peaks_cwt, peak_widths, peak_prominences
from scipy.signal import welch
from scipy.signal import butter, filtfilt
from scipy import signal
from scipy.stats import linregress
import scipy.io as sio
import scipy
import time as tt
#import h5py
import os
import neo
from tqdm import tqdm
import matplotlib as mpl
from neo.io.neomatlabio import NeoMatlabIO
import quantities as pq
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from itertools import cycle, islice
from pandas.plotting import parallel_coordinates
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from scipy import stats, spatial, cluster
from tkinter import filedialog
from tkinter import Tk
import elephant
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d
import elephant.signal_processing as elsig
import elephant.statistics as el #isi, cv, mean_firing_rate, instantaneous_rate, complexity_pdf
from elephant.spike_train_correlation import spike_time_tiling_coefficient
from elephant.spike_train_dissimilarity import van_rossum_dist, victor_purpura_dist
from elephant.spike_train_generation import homogeneous_poisson_process
from scipy.signal import correlate
import datetime
import gc
import re
from numba.typed import List
from numba import njit
#import scipy.io as sio

sns.set()
sns.set_context("notebook", font_scale=2.5, rc={"lines.linewidth": 0.1})
sns.set_style("white")
#mpl.rcParams['lines.linewidth'] = 0.3
#mpl.rcParams['font.size'] = 14
# TO DO

def plotting_palette():
    custom = [(0, 0, 0), (255, 0, 102), (16, 127, 128), (64, 0, 127),
              (107, 102, 255), (102, 204, 254)]
    custom = [tuple(np.array(list(rgb))/255) for rgb in custom]
    custom2 = [sns.light_palette(color, n_colors = 3, reverse=True) for color in custom]
    custom3 = np.array(custom2).reshape([18, 3]).tolist()
    custom4 = custom3[0::3]+custom3[1::3]+custom3[2::3]
    palette = sns.color_palette(custom4)
    #sns.palplot(palette)
    return palette

try:
    palette = plotting_palette()
except:
    print("error setting custom palette")
    palette = sns.color_palette("colorblind")
sns.set_palette(palette)


def file_paths():
    """Create pop up that allows selection of spike file and template file"""
    root = Tk()
    root.filename = filedialog.askopenfilename(initialdir = os.getcwd(), title = "Please select Spike2 file")
    #root.template = filedialog.askopenfilename(initialdir = os.getcwd(), title = "Please select ECG template")

    return root.filename#, root.template


def read_file(file):
    dorsal_root = 0
    spinal_nerve = 0
    phrenic = 0
    pressure = 0
    #reader = neo.io.Spike2IO(filename=file)
    reader = neo.Spike2IO(filename=file)
    data = reader.read(lazy=False)[0]
    try:
        for seg in data.segments[0].analogsignals:

            ch = str(seg.annotations['channel_names']).split(sep="'")[1]

            if ch == 'ppmmHg':
                    pressure = seg
                    #pressure = signal_proc.emg_proc(flex)
            elif (ch == 'Vagus') | (ch == "DR"):
                    dorsal_root = seg
            elif (ch == 'PNA') | (ch == "SN"):
                    spinal_nerve = seg
            elif (ch == 'EMG') | (ch == "PND"):
                    phrenic = seg
    except:
        print("New version of neo appears to remove channel_names and has bundled signals with same sr")
        try:
            for seg in data.segments[0].analogsignals:

                ch = seg.name

                if ch == 'ppmmHg':
                        pressure = seg
                        #pressure = signal_proc.emg_proc(flex)
                elif (ch == 'Channel bundle (PNA,Vagus) ') | (ch == 'Channel bundle (SN,DR) '):
                    #now work with this bundle to tease appart analogsignals
                        DR_sig = seg[:, 0]
                        SN_sig = seg[:, 1]
                        dorsal_root = DR_sig
                        spinal_nerve = SN_sig
                elif (ch == 'EMG') | (ch == "PND"):
                        phrenic = seg
        except:
            print("print unable to parse block")

    return spinal_nerve, dorsal_root, phrenic, pressure

def new_read_file(file, time, condition, Randall = False):
    dorsal_root = 0
    spinal_nerve = 0
    phrenic = 0
    pressure = 0
    randall = 0
    #reader = neo.io.Spike2IO(filename=file)
    reader = neo.Spike2IO(filename=file)
    data = reader.read(lazy=True)[0]
    h, m, s = time.split(':')
    """Time in seconds"""
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

    elif condition == "drug2":
        subtract = -35.001 * pq.s
        add = 85.001 * pq.s

        con_st = t + (-35 * pq.s)
        con_end = t + (-5 * pq.s)
        drug_st = t
        drug_end = t + (30 * pq.s)
        wash_st = t + (50 * pq.s)
        wash_end = t + (80.00 * pq.s)


    elif condition == "noxious":
        subtract = -15.001 * pq.s
        add = 30.001 * pq.s
        con_st = t + (-15 * pq.s)
        con_end = t + (-5 * pq.s)
        drug_st = t
        drug_end = t + (10 * pq.s)
        wash_st = t + (15 * pq.s)
        wash_end = t + (25.00 * pq.s)


    elif condition == "long":
        subtract = -60.001 * pq.s
        add = 160.001 * pq.s
        con_st = t + (-60 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (60 * pq.s)
        wash_st = t + (90 * pq.s)
        wash_end = t + (150.00 * pq.s)

    if condition == "step30":
        subtract = -30.001 * pq.s
        add = 60.001 * pq.s

        con_st = t + (-30 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (30 * pq.s)
        wash_st = t + (30 * pq.s)
        wash_end = t + (60.00 * pq.s)

    if condition == "ramp30":
        subtract = -30.001 * pq.s
        add = 60.001 * pq.s

        con_st = t + (-30 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (30 * pq.s)
        wash_st = t + (30 * pq.s)
        wash_end = t + (60.00 * pq.s)

    if condition == "ramp60":
        subtract = -60.001 * pq.s
        add = 120.001 * pq.s

        con_st = t + (-60 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (60 * pq.s)
        wash_st = t + (60 * pq.s)
        wash_end = t + (120.00 * pq.s)

    if condition == "phasic":
        subtract = -30.001 * pq.s
        add = 60.001 * pq.s

        con_st = t + (-30 * pq.s)
        con_end = t
        drug_st = t
        drug_end = t + (30 * pq.s)
        wash_st = t + (30 * pq.s)
        wash_end = t + (60.00 * pq.s)



    time_slices = {"control": (con_st, con_end),
                     "stim": (drug_st, drug_end),
                      "recover":(wash_st, wash_end)}


    t0, t1 = (t + subtract), (t + add)

    try:
            print("Loading Data")
            for sig in tqdm(data.segments[0].analogsignals):

                ch = sig.name

                if ch == 'ppmmHg':
                        pressure = sig.load(time_slice =(t0, t1))
                        #pressure = pressure.time_slice(pressure.times[0], t1-(5*pq.s)) #why -5 secs?
                        pressure = elsig.butter(pressure, lowpass_freq = 700.00)
                        #pressure = signal_proc.emg_proc(flex)

                elif (ch == 'Channel bundle (PNA,Vagus) ') | (ch == 'Channel bundle (SN,DR) ') | (ch == 'Channel bundle (DR,SN) ') :
                    #now work with this bundle to tease appart analogsignals
                        DR_sig = sig.load(time_slice = (t0, t1), channel_indexes=[1])
                        SN_sig = sig.load(time_slice = (t0, t1), channel_indexes=[0])

                        if DR_sig.name != "DR":
                            DR_sig1 = SN_sig
                            SN_sig = DR_sig
                            DR_sig = DR_sig1

                        #DR_sig = DR_sig.time_slice(DR_sig.times[0], t1-(5*pq.s))
                        #SN_sig = SN_sig.time_slice(SN_sig.times[0], t1-(5*pq.s))

                        #DR_sig = elsig.butter(DR_sig, highpass_freq = 200.00)#
                        #SN_sig = elsig.butter(SN_sig, highpass_freq = 200.00)# lowpass_freq = 700.00)

                        dorsal_root = DR_sig
                        spinal_nerve = SN_sig

                elif (ch == 'EMG') | (ch == "PND"):
                        phrenic = sig.load(time_slice =(t0, t1))
                        #phrenic = phrenic.time_slice(phrenic.times[0], t1-(5*pq.s))
                        phrenic = elsig.butter(phrenic, lowpass_freq = 700.00)

                elif (ch == "Randall") & (Randall==True):
                        randall = sig.load(time_slice =(t0, t1))
                        #randall = randall.time_slice(randall.times[0], t1-(5*pq.s))

    except:
            print("print unable to parse block")


    return dorsal_root, spinal_nerve, phrenic, pressure, time_slices, randall

def join_data(data1, data2):
    """Takes two tuples of data and joins the data in them"""
    dorsal_root = merge_analogsignals(data1[0], data2[0])

    spinal_nerve = merge_analogsignals(data1[1], data2[1])

    phrenic  = merge_analogsignals(data1[2], data2[2])

    pressure  = merge_analogsignals(data1[3], data2[3])

    randall = merge_analogsignals(data1[5], data2[5])

    control = data1[4]
    gaba = data2[4]

    control  = {k+"_nogaba":v for (k, v) in control.items()}
    gaba = {k+"_gaba":v for (k, v) in gaba.items()}
    t_diff = gaba["control_gaba"][0]-control["recover_nogaba"][1]
    gaba = {k:(v1-t_diff, v2-t_diff) for (k, (v1, v2)) in gaba.items()}
    time_slices = {**control, **gaba}
    t_start = control["control_nogaba"][0]
    #time_slices = {k:(v1-t_start, v2-t_start) for (k, (v1, v2)) in time_slices.items()}

    return dorsal_root, spinal_nerve, phrenic, pressure, time_slices, randall

def merge_analogsignals(sig1, sig2):
    """Merges 2 signals by creating new signal with start of first, end of second"""
    raw_signal = np.concatenate([sig1.magnitude.flatten(), sig2.magnitude.flatten()])
    sig2_times = np.array((sig2.times-(sig2.times[0]-sig1.times[-1]))+(1/sig2.sampling_rate))
    sig1_times = np.array(sig1.times)
    times = np.concatenate([sig1_times, sig2_times])*pq.s

    #times = times-times[0]
    name = sig1.name
    sr = sig1.sampling_rate
    t_start = times[0]
    t_stop = times[-1]
    units = sig1.units

    new_sig = neo.AnalogSignal(signal = raw_signal, units = units, t_start = t_start,
                               t_stop = t_stop, sampling_rate = sr, name = name)

    return new_sig



def subset_times(SN, DR, phrenic, pressure, time, condition):#, condition):
    h, m, s = time.split(':')
    seconds = int(h) * 3600 + int(m) * 60 + int(s)
    phrenic_dt  = 1/np.array(phrenic.sampling_rate)
    DR_dt  = 1/np.array(DR.sampling_rate)
    SN_dt  = 1/np.array(SN.sampling_rate)
    pressure_dt = 1/np.array(pressure.sampling_rate)

    phrenic_start = int(seconds/phrenic_dt)
    SN_start = int(seconds/SN_dt)
    DR_start = int(seconds/DR_dt)
    pressure_start = int(seconds/pressure_dt)

    if condition == "drug":
            subtract = 35
            add = 65
    elif condition == "noxious":
        subtract = 15
        add = 25
    elif condition == "long":
        subtract = 60
        add = 150



    phrenic_subset = phrenic[phrenic_start-int(subtract/phrenic_dt) : phrenic_start+int(add/phrenic_dt)]
    SN_subset = SN[SN_start-int(subtract/SN_dt) : SN_start+int(add/SN_dt)]
    DR_subset = DR[DR_start-int(subtract/DR_dt) : DR_start+int(add/DR_dt)]
    pressure_subset = pressure[pressure_start-int(subtract/pressure_dt) : pressure_start+int(add/pressure_dt)]


    return phrenic_subset, SN_subset, DR_subset, pressure_subset

def create_template(subset, invert):
    """Extract a few spikes, average and create waveform template"""

    bin_starts =[0, 3000, 6000, 9000, 12000, 15000]
    if invert == True:
        signal = subset.magnitude.flatten()-(2*subset.magnitude.flatten())
    else:
        signal = subset.magnitude.flatten()
    peaks =[]
    for start in bin_starts:
        found_peaks = find_peaks(signal[start:start+3000],
                          height = np.max(signal[start:start+3000])*0.8)[0]+start

        peaks.append(found_peaks)
    peaks = np.array(peaks).flatten()
    waveforms = []
    if subset.sampling_rate < 9000:
        add = 120 # 30
        subtract = 100 # 50
    else:
        add = 500
        subtract = 450
    for peak in peaks:
        try:
            peak = peak[0]
        except:
            print("create_template_exception caught")
        #print(peak)

        start = peak-subtract
        end = peak+add
        waveform = signal[start:end]
        waveforms.append(waveform)
    wave_df = pd.DataFrame(waveforms).transpose()
    maximum = wave_df.max()
    mean = maximum.mean()
    std = maximum.std()
    thresh = mean+(1*std)
    screen = maximum[maximum<thresh]
    columns = screen.index
    wave_df = wave_df[columns]
    mean_wave = wave_df.mean(axis=1)
    return mean_wave



def get_ecg(subset, template, invert):
    """Will get starts index of ECG"""
    starts = []
    all_starts = []
    if invert == True:
        raw = np.array(subset.magnitude.flatten())-(2*np.array(subset.magnitude.flatten()))
    else:
        raw = np.array(subset.magnitude).flatten()
    print("Extracting ECG")
    for n in tqdm(range(0, len(subset)-len(template), 5)):
    #for n in tqdm(signals[:-len(template)]):
         subset2 = raw[n:n+len(template)]
         #print(np.isnan(subset).any())
         #print(np.isinf(phrenic_subset).any())
         scaling_factor = np.max(subset2)/np.max(template)
         template2 = template*scaling_factor
         try:
             m, c = np.polyfit(subset2, template2, 1)
         except:
             print("Couldnt converge")

         if (m > 0.5) & (scaling_factor>0.5):
             vals = np.arange(n-30, n, 1)
             boolean = np.isin(vals, np.array(starts)).any()
             if boolean==True:
                 pass
             else:
                 starts.append(n)
         elif m > 0.4:
             vals = np.arange(n-30, n, 1)
             boolean = np.isin(vals, np.array(starts)).any()
             if boolean==True:
                 pass
             else:
                 all_starts.append(n)
    try:
        start_times = subset.times[np.array(starts)]
    except:
        print("Probs no phrenic ecg")
        start_times = np.array([])
    return starts, all_starts, start_times


def remove_ecg_artifact(subset, starts, all_starts, template, noise_mask):
    """Removes ecg artifact from phrenic, DR and SN signals
    Now also removes large noise first"""
    #Remove all starts from phrenic
    #global clean_subset, noisy_subset, mask
    subset_dt  = 1/np.array(subset.sampling_rate)

    if subset.sampling_rate < 8000:
         mask = np.array([np.arange(start,start+len(template)+50,1) for start in all_starts], dtype="int64").flatten()

    elif (subset.sampling_rate > 8000) & (subset.sampling_rate < 9000):
         mask = np.array([np.arange(start,start+len(template)+100,1) for start in starts], dtype="int64").flatten() #change back to all_starts if necessary
    elif subset.sampling_rate > 9000:
        print("wider mask")
        mask = np.array([np.arange(start,start+len(template),1) for start in starts], dtype="int64").flatten()
    try:
        #clean_subset = np.array(subset.magnitude).flatten()
        clean_subset = subset*1
        if len(noise_mask)!=0:
            clean_subset[noise_mask] = 0 * subset.units
        clean_subset[mask] = 0 * subset.units
    except:
        print("Index out of bounds error caught")
        #starts = starts[:-1]
        mask = mask[mask<len(clean_subset)]
        if len(noise_mask)!=0:
            noise_mask = noise_mask[noise_mask<len(clean_subset)]
        clean_subset = subset*1
        clean_subset[noise_mask] = 0 * subset.units
        clean_subset[mask] = 0 * subset.units

    noisy_subset = clean_subset*1
    mean = float(np.mean(clean_subset))
    std = float(np.std(clean_subset))
    new_vals = np.random.normal(mean, std, len(mask)) * subset.units
    new_vals = np.reshape(new_vals,(np.shape(new_vals)[0], 1))
    noisy_subset[mask] = new_vals

    return clean_subset, noisy_subset

def remove_noise(phrenic, DR, SN):
    """This finds any peaks that are over 1000 and removes region around them"""
    signals = [phrenic, DR, SN]
    masks = []
    for sig in signals:
        median = np.abs(np.median(sig.magnitude.flatten()))
        std = median/0.6745
        noise_peaks = find_peaks(np.array(sig.magnitude.flatten()), height = median+(std*6))[0].tolist()
        noise_mask = [np.arange(idx-60, idx+60, 1) for idx in noise_peaks]
        masks.append(np.array(noise_mask))

    return masks


def spike_detection(DR_subset, SN_subset, DR_dt, noisy_DR, noisy_SN,
                    time_slices, DR, SN, min_std, max_std, median=True, DR_min_std=None,
                    DR_max_std=None, low_pass_filter=None):
    """Detects spikes in SN and DR channels - assumes DR and SN have same sampling rate"""
    print("Spike Detection")
    if low_pass_filter != None:
        DR_subset = elsig.butter(DR_subset, highpass_freq=50., lowpass_freq=low_pass_filter,
                                 order=4, filter_function='filtfilt')
        SN_subset = elsig.butter(SN_subset, highpass_freq=50., lowpass_freq=low_pass_filter,
                                 order=4, filter_function='filtfilt')
    DR_sub = DR_subset.magnitude.flatten()
    SN_sub = SN_subset.magnitude.flatten()

    #DR_median = np.abs(np.median(noisy_DR.magnitude))
    if median == True:
        DR_median = np.abs(np.median(DR_subset.magnitude))
        DR_std = (1+DR_median)/0.6745
        SN_median = np.abs(np.median(SN_subset.magnitude))
        SN_std  = (1+SN_median)/0.6745
        print("SN_median is {}".format(SN_median))
        print("DR_median is {}".format(DR_median))
    else:
         DR_median = np.abs(np.median(DR_subset.magnitude))
         SN_median = np.abs(np.median(SN_subset.magnitude))
         DR_std = np.std(DR_subset.magnitude)
         SN_std = np.std(SN_subset.magnitude)
         print("SN_median is {}".format(SN_median))
         print("DR_median is {}".format(DR_median))


    if DR_min_std!=None:
        DR_threshold = DR_median+(DR_std*DR_min_std)
    else:
        DR_threshold = DR_median+(DR_std*min_std)

    SN_threshold = SN_median+(SN_std*min_std)
    print("DR threshold is {}".format(DR_threshold))
    print("SN threshold is {}".format(SN_threshold))
    threshold = np.mean([DR_threshold, SN_threshold])
    print("Mean threshold is {}".format(threshold))

    #DR_peaks = find_peaks(DR_sub-(DR_sub*2), height=(DR_median+(DR_std*min_std), DR_median+(DR_std*max_std)),  prominence=DR_median+(DR_std*min_std))[0]
    #SN_peaks = find_peaks(SN_sub-(SN_sub*2), height=(SN_median+(SN_std*min_std), SN_median+(SN_std*max_std)),  prominence=SN_median+(SN_std*min_std))[0]
    if DR_max_std!=None:
        DR_peaks = find_peaks(DR_sub-(DR_sub*2), height=(DR_threshold, DR_median+(DR_std*DR_max_std)),  prominence=DR_threshold)[0]
    else:
        DR_peaks = find_peaks(DR_sub-(DR_sub*2), height=(DR_threshold, DR_median+(DR_std*max_std)),  prominence=DR_threshold)[0]

    SN_peaks = find_peaks(SN_sub-(SN_sub*2), height=(SN_threshold, SN_median+(SN_std*max_std)),  prominence=SN_threshold)[0]


    DR_times = DR_subset.times[DR_peaks]
    SN_times = SN_subset.times[SN_peaks]

    DR_st = neo.SpikeTrain(times=DR_times, units="sec", t_start = DR_subset.t_start, t_stop = DR_subset.t_stop,
                           sampling_rate = DR_subset.sampling_rate, name="DR_st")
    SN_st = neo.SpikeTrain(times=SN_times, units="sec", t_start = SN_subset.t_start, t_stop = SN_subset.t_stop,
                           sampling_rate = SN_subset.sampling_rate, name="SN_st")

    DR_amps = DR_subset.magnitude[DR_peaks]
    SN_amps = SN_subset.magnitude[SN_peaks]

    DR_filters = []
    SN_filters = []

    for key in time_slices.keys():
        st = time_slices[key][0]
        end = time_slices[key][1]

        DR_filter = ((DR_st>st)&(DR_st<end))
        SN_filter = ((SN_st>st)&(SN_st<end))
        DR_filters.append(DR_filter)
        SN_filters.append(SN_filter)

    DR_mask = [any(tup) for tup in zip(*DR_filters)]
    SN_mask = [any(tup) for tup in zip(*SN_filters)]

    DR_peaks = DR_st[DR_mask]
    SN_peaks = SN_st[SN_mask]
    """if len(DR_filters)==2:
        DR_peaks = DR_st[DR_filters[0]|DR_filters[1]]
        SN_peaks = SN_st[SN_filters[0]|SN_filters[1]]
    elif len(DR_filters)==3:
        DR_peaks = DR_st[DR_filters[0]|DR_filters[1]|DR_filters[2]]
        SN_peaks = SN_st[SN_filters[0]|SN_filters[1]|SN_filters[2]]
    else:
        DR_peaks = DR_st[DR_filters[0]]
        SN_peaks = SN_st[SN_filters[0]]"""



    return DR_peaks, DR_amps, SN_peaks, SN_amps, DR_subset, SN_subset
   # return amps and some spiketrains


def old_extract_spikes(clean_channel, peaks):
    """Extract spikes from a channel and return df containing waveforms that have
    been zeroed to baseline"""
    #global df, waveform
    print("Extracting Spikes")
    idxs = int(0.00144  * clean_channel.sampling_rate)
    peak_idx = np.where(np.isin(clean_channel.times, peaks.times)==True)[0]
    starts = peak_idx-idxs
    ends = peak_idx+idxs
    df = pd.DataFrame()
    for n in range(0, len(starts)):
        try:
            waveform = clean_channel.magnitude.flatten()[starts[n]:ends[n]]
            df[n]=waveform
        except:
            print("probably peak at 0")
            waveform = clean_channel.magnitude.flatten()[0:(2*idxs)]
            df[n]=waveform

    df = df.apply(zero_offset)
    return df



def get_start_ends(clean_channel, peaks):
    idxs = int(0.00144  * clean_channel.sampling_rate)
    peak_idx = np.where(np.isin(clean_channel.times, peaks.times)==True)[0]
    starts = (peak_idx-idxs)
    ends = (peak_idx+idxs)
    magnitude = clean_channel.magnitude.flatten()
    return starts, ends, magnitude
@njit
def fast_extract_spikes(starts, ends, magnitude):
    filtered_list = List()
    for n in range(0, len(starts)):
        waveform = magnitude[starts[n]:ends[n]]
        filtered_list.append(waveform)

    return filtered_list

def extract_spikes(clean_channel, peaks):
    print("Extracting Spikes-new version")
    starts, ends, magnitude = get_start_ends(clean_channel, peaks)
    filtered_list = fast_extract_spikes(starts, ends, magnitude)
    waveforms = np.array(filtered_list)
    df = pd.DataFrame(waveforms).transpose()
    df = df.apply(zero_offset)
    return df



def pd_centers(featuresUsed, centers):
	colNames = list(featuresUsed)
	colNames.append('prediction')

	# Zip with a column called 'prediction' (index)
	Z = [np.append(A, index) for index, A in enumerate(centers)]

	# Convert to pandas data frame for plotting
	P = pd.DataFrame(Z, columns=colNames)
	P['prediction'] = P['prediction'].astype(int)
	return P


def clusters(df, AP, n_clusters):
    #sns.clustermap(df, standard_scale=1)
    if AP==True:
        #df= df.apply(zero_offset)
        df = df.dropna(axis=1)
        X=np.array(df)

    else:
        X = StandardScaler().fit_transform(df)
    #need to loop through cluster numbers to identify optimal clusters
    sse={}
    for n in range(1, 20):
        kmeans = KMeans(n_clusters=n, random_state=1)
        model = kmeans.fit(X)
        sse[n]=model.inertia_
    sse_df=pd.DataFrame(sse, index=[0])
    sse_df=sse_df.transpose()
    sse_df.columns=["Value"]
    kmeans = KMeans(n_clusters, random_state=1) #changed from randoms state=1 to none, 3 better, 5 ok
    model = kmeans.fit(X)
    x_df=pd.DataFrame(X, columns=df.columns)
    test=kmeans.fit_predict(X)
    x_df["Clusters"]= test
    centers = model.cluster_centers_
    P = pd_centers(df.columns, centers)
    pca2 = PCA(n_components=2)
    clustered_df = x_df.copy()
    clusters=x_df.pop("Clusters")
    pca2.fit(x_df)
    #PCA on raw x_df with popped clusters column
    pca2_df=pca2.transform(x_df)
    pca2_df=pd.DataFrame(pca2_df)
    pca2_df.columns=["PC1", "PC2"]
    pca2_df["Clusters"]=clusters

    return clustered_df, pca2_df, sse_df

def spectra(array, time_slices):
    fs = array.sampling_rate
    df=pd.DataFrame()
    f=0
    labels = time_slices.keys()
    print("Calulating power spectrum")
    for label in tqdm(labels):
        start = time_slices[label][0]
        end = time_slices[label][1]
        try:
            subset = array.time_slice(start, end)
        except:
            try:
                start = start+(np.float(1/fs)*pq.s)
                end = end-(np.float(1/fs)*pq.s)
                subset = array.time_slice(start, end)
            except:
                print("Spectra failed")


        f, psd =welch(subset.magnitude.flatten(), fs=fs, window="hanning", nperseg=1024, detrend="constant")
        df[label] = psd
    df.index=f
    return df

def calc_freq(events, time_slices):
    """Calculates instant frequency of events """
    #global freq, df
    df = pd.DataFrame()
    labels = time_slices.keys()
    for label in labels:
        #if label == "control":
        start, end = time_slices[label]
        #elif label == "drug":
        #    start, end = time_slices[label]
        #elif label == "wash":
        #    start, end = time_slices[label]
        try:
            freq = el.instantaneous_rate(events.time_slice(start, end),
                                         np.float(1/events.sampling_rate)*pq.s) #events needs to be spiketrain
            time = np.array(freq.times.rescale(pq.s))

            df[label] = pd.Series(freq.magnitude.flatten(), index=time-time[0])
        except:
            print("Couldnt calc instantneous rate")

    return df


def zero_offset(series):
    index0 = series[0]
    series = series - index0
    return series

def wave_clus_export(DR_df, SN_df, DR_peaks, SN_peaks, DR_subset, SN_subset, filepath):
    print("Saving wave_clus export")
    #DR_spike_times = DR_peaks * (1/float(DR_subset.sampling_rate)*1000)
    #SN_spike_times = SN_peaks * (1/float(SN_subset.sampling_rate)*1000)
    DR_spike_times = DR_peaks.times.rescale("ms")
    SN_spike_times = SN_peaks.times.rescale("ms")

    #DR_df = DR_df.apply(zero_offset)
    #SN_df = SN_df.apply(zero_offset)
    DR_matrix = DR_df.transpose().values
    SN_matrix = SN_df.transpose().values

    DR_mat_dict = {"spikes":DR_matrix, "sr":DR_subset.sampling_rate, "index":DR_spike_times}
    SN_mat_dict = {"spikes":SN_matrix, "sr":SN_subset.sampling_rate, "index":SN_spike_times}

    sio.savemat(filepath+"//_DR_matfile.mat", DR_mat_dict)
    sio.savemat(filepath+"//_SN_matfile.mat", SN_mat_dict)



def create_neo_block(DR_sig, DR_subset, SN_sig, SN_subset, DR_st, SN_st, file, time_slices, DR_df, SN_df, randall):
    """Creates neo block containing  1 segments, 2 channels, 2 analogsignals
    3 epochs, """
    print("Creating block")
    bl=neo.Block(name="Block")
    bl.file_origin = file
    seg=neo.Segment(name = "Seg")
    names=["Dorsal root", "Spinal nerve"]
    for n in range(2):
        chx = neo.ChannelIndex(name=names[n], index=[n], channel_id = n, channel_name = names[n][:6])
        bl.channel_indexes.append(chx)

    conditions = list(time_slices.keys())
    t_slices = list(time_slices.values())
    durations = np.array(t_slices)[:,1]-np.array(t_slices)[:,0]

    epc = neo.Epoch(times=np.array(t_slices)[:,0]*pq.s,
                    durations=durations*pq.s,
                     labels=np.array(conditions, dtype="S"), name = "Epoch")

    seg.epochs.append(epc)

    DR_sig = neo.AnalogSignal(DR_sig.magnitude*DR_sig.units, t_start = DR_sig.t_start, t_stop = DR_sig.t_stop,
                              sampling_rate=DR_sig.sampling_rate, name = "DR sig")

    SN_sig = neo.AnalogSignal(SN_sig.magnitude*SN_sig.units, t_start = SN_sig.t_start, t_stop = DR_sig.t_stop,
                              sampling_rate=SN_sig.sampling_rate, name = "SN sig")


    signals = [DR_sig, SN_sig]
    seg.analogsignals.append(DR_sig)
    seg.analogsignals.append(SN_sig)

    if type(randall) == neo.core.analogsignal.AnalogSignal:
        randall = neo.AnalogSignal(randall.magnitude*randall.units,
                                   t_start = randall.t_start, t_stop = randall.t_stop,
                              sampling_rate=randall.sampling_rate, name = "randall sig")
        signals = [DR_sig, SN_sig, randall]
        seg.analogsignals.append(randall)

    DR_idx = np.where(np.isin(DR_sig.times, DR_st.times))[0]
    SN_idx = np.where(np.isin(SN_sig.times, SN_st.times))[0]


    #idxs = int(0.00144  * DR_sig.sampling_rate)
    #DR_waveforms = np.array([DR_sig[idx-idxs:idx+idxs].magnitude for idx in DR_idx])*DR_sig.units
    #SN_waveforms = np.array([SN_sig[idx-idxs:idx+idxs].magnitude for idx in SN_idx])*SN_sig.units

    DR_waveforms = DR_df.to_numpy().transpose()*DR_sig.units
    SN_waveforms = SN_df.to_numpy().transpose()*SN_sig.units

    DR_st.waveforms = DR_waveforms
    SN_st.waveforms = SN_waveforms



    seg.spiketrains.append(DR_st)
    seg.spiketrains.append(SN_st)
    n=0
    for chx in bl.channel_indexes:
        sig = signals[n]
        chx.analogsignals.append(sig)
        sig.channel_index = chx
        n=n+1

    bl.segments.append(seg)
    return bl


def save_block(bl, time, filepath, use):
    """This should save the block in NIXIO format, a standardised neuroscience format"""
    print("Saving block")
    time = ",".join(time.split(":"))

    """Check file exists, if so delete it"""
    path = filepath+"//_"+time+".blk"
    if os.path.exists(path):
        os.remove(path)
    nixfile = neo.NixIO(filepath+"//_"+time+".blk", mode="rw")
    nixfile.write_block(bl, use)
    nixfile.close()


def main_plot(raw, clean, elbow, freq, cluster_df, pca_df, spectra_df, events, filepath, channel_name, time, condition, time_slices):
    print("Plotting")
    palette = plotting_palette()
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 0.3})
    sns.set_style("white")
    fig = plt.figure(figsize=(15, 7))
    dt = 1.0/float(raw.sampling_rate)

    ax1 = fig.add_subplot(421, ylabel = "Raw Trace", xlabel = "datapoint")
    ax1.set_title(channel_name+time)
    ax2 = fig.add_subplot(422, ylabel = "Within Cluster Sum of Squares", xlabel = "n clusters")
    ax3 = fig.add_subplot(423, ylabel = "Clean trace", xlabel = "seconds")

    ax4 = fig.add_subplot(424, ylabel = "Voltage (uV)", xlabel = "datapoint")
    ax5 = fig.add_subplot(425, ylabel = "Spike Freq (Hz)", xlabel = "seconds")
    ax6 = fig.add_subplot(426, ylabel = "Voltage (uV)", xlabel = "datapoint")
    ax7 = fig.add_subplot(427, ylabel = "Power Spectra", xlabel = "frequency")
    ax8 = fig.add_subplot(428,  ylabel = "")

    ax3.plot(clean.times-clean.times[0], clean.magnitude.flatten(), color = "k")
    ymin, ymax = ax3.get_ylim()
    labels  = list(time_slices.keys())
    #colors=sns.color_palette("bright", n_colors=8)
    for n in range(len(labels)):
        label = labels[n]
        start, end = time_slices[label]
        count = len(events.time_slice(start, end))
        x_pos = (1/(len(labels)+1))*(n+1)
        ax3.text(x_pos, 1, "{} peaks = {}".format(label, str(count)), transform=ax3.transAxes,
                 horizontalalignment='center', fontsize=8)
        ax3.plot([start-time_slices[labels[0]][0], end-time_slices[labels[0]][0]],
         [ymax, ymax], color=palette[n], linewidth=4)

    freq.plot(ax=ax5, linewidth =1)
    ax5.set_ylim(bottom=0, top=freq.max().max()*1.4)


    plt.gcf().subplots_adjust(bottom=0.15)
    ax1.plot(raw.times-raw.times[0], raw.magnitude.flatten(), color ="k",
        linewidth=0.3)


    try:
        peak_idx = np.where(np.isin(clean.times, events.times)==True)[0]
        ax3.scatter(events.times-clean.times[0], clean.magnitude.flatten()[peak_idx], c="r")

    except:
        print("Couldnt plot events")


    spectra_df.plot(ax=ax7, linewidth =1, color = palette)
    try:
        """plot waveforms"""
        cluster_df[0:100].plot(ax=ax4, legend=False)
        #subset = cluster_df[cluster_df.columns[:-1]]
        #subset[:100].transpose().plot(ax=ax4, legend=False)
        cluster_df.groupby("Clusters").mean().transpose().plot(ax=ax6, linewidth =1)
        ax6.legend(loc="upper right")
        sns.scatterplot(ax=ax8, x="PC1", y="PC2", data=pca_df, hue="Clusters")
        elbow.plot(ax=ax2, legend=["Deviation"])

    except:
        print("no clustering")
    sns.despine(top=True, right=True, left=False, bottom=False)
    plt.autoscale()
    #fig.tight_layout()
    fig.savefig(filepath+"//_"+channel_name+".svg")

def save(freq, cluster_df, pca_df, spectra_df, time, channel_name, filepath):
    """Save all data to Excel sheets - one wkbk per channel
    Export clean traces to txt file"""
    print("Saving spectra data to excel")
    time = ",".join(time.split(":"))
    with pd.ExcelWriter(filepath+"//_"+channel_name+"_"+time+".xlsx") as writer:

        spectra_df.to_excel(writer, sheet_name='Spectra')
        #try:
        #    freq.to_excel(writer, sheet_name='SpikeFreq')
        #    cluster_df.to_excel(writer, sheet_name='Clusters')
        #    pca_df.to_excel(writer, sheet_name='PCA')
        #except:
        #    pass



def export_txt(clean_phrenic, clean_DR, clean_SN, time, filepath):
    print("Exporting txt file for Spike2")
    time = ",".join(time.split(":"))
    export_df=pd.DataFrame([clean_DR.magnitude.flatten(), clean_SN.magnitude.flatten()]).transpose()
    export_df.to_csv(filepath+"//_"+time+".txt", header=None, index=None, sep=" ", mode="a")



def spike_train_plot(bl):
    """Plots a raster of spike trains-either takes DR and SN spikes or one of them split into associated clusters"""
    print("Plotting spiketrain")
    for seg in bl.segments:
        print("SEG: " + str(seg.file_origin))
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        ax1.set_title(seg.file_origin)
        ax1.set_ylabel('arbitrary units')
        mint = 0 * pq.s
        maxt = np.inf * pq.s
        for i, asig in enumerate(seg.analogsignals):
            times = asig.times.rescale('s').magnitude
            asig = asig.magnitude
            if i == 0:
                label = "DR"
            else:
                 label = "SN"
            ax1.plot(times, asig, label = label)


        trains = [st.rescale('s').magnitude for st in seg.spiketrains]
        #colors = plt.cm.jet(np.linspace(0, 1, len(seg.spiketrains)))
        #ax2.eventplot(trains, colors=colors)
        ax2.eventplot(trains)
        #ax1.scatter(trains[0], np.repeat(-50, len(trains[0])), c="r")
        lim = ax1.get_xlim()
        ax2.set_xlim(lim)
        ax1.legend()



def burst_detection(clean_phrenic, phrenic_subset, phrenic_dt, time):
    """Integrates noisy phrenic to define start and stop of phrenic bursts"""
    print("Detecting phrenic bursts")
    ##get peaks, cluster close ISIs somehow, add possion surpise if this is not good enough
    #peaks= find_peaks(noisy_phrenic, height=np.mean(noisy_phrenic)+(np.std(noisy_phrenic))*1,
                     # prominence=np.mean(noisy_phrenic)+(np.std(noisy_phrenic))*1)[0]

    peaks= find_peaks(clean_phrenic, height=np.mean(clean_phrenic)+(np.std(clean_phrenic))*1,
                      prominence=np.mean(clean_phrenic)+(np.std(clean_phrenic))*1)[0]


    bandwidth = 1000
    kde = gaussian_kde(peaks, bw_method=bandwidth / peaks.std(ddof=1))
    x=np.arange(len(clean_phrenic))
    test=kde.evaluate(x)
    thresh = test.mean()+(test.std())

    kde_peaks, kde_params = find_peaks(test, height = np.median(test), distance =3000, width=500, rel_height=0.5)
    plt.figure()
    plt.plot(test)
    plt.scatter(kde_peaks, test[kde_peaks], c="r")
    plt.figure()
    #plt.plot(noisy_phrenic)
    plt.plot(clean_phrenic)
    plt.scatter(kde_peaks, test[kde_peaks], c="r")
    starts = kde_params["left_ips"]
    ends = kde_params["right_ips"]

    for start in starts:
        plt.plot([start, start], [0, 200], c="b")
    for end in ends:
        plt.plot([end, end], [0, 200], c="r")

    # find burst duration and burst rate (per minute), (area under curve too?")
    ten = int(10/phrenic_dt)
    fifteen = int(15/phrenic_dt)
    twenty5 = int(25/phrenic_dt)
    thirty = int(30/phrenic_dt)

    control_bursts = [(starts, ends) for starts, ends in zip(starts, ends) if (starts>ten)&(starts<fifteen)]
    drug_bursts = [(starts, ends) for starts, ends in zip(starts, ends) if (starts>fifteen)&(starts<twenty5)]
    wash_bursts = [(starts, ends) for starts, ends in zip(starts, ends) if (starts>thirty)]

    control_rate = len(control_bursts)*6
    drug_rate = len(drug_bursts)*6
    wash_rate =len(wash_bursts)*6

    control_durations = [(burst[1]-burst[0])*phrenic_dt for burst in control_bursts]
    drug_durations = [(burst[1]-burst[0])*phrenic_dt for burst in drug_bursts]
    wash_durations = [(burst[1]-burst[0])*phrenic_dt for burst in wash_bursts]

    burst_rate_df = pd.DataFrame([control_rate, drug_rate, wash_rate]).transpose()
    duration_df = pd.DataFrame([control_durations, drug_durations, wash_durations]).transpose()
    columns = ["control", "drug", "wash"]
    burst_rate_df.columns = columns
    duration_df.columns = columns
    time = ",".join(time.split(":"))
    with pd.ExcelWriter("./"+"burst_analysis_"+time+".xlsx") as writer:
        burst_rate_df.to_excel(writer, sheet_name='Burst_rate')
        duration_df.to_excel(writer, sheet_name='Burst_durations')

    return burst_rate_df, duration_df



"""Post analysis working with bl structure"""
def wave_clus_import(mat_file):
    mat = sio.loadmat(mat_file)
    mat["cluster_class"]
    clusters = pd.DataFrame(mat["cluster_class"]).iloc[:,0]
    return clusters

def read_block(nixfile):
    """Reads nixfile and returns bl, st1, st2"""
    print("Reading block")
    reader = neo.NixIO(nixfile)
    bl=reader.read(lazy=False)[0]
    st1 = bl.segments[0].spiketrains[1]
    st2 = bl.segments[0].spiketrains[0]
    SN = bl.segments[0].analogsignals[1]
    DR = bl.segments[0].analogsignals[0]

    try:
        randall = bl.segments[0].analogsignals[2]
    except:
        print("no randall trace")
        randall = None

    return bl, st1, st2, SN, DR, randall

def shift_time(series, mean_lat):
    #print(series)
    random_intervals = np.abs(np.random.normal(loc= mean_lat, scale = mean_lat/4, size=series.shape[0]))
    max_interval = random_intervals.max()
    series = series + random_intervals
    return series, max_interval



def spike_latency(st1, st2, sim, slow_c):
       """
        Create a matrix of all time differences between each spike, rows are st1, columns are st2
        The minimum value of the row is the closest st2 spike to the spike that is the row in st1

       """
       st1_matrix = np.matrix([st1]*len(st2))
       st1_matrix = st1_matrix.transpose()
       st2_matrix = np.matrix([st2]*len(st1))
       result = st2_matrix-st1_matrix

       df = pd.DataFrame(result)
       subset =df[(df>=0)&(df<slow_c)]
       df = df[df>=0]
       if sim ==False:
           fig1=plt.figure(figsize = (15, 7))
           ax1=fig1.add_subplot(211)
           ax2=fig1.add_subplot(212)
           sns.heatmap(np.abs(result), ax=ax1)
            #this subsets based on 10 ms need to add slowest latency here
           for ax in [ax1, ax2]:
                ax.set_ylabel("Spinal Nerve Spiketrain")
                ax.set_xlabel("Dorsal Root Spiketrain")
           sns.heatmap(subset, ax=ax2)
           fig1.tight_layout()
       else:
           fig1 = None
       #fig1.savefig("./latency_exmaple.png")

       return df, fig1

def waveform_correlation(st1, st2, sim):
    """loops through st2_wv_df and performs correlation with each and every columns in st1_wv_df"""
    st1_wv = st1.waveforms
    st2_wv = st2.waveforms

    st1_df = pd.DataFrame(np.reshape(st1_wv, np.shape(st1_wv)[:2])).transpose()
    st2_df = pd.DataFrame(np.reshape(st2_wv, np.shape(st2_wv)[:2])).transpose()

    #spit out measurements of these like width and amplitude

    corr_df = pd.DataFrame()
    subset2 = corr_df[corr_df>0.9]
    #print("Calculating waveform correlation")
    #for col in tqdm(st1_df.columns): #loop through each central spike waveform
    #    corr = st2_df.corrwith(st1_df[col])
    #    corr_df[col] = corr #each row in corr df represents a perip spike and each row is central
    if sim == False:
        #corr_df = corr_df.transpose()
        fig1=plt.figure(figsize = (15, 7))
        ax1=fig1.add_subplot(211)
        ax2=fig1.add_subplot(212)
        for ax in [ax1, ax2]:
            ax.set_ylabel("Spinal Nerve Spiketrain")
            ax.set_xlabel("Dorsal Root Spiketrain")
        #corr_df=corr_df.transpose()
        #sns.heatmap(corr_df, ax=ax1)

        #sns.heatmap(subset2, ax=ax2)
        fig1.tight_layout()
    else:
        fig1= None
    #fig1.savefig("./waveform_correlation_examples.png")
    return subset2, st1_df, st2_df, fig1


def create_units(st, st_wv, clusters_dict, bl, chx_index):
    """Quick way to make spiketrains according to cluster
       Use clusters dict for spinal units and dorsal dict for dorsal units
       Chx_index[0] is dorsal root, chx_index[1] is spinal """
    df = pd.DataFrame(st.times, columns=["Times"])
    df["Clusters"]= df.index.map(clusters_dict)
    try:
        df["Clusters"] = df["Clusters"].astype("int64")
    except:
        print("couldnt convert to int")
    st_wv_T = st_wv.transpose()
    st_wv_T["Clusters"] = st_wv_T.index.map(clusters_dict)


    for n in np.sort(df.Clusters.unique()):
        try:
            if np.isnan(n):
                pass
            else:
                subset = df[df.Clusters==n]
                times = subset.Times.to_list() * pq.s
                unit_st = neo.SpikeTrain(times, units="sec",
                                         t_start = st.t_start,
                                         t_stop = st.t_stop,
                                         sampling_rate = st.sampling_rate)

                wv_subset = st_wv_T[st_wv_T.Clusters==n]
                wvs = np.array(wv_subset.iloc[:, :-1].values.tolist())*pq.uV # to drop clusters column
                unit_st.waveforms = wvs
                unit = neo.Unit(name = str(int(n)))
                unit.spiketrains.append(unit_st)

                bl.channel_indexes[chx_index].units.append(unit)
        except:
            print("create units error")
    return bl

def get_epoch(bl):
    starts= bl.segments[0].epochs[0].times
    durations = bl.segments[0].epochs[0].durations
    ends = starts + durations
    conditions = bl.segments[0].epochs[0].labels.astype("U")
    t_slices = {}
    t = zip(starts, ends, conditions)
    for start, end, condition in t:
        t_slices[condition] = (start, end)
    return t_slices

def match(st1, st2, cluster_labels, latency, new_method, new_new_method, slow_c, t_slices, distance):
    """This function matches the spikes based on latency"""

    st1_isi = el.isi(st1)
    st2_isi = el.isi(st2)
    min_idx = latency.idxmin(axis=0) #this seems to find the min row index (st1) for each column (st2)
    min_latency = latency.min(axis=0)
    print("Min_idx length = {}".format(min_idx.shape[0]))
    na_min_idx = min_idx.dropna()
    print("New min_idx length = {}".format(na_min_idx.shape[0]))

    latency_array = np.array(latency)
    latency_array[np.array(na_min_idx, dtype="int64"), np.array(na_min_idx.index, dtype="int64")] = np.nan
    second_latency = pd.DataFrame(latency_array)
    second_min_idx = second_latency.idxmin(axis=0)
    second_min_lat = second_latency.min(axis=0)

    df = pd.DataFrame()

    df["st1"] = st1.times
    df["st1_isi"] = pd.Series(st1_isi)
    df["st1_label"] = cluster_labels
    ##How to map st2 and st2_isi here
    st2_idx = np.arange(len(st2.times))
    df2 = pd.DataFrame([st2.times, st2_isi, min_idx, min_latency, second_min_idx, second_min_lat, st2_idx]).transpose()

    df2.columns = ["st2", "st2_isi", "min_idx", "min_latency", "2_min_idx", "2_min_lat", "st2_idx"]
    assigned = df2[df2.min_idx.notnull()]
    unassigned = df2[df2.min_idx.isnull()]
    new_assigned = assigned.drop_duplicates("min_idx", keep="first")
    new_unassigned = pd.concat([unassigned, assigned[assigned.min_idx.duplicated("first")]], axis=0) #last seems to work
    new_assigned.index = new_assigned.min_idx

    df = pd.concat([df, new_assigned], axis=1)

    ##outliers
    if new_method==True:
        """Pick out ISI that are faster than the slowest conduction and that have latency greater than ISI"""
        isi_subset = df[(df.st1_isi<slow_c)&(df.st2.notnull())&(df.min_latency>df.st1_isi)]
        new_idx = df.min_idx.copy()
        new_idx.loc[isi_subset.index] = isi_subset["2_min_idx"]
        df["new_idx"] = new_idx
        df["st1_label"] = df.new_idx.map(df.st1_label)


        if new_new_method==True:
            new_subset=df[(df.st1_isi<slow_c)&(df.st2.notnull())&
                          (df.min_latency>df.st2_isi)&
                          (df.st2_isi<df.st1_isi)&
                          (df["2_min_lat"]>df.st1_isi)]
            new_new_idx = new_idx.copy()
            new_new_idx.loc[new_subset.index] = new_subset["2_min_idx"]
            df["new_new_idx"] = new_new_idx
            df["new_new_label"] = df.new_new_idx.map(df.st1_label) #scarp the new_new_label stuff i think

    df.min_latency = (df.min_latency * 1e3).astype("float64")
    df["Condition"]  = [np.nan]*df.shape[0]
    for n in t_slices.keys():
        start, stop = t_slices[n]
        start, stop  = float(start), float(stop)
        df.loc[(df.st1>start)&(df.st1<stop), "Condition"] = str(n)
    df["Fibre"]  = [np.nan]*df.shape[0]
    fast = (distance/10)*1e3
    slow =(distance/1.2)*1e3
    fibre_types = {"alpha/beta": (fast, 0), "delta":(slow, fast), "c":(slow*100, slow)}

    for n in fibre_types:
        start, stop = fibre_types[n]
        #start, stop  = float(start), float(stop)
        df.loc[(df.min_latency<start)&(df.min_latency>stop), "Fibre"] = str(n)

    #fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, figsize = (15, 7))atc
    fig, axes = plt.subplots(nrows=2,ncols =2, figsize = (15, 7))
    (ax1, ax2), (ax3, ax4) = axes
    lineoffsets1 = [1, 0]
    colors = sns.color_palette("bright", n_colors = 20)
    print("Plotting events")
    for n in tqdm(range(len(df.st1_label.unique()))):

        label = df.st1_label.unique()[n]
        subset = df[df.st1_label==label]
        ax1.eventplot([subset.st1, subset.st2], color=[colors[n], colors[n]], lineoffsets =lineoffsets1, linewidth=1, label=label)

    ax1.eventplot(new_unassigned.st2, color = "k", lineoffsets=lineoffsets1[1], linewidth=0.5)
    #ax1.legend(df.st1_label.unique())
    #st1_counts = df.st1_label.value_counts()
    st1_counts = df.groupby("Condition").st1_label.value_counts()
    #st2_counts = df[df.st2_idx.notnull()].st1_label.value_counts()
    st2_counts = df[df.st2_idx.notnull()].groupby("Condition").st1_label.value_counts()
    match_counts = pd.concat([st1_counts, st2_counts], axis=1)
    match_counts.columns = ["st1_units", "st2_units"]
    match_counts["Percentage_Remaining"] = (match_counts.st2_units/match_counts.st1_units)*100
    match_counts["Percentage_Change"] = ((match_counts.st2_units-match_counts.st1_units)/match_counts.st1_units)*100
    match_counts.plot(y="Percentage_Change", kind="bar", ax=ax2, color = palette, legend=False)
    ax2.set_ylabel("Percentage Change")

    ##add latency kde to plot
    na_subset = df[df.st2_idx.notnull()]
    for n in tqdm(range(len(na_subset.st1_label.unique()))):
        label = na_subset.st1_label.unique()[n]
        subset = na_subset[(na_subset.st1_label == label)]
        sns.distplot(subset.min_latency, hist=False, color= palette[n], kde_kws={"shade":True}, ax=ax3, label=label)
    ax3.set_ylabel("Kernel Density")
    ax3.set_xlabel("Latency")
    sns.barplot(data=na_subset, x="st1_label", y="min_latency", ax=ax4)
    ax4.set_ylabel("Latency (ms)")
    ax4.set_xlabel("Cluster")
    fig.tight_layout()

    return df, unassigned, fig, match_counts

def new_correlation(wv, DR_slice):
    idx = (int(0.00144  * DR_slice.sampling_rate))
    window = 2*idx
    slice = DR_slice.magnitude.flatten()[:-window]
    range_idx = np.arange(len(slice))
    df = pd.DataFrame([slice[n:n+window] for n in range_idx]).transpose()
    df = df.dropna(axis=1)
    df["spike"] = pd.Series(wv)
    #df = df/np.sqrt((df**2).sum())
    result  = df.iloc[:, :-1].corrwith(df.spike)
    try:
        min_idx = result[result>0.8].index[0] # used greater than 0.8 here
        min_latency = pd.Series(float(min_idx/DR_slice.sampling_rate)*1e3)
        min_idx_corr = pd.Series(result[result>0.8].iloc[0])
    except:
        print("new corr failed")
        min_latency = pd.Series(np.nan)
        min_idx_corr = pd.Series(np.nan)
        #min_idx = result[result>0.8].index[0] # used greater than 0.8 here
        #min_latency = float(min_idx/DR_slice.sampling_rate)
        #min_idx_corr = result[result>0.8].iloc[0]
    return min_latency, min_idx_corr


def cross_correlation(st1, SN, DR, cluster_labels, slow_c):
    """This function computes the cross correlation of SN with DR for 100ms for
    every spike in SN"""
    if slow_c >2e-3:
        duration = ((slow_c)*pq.s)/50
    else:
        duration = (slow_c)*pq.s
    cross_corr = pd.DataFrame()
    print("Computing cross correlation")
    for time in tqdm(st1.times):
        if (time+duration < SN.t_stop)&((time-(1e-3*pq.s))>SN.t_start):

            SN_slice = SN.time_slice(time-(1e-3*pq.s), time+duration)
            DR_slice = DR.time_slice(time-(1e-3*pq.s), time+duration)
        else:
            SN_slice = SN.time_slice(time-(1e-3*pq.s), SN.t_stop-(1e-3*pq.s))
            DR_slice = DR.time_slice(time-(1e-3*pq.s), DR.t_stop-(1e-3*pq.s))


        wv = st1[st1.times==time].waveforms.flatten()
        #min_latency, min_idx_corr = new_correlation(wv, DR_slice)

        df = pd.DataFrame([SN_slice.magnitude.flatten(), DR_slice.magnitude.flatten()]).transpose()

        df = df/np.sqrt((df**2).sum())
        df.index = SN_slice.times
        x1 = df.iloc[:, 0]
        x2 = df.iloc[:, 1]

        corr = correlate(x1, x2, mode="full")
        index = pd.Series(df.index)-df.index[0]
        new_index = pd.concat([-index[::-1], index])
        final_corr = pd.Series(corr)

        if index.shape[0]%2 > 0:
            final_corr.index = new_index[1:]
        else:
            final_corr.index = new_index[1:]
        #fig, (ax1, ax2) = plt.subplots(nrows =2)
        #final_corr.plot(ax=ax2)
        #df.plot(ax=ax1)
        #ax1.scatter(time, df[df.index>=time].iloc[0,0])
        final_corr = final_corr.to_frame(str(float(time)))
        max_idx_val = pd.concat([final_corr.idxmax(), final_corr.max()], axis=1)
                                 #min_latency, min_idx_corr],
        cross_corr = pd.concat([cross_corr, max_idx_val], axis=0)
    cross_corr.reset_index(inplace=True)
    #cross_corr.columns = ["st1", "time_lag", "Corr", "min_latency", "new_corr"]
    cross_corr.columns = ["st1", "time_lag", "Corr"]
    cross_corr["st1_label"] = cluster_labels
    cross_corr.time_lag = cross_corr.time_lag*1e3

    """Produce plot describing and bin corr and time lag"""
    #palette = sns.color_palette("bright", n_colors = 8)
    fig, axes = plt.subplots(nrows =2, ncols = 2, figsize = (15, 7))
    (ax1, ax2), (ax3, ax4) = axes
    print("Plotting kernel density estimates")
    for n in tqdm(range(len(cross_corr.st1_label.unique()))):
        try:
            label = cross_corr.st1_label.unique()[n]
            subset = cross_corr[(cross_corr.st1_label == label)]
            sns.distplot(subset.Corr, hist=False, color= palette[n], kde_kws={"shade":True}, ax=ax1, label=label)
            sns.distplot(subset[subset.Corr>0.85].time_lag, hist=False, color= palette[n], kde_kws={"shade":True}, ax=ax2, label=label)
        except:
            print("Couldnt plot KDE")
    fig.tight_layout()

    count_total =cross_corr.groupby("st1_label").Corr.count()
    count_80 = (cross_corr[(cross_corr.Corr>0.8)&(cross_corr.Corr<0.9)].groupby("st1_label").Corr.count()/count_total)*100
    count_90 = (cross_corr[cross_corr.Corr>=0.9].groupby("st1_label").Corr.count()/count_total)*100
    count_rest = (cross_corr[cross_corr.Corr<0.8].groupby("st1_label").Corr.count()/count_total)*100
    count_df = pd.concat([count_total, count_rest, count_80, count_90], axis=1)
    count_df.columns = ["total", "<0.8 (Low Correlation)", "0.8-0.9", ">0.9 (High correlation)"]
    count_df[count_df.columns[1:]].plot(kind="bar", stacked=True, color=palette, ax=ax3)
    try:
        sns.barplot(data = cross_corr[cross_corr.Corr>0.9], x = "st1_label", y= "time_lag", ci="sd", ax=ax4)
    except:
        print("Couldnt plot cross corr bar plot")
    ax4.set_ylabel("Time Lag (ms)")
    ax4.set_xlabel("Unit")
    ax3.set_xlabel("Unit")
    ax1.set_ylabel("Kernel Density")
    ax2.set_ylabel("Kernel Density")
    fig.tight_layout()
    #savefig
    try:
        time_lag_desc = cross_corr[cross_corr.Corr>=0.9].groupby("st1_label").time_lag.describe()
    except:
        try:
            time_lag_desc = cross_corr[cross_corr.Corr>=0.8].groupby("st1_label").time_lag.describe()
        except:
            print("no time lag desc")
            time_lag_desc = pd.DataFrame()
    return cross_corr, time_lag_desc, count_df, fig

def rel_freq(series):
    max_val = series.max()
    if max_val/10 >1:
        max_val = int(np.ceil(series.max()/100))*100
    else:
        max_val = max_val
    #print(max_val)
    bin_width = max_val/40
    bins = np.arange(0, max_val, bin_width)

    hist, edges = np.histogram(series, bins)
    rel_freq = hist/float(hist.sum())
    return pd.Series(rel_freq, index = bins[:-1])

def spike_loop(distance):
    """Loops spike_analysis"""
    folder = filedialog.askdirectory(initialdir = os.getcwd(), title = "Please select folder containing result folder")
    sub_dirs = os.listdir(folder)
    failures = []
    for sub_dir in tqdm(sub_dirs):
        files = os.listdir(os.path.join(folder, sub_dir))
        try:
            matfile = [file for file in files if file =="times__SN_matfile.mat"][0]
            matfile_path = os.path.join(folder, sub_dir, matfile)
            nixfile = [file for file in files if ("blk" in file)&("NEW" not in file)][0]
            nixfile_path = os.path.join(folder, sub_dir, nixfile)

            spike_analysis(matfile_path, nixfile_path, distance)
        except:
            print("failed")
            failures.append(sub_dir)
    print("These timepoint folders failed {}".format(failures))





def spike_analysis(matfile, nixfile, distance):
    """Select wave_clus matlab file output
    Select nix file"""
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 0.3})
    sns.set_style("white")
    try:
        palette = plotting_palette()
    except:
        print("error setting custom palette")
        palette = sns.color_palette("colorblind")
    if (matfile == None)&(nixfile ==None):
        matfile = filedialog.askopenfilename(initialdir = os.getcwd(), title = "Please select waveclus matlab file")

        nixfile = filedialog.askopenfilename(initialdir = os.getcwd(), title = "Please select nix file")

    filepath = os.path.dirname(matfile)
    cluster_labels = wave_clus_import(matfile)
    cluster_labels = cluster_labels.astype("int64")
    bl, st1, st2, SN, DR, randall = read_block(nixfile)
    """if randall in bl - might have to downsample unit rates or calculate them at
    randall freq"""

    t_slices = get_epoch(bl)
    dt = float(1/st1.sampling_rate)
    slow_c = distance/0.1 #0.1 m/s is slowest expected c fibre velocity
    latency, latency_fig = spike_latency(st1, st2, False, slow_c)

    latency_subset = latency[(latency>=0)&(latency<slow_c)]
    match_df, unassigned_df, match_fig, match_counts = match(st1, st2, cluster_labels, latency_subset, False, False, slow_c, t_slices, distance)
    cross_corr, time_lag_desc, count_df, corr_fig = cross_correlation(st1, SN, DR, cluster_labels, slow_c)

    waveform_corr, st1_df, st2_df, waveform_corr_fig = waveform_correlation(st1, st2, False)
    clusters_dict = cluster_labels.to_dict()
    match_subset = match_df.loc[match_df.min_idx.notnull(), ["st1_label", "st2_idx"]]
    match_subset.set_index("st2_idx", inplace=True)
    try:
        match_subset = match_subset.astype("int64")
    except:
        print("couldnt convert match_subset to int-could be errors in dorsal dict")

    dorsal_dict = match_subset.st1_label.to_dict()

    #st1_melt, st2_melt, dorsal_dict, sscore_df = spike_similarity_score(waveform_corr, spike_latencies, st1_df, st2_df, clusters_dict, dt)
    #st1_df.index = np.arange(len(st1_df.index))/SN.sampling_rate
    #st2_df.index = np.arange(len(st2_df.index))/DR.sampling_rate

    st1_feat = get_spike_features(st1, st1_df, clusters_dict)
    st2_feat = get_spike_features(st2, st2_df, dorsal_dict)
    st1_mean_feat = st1_feat.groupby("Clusters")["Width", "Prominence", "Peak_height"].describe()
    st2_mean_feat = st2_feat.groupby("Clusters")["Width", "Prominence", "Peak_height"].describe()


    rel_st1 = st1_feat[["Width", "Prominence", "Peak_height"]].apply(rel_freq)
    rel_st2 = st2_feat[["Width", "Prominence", "Peak_height"]].apply(rel_freq)

    cumsum_st1 = rel_st1.cumsum(axis=0)
    cumsum_st2 = rel_st2.cumsum(axis=0)

    print("Measuring waveforms")
    shape_fig = plt.figure(figsize = (15, 7))
    for col in st1_feat.columns[0:3]:
        n=int(np.where(st1_feat.columns == col)[0])+1
        ax1=shape_fig.add_subplot(2, 3, n)
        ax2 = shape_fig.add_subplot(2, 3, n+3)
        sns.barplot(x="Clusters", y=col, data = st1_feat, ax=ax1)
        sns.barplot(x="Clusters", y=col, data = st2_feat, ax=ax2)
        ax1.set_ylim(bottom = 0)
        ax2.set_ylim(bottom = 0)
    shape_fig.tight_layout()

    new_bl = create_units(st1, st1_df, clusters_dict, bl, 1) # adds spinal nerve units to block
    new_bl = create_units(st2, st2_df, dorsal_dict, new_bl, 0) # adds dorsal root units to block


    raster_fig = plt.figure(figsize = (15, 7))
    raster_ax1 = raster_fig.add_subplot(211)
    raster_ax2 = raster_fig.add_subplot(212)
    isi_fig, isi_axes = plt.subplots(nrows = 2, ncols = len(st1_feat.Clusters.unique()),
                                     figsize = (15, 7))
    SN_axes = isi_axes[0]
    DR_axes = isi_axes[1]
    SN_sts, SN_isis = unit_plot(new_bl, 1, raster_ax1, SN_axes, 0, palette)
    DR_sts, DR_isis = unit_plot(new_bl, 0, raster_ax2, DR_axes, len(SN_sts), palette)
    SN_isis_desc = SN_isis[SN_isis<3e-3].describe()
    DR_isis_desc = DR_isis[DR_isis<3e-3].describe()
    isi_fig.tight_layout()
    raster_fig.tight_layout()
    print("Getting spiketrain stats")
    global_stats, global_rates = get_statistics(st1, st2, t_slices) # full_df, control_df, test_df, recover_df
    unit_stats, unit_rates= unit_statistics(SN_sts, DR_sts, t_slices)
    SN_units_df = pd.DataFrame()
    DR_units_df = pd.DataFrame()
    for key in unit_rates.keys():
        try:
            SN_rate = pd.Series(unit_rates[key]["Full"][0].magnitude.flatten(),
            index = unit_rates[key]["Full"][0].times)
            DR_rate = pd.Series(unit_rates[key]["Full"][1].magnitude.flatten(),
            index = unit_rates[key]["Full"][1].times)
            SN_units_df[key] = SN_rate
            DR_units_df[key] = DR_rate
        except:
            print("probs no rate calculated")

    rate_plot(st1, st2, global_rates, filepath)

    global_rates_values = list(global_rates.values())
    global_rates_keys = list(global_rates.keys())
    global_rates = [y.magnitude.flatten() for z in global_rates_values for y in z]

    rate_df = pd.DataFrame(global_rates).transpose()
    #if rate_df.shape[0]>1048576:
    rate_df = rate_df.dropna()
    columns = [[x+"_SN", x+"_DR"] for x in global_rates_keys]
    columns = [y for z in columns for y in z]
    rate_df.columns = columns
    """Make and save figures"""
    """SN spikes"""

    

    st1_df["Index"] = ((st1_df.index)/st1.sampling_rate)*1e3
    st1_melt = st1_df.melt(id_vars = "Index")
    st1_melt["Clusters"] = st1_melt.variable.map(clusters_dict)

    st2_df["Index"] = (st2_df.index/st2.sampling_rate)*1e3
    st2_melt = st2_df.melt(id_vars = "Index")
    st2_melt["Clusters"] = st2_melt.variable.map(dorsal_dict)
    st2_melt.loc[st2_melt.Clusters.isnull(), "Clusters"] = "Unassigned"
    #palette = sns.color_palette("bright", n_colors = 20)
    no_units = len(st1_melt.Clusters.dropna().unique())
    col_no = round(no_units/2)
    g=sns.relplot(x="Index", y="value", data=st1_melt, hue="Clusters", kind="line", col="Clusters", col_wrap=col_no,
                palette = palette[:len(st1_melt.Clusters.dropna().unique())], units = "variable", estimator=None)
    #average plots
    f=sns.relplot(x="Index", y="value", data=st1_melt, hue="Clusters", kind="line", col="Clusters", col_wrap=col_no,
                palette = palette[:len(st1_melt.Clusters.dropna().unique())])


    e=sns.relplot(x="Index", y="value", data=st2_melt, hue="Clusters", kind="line", col="Clusters", col_wrap=col_no,
                palette = palette[:len(st2_melt.Clusters.unique())], units = "variable", estimator=None)
    #average plots
    d=sns.relplot(x="Index", y="value", data=st2_melt, hue="Clusters", kind="line", col="Clusters", col_wrap=col_no,
                palette = palette[:len(st2_melt.Clusters.unique())])

    for facet in [g, f, e, d]:
        facet.set(xlabel ="Time (ms)", ylabel ="Voltage (uV)")

    print("Saving figures")

    figures = [g.fig, f.fig, e.fig, d.fig, isi_fig, raster_fig, shape_fig, waveform_corr_fig, latency_fig, match_fig, corr_fig]
    try:
        [fig.suptitle(str(datetime.timedelta(seconds = float(t_slices["drug"][0])))) for fig in figures]
    except:
        print("no drug")
    g.savefig(filepath+"//SN_spikes_ind.svg")
    f.savefig(filepath+"//SN_spikes_mean.svg")
    e.savefig(filepath+"//DR_spikes_ind.svg")
    d.savefig(filepath+"//DR_spikes_mean.svg")

    isi_fig.savefig(filepath+"//isi.svg")
    raster_fig.savefig(filepath+"//raster_fig.svg")
    shape_fig.savefig(filepath+"//shape_fig.svg")
    match_fig.savefig(filepath+"//match_fig.svg")
    corr_fig.savefig(filepath+"//corr_fig.svg")
    #waveform_corr_fig.savefig(filepath+"//waveform_corr_fig.svg")
    #latency_fig.savefig(filepath+"//latency_fig.svg")

    g.savefig(filepath+"//SN_spikes_ind.png")
    f.savefig(filepath+"//SN_spikes_mean.png")
    e.savefig(filepath+"//DR_spikes_ind.png")
    d.savefig(filepath+"//DR_spikes_mean.png")

    isi_fig.savefig(filepath+"//isi.png")
    raster_fig.savefig(filepath+"//raster_fig.png")
    shape_fig.savefig(filepath+"//shape_fig.png")
    match_fig.savefig(filepath+"//match_fig.png")
    corr_fig.savefig(filepath+"//corr_fig.png")
    waveform_corr_fig.savefig(filepath+"//waveform_corr_fig.png")
    latency_fig.savefig(filepath+"//latency_fig.png")

    [plt.close(fig) for fig in figures[:-2]]
    [fig.clf() for fig in figures[:-2]]

    """Summarising latencies"""
    latency_desc = match_df.groupby("st1_label").min_latency.describe()
    fast = (distance/10)*1e3
    slow =(distance/1.2)*1e3
    fibre_types = {"alpha/beta": (fast, 0), "delta":(slow, fast), "c":(slow*100, slow)}
    latency_desc["Fibre"] = [np.nan]*latency_desc.shape[0]

    for n in fibre_types:
        start, stop = fibre_types[n]
        #start, stop  = float(start), float(stop)
        latency_desc.loc[(latency_desc["mean"]<start)&(latency_desc["mean"]>stop), "Fibre"] = str(n)

    fibre_desc = match_df.groupby("st1_label").Fibre.describe()


    """Write dataframes to sheets of an excel file"""
    print("Writing data to excel")
    with pd.ExcelWriter(filepath+"//block_analysis.xlsx") as writer:
        st1_feat.to_excel(writer, sheet_name='SN_features')
        st2_feat.to_excel(writer, sheet_name='DR_features')

        st1_mean_feat.to_excel(writer, sheet_name='SN_mean_features')
        st2_mean_feat.to_excel(writer, sheet_name='DR_mean_features')

        st1_df.to_excel(writer, sheet_name='SN_waveforms')
        st2_df.to_excel(writer, sheet_name="DR_waveforms")

        global_stats.to_excel(writer, sheet_name="Global_spike_stats")
        unit_stats.to_excel(writer, sheet_name="Unit_spike_stats")

        SN_isis.to_excel(writer, sheet_name="SN_isis")
        DR_isis.to_excel(writer, sheet_name="DR_isis")

        rate_df.to_excel(writer, sheet_name= "inst firing rates")
        match_df.to_excel(writer, sheet_name= "match_df")
        unassigned_df.to_excel(writer, sheet_name= "unassigned_df")
        match_counts.to_excel(writer, sheet_name= "match_counts")
        latency_desc.to_excel(writer, sheet_name = "latency_desc")
        fibre_desc.to_excel(writer, sheet_name = "fibre_desc")

        cross_corr.to_excel(writer, sheet_name= "cross_correlation")
        time_lag_desc.to_excel(writer, sheet_name= "time_lag_summary")
        count_df.to_excel(writer, sheet_name= "cross_corr_counts")

        rel_st1.to_excel(writer, sheet_name = "rel_st1")
        rel_st2.to_excel(writer, sheet_name = "rel_st2")
        cumsum_st1.to_excel(writer, sheet_name = "cumsum_st1")
        cumsum_st2.to_excel(writer, sheet_name = "cumsum_st2")

        SN_isis_desc.to_excel(writer, sheet_name = "SN_isi<3ms")
        DR_isis_desc.to_excel(writer, sheet_name = "DR_isi<3ms")
        try:
            SN_units_df.to_excel(writer, sheet_name = "SN_unit_rates")
            DR_units_df.to_excel(writer, sheet_name = "DR_unit_rates")
        except:
            print("Couldn't save rates to excel")
            #save as csv

    """Save new block"""
     #work on this
    new_bl.annotations = {}
    new_bl.annotations["rec_datetime"] = new_bl.rec_datetime.strftime("%m/%d/%Y, %H:%M:%S")
    new_bl.segments[0].annotations = {}

    for n, seg in enumerate(new_bl.segments):
        seg.annotations["rec_datetime"] = seg.rec_datetime.strftime("%m/%d/%Y, %H:%M:%S")
        seg.rec_datetime = None

    new_bl.rec_datetime = None
    new_bl.segments[0].analogsignals[0].annotations = {k:v for k,v in new_bl.segments[0].analogsignals[0].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}
    new_bl.segments[0].analogsignals[1].annotations = {k:v for k,v in new_bl.segments[0].analogsignals[1].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}

    try:
        new_bl.segments[0].analogsignals[2].annotations = {k:v for k,v in new_bl.segments[0].analogsignals[1].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}
    except:
        print("Couldnt add randall")

    new_bl.segments[0].epochs[0].annotations = {k:v for k,v in new_bl.segments[0].epochs[0].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}

    new_bl.segments[0].spiketrains[0].annotations = {k:v for k,v in new_bl.segments[0].spiketrains[0].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}
    new_bl.segments[0].spiketrains[1].annotations = {k:v for k,v in new_bl.segments[0].spiketrains[1].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}

    new_bl.channel_indexes[0].annotations = {k:v for k,v in new_bl.channel_indexes[0].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}
    new_bl.channel_indexes[1].annotations = {k:v for k,v in new_bl.channel_indexes[1].annotations.items()
                                                        if (k!="nix_name")&(k!="neo_name")}
    print("Saving block")
    save_block(new_bl, "NEW", filepath, True) # seems to save it as two blocks
    print("Finished spike analysis")
    gc.collect()

def get_spike_features(st, st_df, clusters_dict):
    """measure width and peak amplitude"""
    feature_df = st_df.apply(get_features).transpose()
    feature_df.columns = ["Prominence", "Peak_height", "Width"]
    feature_df.Width = (feature_df.Width / st.sampling_rate)*1e3
    feature_df["Clusters"]= feature_df.index.map(clusters_dict)
    feature_df["Time (s)"] = np.array(st.times)

    return feature_df

def get_features(series):
    """Calculates FWHM, prominence and amplitude"""
    series = series - (2*series)
    #peaks = find_peaks(series, height =10, prominence=10, distance = 25, width=1, rel_height=0.5)
    #prominence = peaks[1]["prominences"]
    #peak_height = peaks[1]["peak_heights"]
    #width = peaks[1]["widths"]

    #if  prominence == 0 :
    peak = series[int(series.shape[0]/2)]
    peak_idx = int(series.shape[0]/2)
    width, wh, li, ri = peak_widths(series, [peak_idx], rel_height = 0.5)
    prominence, lb, rb = peak_prominences(series, [peak_idx])

    return pd.Series((prominence[0], peak, width[0]))



def unit_plot(bl, chx_index, ax, isi_axes, ax_num, palette):


    print("Calculating ISI")
    unit_no = len(bl.channel_indexes[chx_index].units)
    print(unit_no)
    labels = []
    sts=[]
    isis = []
    for n in range(unit_no):

        label = bl.channel_indexes[chx_index].units[n].name
        #print(label)
        if int(label) != 0:
            try:
                st = bl.channel_indexes[chx_index].units[n].spiketrains[0]
                st.name = label
                isi = el.isi(st)
                sts.append(st)
                isis.append(isi)
                labels.append(label)
            except:
                print("no spiketrain")

    #colors = sns.color_palette("bright", n_colors = 20)
    colors=palette
    #print(len(sts))
    #labels = np.arange(len(sts)).tolist()

    #labels = [int(label) for label in labels]
    #lineoffsets = [int(re.search(r'\d+', label).group()) for label in labels]

    lineoffsets = [int(label) for label in labels]
    colors = np.array(colors)[np.array(lineoffsets)].tolist()
    #print(lineoffsets)
    ax.eventplot(sts, colors = colors, linewidth = 1,
                 lineoffsets = lineoffsets)

    ax.legend(labels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Unit")
    ax.set_ylim(bottom = 0.5, top=unit_no)
    if chx_index == 1:
        ax.set_title("Spinal Nerve")
        row_idx = 0
    else:
        ax.set_title("Dorsal Root")
        row_idx = 1
    #try:

     #isis is shorter in second nerve
    for idx, isi in enumerate(isis):
        #print(idx+ax_num+1)
        #isi_ax = isi_fig.add_subplot(2,len(isis), idx+ax_num+1)
        isi_ax = isi_axes[idx]
        #print(isi_ax)
        isi_ax.hist(isi.rescale(pq.ms), bins = 200, range = (0, 3000),
                    color = colors[idx], linewidth = 1)
        #isi_ax.set_ylabel("Density")
        #isi_ax.set_xlabel("Time (ms)")
    #except:
    #    print("couldnt do ISI")

    isis = [np.array(isi) for isi in isis]
    isis_df = pd.DataFrame(isis).transpose()

    return sts, isis_df

def get_stats(st1, st2):
    sttc = spike_time_tiling_coefficient(st1, st2)
    vrd = van_rossum_dist([st1, st2])[0]

    st1_mean_rate = float(el.mean_firing_rate(st1))
    st2_mean_rate = float(el.mean_firing_rate(st2))

    st1_cv = el.cv(el.isi(st1))
    st2_cv = el.cv(el.isi(st2))
    try:
        st1_cv2 = el.cv2(st1.times)
        st2_cv2 = el.cv2(st2.times)
    except:
        print("Couldnt calc cv2")
        st1_cv2 = None
        st2_cv2 = None

    #try:
    st1_inst_rate = el.instantaneous_rate(st1, 0.5*pq.s) #change sr from np.float(1/st1.sampling_rate)
    st2_inst_rate = el.instantaneous_rate(st2, 0.5*pq.s)
    st1_inst_rate.t_start = st1.t_start
    st2_inst_rate.t_start = st2.t_start
    #st1_inst_rate.times = st1_inst_rate.times.rescale(pq.s)-st1_inst_rate.t_start
    #st2_inst_rate.times = st2_inst_rate.times.rescale(pq.s)-st2_inst_rate.t_start
    #st1_inst_rate.t_stop = st1.t_stop
    #st2_inst_rate.t_stop = st2.t_stop

    #except:
    #print("couldnt calc inst rate")
    #st1_inst_rate, st2_inst_rate = (None, None)
    inst_rates = [st1_inst_rate, st2_inst_rate]

    data  = [[sttc, sttc], vrd, [st1_mean_rate, st2_mean_rate],
             [st1_cv, st2_cv], [st1_cv2, st2_cv2]]
    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ["STTC", "VRD", "MeanFiringRate", "CV1", "CV2"]
    return df, inst_rates

def get_statistics(st1, st2, t_slices):
    """Global stats of spinal nerve with dorsal root-return global stats df"""
    dfs = {}
    rates = {}
    labels = t_slices.keys()
    for n in range(len(labels)):
        label = list(labels)[n]
        start, end = t_slices[label]
        new_st1 = st1.time_slice(start, end)
        new_st2 = st2.time_slice(start, end)
        try:
            df, rate = get_stats(new_st1, new_st2)
            dfs[label] = df
            rates[label] = rate
        except:
            print("Couldnt get condition rates")

    full_df, full_rates = get_stats(st1, st2)
    dfs["Full"] = full_df
    rates["Full"] = full_rates


    #con_st = t_slices["control"][0]
    #con_end = t_slices["control"][1]
    #drug_st = t_slices["drug"][0]
    #drug_end = t_slices["drug"][1]
    #wash_st = t_slices["wash"][0]
    #wash_end = t_slices["wash"][1]

    #st1_control = st1.time_slice(con_st, con_end)
    #st2_control = st2.time_slice(con_st, con_end)

    #st1_test = st1.time_slice(drug_st, drug_end)
    #st2_test = st2.time_slice(drug_st, drug_end)

    #st1_recover = st1.time_slice(wash_st, wash_end)
    #st2_recover = st2.time_slice(wash_st, wash_end)

    #control_df, control_rates = get_stats(st1_control, st2_control)
    #test_df, test_rates = get_stats(st1_test, st2_test)
    #recover_df, recover_rates = get_stats(st1_recover, st2_recover)

    #dfs = [full_df, control_df, test_df, recover_df]
    #epochs = ["Full", "Control", "Test", "Recover"]
    epochs = list(dfs.keys())
    for n in range(len(epochs)):
        epoch = epochs[n]
        df = dfs[epoch]
        df["Epoch"] = [epochs[n]]*df.shape[0]
        df["Spiketrain"] = ["Spinal", "Dorsal"]



    final_df = pd.concat(dfs)
    #rates = [full_rates, control_rates, test_rates, recover_rates] # this is [st1, st2] inside []
    return final_df, rates


def unit_statistics(sts1, sts2, t_slices):
    """Unit stats of spinal nerve and dorsal root"""
    units = []
    units_rates = {}
    for n in range(len(sts1)):
        try:
            label = sts1[n].name
            unit_df, unit_rates = get_statistics(sts1[n], sts2[n], t_slices)
            unit_df["Unit"] = [label]*unit_df.shape[0]
            units.append(unit_df)
            units_rates[label] = unit_rates
        except:
            unit_df, unit_rates = get_statistics(sts1[n], sts1[n], t_slices)
            unit_df = unit_df[unit_df.Spiketrain=="Spinal"]
            unit_df["Unit"] = [label]*unit_df.shape[0]
            units.append(unit_df)
            units_rates[label] = unit_rates
            print ("unit had no spikes")

    final = pd.concat(units)
    return final, units_rates

def rate_plot(st1, st2, rates, filename):
    #rate1 = rates[0][0]
    #rate2 = rates[0][1]
    rate1 = rates["Full"][0]
    rate2 = rates["Full"][1]
    fig1=plt.figure(figsize = (15, 7))
    ax1=fig1.add_subplot(211)
    ax2=fig1.add_subplot(212)
    ax1.plot(rate1.times.rescale(pq.s), rate1.magnitude, color="k")
    ax1.plot(rate2.times.rescale(pq.s), rate2.magnitude, color="r")
    ax1.set_ylabel("Frequency")
    ax1.set_xlabel("Time (s)")
    ax2.eventplot([st2, st1], colors = ["r", "k"])
    ax2.set_ylabel("Spiketrain")
    ax2.set_xlabel("Time (s)")
    fig1.tight_layout()
    fig1.savefig(filename+"//rate_plot.svg")
    fig1.savefig(filename+"//rate_plot.png")
    plt.close(fig1)

def create_st2(st1, alpha_vel, beta_vel, delta_vel, c_vel, delete, jitter, new_method, new_new_method, distance, poisson):
    """This function will take the SN spiketrain and create an artifical spiketrain.
    Each spike time in the new train will the be randomly shifted within the range of
    the conduction velocities (and randomly low pass filtered to compare between pure time
    and also shape). This way an algorithm can be evaluated to produce the right match
    Conduction velocites over a 2mm distance :
        alpha = 100 m/s = 2 0us
        beta = 50 m/s = 30 us
        delta = 20 m/s = 100 us
        c = 1 m/s = 2 ms
    A range between 20 us and 2 ms
    Maybe need to define clusters specific to those speeds every 4th is this cluster
    starting from 1 then 2 then 3 then 4.

    Usage create_test_st(st1, 100, 50, 20, 1)
    """
    #alpha_vel, beta_vel, delta_vel, c_vel = [100, 50, 20, 2]
    alpha_lat, beta_lat, delta_lat, c_lat = (distance/np.array([alpha_vel, beta_vel, delta_vel, c_vel]))*pq.s
    global min_idx, df, alpha1, labels, min_latency, min_df, new_st, matches, alpha2
    new_st = st1

    df = pd.DataFrame({"st1": pd.Series(new_st.times), "st2": pd.Series(new_st.times)}, dtype = "float64")

    alpha1 = df.iloc[::4, 0]
    beta1 = df.iloc[1::4, 0]
    delta1 = df.iloc[2::4, 0]
    c1 = df.iloc[3::4, 0]

    labels = {"alpha":alpha1, "beta":beta1, "delta":delta1, "c":c1}
    series_label = pd.Series()
    for k, v in labels.items():
        label_series = pd.Series([k]*v.shape[0], index=v.index)
        series_label = pd.concat([series_label, label_series])

    df["st1_label"] = series_label
    """Jitter and reorder entire df after initial label has been assigned"""
    if jitter == True:
        df.st1_label = np.random.permutation(df.st1_label)

    alpha1 = df.loc[df.st1_label == "alpha", "st2"]
    beta1 = df.loc[df.st1_label == "beta", "st2"]
    delta1 = df.loc[df.st1_label == "delta", "st2"]
    c1 = df.loc[df.st1_label == "c", "st2"]

    alpha2 = df.loc[df.st1_label == "alpha", "st2"]
    beta2 = df.loc[df.st1_label == "beta", "st2"]
    delta2 = df.loc[df.st1_label == "delta", "st2"]
    c2 = df.loc[df.st1_label == "c", "st2"]

    df.loc[df.st1_label == "alpha", "st2"], slow_a  = shift_time(alpha2, alpha_lat)
    df.loc[df.st1_label == "beta", "st2"], slow_b = shift_time(beta2, beta_lat)
    df.loc[df.st1_label == "delta", "st2"], slow_d = shift_time(delta2, delta_lat)
    df.loc[df.st1_label == "c", "st2"], slow_c = shift_time(c2, c_lat)

    df["st1_isi"] = pd.Series(el.isi(df.iloc[:, 0]))
    df["st2_isi"] = pd.Series(el.isi(df.iloc[:, 1]))

    try:
        new_st = neo.SpikeTrain(df.st2*pq.s, t_stop=new_st.t_stop, units = new_st.units, sampling_rate = new_st.sampling_rate)
    except:
        print("st2 > t_stop")
        #df.st2 = df.st2[df.st2<float(new_st.t_stop)]
        #new_st = neo.SpikeTrain(df.st2*pq.s, t_stop=new_st.t_stop, units = new_st.units, sampling_rate = new_st.sampling_rate)
        new_st = neo.SpikeTrain(df.st2*pq.s, t_stop=(df.st2.iloc[-1]+0.5)*pq.s, units = new_st.units, sampling_rate = new_st.sampling_rate)


    if delete>0:
        random_indexes = np.random.randint(0, df.st2.shape[0], int(len(new_st)*delete))
        new_st=np.delete(new_st, random_indexes)
        df.iloc[random_indexes, 1] = np.nan

    latency, fig = spike_latency(st1, new_st, True, slow_c)
    latency = latency[latency<slow_c]
    min_idx = latency.idxmin(axis=0) #this seems to find the min row index (st1) for each column (st2)
    min_latency = latency.min(axis=0)

    print("Min_idx length = {}".format(min_idx.shape[0]))
    na_min_idx = min_idx.dropna()
    print("New min_idx length = {}".format(na_min_idx.shape[0]))

    latency_array = np.array(latency)
    latency_array[np.array(na_min_idx, dtype="int64"), np.array(na_min_idx.index, dtype="int64")] = np.nan


    second_latency = pd.DataFrame(latency_array)
    second_min_idx = second_latency.idxmin(axis=0)
    second_min_lat = second_latency.min(axis=0)

    min_idx.index = df.st2.dropna().index
    min_latency.index = df.st2.dropna().index
    second_min_idx.index = df.st2.dropna().index
    second_min_lat.index = df.st2.dropna().index

    df["min_idx"] = min_idx
    df["min_latency"] = min_latency
    df["2_min_idx"] = second_min_idx
    df["2_min_lat"] = second_min_lat

    df["st2_label"] = df.min_idx.map(df.st1_label)
    na_df = df[df.st2.notnull()]
    matches = na_df.st1_label == na_df.st2_label
    df["matches"] = matches
    failures = df[df.matches==False].st1_label.value_counts()
    failures = (failures/df.st1_label.value_counts()) * 100.

    if new_method==True:
        """Pick out ISI that are faster than the slowest conduction and that have latency greater than ISI"""
        isi_subset = df[(df.st1_isi<slow_c)&(df.st2.notnull())&(df.min_latency>df.st1_isi)]
        new_idx = min_idx.copy()
        new_idx.loc[isi_subset.index] = isi_subset["2_min_idx"]
        df["new_idx"] = new_idx
        df["new_label"] = df.new_idx.map(df.st1_label)
        na_df = df[df.st2.notnull()]
        new_matches = na_df.st1_label == na_df.new_label
        df["new_matches"] = new_matches
        matches = new_matches
        failures = df[df.new_matches==False].st1_label.value_counts()
        failures = (failures/df.st1_label.value_counts()) * 100.

        if new_new_method==True:
            new_subset=df[(df.st1_isi<slow_c)&(df.st2.notnull())&
                          (df.min_latency>df.st2_isi)&
                          (df.st2_isi<df.st1_isi)&
                          (df["2_min_lat"]>df.st1_isi)]
            new_new_idx = new_idx.copy()
            new_new_idx.loc[new_subset.index] = new_subset["2_min_idx"]
            df["new_new_idx"] = new_new_idx
            df["new_new_label"] = df.new_new_idx.map(df.st1_label)
            na_df = df[df.st2.notnull()]
            new_new_matches = na_df.st1_label == na_df.new_new_label
            df["new_new_matches"] = new_new_matches
            matches = new_new_matches
            failures = df[df.new_new_matches==False].st1_label.value_counts()
            failures = (failures/df.st1_label.value_counts()) * 100.
    #check if latency 1 is quicker than latency 2 +isi
    # if lat 1 is > isi then assign 1st to its second min index/alternatively make 1=1 and 2=2
    # if lat1 < ISI leave
    # if both lat are slower than ISI then 1=1, 2=2

    try:
        true, false = matches.value_counts()

    except:
        true = matches.value_counts()

    success_rate = float((true/matches.shape[0]) * 100)
    print("success rate is {}".format(success_rate))

    """Which fibre types fail?"""

    try:
        st1_wv = st1.waveforms
        st2_wv = st1.waveforms
        st1_df = pd.DataFrame(np.reshape(st1_wv, np.shape(st1_wv)[:2])).transpose()
        st2_df = pd.DataFrame(np.reshape(st1_wv, np.shape(st1_wv)[:2])).transpose()
        #st1_feat = get_spike_features(st1, st1_df, cluster_dict)
        #st2_feat = get_spike_features(new_st, st2_df, dorsal_dict)
    except:
        st1_df=pd.DataFrame()
        st2_df=pd.DataFrame()
        print("no waveforms")
    if poisson == False:
        bl2=0
        bl2 = neo.Block()
        ch2 = neo.ChannelIndex(1)
        ch1 = neo.ChannelIndex(0)
        bl2.channel_indexes.append(ch1)
        bl2.channel_indexes.append(ch2)

        new_bl2=0
        new_bl2 = create_units(st1, st1_df, df.st1_label.to_dict(), bl2, 1) # adds spinal nerve units to block
        new_bl2 = create_units(new_st, st2_df, df.st2_label.to_dict(), new_bl2, 0) # adds dorsal root units to block
        new_bl2.channel_indexes[0].units = new_bl2.channel_indexes[0].units[0:4]
        new_bl2.channel_indexes[1].units = new_bl2.channel_indexes[0].units[0:4]

        raster_fig = plt.figure(figsize = (15, 7))                                              ##UNDO these unhighlights during real tests
        raster_ax1 = raster_fig.add_subplot(211)
        raster_ax2 = raster_fig.add_subplot(212)
        isi_fig = plt.figure()
        SN_sts, SN_isis = unit_plot(new_bl2, 1, raster_ax1, isi_fig, 0)
        DR_sts, DR_isis = unit_plot(new_bl2, 0, raster_ax2, isi_fig, 4)
        isi_fig.tight_layout()


    return success_rate, failures
    # add a score for percentage that were correct, jitter fibre types, randomly delete some too, need a definitive map of identity
    # make number of spikes deleted a percentage of firing rate, add case by case analysis to try to improve score

def create_simulation(st, freq, alpha_vel, beta_vel, delta_vel, c_vel, new_method, new_new_method, distance, poisson):
    """This functions creates a poisson spike train and uses create st2 to test
    for different firing frequencies and conduction velocities and different deletions"""
    global success_df
    if poisson == True:
        """Create fake spiketrain"""

        st = homogeneous_poisson_process(rate = freq*pq.Hz, t_start = 0*pq.s, t_stop = 10*pq.s)

    else:
        st = st
    success_rates = []

    failure_rates = []
    """Vary percentage of deletion. Ensure jitter of spike labels"""
    #success_rate = create_st2(st, alpha_vel, beta_vel, delta_vel, c_vel, 0, False)
    percentages = np.arange(0, 0.5, 0.1)
    for percentage in percentages:
        success_rate, failures = create_st2(st, alpha_vel, beta_vel, delta_vel, c_vel, percentage,
                                            True, new_method, new_new_method, distance, poisson)
        success_rates.append(success_rate)

        failure_rates.append(failures)

        plt.close("all")

    success_df = pd.DataFrame(success_rates, index = percentages)

    velocities = np.tile(np.array([freq, alpha_vel, beta_vel ,delta_vel, c_vel]), 5).reshape(5, 5)
    velocity_df = pd.DataFrame(velocities, index = success_df.index)

    success_df = pd.concat([success_df, velocity_df], axis=1)
    success_df.reset_index(inplace=True)
    success_df.columns = ["deleted", "accuracy", "freq", "alpha", "beta", "delta", "c"]

    failure_df = pd.DataFrame(failure_rates, index=percentages)


    failure2_df = pd.concat([failure_df, velocity_df], axis=1)
    failure2_df.reset_index(inplace=True)
    failure2_df.columns = ["deleted", *failure_df.columns ,"freq", "alpha_vel", "beta_vel", "delta_vel", "c_vel"]


    return success_df, failure2_df

def io(folder):
    file = pd.read_excel(os.path.join(folder,"block_analysis.xlsx"), sheet_name=None)
    return file

def facet_plot(data, var, path):
    #palette = sns.color_palette("bright", n_colors= 20 )
    if (var == "STTC") | (var == "VRD"):
        data = data[data.Spiketrain=="Dorsal"]
        grid = sns.FacetGrid(data=data, col="Unit",
                             col_wrap=5, height=1.5, palette=palette)
    else:
        grid = sns.FacetGrid(data=data, col="Unit", hue="Spiketrain",
                             col_wrap=5, height=1.5, palette=palette)
    grid.map(plt.plot, "Epoch", var, marker="o")
    grid.add_legend()
    #grid.fig.tight_layout()
    grid.savefig(path+"//"+var+".svg")


def scrape_block_analysis():
    directory = filedialog.askdirectory()
    subdirs = os.listdir(directory)
    subdirs = [os.path.join(directory, subdir) for subdir in subdirs]
    core_df = pd.DataFrame()
    core_global = pd.DataFrame()
    for n, subdir in enumerate(subdirs):
        file = io(subdir)

        """Spinal nerve features"""
        print("Getting spinal nerve features")
        SN_feats = file["SN_features"].copy()
        SN_feats.set_index(SN_feats.columns[0], inplace=True)
        SN_feats_mean = SN_feats.groupby("Clusters")["Prominence", "Peak_height", "Width"].mean()
        SN_feats_mean["Sample"] = [n+1]*SN_feats_mean.shape[0]
        concat = SN_feats_mean

        """Get mean unit latency"""
        print("Getting mea unit latency")
        latency = file["latency_desc"].copy()
        spike_type_map = latency["Fibre"].to_dict()

        latency_mean = latency["mean"]
        concat["latency"] = latency_mean

        """Extract % filtered by condition"""
        print("Extracting % filtered by condition")
        filtered = file["match_counts"].copy()
        filtered.set_index(filtered.iloc[:,0].fillna(method="ffill"), inplace=True)
        columns = filtered.index.unique()
        filt_pivot = filtered.pivot(columns = "st1_label", values="Percentage_Change").transpose()
        filt_pivot.columns = columns
        concat = pd.concat([concat, filt_pivot], axis=1)

        """Average unit stats"""
        print("Averaging unit stats")
        unit_stats = file["Unit_spike_stats"].copy()
        unit_stats.set_index(unit_stats.columns[0], inplace=True)
        unit_stats=unit_stats[(unit_stats.Epoch!="Full")&(unit_stats.Unit!=0)]
        facet_plot(unit_stats, "MeanFiringRate", subdir)
        facet_plot(unit_stats, "STTC", subdir)
        facet_plot(unit_stats, "VRD", subdir)
        facet_plot(unit_stats, "CV2", subdir)

        #"""Calculate freq change between control and capsaicin"""

        #unit_freq_pivot =unit_stats[(unit_stats.Spiketrain=="Spinal")].pivot(columns = "Epoch", index="Unit", values = "MeanFiringRate")
        #unit_freq_pivot.fillna(0, inplace=True)
        #percent_change = unit_freq_pivot.cap-unit_freq_pivot.control

        #"""Define threshold-using 2Hz"""
        #thresh_change = percent_change>2
        #cap_sens_dict = thresh_change.to_dict()
        #concat["firingchange"] = percent_change
        #concat["Sens"] = concat.index.map(cap_sens_dict)
        concat["Fibre"] = concat.index.map(spike_type_map)
        concat = concat.iloc[1:, :]
        core_df = pd.concat([core_df, concat])

        global_stats = file["Global_spike_stats"]
        global_stats = global_stats[(global_stats.Epoch!="Full")]
        global_stats["Sample"] = [n+1]*global_stats.shape[0]
        core_global = pd.concat([core_global, global_stats])

    #sns.catplot(kind="point", data=core_global, x="Epoch", y="MeanFiringRate",
                   # hue ="Spiketrain", ci=68, linewidth = 1.5, palette=palette)
    print("new")
    core_df.reset_index(inplace=True)
    core_df.to_csv(directory+"//block_summary1.csv")
    core_global.to_csv(directory+"//block_summary2.csv")
    return core_df, core_global

#final_df, final_failure_df, failures, alpha, beta, delta, c=multiple_simulation(True, True, 2e-3, False)
"""Junk"""################################################################################



def make_null(clean_phrenic, clean_DR, clean_SN, phrenic, DR, SN):
    """Sets 5 seconds between control and drug and between drug and wash to null"""
    phrenic_dt  = 1/np.array(phrenic.sampling_rate)
    DR_dt  = 1/np.array(DR.sampling_rate)
    SN_dt  = 1/np.array(SN.sampling_rate)

    dts = [phrenic_dt, DR_dt, SN_dt]
    sigs = [clean_phrenic, clean_DR, clean_SN]
    for dt in dts:
        index = dts.index(dt)
        first_5_start, first_5_end = (int(10/dt), int(15/dt))
        print("start and end is {} and {}".format(first_5_start, first_5_end))
        second_5_start, second_5_end  = (int(25/dt), int(30/dt))
        print("start and end is {} and {}".format(second_5_start, second_5_end))
        mask = np.append(np.arange(first_5_start, first_5_end, 1), np.arange(second_5_start, second_5_end, 1))
        sig = sigs[index]
        sig[mask[:-1]] = np.nan
        sigs[index] = sig

    return sigs

def analyse_ecg():
    """Analyses ecg frequency before during and after stimulus given at time 10 seconds
    Could extend to analysis of PQRST if necessary-useful if given KCNQ blockers?"""
    pass

def remove_artifacts(phrenic_subset, DR_subset, SN_subset, phrenic_mask, DR_mask, SN_mask):
    signals = [phrenic_subset, DR_subset, SN_subset]
    masks = [phrenic_mask, DR_mask, SN_mask]
    cleaned=[]
    noisy = []
    for sig in signals:
       try:
            index = signals.index(sig)
            mask = masks[index]
            print(mask)
            sig = sig.magnitude.flatten()
            clean = sig.copy()
            noisy_subset = sig.copy()
            mean = np.mean(sig)
            std = np.std(sig)
            new_vals = np.random.normal(mean, std, mask.shape)
            noisy_subset[mask] = new_vals

            clean[mask]=0

       except:
            print("not cleaned")
            noisy_subset = sig
            clean = sig
       noisy.append(noisy_subset)
       cleaned.append(clean)


    return (cleaned, noisy)

def eplot(st1, st2):
    fig1=plt.figure()
    ax1=fig1.add_subplot(111)
    ax1.eventplot([st2, st1], colors = ["r", "k"])
    fig1.savefig("./example raster.svg")


def spike_analysis2(bl):
    """Create sample spike trains, work with these to create a solution to nearest neighbour,
    also look at waveforms, produce output of ISI, psd"""
    from elephant.spike_train_generation import homogeneous_poisson_process
    from elephant.statistics import isi, cv, mean_firing_rate, instantaneous_rate, complexity_pdf
    from elephant.spike_train_correlation import spike_time_tiling_coefficient
    from elephant.spike_train_dissimilarity import van_rossum_dist, victor_purpura_dist
    #st = homogeneous_poisson_process(rate=10.0*pq.Hz, t_start = 0.0*pq.s, t_stop=100*pq.s)
    #rand = np.random.choice(st, 300)
    #st2 = st[~np.isin(st, rand)]
    #st2 = st2 +(0.005*pq.s)
    #st2.t_stop = st.t_stop
    #st2.t_start = st.t_start

    st1 = bl.segments[0].spiketrains[1]
    st2 = bl.segments[0].spiketrains[0]


    eplot(st1, st2)
    spike_time_tiling_coefficient(st1, st2)
    van_rossum_dist([st1, st2])
    plt.figure()
    st1_isi = isi(st1)
    st2_isi = isi(st2)
    plt.scatter(st1_isi[:-1], st1_isi[1:], c="k")
    plt.scatter(st2_isi[:-1], st2_isi[1:], c="r")
    plt.figure()
    plt.hist(np.array(st1_isi.rescale(pq.ms)), bins = 200, range = (0, 5000))
    plt.figure()
    plt.hist(np.array(st2_isi.rescale(pq.ms)), bins = 200, range = (0, 5000))
    #plot isis against each other
    """get widths"""
    st1_cv = cv(st1_isi)
    st2_cv = cv(st2_isi)


def export_matlab(clean_DR, DR_subset, clean_SN, SN_subset):

    subset_size = len(clean_DR)/3
    windows = np.arange(subset_size, len(clean_DR)+1, subset_size, dtype="int64")

    for win in tqdm(windows):
        bl=neo.Block()
        seg=neo.Segment()
        index = np.where(windows==win)
        clean_DR_2 = clean_DR[win-int(subset_size) : win]
        clean_SN_2 = clean_SN[win-int(subset_size) : win]

        DR_sig = neo.AnalogSignal(clean_DR_2*(pq.uV), t_stop = DR_subset.t_start, sampling_rate=DR_subset.sampling_rate)
        SN_sig = neo.AnalogSignal(clean_SN_2*(pq.uV), t_stop = SN_subset.t_start, sampling_rate=SN_subset.sampling_rate)
        seg.analogsignals.append(DR_sig)
        seg.analogsignals.append(SN_sig)
        bl.segments.append(seg)
        #from neo.io.neomatlabio import NeoMatlabIO
        w = NeoMatlabIO(filename='./block'+str((index[0][0]*10)+10)+'.mat')
        w.write_block(bl)

    full_block = neo.Block()
    full_seg = neo.Segment()
    full_DR_sig = neo.AnalogSignal(clean_DR*(pq.uV), t_stop = DR_subset.t_start, sampling_rate=DR_subset.sampling_rate)
    full_SN_sig = neo.AnalogSignal(clean_SN*(pq.uV), t_stop = SN_subset.t_start, sampling_rate=SN_subset.sampling_rate)
    full_seg.analogsignals.append(full_DR_sig)
    full_seg.analogsignals.append(full_SN_sig)
    full_block.segments.append(full_seg)
    w2 = NeoMatlabIO(filename='./fullblock'+'.mat')
    w2.write_block(full_block)



def neo_plot_check(bl):
    """This plots two subplots, one for DR and one for SN. It also plots
    the trains associated with each to check that the timing is correct"""
    for seg in bl.segments:
        fig= plt.figure()
        sr = seg.analogsignals[0].sampling_rate

        for n in range(0, 2):
            asig = seg.analogsignals[n]
            times = asig.times.rescale('s').magnitude
            asig = asig.magnitude
            ax1 = fig.add_subplot(2, 1, n+1)
            ax1.plot(times, asig)

            st = seg.spiketrains[n]
            train = st.rescale('s').magnitude
            spike_indices = (train/np.float((1/sr))).astype("int64")

            ax1.scatter(train, asig[spike_indices], c="r")
"""Spectral decmposition with ICA to remove ecg artifacts"""


def singular_spectral(array, x=130, groups =12):
    """Need to bin signal as crashes computer"""
    from pyts.decomposition import SingularSpectrumAnalysis
    array = np.array([array], dtype="float64")
    ssa = SingularSpectrumAnalysis(window_size=x, groups = groups)
    X_ssa = ssa.fit_transform(array)

    plt.figure()
    index = 0
    for n in X_ssa:
        plt.subplot(groups, 1, index+1)
        plt.plot(np.transpose(n))
        index= index+1

    from sklearn.decomposition import FastICA
    transformer = FastICA(n_components=groups)
    X_transformed = transformer.fit_transform(np.transpose(X_ssa))
    plt.figure()
    for n in range(X_transformed.shape[1]):
        plt.subplot(X_transformed.shape[1], 1, n+1)
        plt.plot(X_transformed[:, n])

    return X_transformed

def latest_ss(test):
    decomp =  elephant.signal_processing.wavelet_transform(test, [10, 50,
                                                        100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000])
    decomp2 = decomp.reshape(50000, 13)
    plt.figure()
    index = 0
    for n in range(decomp2.shape[1]):
        plt.subplot(decomp2.shape[1], 1, n+1)
        plt.plot(decomp2[:, n])
        index= index+1
    from sklearn.decomposition import FastICA
    transformer = FastICA(n_components=13)
    X_transformed = transformer.fit_transform(np.abs(decomp2))
    plt.figure()
    for n in range(X_transformed.shape[1]):
        plt.subplot(X_transformed.shape[1], 1, n+1)
        plt.plot(X_transformed[:, n])
    #y = elephant.signal_processing.butter(DR_subset, highpass_freq = 300)
    pass

def ss(array):

    L = 130
    N = len(array)
    K=N-L+1
    test=[]
    for n in range(K):
        subset = array[n:n+L]
        test.append(subset)
    test2 = np.transpose(np.matrix(test))
    #Xt=np.cov(test2)
    #print(np.shape(test2))
    #print(np.shape(Xt))
    U, s, V = scipy.linalg.svd(test2[:, 0:10000])



def test_plot(phrenic, DR, SN):
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax1.plot(phrenic)
    ax2.plot(DR)
    ax3.plot(SN)



def spike_similarity_score(waveform_corr, spike_latencies, st1_df, st2_df, clusters_dict, dt):
    """add inverted waveform corr and spike latency together"""
    subset3  = 1-waveform_corr
    ss_score = subset3+spike_latencies
    sscore_df = pd.concat([ss_score.idxmin(), ss_score.min()], axis=1)
    sscore_df.columns = ["Index", "Min"]
    sscore_df = sscore_df.sort_values(by=["Min"])
    not_dup = sscore_df.drop_duplicates("Index", keep="first")
    not_dup=not_dup.dropna()
    not_dup["Clusters"]=not_dup.Index.map(clusters_dict)
    dorsal_dict = not_dup.Clusters.to_dict()

    st2_melt = st2_df.melt()
    st2_melt["Clusters"] = st2_melt.variable.map(dorsal_dict)
    st2_melt["Index"] = np.tile(np.arange(0, 24), st2_df.shape[1])
    st2_melt.Index =  st2_melt.Index*(dt)    # change to global dt

    st1_melt = st1_df.melt()
    st1_melt["Clusters"] = st1_melt.variable.map(clusters_dict)
    st1_melt["Index"] = np.tile(np.arange(0, 24), st1_df.shape[1])
    st1_melt.Index =  st1_melt.Index*(dt)

    return st1_melt, st2_melt, dorsal_dict, sscore_df
