import numpy as np
from scipy.io import wavfile
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import matplotlib.mlab as mlab
from matplotlib import pyplot as plt

from IPython.display import Audio

import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')
import readligo as rl


# function to whiten data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    freqs1 = np.linspace(0, 2048., Nt//2+1)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    norm = 1./np.sqrt(1./(dt*2))
    white_hf = hf / np.sqrt(interp_psd(freqs)) * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht


# function to keep the data within integer limits, and write to wavfile:
def write_wavfile(filename,fs,data):
    filename = '../audio/' + filename
    # open(filename, 'w').close()
    d = np.int16(data/np.max(np.abs(data)) * 32767 * 0.9)
    wavfile.write(filename,int(fs), d)


# function that shifts frequency of a band-passed signal
def reqshift(data,fshift=100,sample_rate=4096):
    """Frequency shift the signal by constant
    """
    x = np.fft.rfft(data)
    T = len(data)/float(sample_rate)
    df = 1.0/T
    nbins = int(fshift/df)
    # print T,df,nbins,x.real.shape
    y = np.roll(x.real,nbins) + 1j*np.roll(x.imag,nbins)
    y[0:nbins]=0.
    z = np.fft.irfft(y)
    return z


# funtion for creating a plot
def create_plot(x, y, pcolor, title='', plot_label='', xlabel=None, ylabel=None, xlim=None, ylim=None, \
 grid=True, legend_loc='upper left', loglog=False):
    if (loglog):
        plt.loglog(x, y, pcolor, label=plot_label)
    else:
        plt.plot(x, y, pcolor, label=plot_label)
    if ylim:
        plt.ylim(ylim)
    if xlim:
        plt.xlim(xlim)
    if grid:
        plt.grid('on')
    if ylabel:
        plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    plt.legend(loc=legend_loc)
    if title != '':
        plt.title(plot_label)


# test 1
def test_write_wavfile():
    _1, time, _2 = rl.loaddata('../data/H-H1_LOSC_4_V2-1126259446-32.hdf5', 'H1')
    fs = 4096

    f_template = h5py.File('../data/GW150914_4_template.hdf5', "r")
    template_p, _ = f_template["template"][...]
    indxt = np.where((time >= (time[0]+16-2)) & (time < (time[0]+16+2)))
    write_wavfile("test_write_wavfile.wav", int(fs), template_p[indxt])
    d = wavfile.read("../audio/test_write_wavfile.wav") # if there is no error, test is passed



# test 2
def test_whiten():
    strain_H1, time, chan_dict_H1 = rl.loaddata('../data/H-H1_LOSC_4_V2-1126259446-32.hdf5', 'H1')
    dt = time[1] - time[0]
    fs = 4096
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = 4*fs)
    # fband = [43.0,300.0]
    # psd_H1 = interp1d(freqs, Pxx_H1)
    # strain_H1_whiten = whiten(strain_H1, psd_H1, dt)
    # bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    # normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
    # strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten) / normalization


    Pxx = (1.e-22*(18./(0.1+freqs))**2)**2+0.7e-23**2+((freqs/2000.)*4.e-23)**2
    psd_smooth = interp1d(freqs, Pxx)
    f_template = h5py.File('../data/GW150914_4_template.hdf5', "r")
    template_p, template_c = f_template["template"][...]
    template_p_smooth = whiten(template_p, psd_smooth, dt)
    indxt = np.where((time >= (time[0]+16-2)) & (time < (time[0]+16+2)))
    write_wavfile("test_template_whiten.wav", int(fs), template_p_smooth[indxt])
    Audio("../audio/test_template_whiten.wav") # if there is no error, test is passed


# test 3
def test_reqshift():
    strain_H1, time, chan_dict_H1 = rl.loaddata('../data/H-H1_LOSC_4_V2-1126259446-32.hdf5', 'H1')
    dt = time[1] - time[0]
    fs = 4096
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = 4*fs)
    Pxx = (1.e-22*(18./(0.1+freqs))**2)**2+0.7e-23**2+((freqs/2000.)*4.e-23)**2
    psd_smooth = interp1d(freqs, Pxx)
    f_template = h5py.File('../data/GW150914_4_template.hdf5', "r")
    template_p, template_c = f_template["template"][...]
    template_p_smooth = whiten(template_p, psd_smooth, dt)

    fshift = 400.
    tevent = 1126259462.44
    deltat = 5
    template_p_shifted = reqshift(template_p_smooth,fshift=fshift,sample_rate=fs)
    indxt = np.where((time >= (time[0]+16-2)) & (time < (time[0]+16+2)))
    write_wavfile("test_template_shifted.wav",int(fs), template_p_shifted[indxt])
    Audio("../audio/test_template_shifted.wav") # if there is no error, test is passed


# test 4
def test_create_plot():
    plt.figure(figsize=(10,8))
    x = np.arange(1, 20)
    y = np.array([el * 2 for el in x])
    create_plot(x, y, 'r', title='test', plot_label='test plot', xlabel='X', ylabel='Y')
    plt.savefig('../figures/test_plot.png')
    # if there is no error, test is passed

