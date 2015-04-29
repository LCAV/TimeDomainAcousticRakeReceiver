
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz
from scipy.io import wavfile
from scipy.signal import resample,fftconvolve

import pyroomacoustics as pra

# simulation parameters
bf_names = ('Rake-MaxSINR', 'Rake-Perceptual', 'Rake-MVDR')
bf_designs = (pra.Beamformer.rakeMaxSINRFilters, 
              pra.Beamformer.rakePerceptualFilters, 
              pra.Beamformer.rakeMVDRFilters)
max_source = 10
loops = 10000
loops = 10
SINR = np.zeros((max_source, len(bf_designs), loops))

# Beam pattern figure properties
freq=[800, 1600]
figsize=(1.88,2.24)
xlim=[-4,8]
ylim=[-4.9,9.4]

# Some simulation parameters
Fs = 8000
t0 = 1./(Fs*np.pi*1e-2)  # starting time function of sinc decay in RIR response
absorption = 0.90
max_order_sim = 3
sigma2_n = 1e-15
SNRdB = 10.

# Room 1 : Shoe box
room_dim = np.array([4, 6])

# the good source is fixed for all 
good_source = [1, 4.5]       # good source
normal_interferer = [2.8, 4.3]   # interferer
hard_interferer = [1.5, 3]   # interferer in direct path
#normal_interferer = hard_interferer

# microphone array design parameters
mic1 = [2, 1.5]         # position
M = 8                    # number of microphones
d = 0.08                # distance between microphones
phi = 0.                # angle from horizontal
max_order_design = 1    # maximum image generation used in design
shape = 'Linear'        # array shape
Lg_t = 0.03             # Filter size in seconds
Lg = np.ceil(Lg_t*Fs)   # Filter size in samples
delay = 0.02            # delay of beamformer

# define the FFT length
N = 1024

# create a microphone array
if shape is 'Circular':
    R = pra.circular2DArray(mic1, M, phi, d*M/(2*np.pi)) 
else:
    R = pra.linear2DArray(mic1, M, phi, d) 
center = R.mean(axis=1, keepdims=True)
mics = pra.Beamformer(R, Fs, N=N, Lg=Lg)

# create the room with sources and mics
room1 = pra.Room.shoeBox2D(
    [0,0],
    room_dim,
    Fs,
    t0 = t0,
    max_order=max_order_sim,
    absorption=absorption,
    sigma2_awgn=sigma2_n)

for n in np.arange(loops):

    # Add new source and interferer to the room
    source = np.random.random(2)*room_dim
    room1.addSource(source)
    room1.sources[0].setOrdering('nearest', ref_point=center)

    interferer = np.random.random(2)*room_dim
    room1.addSource(interferer)
    room1.sources[1].setOrdering('nearest', ref_point=center)

    # compute noise power for 20dB SNR
    dc = np.sqrt(((source[:,np.newaxis] - center)**2).sum())
    sigma2_n = 10.**(-SNRdB/10.)/(4.*np.pi*dc)**2

    for i in np.arange(max_source):

        # Select their nearest image sources (from the array center)
        good_sources = room1.sources[0][:i+1]
        bad_sources = room1.sources[1][:i+1]

        # compute SNR for all three beamformers
        mics = pra.Beamformer(R, Fs, N=N, Lg=Lg)
        for t in np.arange(len(bf_designs)):
            SINR[i,t,n] = bf_designs[t](mics, good_sources, bad_sources, sigma2_n*np.eye(mics.Lg*mics.M), delay=delay)

    # remove the source and interferer from the room
    room1.sources.pop()
    room1.sources.pop()

np.save('data/SINR_data.npy', SINR)
