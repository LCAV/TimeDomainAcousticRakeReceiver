
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import resample

import pyroomacoustics as pra
import TDBeamformers as tdb

# Beam pattern figure properties
freq=[800, 1600]
figsize=(1.88,2.24)
xlim=[-4,8]
ylim=[-4.9,9.4]

# Some simulation parameters
Fs = 8000
t0 = 1./(Fs*np.pi*1e-2)  # starting time function of sinc decay in RIR response
absorption = 0.90
max_order_sim = 10
sigma2_n = 1e-7

# Room 1 : Shoe box
room_dim = [4, 6]

# the good source is fixed for all 
good_source = [1, 4.5]       # good source
normal_interferer = [2.8, 4.3]   # interferer
hard_interferer = [1.5, 3]   # interferer in direct path

# microphone array design parameters
mic1 = [2, 1.5]         # position
M = 8                    # number of microphones
d = 0.08                # distance between microphones
phi = 0.                # angle from horizontal
max_order_design = 1    # maximum image generation used in design
shape = 'Linear'        # array shape
Lg_t = 0.050            # Filter size in seconds
Lg = np.ceil(Lg_t*Fs)   # Filter size in samples
delay = 0.020           # beamformer delay

# create a microphone array
if shape is 'Circular':
    R = pra.circular2DArray(mic1, M, phi, d*M/(2*np.pi)) 
else:
    R = pra.linear2DArray(mic1, M, phi, d) 

# define FFT length for beamformer analysis
N = 1024

# create the room with sources and mics
room1 = pra.Room.shoeBox2D(
    [0,0],
    room_dim,
    Fs,
    t0 = t0,
    max_order=max_order_sim,
    absorption=absorption,
    sigma2_awgn=sigma2_n)

# add source and interferer
room1.addSource(good_source)
room1.addSource(normal_interferer)

''' 
SCENARIO 1
Only one source of interest and one interferer (easy)
RakeMVDR
'''
print 'Scenario1...'

# Compute the beamforming weights depending on room geometry
good_sources = room1.sources[0].getImages(max_order=max_order_design)
bad_sources = room1.sources[1].getImages(max_order=max_order_design)
mics = tdb.RakeMVDR_TD(R, Fs, N, Lg=Lg)
mics.computeWeights(good_sources, bad_sources, sigma2_n*np.eye(mics.Lg*mics.M), delay=delay)
mics.weightsFromFilters()

room1.addMicrophoneArray(mics)

# plot the room and beamformer
fig, ax = room1.plot(img_order=np.minimum(room1.max_order, 1), 
        freq=freq,
        mic_marker_size=2,
        figsize=figsize, no_axis=True,
        xlim=xlim, ylim=ylim,
        autoscale_on=False)
fig.savefig('figures/scenario_RakeMVDR.pdf',
            facecolor=fig.get_facecolor(), edgecolor='none')

'''
SCENARIO 2
One source or interest and one interefer (easy)
RakePerceptual
'''
print 'Scenario2...'

# Compute the beamforming weights depending on room geometry
mics = tdb.RakePerceptual_TD(R, Fs, N, Lg=Lg)
mics.computeWeights(good_sources, bad_sources, sigma2_n*np.eye(mics.Lg*mics.M), delay=delay)
mics.weightsFromFilters()

room1.addMicrophoneArray(mics)

# plot the room and beamformer
fig, ax = room1.plot(img_order=np.minimum(room1.max_order, 1), 
        freq=freq,
        mic_marker_size=2,
        figsize=figsize, no_axis=True,
        xlim=xlim, ylim=ylim,
        autoscale_on=False)
fig.savefig('figures/scenario_RakePerceptual.pdf',
            facecolor=fig.get_facecolor(), edgecolor='none')

'''
SCENARIO 3
One source or interest and one interefer (hard)
RakePerceptual
'''
print 'Scenario3...'

room1.sources.pop()
room1.addSource(hard_interferer)

bad_sources = room1.sources[1].getImages(max_order=max_order_design)

# Compute the beamforming weights depending on room geometry
mics = tdb.RakePerceptual_TD(R, Fs, N, Lg=Lg)
mics.computeWeights(good_sources, bad_sources, sigma2_n*np.eye(mics.Lg*mics.M), delay=delay)
mics.weightsFromFilters()

room1.addMicrophoneArray(mics)

# plot the room and beamformer
fig, ax = room1.plot(img_order=np.minimum(room1.max_order, 1), 
        freq=freq,
        mic_marker_size=2,
        figsize=figsize, no_axis=True,
        xlim=xlim, ylim=ylim,
        autoscale_on=False)
fig.savefig('figures/scenario_RakePerceptual_hard.pdf',
            facecolor=fig.get_facecolor(), edgecolor='none')

