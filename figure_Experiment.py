


import numpy as np
from scipy.io import wavfile
from scipy.signal import resample, decimate
from os import getpid

import matplotlib.pyplot as plt

import pyroomacoustics as pra

data_folder = 'Experiment_Data/'

Fs = 44100.
N = 1024
Lg = int(0.100*Fs) # 350 ms long filter
sigma2_n = 1e-6
c = 345.35589778869274

speech_sample1 = 'samples/fq_sample1.wav'
speech_sample2 = 'samples/fq_sample2.wav'

# read in microphones and sources locations
f = open(data_folder + 'RIRs.positions', 'r')
lines = f.readlines()
f.close()

# count # of mics and speakers
n_src = len([l for l in lines if l.split(' ')[0] == 's'])
n_mic = len([l for l in lines if l.split(' ')[0] == 'm'])

# reflection coefficients from the walls (hand-waving)
reflection = {'ground':0.7, 'south':0.7, 'west':0.7, 'north':0.7, 'east':0.7, 'ceilling':0.4}

# read in the speakers and mics positions
sources = np.zeros((3,n_src))
mics = np.zeros((3,n_mic))
room_dim = np.zeros(3)

for line in lines:
    l = line.split(' ')
    if l[0] == 'room':
        room_dim = np.array([float(l[1]), float(l[2]), float(l[3])])
    elif l[0] == 'm':
        mics[:,int(l[1])-1] = np.array([float(l[2]), float(l[3]), float(l[4])])
    elif l[0] == 's':
        sources[:,int(l[1])-1] = np.array([float(l[2]), float(l[3]), float(l[4])])
    else:
        continue


# Create the room
room = pra.ShoeBox3D(np.zeros(3), room_dim, Fs, 
        max_order=1, 
        absorption=reflection,
        sigma2_awgn=sigma2_n)

# Create the beamformer
bf = pra.Beamformer(mics, Fs, N=N, Lg=Lg)
room.addMicrophoneArray(bf)

# create a single reference mic at position of microphone 4
ref_mic_n = 4
ref_mic = pra.MicrophoneArray(bf.R[:,ref_mic_n,np.newaxis], Fs)

# Read in all the RIRs
f = open(data_folder + 'RIRs.files', 'r')
lines = f.readlines()
f.close()

# build empty rir array for the room
RIRs = [[] for r in range(n_mic)]

for line in lines:
    l = line.split(' ')
    r = int(l[0]) - 1  # microphone number
    s = int(l[1]) - 1  # speaker number

    # read wav file
    rir_fs,rir = wavfile.read(data_folder + l[3][:-1])

    # resample at new sampling rate if necessary
    if rir_fs != Fs:
        rir = resample(rir, np.ceil(rir.shape[0]/float(rir_fs)*Fs))

    RIRs[r].insert(s, rir)

# receptacle arrays
pesq_input = np.zeros(2)
pesq_bf = np.zeros(2)

# since we run multiple thread, we need to uniquely identify filenames
#pid = str(getpid())
pid = '0'

file_ref  = 'output_samples/fqref' + pid + '.wav'
file_suffix = '-' + pid + '.wav'
files_bf = 'output_samples/fq' + str(1) + file_suffix
file_raw  = 'output_samples/fqraw' + pid + '.wav'

# index of good and bad sources
good = 5
bad =  0

# Read the two speech samples used
rate, good_signal = wavfile.read(speech_sample1)
good_signal = np.array(good_signal, dtype=float)
good_signal = pra.normalize(good_signal)
good_signal = pra.highpass(good_signal, rate)
good_len = good_signal.shape[0]/float(Fs)

rate, bad_signal = wavfile.read(speech_sample2)
bad_signal = np.array(bad_signal, dtype=float)
bad_signal = pra.normalize(bad_signal)
bad_signal = pra.highpass(bad_signal, rate)
bad_len = bad_signal.shape[0]/float(Fs)

# variance of good signal
good_sigma2 = np.mean(good_signal**2)

# normalize interference signal to have equal power with desired signal
bad_signal *= good_sigma2/np.mean(bad_signal**2)

# pick good source position at random
good_distance = np.linalg.norm(bf.center[:,0] - np.array(sources[:,good]))

# pick bad source position at random
bad_distance = np.linalg.norm(bf.center[:,0] - np.array(sources[:,bad]))

if good_len > bad_len:
    good_delay = 0
    bad_delay = (good_len - bad_len)/2.
else:
    bad_delay = 0
    good_delay = (bad_len - good_len)/2.


# create the reference room for freespace, noisless, no interference simulation
ref_room = pra.ShoeBox3D(
    [0,0,0],
    room_dim,
    Fs,
    max_order=0)
ref_room.addSource(sources[:,good], signal=good_signal, delay=good_delay)
ref_room.addMicrophoneArray(ref_mic)
ref_room.compute_RIR()
ref_room.simulate()
reference = pra.highpass(ref_mic.signals[0], Fs)
reference_n = pra.normalize(reference)

# save the reference desired signal
wavfile.write(file_ref, Fs, pra.to_16b(reference_n))

# add the sources to the 'real' room
room.addSource(sources[:,good], 
        signal=good_signal, delay=good_delay)
room.addSource(sources[:,bad], 
        signal=bad_signal, delay=bad_delay)

# add the measured RIR to the room
for r in range(n_mic):
    room.rir.append([])
    room.rir[r].append(RIRs[r][good])
    room.rir[r].append(RIRs[r][bad])

# compute the input signal to the microphones
room.simulate()

# save degraded signal at reference microphone
raw = bf.signals[ref_mic_n]
raw_n = pra.normalize(pra.highpass(raw, Fs))
wavfile.write(file_raw, Fs, pra.to_16b(raw_n))

'''
pesq_input = pra.pesq(file_ref, file_raw, Fs=Fs)
print 'Input PESQ: ',pesq_input

for src in room.sources:
    src.setOrdering('nearest', ref_point=bf.center)

for i in range(7):

    good_img = room.sources[0][:i+1]
    bad_img = room.sources[1][:i+1]

    bf.rakePerceptualFilters(good_img, bad_img, sigma2_n*np.eye(n_mic*Lg))

    # run beamformer
    output = bf.process()
    output = pra.normalize(pra.highpass(output, Fs))
    output = pra.time_align(reference_n, output)

    # save files for PESQ evaluation
    wavfile.write(files_bf, Fs, pra.to_16b(output))

    # compute PESQ
    pesq_bf = pra.pesq(file_ref, files_bf, Fs=Fs)

    print str(i) + ' Image sources: ',pesq_bf
'''


plt.figure()
for m in range(n_mic):

    rir_sim = room.sources[0].getRIR(mics[:,m], Fs)
    plt.subplot(3,3,m+1)
    plt.plot(room.rir[m][0][:rir_sim.shape[0]])
    plt.plot(rir_sim)

plt.show()
