
import numpy as np
from scipy.io import wavfile
from scipy.signal import resample, decimate
from os import getpid

import matplotlib.pyplot as plt

import pyroomacoustics as pra

data_folder = 'Experiment_Data/'
#data_folder = 'Raven/'

# this is the speed of sound
# it is now set separately in the pyroomacoustics package.
# We will need a way to set it package wide.
c = 344.5

Fs = 8000.
N = 1024
Lg = int(0.150*Fs) # 350 ms long filter
delay_bf = 0.04
sigma2_n = 1e-6

speech_sample1 = 'samples/fq_sample1_8000.wav'
speech_sample2 = 'samples/fq_sample2_8000.wav'

# read in microphones and sources locations
f = open(data_folder + 'RIRs.positions', 'r')
lines = f.readlines()
f.close()

# count # of mics and speakers
n_src = len([l for l in lines if l.split(' ')[0] == 's'])
n_mic = len([l for l in lines if l.split(' ')[0] == 'm'])

# reflection coefficients from the walls (hand-waving)
reflection = {'ground':0.8, 'south':0.8, 'west':0.8, 'north':0.8, 'east':0.8, 'ceilling':0.5}

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
    rir = np.array(rir, dtype='float64')

    if rir_fs != Fs:
        import scikits.samplerate as sr
        rir = sr.resample(rir, Fs/float(rir_fs), 'sinc_best')

        # the factor 2 was empirically determined to be necessary to get
        # amplitude of RIR in the correct ballpark.
        rir *= 2.

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
good_signal = np.array(good_signal, dtype='float64')
good_signal = pra.normalize(good_signal)
good_signal = pra.highpass(good_signal, rate)
good_len = good_signal.shape[0]/float(Fs)

rate, bad_signal = wavfile.read(speech_sample2)
bad_signal = np.array(bad_signal, dtype='float64')
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
#wavfile.write(file_ref, Fs, pra.to_16b(reference_n))

new_ref = good_signal.copy()
new_ref = pra.normalize(pra.highpass(new_ref, Fs))
wavfile.write(file_ref, Fs, pra.to_16b(new_ref))

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

pesq_input = pra.pesq(file_ref, file_raw, Fs=Fs)
print 'Input PESQ: ',pesq_input

for src in room.sources:
    src.setOrdering('nearest', ref_point=bf.center)

for i in range(7):

    good_img = room.sources[0][:i+1]
    bad_img = room.sources[1][:i+1]

    bf.rakePerceptualFilters(good_img, bad_img, sigma2_n*np.eye(n_mic*Lg), delay=delay_bf)

    # run beamformer
    output = bf.process()
    output = pra.normalize(pra.highpass(output, Fs))
    output = pra.time_align(reference_n, output)

    files_bf = 'output_samples/fq' + str(i) + file_suffix

    # save files for PESQ evaluation
    wavfile.write(files_bf, Fs, pra.to_16b(output))

    # compute PESQ
    pesq_bf = pra.pesq(file_ref, files_bf, Fs=Fs)

    print str(i) + ' Image sources: ',pesq_bf


plt.figure()
for m in range(n_mic):

    rir_sim = room.sources[0].getRIR(mics[:,m], Fs)
    plt.subplot(3,3,m+1)
    plt.plot(room.rir[m][0][:rir_sim.shape[0]])
    plt.plot(rir_sim)

plt.show()
'''

if __name__ == '__main__':

    import numpy as np
    import sys
    import time

    if len(sys.argv) == 3 and sys.argv[1] == '-s':
        parallel = False
        Loops = int(sys.argv[2])
    elif len(sys.argv) == 2:
        parallel = True
        Loops = int(sys.argv[1])
    else:
        print 'Usage: ipython figure_quality_sim.py -- [-s] <loop_number>'
        print '       -s: Serial loop, no parallelism used.'
        sys.exit(0)

    # we restrict sources to be in a square 1m away from every wall and from the array
    bbox_size = np.array([[2.,2.5]])
    bbox_origin = np.array([[1.,2.5]])

    # draw all target and interferer at random
    good_source = np.random.random((Loops,2))*bbox_size + bbox_origin
    bad_source = np.random.random((Loops,2))*bbox_size + bbox_origin

    # start timing simulation
    start = time.time()

    if parallel is True:
        # Launch many workers!
        from IPython import parallel

        # setup parallel computation env
        c = parallel.Client()
        print c.ids
        c.blocks = True
        view = c.load_balanced_view()

        out = view.map_sync(perceptual_quality_evaluation, good_source, bad_source)

    else:
        # Just one boring loop...
        out = []
        for i in xrange(Loops):
            out.append(perceptual_quality_evaluation(good_source[i,:], bad_source[i,:]))

    # How long was this ?
    ellapsed = time.time() - start

    # how long was this ?
    print('Time ellapsed: ' + str(ellapsed))

    # recover all the data
    pesq_input = np.array([o[0] for o in out])
    pesq_trinicon = np.array([o[1] for o in out])
    pesq_bf = np.array([o[2] for o in out])
    isinr = np.array([o[3] for o in out])
    osinr_trinicon = np.array([o[4] for o in out])
    osinr_bf = np.array([o[5] for o in out])

    # save the simulation results to file
    filename = 'sim_data/quality_' + time.strftime('%Y%m%d-%H%M%S') + '.npz'
    np.savez_compressed(filename, good_source=good_source, bad_source=bad_source,
            isinr=isinr, osinr_bf=osinr_bf, osinr_trinicon=osinr_trinicon,
            pesq_bf=pesq_bf, pesq_input=pesq_input, pesq_trinicon=pesq_trinicon)

'''
