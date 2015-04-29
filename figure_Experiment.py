

def perceptual_quality_evaluation(room_dim, mics, good_pos, good_index, bad_pos, bad_index, rir_location):
    print 'start'

    import numpy as np
    from scipy.io import wavfile
    from os import getpid

    import pyroomacoustics as pra

    # number of sources to  consider
    n_sources = np.arange(1,8)
    S = n_sources.shape[0]

    # number of mics
    n_mic = mics.shape[1]

    # this is the speed of sound
    # it is now set separately in the pyroomacoustics package.
    # We will need a way to set it package wide.
    c = 344.5

    Fs = 8000.
    N = 1024
    Lg = int(0.03*Fs) # 350 ms long filter
    delay_bf = 0.02
    sigma2_n = 1e-6

    # reflection coefficients from the walls (hand-waving)
    reflection = {'ground':0.8, 'south':0.8, 'west':0.8, 'north':0.8, 'east':0.8, 'ceilling':0.5}

    speech_sample1 = 'samples/fq_sample1_8000.wav'
    speech_sample2 = 'samples/fq_sample2_8000.wav'

    # Create the room
    room = pra.ShoeBox3D(np.zeros(3), room_dim, Fs, 
            max_order=1, 
            absorption=reflection,
            sigma2_awgn=sigma2_n)

    # Create the beamformer
    bf = pra.Beamformer(mics, Fs, N=N, Lg=Lg)
    room.addMicrophoneArray(bf)

    # data receptacles
    beamformer_names = ['Rake Perceptual',
                        'Rake MVDR']
    bf_weights_fun   = [bf.rakePerceptualFilters,
                        bf.rakeMVDRFilters]
    bf_fnames = ['1','2']
    NBF = len(beamformer_names)

    # receptacle arrays
    pesq_input = np.zeros(2)
    pesq_bf = np.zeros((2,NBF,S))

    # create a single reference mic at position of microphone 4
    ref_mic_n = 4
    ref_mic = pra.MicrophoneArray(bf.R[:,ref_mic_n,np.newaxis], Fs)

    # since we run multiple thread, we need to uniquely identify filenames
    pid = str(getpid())

    file_ref  = 'output_samples/fqref' + pid + '.wav'
    file_suffix = '-' + pid + '.wav'
    files_bf = ['output_samples/fq' + str(i+1) + file_suffix for i in xrange(NBF)]
    file_raw  = 'output_samples/fqraw' + pid + '.wav'

    # index of good and bad sources
    good = good_index
    bad =  bad_index

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
    good_distance = np.linalg.norm(bf.center[:,0] - good_pos)

    # pick bad source position at random
    bad_distance = np.linalg.norm(bf.center[:,0] - bad_pos)

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
    ref_room.addSource(good_pos, signal=good_signal, delay=good_delay)
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
    room.addSource(good_pos, signal=good_signal, delay=good_delay)
    room.addSource(bad_pos, signal=bad_signal, delay=bad_delay)

    # read in the RIR from file
    for r in range(n_mic):
        for s in [good_index, bad_index]:

            # read wav file
            fname_rir = rir_location % (r+1,s+1)
            rir_fs,rir = wavfile.read(fname_rir)
            rir = np.array(rir, dtype='float64')

            if rir_fs != Fs:
                raise NameError('The RIR and the signals do not have the same sampling rate.')
                '''
                import scikits.samplerate as sr
                rir = sr.resample(rir, Fs/float(rir_fs), 'sinc_best')

                # the factor 2 was empirically determined to be necessary to get
                # amplitude of RIR in the correct ballpark.
                rir *= 2.
                '''

            room.rir.append([])
            room.rir[r].append(rir)

    # compute the input signal to the microphones
    room.simulate()

    # save degraded signal at reference microphone
    raw = bf.signals[ref_mic_n]
    raw_n = pra.normalize(pra.highpass(raw, Fs))
    wavfile.write(file_raw, Fs, pra.to_16b(raw_n))

    pesq_input = pra.pesq(file_ref, file_raw, Fs=Fs)

    for src in room.sources:
        src.setOrdering('strongest', ref_point=bf.center)

    for k,s in enumerate(n_sources):

        good_img = room.sources[0][:s]
        bad_img = room.sources[1][:s]

        for i, bfr in enumerate(beamformer_names):

            bf_weights_fun[i](good_img, bad_img, sigma2_n*np.eye(n_mic*Lg), delay=delay_bf)

            # run beamformer
            output = bf.process()
            output = pra.normalize(pra.highpass(output, Fs))
            output = pra.time_align(reference_n, output)

            # save files for PESQ evaluation
            wavfile.write(files_bf[i], Fs, pra.to_16b(output))

            # compute PESQ
            x = pra.pesq(file_ref, files_bf[i], Fs=Fs)
            pesq_bf[:,i,k] = pra.pesq(file_ref, files_bf[i], Fs=Fs).T

    ''' This is how you can compare the true RIRs with the image src model generated one
    plt.figure()
    for m in range(n_mic):

        rir_sim = room.sources[0].getRIR(mics[:,m], Fs)
        plt.subplot(3,3,m+1)
        plt.plot(room.rir[m][0][:rir_sim.shape[0]])
        plt.plot(rir_sim)

    plt.show()
    '''

    print 'Finished'

    return pesq_input, pesq_bf

if __name__ == '__main__':

    import numpy as np
    import sys
    import time

    # Process arguments
    ###################

    if len(sys.argv) == 3 and sys.argv[1] == '-s':
        parallel = False
        Loops = int(sys.argv[2])
    elif len(sys.argv) == 2:
        parallel = True
        Loops = int(sys.argv[1])
    else:
        print 'Usage: ipython %s -- [-s] <loop_number>' % (sys.argv[0])
        print '       -s: Serial loop, no parallelism used.'
        sys.exit(0)

    # This is the location of simulated/measured RIRs
    #################################################

    rir_type = 'measured'
    #rir_type = 'simulated'

    if rir_type == 'measured':
        data_folder = 'Experiment_Data/'
    elif rir_type == 'simulated':
        data_folder = 'Raven/'
    else:
        raise NameError('RIR type needs to be ''measured'' or ''simulated''.')


    # read in room size, microphones and sources locations
    ######################################################

    f = open(data_folder + 'RIRs.positions', 'r')
    lines = f.readlines()
    f.close()

    # count # of mics and speakers
    n_src = len([l for l in lines if l.split(' ')[0] == 's'])
    n_mic = len([l for l in lines if l.split(' ')[0] == 'm'])

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


    # PREPARE PARAMETERS
    ####################

    if rir_type == 'measured':

        Loops = np.minimum(n_src*(n_src-1), Loops)
        
        good_src_index = np.zeros(Loops)
        bad_src_index = np.zeros(Loops)
        good_source = np.zeros((sources.shape[0], Loops))
        bad_source = np.zeros((sources.shape[0], Loops))
        from itertools import product
        i = 0
        for s1,s2 in product(np.arange(n_src),np.arange(n_src)):

            if i == Loops:
                break

            if s1 == s2:
                continue
            else:
                good_src_index[i] = s1
                good_source[:,i] = sources[:,s1]
                bad_src_index[i] = s2
                bad_source[:,i] = sources[:,s2]
                i += 1


    elif rir_type == 'simulated':

        Loops = np.minimum(n_src/2, Loops)

        good_src_index = np.arange(n_src/2)
        bad_src_index = np.arange(n_src/2,n_src)
        good_source = sources[:,:n_src/2]
        bad_source = sources[:,n_src/2:]

    # room dimension and location RIR files really just need to
    # to be repeated in an array
    rir_location_arr = [data_folder + '/RIRs/rir_%d_%d.wav' for i in range(Loops)]
    room_dim_arr = [room_dim for i in range(Loops)]
    mics_arr = [mics for i in range(Loops)]

    # SIMULATION
    ############

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

        out = view.map_sync(
                perceptual_quality_evaluation, 
                room_dim_arr, 
                mics_arr,
                good_source.T, good_src_index, 
                bad_source.T, bad_src_index, 
                rir_location_arr
                )

    else:
        # Just one boring loop...
        out = []
        for i in xrange(Loops):
            out.append(perceptual_quality_evaluation(
                room_dim_arr[i], 
                mics_arr[i],
                good_source[:,i], good_src_index[i], 
                bad_source[:,i], bad_src_index[i],
                rir_location_arr[i]
                ))

    # How long was this ?
    ellapsed = time.time() - start

    # how long was this ?
    print('Time ellapsed: ' + str(ellapsed))

    # recover all the data
    pesq_input = np.array([o[0] for o in out])
    pesq_bf = np.array([o[1] for o in out])

    # save the simulation results to file
    filename = 'sim_data/quality_' + rir_type + '_rir_' + time.strftime('%Y%m%d-%H%M%S') + '.npz'
    np.savez_compressed(filename, good_source=good_source, bad_source=bad_source,
            good_src_index=good_src_index, bad_src_index=bad_src_index,
            pesq_bf=pesq_bf, pesq_input=pesq_input)

