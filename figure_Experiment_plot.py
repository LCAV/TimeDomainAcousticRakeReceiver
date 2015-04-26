import numpy as np
import matplotlib.pyplot as plt

import sys
import os
import fnmatch

import pyroomacoustics as pra

max_sources = 7
sim_data_dir = './sim_data/'

beamformer_names = ['Rake Perceptual',
                    'Rake MVDR',]
bf_dict = dict(zip(beamformer_names, 
               range(len(beamformer_names))))
NBF = len(beamformer_names)

loops = 0

if len(sys.argv) == 0:
    # if no argument is specified, use all available files
    name_pattern = './sim_data/quality_2015*.npz'
    files = [file for file in os.listdir(sim_data_dir) if fnmatch.fnmatch(file, name_pattern)]
else:
    files = sys.argv[1:]

# Empty data containers
good_source = np.zeros((3,0))
bad_source = np.zeros((3,0))
ipesq = np.zeros((0,2))
opesq_bf = np.zeros((0,2,NBF,max_sources))

# Read in all the data
for fname in files:
    print 'Loading from',fname

    a = np.load(fname)

    good_source = np.concatenate((good_source, a['good_source']), axis=1)
    bad_source = np.concatenate((bad_source, a['bad_source']), axis=1)

    ipesq = np.concatenate((ipesq,a['pesq_input'][:,:,0]), axis=0)
    opesq_bf = np.concatenate((opesq_bf,a['pesq_bf']), axis=0)

loops = good_source.shape[1]

print 'Number of loops:',loops
print 'Median input Raw MOS',np.median(ipesq[:,0])
print 'Median input MOS LQO',np.median(ipesq[:,1])

def nice_plot(x, ylabel, bf_order=None):
    '''
    Define a function to plot consistently the data
    '''

    if bf_order is None:
        bf_order = beamformer_names

    ax1 = plt.gca()

    newmap = plt.get_cmap('gist_heat')
    from itertools import cycle

    # totally a hack to get the same line styles as Fig6/7
    lines = ['-D','-v','->','-s','-o']
    linecycler = cycle(lines)

    # totally a hack to get the same line styles as Fig6/7
    map1 = [newmap( k ) for k in np.linspace(0.25,0.9,5)]
    map2 = [map1[3],map1[2],map1[4],map1[0],map1[1]]

    ax1.set_color_cycle(map2)

    # no clipping of the beautiful markers
    plt.setp(ax1,'clip_on',False)

    for bf in bf_order:
        i = bf_dict[bf]
        med, ci = pra.median(x, axis=0)
        p, = plt.plot(range(0, max_sources), 
                #np.median(x[:,i,:], axis=0),
                med[i,:],
            next(linecycler),
            linewidth=1,
            markersize=4,
            markeredgewidth=.5,
            clip_on=False)

        # confidence interval for the median
        plt.fill_between(range(0, max_sources),
                med[i,:]+ci[0,i,:], med[i,:]+ci[1,i,:],
                #color='grey',
                color=p.get_color(),
                linewidth=0.05,
                edgecolor='k',
                alpha=0.3)

        # Hide right and top axes
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_position(('outward', 10))
        ax1.spines['left'].set_position(('outward', 15))
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')

        # Make ticks nicer
        ax1.xaxis.set_tick_params(width=.3, length=3)
        ax1.yaxis.set_tick_params(width=.3, length=3)

        # Make axis lines thinner
        for axis in ['bottom','left']:
          ax1.spines[axis].set_linewidth(0.3)

        # Set ticks fontsize
        plt.xticks(size=9)
        plt.yticks(size=9)

        # Set labels
        plt.xlabel(r'Number of images $K$', fontsize=10)
        plt.ylabel(ylabel, fontsize=10)

        plt.legend(bf_order, fontsize=7, loc='upper left', frameon=False, labelspacing=0)

# Here we plot the figure used in the paper (Fig. 10)
plt.figure(figsize=(4,3))
nice_plot(opesq_bf[:,0,:,:], 'PESQ [MOS]', 
        bf_order=['Rake Perceptual','Rake MVDR'])

# plot input SNR
med, ci = pra.median(ipesq[:,0])

o = np.ones(max_sources)
p, = plt.plot(np.arange(max_sources), np.median(ipesq[:,0])*o)
plt.text(5.0, 1.28, 'Input PESQ', fontsize=7)

'''
# confidence interval for the median
plt.fill_between(range(0, max_sources),
        (med+ci[0])*o, (med+ci[1])*o,
        #color='grey',
        color=p.get_color(),
        linewidth=0.05,
        edgecolor='k',
        alpha=0.3)
'''


plt.tight_layout(pad=0.2)
plt.savefig('figures/pesq_measured_rir.pdf')
plt.show()

