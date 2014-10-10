
import numpy as np
import matplotlib.pyplot as plt

import pyroomacoustics as pra

beamformer_names = ['Rake MaxSINR', 'Rake Perceptual', 'Rake MVDR']
SINR = np.load('data/SINR_data.npy')
# uncomment the following line to use the simulated data used in the paper
#SINR = np.load('data/SINR_data_Lg30ms_d20ms_20141008.npy')

max_K, n_bf, n_monte_carlo = SINR.shape
SINR_med = np.percentile(SINR, [50, 5, 95], axis=-1)

SINR_gain_5sources = pra.dB(SINR_med[0,5,:]) - pra.dB(SINR_med[0,0,:])
print 'SNR gain of using 5 sources instead of one:'
print 'Rake MaxSINR: %.2f dB' % SINR_gain_5sources[0]
print 'Rake Perceptual: %.2f dB' % SINR_gain_5sources[1]
print 'Rake MVDR: %.2f dB' % SINR_gain_5sources[2]

#---------------------------------------------------------------------
# Export the SNR figure
#---------------------------------------------------------------------

plt.figure(figsize=(4, 3))

newmap = plt.get_cmap('gist_heat')
ax1 = plt.gca()
ax1.set_color_cycle([newmap( k ) for k in np.linspace(0.25,0.9,len(beamformer_names))])

from itertools import cycle
lines = ['-s','-o','-v','-D','->']
linecycler = cycle(lines)

k_axis = np.arange(0,max_K, dtype=float)
k_axis[0] += 0.1
k_axis[-1] -= 0.1
print k_axis

for i in np.arange(n_bf):
    p, = plt.plot(k_axis,
                  pra.dB(SINR_med[0,:,i]),
                  next(linecycler),
                  linewidth=1,
                  markersize=4,
                  markeredgewidth=.5)

plt.fill_between(k_axis,
                 pra.dB(SINR_med[1,:,0]), pra.dB(SINR_med[2,:,0]),
                 color='grey',
                 linewidth=0.3,
                 edgecolor='k',
                 alpha=0.7)

# Hide right and top axes
ax1.set_xlim(0,9)
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
plt.ylabel('Output SINR [dB]', fontsize=10)
plt.tight_layout()


plt.legend(beamformer_names, fontsize=7, loc='upper left', frameon=False, labelspacing=0)

plt.savefig('figures/SINR_vs_K.pdf')

plt.close()

