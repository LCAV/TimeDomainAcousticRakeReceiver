Time Domain Acoustic Rake Receiver
==================================

This repository contains all the code to reproduce the results of the paper
[*Raking Echoes in the Time Domain*](https://infoscience.epfl.ch/record/202223).

Using the simple python room acoustics framework created for the 
[Acoustic Rake Receiver](https://github.com/LCAV/AcousticRakeReceiver),
we demonstrate two time domain beamformer designs that use echoes constructively
to improve the signal to interference and noise ratio.
These designs are explained in details in the paper *Raking Echoes in the Time Domain*.
All the figures of the paper can be recreated by calling
simple scripts leveraging this framework. In addition to the results of
the paper, we include spectrograms and samples of speech samples.
We strongly hope that this code will be useful beyond the scope of this paper
and plan to develop it into a standalone python package in the future.

We are available for any question or request relating to either the code
or the theory behind it. Just ask!


Abstract
--------

The geometry of room acoustics is such that the reverberant signal can be seen
as a single signal emitted from multiple locations. In analogy with the rake
receiver from wireless communications, we propose several beamforming
strategies that exploit, rather than suppress, this extra temporal diversity.
Unlike earlier work in the frequency domain, time domain designs allow to
control important timing parameters. Particularly important is the ability to
control perceptually relevant parameters such as the amount of early echoes or
the length of the beamformer response.  

Relying on the knowledge of the position of sources, we derive several optimal
beamformers. Leveraging perceptual cues, we stray from distortionless design
and improve interference and noise reduction without degrading the perceptual
quality. The designs are validated through simulation. Using early echoes is
shown to strictly improve the signal to interference and noise ratio. Code and
speech samples are available online in this repository.


Authors
-------

Robin Scheibler, Ivan Dokmanić, and Martin Vetterli are with 
Laboratory for Audiovisual Communications ([LCAV](http://lcav.epfl.ch)) at 
[EPFL](http://www.epfl.ch).

<img src="http://lcav.epfl.ch/files/content/sites/lcav/files/images/Home/LCAV_anim_200.gif">

#### Contact

[Robin Scheibler](mailto:robin[dot]scheibler[at]epfl[dot]ch) <br>
EPFL-IC-LCAV <br>
BC Building <br>
Station 14 <br>
1015 Lausanne


Recreate the figures and sound samples
--------------------------------------

In a UNIX terminal, run the following script.

    ./make_all_figures.sh

Alternatively, type in the following commands in an ipython shell.

    run figure_beampatterns.py
    run figure_SINR_sim.py
    run figure_SINR_plot.py

The figures and sound samples generated are collected in `figures` and
`output_samples`, respectively.

The SINR simulation can be very long. For that reason, we split
the simulation loop and the plotting routine into two files.
We recommend to first modify the value of the `loops` variable
in `figure_SINR_sim.py` to 10 and time the simulation to have an
estimate of the total time required to run the simulatioin.
The output of the simulation is saved in `data/SINR_data.npy`.
The data can be plotted using `figure_SINR_plot.py`.


Simulation Data
---------------

The simulation data from the paper is available in `data/SINR_data_Lg30ms_d20ms_20141008.npy`
and can be processed with `figure_SINR_plot.py`.


Sound Samples and Spectrograms
------------------------------

While it was omitted in the paper for space reasons, we provide
here simulation results for speech samples.

* [sample1](https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/output_samples/input.wav) 
  Simulated microphone input signal.
* [sample2](https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/output_samples/output_DirectMVDR.wav) 
  Output of MVDR using the direct sound only.
* [sample3](https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/output_samples/output_RakeMVDR.wav) 
  Output of Rake MVDR using the direct sound and 1st order echoes.
* [sample4](https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/output_samples/output_DirectPerceptual.wav) 
  Output of Perceptual beamformer using the direct sound only.
* [sample5](https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/output_samples/output_RakePerceptual.wav) 
  Output of Rake Perceptual using the direct sound and 1st order echoes.

The spectrogram of all samples as well as the desired sound are provided.

<img src="https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/figures/spectrograms.png">

The samples and the spectrograms are created by the script `figure_spectrogam_and_samples.py`.


Extra Scripts
-------------

We also include extra scripts that let us play with the different beamformers.
The script names should be self-explanatory.

    RakeLstSqFromFD.py
    RakeMVDR.py
    RakeMaxSINR.py
    RakeMaxUDR.py
    RakeOF.py
    RakePerceptual.py


Dependencies
------------

* A working distribution of [Python 2.7](https://www.python.org/downloads/).
* The code relies heavily on [Numpy](http://www.numpy.org/),
  [Scipy](http://www.scipy.org/), and [matplotlib](http://matplotlib.org).
* We use the distribution [anaconda](https://store.continuum.io/cshop/anaconda/) to simplify the setup
  of the environment.
* Computations are very heavy and we use the
  [MKL](https://store.continuum.io/cshop/mkl-optimizations/) extension of
  Anaconda to speed things up. There is a [free license](https://store.continuum.io/cshop/academicanaconda) for academics.


Systems Tested
--------------

###Linux

| Machine | IBM System X iDataPlex dx360 M3 |
|---------|---------------------------------|
| System  | Ubuntu 12.04.5                  |
| CPU     | 2x X5675                        |
| RAM     | 96 GB                           |

    System Info:
    Linux 3.11.0-26-generic #45~precise1-Ubuntu SMP Tue Jul 15 04:02:35 UTC 2014 x86_64 x86_64 x86_64 GNU/Linux

    Python Info:
    Python 2.7.8 :: Anaconda 2.0.1 (64-bit)

    Python Packages Info:
    accelerate                1.7.0               np19py27_p0  
    anaconda                  2.0.1                np18py27_0  
    ipython                   2.1.0                    py27_2  
    ipython-notebook          2.1.0                    py27_0  
    ipython-qtconsole         2.1.0                    py27_0  
    matplotlib                1.3.1                np18py27_1  
    mkl                       11.1                np19py27_p3  
    mkl-rt                    11.1                         p0  
    mkl-service               1.0.0                   py27_p1  
    mklfft                    1.0                 np19py27_p0  
    numpy                     1.9.0                   py27_p0  [mkl]
    scipy                     0.14.0              np19py27_p0  [mkl]

###OS X

| Machine | MacBook Pro Retina 15-inch, Early 2013 |
|---------|----------------------------------------|
| System  | OS X Maverick 10.9.4                   |
| CPU     | Intel Core i7                          |
| RAM     | 16 GB                                  |

    System Info:
    Darwin 13.3.0 Darwin Kernel Version 13.3.0: Tue Jun  3 21:27:35 PDT 2014; root:xnu-2422.110.17~1/RELEASE_X86_64 x86_64

    Python Info:
    Python 2.7.8 :: Anaconda 2.0.1 (x86_64)

    Python Packages Info:
    accelerate                1.7.0               np19py27_p0  
    anaconda                  2.0.1                np18py27_0  
    ipython                   2.1.0                    py27_2  
    ipython-notebook          2.1.0                    py27_0  
    ipython-qtconsole         2.1.0                    py27_0  
    matplotlib                1.3.1                np18py27_1  
    mkl                       11.1                np19py27_p3  
    mkl-rt                    11.1                         p0  
    mkl-service               1.0.0                   py27_p1  
    mklfft                    1.0                 np19py27_p0  
    numpy                     1.9.0                   py27_p0  [mkl]
    scipy                     0.14.0              np19py27_p0  [mkl]


License
-------

Copyright (c) 2014, Robin Scheibler, Ivan Dokmanić, Martin Vetterli

This code is free to reuse for non-commercial purpose such as academic or
educational. For any other use, please contact the authors.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Time Domain Acoustic Rake Receiver</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://lcav.epfl.ch" property="cc:attributionName" rel="cc:attributionURL">Robin Scheibler, Ivan Dokmanić, Martin Vetterli</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/LCAV/AcousticRakeReceiver" rel="dct:source">https://github.com/LCAV/AcousticRakeReceiver</a>.

