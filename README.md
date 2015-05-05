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


Measured Room Impulse Responses
-------------------------------

Beyond the simulation results presented in the paper, the raking beamformers
were validated with room impulse responses (RIR) measured in a classroom at
EPFL. The RIR are stored as wav files in `BC329_RIR_8kHz/RIRs` and the
methodology for the measurements is outlined in `BC329_RIR_8kHz/README.md`.

The RIR were measured for a static linear array of 8 microphones spaced by 8
cm. RIRs were measured for 12 different source positions. In the validation all
132 combinations of desired/interfering sources are tested. Sound propagation
is simulated by convolving speech samples with the RIR. 

Using the measured room dimension and locations of sources, as well as
approximate reflection coefficients for the walls, both Rake MVDR and Rake
Perceptual beamformers are computed. The output of the beamformers is evaluated
with respect to the clean desired speech signal using the perceptual evaluation
of speech quality (PESQ) metric.

Even though the calibration of room characteristics, microphones, and speakers locations
was only approximate, we observe the same kind of thread than in the pure SINR simulation.
That is speech quality is improved by using extra reflections in the room.

<img src="https://raw.githubusercontent.com/LCAV/TimeDomainAcousticRakeReceiver/master/figures/pesq_measured_rir.png">

The simulation can be run as `ipython figure_Experiment.py 132` and the output of the
simulation can be plotted using `ipython figure_Experiment_plot.py` (with optional
simulation output filename argument).


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
* An implementation of PESQ is needed to run `figure_Experiment.py`. In a pinch, the reference
  implementation of the ITU can be used. To compile it follow the instructions.

### PESQ Tool

Download the [source files](http://www.itu.int/rec/T-REC-P.862-200511-I!Amd2/en) 
of the ITU P.862 compliance tool from the ITU website.  Follow the instructions
to compile the tool and copy the executable produced in the `bin` directory.

#### Unix compilation (Linux/Mac OS X)

Execute the following sequence of commands to get to the source code.

    mkdir PESQ
    cd PESQ
    wget 'https://www.itu.int/rec/dologin_pub.asp?lang=e&id=T-REC-P.862-200511-I!Amd2!SOFT-ZST-E&type=items'
    unzip dologin_pub.asp\?lang\=e\&id\=T-REC-P.862-200511-I\!Amd2\!SOFT-ZST-E\&type\=items
    cd Software
    unzip 'P862_annex_A_2005_CD  wav final.zip'
    cd P862_annex_A_2005_CD/source/

In the `Software/P862_annex_A_2005_CD/source/` directory, create a file called `Makefile` and copy
the following into it.

    CC=gcc
    CFLAGS=-O2

    OBJS=dsp.o pesqdsp.o pesqio.o pesqmod.o pesqmain.o
    DEPS=dsp.h pesq.h pesqpar.h

    %.o: %.c $(DEPS)
      $(CC) -c -o $@ $< $(CFLAGS)

    pesq: $(OBJS)
      $(CC) -o $@ $^ $(CFLAGS)

    .PHONY : clean
    clean :
      -rm pesq $(OBJS)

Execute compilation by typing this.

    make pesq

Finally move the `pesq` binary to `<repo_root>/bin/`.

Notes:
* The files input to the pesq utility must be 16 bit PCM wav files.
* File names longer than 14 characters (suffix included) cause the utility to
  crash with the message `Abort trap(6)` or similar.

#### Windows compilation

1. Open visual studio, create a new project from existing files and select the directory
  containing the source code of PESQ (`Software\P862_annex_A_2005_CD\source\`).

          FILE -> New -> Project From Existing Code...

2. Select `Visual C++` from the dropdown menu, then next.
    * *Project file location* : directory containing source code of pesq (`Software\P862_annex_A_2005_CD\source\`).
    * *Project Name* : pesq
    * Then next.
    * As *project type*, select `Console application` project.
    * Then finish.

3. Go to

          BUILD -> Configuration Manager...

    and change active solution configuration from `Debug` to `Release`. Then Close.

4. Then 

          BUILD -> Build Solution

5. Copy the executable `Release\pesq.exe` to the bin folder.

*(tested with Microsoft Windows Server 2012)*


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

