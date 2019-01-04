# This code creates a spectrogram of the data from a given channel, 
# then estimates the fringe frequency of scattered light based on the 
# top-mass motion of LIGO suspensions and overlays that on the spectrogram.
# Josh Smith 2015
# 
# To do: 
# Input argument handling
# Use timeseries dict to get all channels at once
# Make plots using Q-transform rather than spectrogram to better visualize 
# the fringes. 
# Load position channels in as a text file. 

# tell matplotlib to use a non-interactive backend for generating png images
from matplotlib import use
use('agg')
from gwpy.timeseries import TimeSeries
from gwpy.time import from_gps
from numpy import sqrt, array, arange, nonzero, argmax, log10
from pylab import savefig, specgram
import sys, string 
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
import matplotlib.mlab as mlab
from os import makedirs, path
import string 
import argparse
from gwpy.time import to_gps

__author__ = 'Joshua Smith'
__email__ = 'joshua.smith@ligo.org'
__version__ = '0.0.1'
prog = 'ScatMon3'

# command line options
parser = argparse.ArgumentParser(description=__doc__, prog=prog)

parser.add_argument('--start', type=to_gps, action='store', default=1230524574,
                    help='Start time GPS or yyyy-mm-dd hh:mm:ss.mmm ')
parser.add_argument('--end', type=to_gps, action='store', default=1230524584,
                    help='End time GPS or yyyy-mm-dd hh:mm:ss.mmm ')
parser.add_argument('-i', '--ifo', default='H1')
parser.add_argument('-p', '--primary', default='ASC-Y_TR_B_PIT_OUT_DQ')

args = parser.parse_args()

start_time = args.start
end_time = args.end
ifo = args.ifo
witness_base = args.primary

witness_chan = ifo + ':' + witness_base
print witness_chan

# load timeseries for witness channel
witness=TimeSeries.fetch(witness_chan, start_time, end_time, verbose=True)
if witness_base=="GDS-CALIB_STRAIN":
	witness=witness.highpass(20,gpass=3) # highpass the witness data

# Calculate DARM spectrogram 
secsPerFFT = .75 # Hz
overlap = 0.9 # fractional overlap
Fs = witness.sample_rate.value
NFFT = int(round(Fs*secsPerFFT))
noverlap = int(round(overlap*NFFT))
Pxx, freq, t, im = specgram(witness.value,NFFT=NFFT,Fs=witness.sample_rate.value,noverlap=noverlap,scale_by_freq='magnitude',detrend=mlab.detrend_linear,window=mlab.window_hanning)

### Known SUS displacement w/r/t sus frame channels (in microns)
position_chans = [\
'SUS-BS_M1_DAMP_L_IN1_DQ',\
'SUS-ETMX_M0_DAMP_L_IN1_DQ',\
'SUS-ETMX_R0_DAMP_L_IN1_DQ',\
'SUS-ETMY_M0_DAMP_L_IN1_DQ',\
'SUS-ETMY_R0_DAMP_L_IN1_DQ',\
'SUS-IM1_M1_DAMP_L_IN1_DQ',\
'SUS-IM2_M1_DAMP_L_IN1_DQ',\
'SUS-IM3_M1_DAMP_L_IN1_DQ',\
'SUS-IM4_M1_DAMP_L_IN1_DQ',\
'SUS-ITMX_M0_DAMP_L_IN1_DQ',\
'SUS-ITMX_R0_DAMP_L_IN1_DQ',\
'SUS-ITMY_M0_DAMP_L_IN1_DQ',\
'SUS-ITMY_R0_DAMP_L_IN1_DQ',\
'SUS-MC1_M1_DAMP_L_IN1_DQ',\
'SUS-MC2_M1_DAMP_L_IN1_DQ',\
'SUS-MC3_M1_DAMP_L_IN1_DQ',\
'SUS-OM1_M1_DAMP_L_IN1_DQ',\
'SUS-OM2_M1_DAMP_L_IN1_DQ',\
'SUS-OM3_M1_DAMP_L_IN1_DQ',\
'SUS-OM1_M1_DAMP_L_IN1_DQ',\
'SUS-OM2_M1_DAMP_L_IN1_DQ',\
'SUS-OM3_M1_DAMP_L_IN1_DQ',\
'SUS-OMC_M1_DAMP_L_IN1_DQ',\
'SUS-PR2_M1_DAMP_L_IN1_DQ',\
'SUS-PR3_M1_DAMP_L_IN1_DQ',\
'SUS-PRM_M1_DAMP_L_IN1_DQ',\
'SUS-RM1_M1_DAMP_L_IN1_DQ',\
'SUS-RM2_M1_DAMP_L_IN1_DQ',\
'SUS-SR2_M1_DAMP_L_IN1_DQ',\
'SUS-SR3_M1_DAMP_L_IN1_DQ',\
'SUS-SRM_M1_DAMP_L_IN1_DQ',\
'SUS-TMSX_M1_DAMP_L_IN1_DQ',\
'SUS-TMSY_M1_DAMP_L_IN1_DQ',\
'OMC-PZT2_MON_DC_OUT_DQ',\
]

### get data for each position channel, calculate fringe freq
for channel in position_chans:
	position_chan = '%s:%s' % (ifo,channel)
	print position_chan
	position=TimeSeries.fetch(position_chan, start_time, end_time, verbose=True)
	position=position.detrend()
	position=position.lowpass(5,gpass=3) # lowpass position channel
	velocity = (position[1:]-position[:-1])*(position.sample_rate.value)
	times = (arange(len(velocity.value))+0.5)/(position.sample_rate.value)
	fudge = 1.0
	scatf1 = abs(fudge*2.0*velocity.value/1.064) # single bounce scatter frequency Virgo Scatter eqn 3
	scatf2 = 2.0*scatf1
	scatf3 = 3.0*scatf1
	scatf4 = 4.0*scatf1
	scatf8 = 8.0*scatf1
	scatf13 = 13.0*scatf1
	# print max scatter values and their times
	print 'max scatter f1 = ' + str(max(scatf1)) + ' Hz'
	tofmax = times[argmax(scatf2)]
	tofmaxgps = tofmax + start_time
	print 'time of max f2 = ' + str(tofmax) + ' s, GPS=' + str(tofmaxgps)

	fig = plt.figure(figsize=(12,12))
	ax1 = fig.add_subplot(211)
	# Plot Spectrogram
        im1 = NonUniformImage(ax1, interpolation='bilinear',extent=(min(t),max(t),10,55),cmap='jet')
        im1.set_data(t,freq,20.0*log10(Pxx))
        if witness_base=="GDS-CALIB_STRAIN":
		print "setting color limits for STRAIN"
		im1.set_clim(-1000,-800)
        #elif witness_base=="ASC-AS_B_RF45_Q_YAW_OUT_DQ" or witness_base=="ASC-AS_B_RF36_Q_PIT_OUT_DQ" or witness_base=="ASC-AS_A_RF45_Q_PIT_OUT_DQ" or witness_base=="LSC-MICH_IN1_DQ":
	#	im1.set_clim(-200,20)
	#elif witness_base == "OMC-LSC_SERVO_OUT_DQ":
	im1.set_clim(-400,-220)
	ax1.images.append(im1)
        cbar1 = fig.colorbar(im1)
        #cbar1.set_clim(-540,-240)
	
	# plot fringe prediction timeseries
	#ax1.plot(times,scatf5, c='blue', linewidth='0.2', label='f5')
	ax1.plot(times,scatf13, c='purple', linewidth='0.6', label='f13')
	ax1.plot(times,scatf8, c='green', linewidth='0.6', label='f8')
	ax1.plot(times,scatf4, c='blue', linewidth='0.6', label='f4')
	ax1.plot(times,scatf2, c='black', linewidth='0.6', label='f2')
	#if plotspec < 1 and dur <= 3600:
        #	ax1.plot(highscattimes, highscatf2, c='r', label='veto', marker='.',linestyle='')
	ax1.legend( loc='upper right' )
	title = 'Specgram of ' + string.replace(witness_chan,'_','\_') + '\nScattering fringe frequency for ' + string.replace(position_chan,'_','\_')
       	#title = 'Scattering fringe frequency for ' + string.replace(position_chan,'_','\_')
	# plt.title(title,fontsize=16)
        plt.yscale('linear')
	plt.ylabel('Frequency [Hz]',fontsize=16)
	xlab = 'Time [seconds] from ' + str(args.start)
        plt.xlim(0,max(times))
        plt.ylim(0,35)
	plt.title(title)

	ax2 = fig.add_subplot(212)
	#print(len(position.data))
	#print(position.data)
	ax2.plot(times,position[1:])
	plt.xlabel(xlab,fontsize=16)
	plt.xlim(0,max(times))
	plt.ylabel('Position [um]',fontsize=16)

	# save figure
	dir = str(start_time) + '-' + str(end_time-start_time) + '/'
	if not path.exists(dir):
		makedirs(dir)
	filename = dir + ifo + '-Fscatter-' + channel + '-' + str(start_time) + '-' + str(end_time-start_time) + '.png'
	savefig(filename)
	plt.close(fig)

quit()
# EOF
