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

# read input argumentsi (update this to use argparse)
if len(sys.argv) == 6:
	ifo = sys.argv[1] # H1 or L1
	start_time = int(sys.argv[2]) # GPS start time
	startutc = str(from_gps(start_time)) # get UTC time
	dur = int(sys.argv[3]) # duration in seconds
	thresh = int(sys.argv[4]) # threshold frequency [Hz] (add as input argument later)
	plotspec = int(sys.argv[5]) # 0 for no spectrogram, 1 for spectrogram
	end_time = int(start_time) + int(dur)
else:
	print "Usage: python scatMon.py ifo start_time dur thresh plospec"
	sys.exit(2)

### Set the channel that witnesses fringes (to plot specgram for) 
#witness_base = "LSC-DARM_IN1_DQ"
#witness_base = "GDS-CALIB_STRAIN"
#witness_base = "LSC-MICH_IN1_DQ"
#witness_base = '%s:ASC-Y_TR_A_NSUM_OUT_DQ' % ifo
#witness_base = 'LSC-SRCL_IN1_DQ'
#witness_base = "LSC-REFL_A_RF9_Q_ERR_DQ"
#witness_base = "ASC-AS_A_RF45_Q_PIT_OUT_DQ" 
#witness_base = "ASC-AS_B_RF36_Q_PIT_OUT_DQ"
#witness_base = "OMC-LSC_SERVO_OUT_DQ"
witness_base = "ASC-Y_TR_B_PIT_OUT_DQ"

if plotspec==1:
	witness_chan = ifo + ':' + witness_base
	print witness_chan
	# load timeseries for witness channel
        witness=TimeSeries.fetch(witness_chan, start_time, start_time+dur, verbose=True)
	if witness_base=="GDS-CALIB_STRAIN":
		witness=witness.highpass(20,gpass=3) # highpass the witness data
        elif witness_base=="LSC-DARM_IN1_DQ" or witness_base=="ASC-AS_B_RF45_I_PIT_OUT_DQ" or witness_base=="ASC-AS_B_RF36_Q_PIT_OUT_DQ" or witness_base=="LSC-MICH_IN1_DQ" or witness_base=="ASC-AS_A_RF45_Q_PIT_OUT_DQ":
		witness=witness.highpass(15,gpass=3) # highpass the witness data
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
	position=TimeSeries.fetch(position_chan, start_time, start_time+dur, verbose=True)
	position=position.detrend()
	position=position.lowpass(5,gpass=3) # lowpass position channel
	velocity = (position[1:]-position[:-1])*(position.sample_rate.value)
	times = (arange(len(velocity.value))+0.5)/(position.sample_rate.value)
	fudge = 1.0
	scatf1 = abs(fudge*2.0*velocity.value/1.064) # single bounce scatter frequency Virgo Scatter eqn 3
	scatf2 = 2.0*scatf1
	scatf3 = 3.0*scatf1
	scatf4 = 4.0*scatf1
	scatf5 = 5.0*scatf1
	# print max scatter values and their times
	print 'max scatter f1 = ' + str(max(scatf1)) + ' Hz'
	tofmax = times[argmax(scatf2)]
	tofmaxgps = tofmax + start_time
	print 'time of max f2 = ' + str(tofmax) + ' s, GPS=' + str(tofmaxgps)

	fig = plt.figure(figsize=(12,12))
	ax1 = fig.add_subplot(211)
	# Plot Spectrogram
	if plotspec==1:
        	im1 = NonUniformImage(ax1, interpolation='bilinear',extent=(min(t),max(t),10,55),cmap='jet')
        	im1.set_data(t,freq,20.0*log10(Pxx))
        	if witness_base=="GDS-CALIB_STRAIN":
			print "setting color limits for STRAIN"
			im1.set_clim(-1000,-800)
        	elif witness_base=="ASC-AS_B_RF45_Q_YAW_OUT_DQ" or witness_base=="ASC-AS_B_RF36_Q_PIT_OUT_DQ" or witness_base=="ASC-AS_A_RF45_Q_PIT_OUT_DQ" or witness_base=="LSC-MICH_IN1_DQ":
			im1.set_clim(-200,20)
		elif witness_base == "OMC-LSC_SERVO_OUT_DQ":
			im1.set_clim(-240,-85)
		ax1.images.append(im1)
        	#cbar1 = fig.colorbar(im1)
        	#cbar1.set_clim(-120,-40)
	
	# plot fringe prediction timeseries
	#ax1.plot(times,scatf5, c='blue', linewidth='0.2', label='f5')
	ax1.plot(times,scatf4, c='purple', linewidth='0.4', label='f4')
	ax1.plot(times,scatf3, c='green', linewidth='0.4', label='f3')
	ax1.plot(times,scatf2, c='blue', linewidth='0.4', label='f2')
	ax1.plot(times,scatf1, c='black', linewidth='0.4', label='f1')
	#if plotspec < 1 and dur <= 3600:
        #	ax1.plot(highscattimes, highscatf2, c='r', label='veto', marker='.',linestyle='')
	ax1.legend( loc='upper right' )
	if plotspec==1:
		title = 'Specgram of ' + string.replace(witness_chan,'_','\_') + '\nScattering fringe frequency for ' + string.replace(position_chan,'_','\_')
	else: 
       		title = 'Scattering fringe frequency for ' + string.replace(position_chan,'_','\_')
	# plt.title(title,fontsize=16)
        plt.yscale('linear')
	plt.ylabel('Frequency [Hz]',fontsize=16)
	xlab = 'Time [seconds] from ' + str(startutc) + ' UTC'
        plt.xlim(0,max(times))
        plt.ylim(0,50)
	ax2 = fig.add_subplot(212)
	#print(len(position.data))
	#print(position.data)
	ax2.plot(times,position[1:])
	plt.xlabel(xlab,fontsize=16)
	plt.xlim(0,max(times))
	plt.ylabel('Position [um]',fontsize=16)

	# save figure
	dir = str(start_time) + '-' + str(dur) + '/'
	if not path.exists(dir):
		makedirs(dir)
	filename = dir + ifo + '-Fscatter-' + channel + '-' + str(start_time) + '-' + str(dur) + '.png'
	savefig(filename)
	plt.close(fig)

quit()

############### Below here is old code for making DQ flags, needs to be cleaned up

if plotspec < 1:
	# determine times when f2 above some threshold
	idxhighscat = argwhere(scatf2>=thresh)
	highscatf2 = scatf2[idxhighscat]
	highscattimes = times[idxhighscat]
	highscattimesgps = highscattimes+start_time

	# save text file with values above threshold [GPS f2 index]
	outdata = hstack((highscattimesgps,highscatf2,idxhighscat))
	savetxt('%s-ALL-TIMES-SCATTER-GT%dHZ-%d-%d.txt' % (ifo,thresh,start_time,dur),outdata,fmt='%f %f %i')

	# save segments XML file with segments (based off code from Duncan Macleod)
	from math import (floor, ceil)
	from gwpy.segments import (Segment, DataQualityFlag)

	flag = '%s:DCH-SCATTERED_LIGHT_GT%dHZ:1' % (ifo,thresh) 
	flag = DataQualityFlag(flag)
	segs = []
	append = segs.append

	for gps in highscattimesgps:
    		if len(segs) and gps in segs[-1]:
        		continue
    		seg = Segment(floor(gps), ceil(gps))
    		append(seg)

	flag.active = segs
	flag.known = [Segment(start_time, end_time)]
	flag.coalesce()
	flag.write('%s-%s_%d-%d-%d.xml.gz' % (flag.ifo, flag.tag.replace('-', '_'), flag.version,start_time, dur))

#EOF
