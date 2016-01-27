from gwpy.timeseries import TimeSeries
from gwpy.time import from_gps
import matplotlib.mlab as mlab
from pylab import savefig, specgram
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from os import makedirs, path
from numpy import sqrt, array, arange, amax, amin, nonzero, argmax, log10, loadtxt
from scipy.signal import iirfilter

# Set input variables
start_time = 1133737217 # GPS start time
startutc = str(from_gps(start_time)) # get UTC time
dur = 59*60
filter_pad = 3 # Extra data at start to avoid filter transients
ifo = 'L1'
#darmbase = ':OMC-DCPD_SUM_OUT_DQ' 
darmbase = ':GDS-CALIB_STRAIN'
darmchan = ifo + darmbase
rangechan = ifo + ':DMT-SNSH_EFFECTIVE_RANGE_MPC.mean'
fnm = 'FMCLVE_channels.txt'

# load channel names from text file
# Note: to generate this file, ran: nds_query -l -n nds.ligo-la.caltech.edu -t raw LVE-* > LVE_channels.txt
lines = loadtxt(fnm,
                   dtype='str',
                   usecols=[0],skiprows=2)

# get BNS range data
#if dur>=60*60:
#else:
range=TimeSeries.fetch(rangechan, start_time, start_time+dur, verbose=True,host='nds.ligo.caltech.edu')
trange = arange(0.0,len(range.value))/range.sample_rate.value

# Get DARM data
darm=TimeSeries.fetch(darmchan, start_time-filter_pad, start_time+dur, verbose=True,host='nds.ligo.caltech.edu')
darm=darm.detrend()
darm=darm.highpass(20)
# build notch filter for 60Hz line
Norder=2
zpk = iirfilter(Norder, [59.8/(darm.sample_rate.value/2.0), 60.2/(darm.sample_rate.value/2.0)], btype='bandstop', ftype='butter',  output='zpk')
darm = darm.filter(*zpk)

# Calculate DARM BLRMS 
stride=30 # seconds
flower=75 # Hz
fupper=95 # Hz
dbpdata = darm.bandpass(flower,fupper)
dbpdata = dbpdata[int(filter_pad*darm.sample_rate.value):]
dbrms = dbpdata.rms(stride)
#dbrms = dbrms.detrend()
dtbrms = arange(0.0,len(dbrms))*stride

# create scaled versions of data to compare to each other
srange=range.detrend()
rangemax = amax(srange.value)
rangemin = amin(srange.value)
sdbrms=dbrms.detrend()
dbrmsmax=amax(sdbrms.value)
sdbrms=sdbrms*(-rangemax/dbrmsmax) # scale BRMS to occupy same range as Range


for line in lines:
    	motionchan = line.replace('b','')
    	print(motionchan)
	basechan=motionchan.replace(':','-')
	try:
		motion=TimeSeries.fetch(motionchan, start_time-filter_pad, start_time+dur, verbose=True,host='nds.ligo-la.caltech.edu')
		pass
	except:
		continue
	#motion=motion.detrend()
	#motion=motion.lowpass(0.5)
	#motion = motion[int(filter_pad*motion.sample_rate.value):]
	tmotion = arange(0.0,len(motion.value))/motion.sample_rate.value
	
	# Calculate motion BLRMS 
	#stride=0.25 # seconds
	#flower=75 # Hz
	#fupper=95 # Hz
	#bpdata = motion.bandpass(flower,fupper)
	#bpdata = bpdata[int(filter_pad*motion.sample_rate.value):]
	#brms = bpdata.rms(stride)
	#brms = brms.detrend()
	#tbrms = arange(0.0,len(brms))*stride

	fig = plt.figure(figsize=(12,12))
	figtitle = "BLRMS: " + str(flower) + " -- " + str(fupper) + " Hz, " + str(stride) +  "s stride"
	ax2 = fig.add_subplot(311)
	ax2.plot(dtbrms,dbrms.value,label=darmchan)
	plt.legend()
	plt.ylabel("BLRMS",fontsize=16)
	plt.title(figtitle)
	plt.xlim(0,dur)
	ax3 = fig.add_subplot(312)
	ax3.plot(tmotion,motion.value,label=motionchan)
	plt.legend()
	plt.ylabel("Counts",fontsize=16)
	xlab = 'Time [seconds] from ' + str(startutc) + ' UTC'
	plt.xlim(0,dur)
	ax4 = fig.add_subplot(313)
	ax4.plot(trange,range.value,label=rangechan)
	plt.legend()
	plt.xlabel(xlab,fontsize=16)
	plt.ylabel('BNS Range [Mpc]')
	plt.xlim(0,dur)	
	# save figure
	dir = 'RangeCorrelation' + str(start_time) + '-' + str(dur) + '/'
	if not path.exists(dir):
		makedirs(dir)
	filename = dir + ifo + basechan + str(start_time) + '-' + str(dur) + '.png'
	savefig(filename)
	plt.close(fig)
	
	# make combined plot
	smotion=motion.detrend()
	motionmax=amax(smotion.value)
	if motionmax==0: # don't scale or plot if max for channel is identically zero
		continue
	else:
		smotion=smotion*(rangemax/motionmax) # scale motion channel to match range

	fig = plt.figure(figsize=(12,12))
	ax1 = fig.add_subplot(111)
	ax1.plot(dtbrms,sdbrms.value,label=darmchan)
	ax1.plot(tmotion,smotion.value,label=motionchan)
	ax1.plot(trange,srange.value,label=rangechan)
	plt.legend()
	plt.xlabel(xlab,fontsize=16)
	filename = dir + ifo + basechan + '-combined-' + str(start_time) + '-' + str(dur) + '.png'
        savefig(filename)
        plt.close(fig)
