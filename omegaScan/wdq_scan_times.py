# Scan some times using wdq
# Reads a text file with a column of times
# Sets the condor accounting_group correctly for omega scans
# nohup condor run's a wdqER omega scan for each time

import os
import numpy as np

ifo='H1'
filen = 'nikhil_tomte.txt';
times = np.loadtxt(filen,float)
#print times

wdq = '/home/detchar/bin/wdqER'
options = '-a accounting_group=ligo.prod.o1.detchar.user_req.omegascan'
outdir = '/home/jrsmith/public_html/omegaScans/tomtes/'

for time in times:
	stime =  "%f" % (time) # convert times to a string
	command = "nohup condor_run " + wdq + " " + options + " " + ifo \
	+ " " + stime + " " + outdir + stime + "/" \
	+ " < /dev/null &>/dev/null &"

	print "Running wdq for ifo=" + ifo + " GPS=" + stime
	#print command
	os.system(command)
