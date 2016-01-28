Readme for correlate-range-channels.py

run using:  python correlate-range-channels.py

the script as written produces this output: https://ldas-jobs.ligo-la.caltech.edu/~jrsmith/detchar/O1/range_variations/RangeCorrelation1133737217-3540/

Those are the current time and durations hardcoded in the script. 

It makes two outputs per channel one like this with strain BLRM, channel, then range: 
https://ldas-jobs.ligo-la.caltech.edu/~jrsmith/detchar/O1/range_variations/RangeCorrelation1133737217-3540/L1L0-FMC-EX_AHU_FAN2_DISCH_RH1133737217-3540.png

one like this where it tries to scale them all to look for a match: 
https://ldas-jobs.ligo-la.caltech.edu/~jrsmith/detchar/O1/range_variations/RangeCorrelation1133737217-3540/L1L0-FMC-EX_AHU_FAN2_DISCH_TEMP-combined-1133737217-3540.png

This script is a purely by eye match. But it does find some interesting things from time to time. When the script is done running, I recommend tarring up the whole directory and downloading it and opening all images to quickly flip through them. 
