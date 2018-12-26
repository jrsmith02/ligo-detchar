# This is rough code to respond to the request for an ALS glitch finder in:
# https://alog.ligo-wa.caltech.edu/aLOG/index.php?callRep=45924
# Josh Smith and Joe Areeda

# Pseudo Code:
# 1) Load the segment LOCK ALS ARMS 
# 2) load raw data for H1:ALS-C_TRX_A_LF_OUT_DQ and TRY 
# 3) Calculate the timeseries standard deviation and identify the outliers
# that exceed N standard deviations, where N is user-specified
# 4) Plot both timeseries (blue and green) and then highlight what outliers were identified in red timeseries overlay fragments. 
# 5) Tune Nsigma until the correct glitches are being found without false alarms.)

import time
start_of_run = time.time()  #nopep8

import numpy
import argparse
import os
import logging
import warnings
from collections import OrderedDict
from scipy.interpolate import UnivariateSpline


from gwpy.timeseries import TimeSeries, TimeSeriesDict
from gwpy.time import from_gps
from gwpy.plot import Plot
from matplotlib import use
use('agg')  # nopep8

from gwpy.segments import DataQualityFlag
from gwpy.segments import SegmentList
from gwpy.segments import Segment
from gwpy.time import to_gps

__author__ = 'Joshua Smith'
__email__ = 'joshua.smith@ligo.org'
__version__ = '0.0.1'
prog = 'ALSglitchFind'


def save_figure(fig, pngfile, **kwargs):
    try:
        fig.save(pngfile, **kwargs)
    except (RuntimeError, IOError, IndexError):
        try:
            fig.save(pngfile, **kwargs)
        except (RuntimeError, IOError, IndexError) as e:
            warnings.warn('Error saving {0}: {1}'.format(pngfile, str(e)))
            return
    fig.close()
    logger.info('Wrote {:s}'.format(pngfile))
    return pngfile

def unique(list_):
    """Returns a unique version of the input list preserving order

    Examples
    --------
    >>> unique(['b', 'c', 'a', 'a', 'd', 'e', 'd', 'a'])
    ['b', 'c', 'a', 'd', 'e']
    """
    return list(OrderedDict.fromkeys(list_).keys())

# command line options
parser = argparse.ArgumentParser(description=__doc__, prog=prog)

parser.add_argument('--start', type=to_gps, action='store', default=1228730598,
                    help='Start time GPS or yyyy-mm-dd hh:mm:ss.mmm ')
parser.add_argument('--end', type=to_gps, action='store', default=1228731098,
                    help='End time GPS or yyyy-mm-dd hh:mm:ss.mmm ')
parser.add_argument('-j', '--nproc', type=int, default=4,
                    help='number of processes for reading')
parser.add_argument('-o', '--outdir', default='./')
parser.add_argument('-i', '--ifo', default='H1')
parser.add_argument('-N', '--sigma', type=float, default=2.0,
                    help='number of sigma to use for identifying outliers')

parser.add_argument('-v', '--verbose', action='count', default=0)

args = parser.parse_args()

verbosity = args.verbose
logging.basicConfig()
logger = logging.getLogger(prog)
logger.setLevel(logging.DEBUG)

if verbosity < 1:
    logger.setLevel(logging.CRITICAL)
elif verbosity < 2:
    logger.setLevel(logging.INFO)
else:
    logger.setLevel(logging.DEBUG)

logger.info('Args: {:s}'.format(args))
outdir = args.outdir
if not os.path.exists(outdir):
    os.makedirs(outdir)

start = args.start
end = args.end
ifo = args.ifo
nproc = args.nproc
N = args.sigma

xchan = '%s:ALS-C_TRX_A_LF_OUT_DQ' % ifo
ychan = '%s:ALS-C_TRY_A_LF_OUT_DQ' % ifo
xgrd = '%s:GRD-ALS_XARM_STATE_N' % ifo
ygrd = '%s:GRD-ALS_YARM_STATE_N' % ifo

# 1) Load the segment LOCK ALS ARMS
grd_tsd =  TimeSeriesDict.get([xgrd, ygrd], start, end, verbose=(verbosity>1),
                            nproc=nproc)

xflag = grd_tsd[xgrd] > 20 * grd_tsd[xgrd].unit
xsegs = xflag.to_dqflag()
xsegs.name = 'X-Guardian state > 20'

yflag = grd_tsd[ygrd] > 20 * grd_tsd[ygrd].unit
ysegs = yflag.to_dqflag()
ysegs.name = 'Y-Guardian state > 20'

# 2) load raw data for H1:ALS-C_TRX_A_LF_OUT_DQ and TRY

gpsstub = '%d-%d' % (start, end-start)
tsd = TimeSeriesDict.get([xchan, ychan], start, end, verbose=(verbosity>1),
                         nproc=nproc)
xts = tsd[xchan]
yts = tsd[ychan]

# 3) Calculate the timeseries standard deviation and identify the outliers
# that exceed N standard deviations, where N is user-specified
def find_outliers(ts, N):
    ts = ts.value  # strip out Quantity extras
    return numpy.nonzero(abs(ts - numpy.mean(ts)) > N*numpy.std(ts))[0]

xoutliers = find_outliers(xts, N)
logger.info('x outliers: {}'.format(xoutliers[1:5]))
xglitches=xts[xoutliers]

youtliers = find_outliers(yts, N)
logger.info('x outliers: {}'.format(youtliers[1:5]))
yglitches=xts[youtliers]
#print xoutliers[1:10]
xglitchtimes=xts[xoutliers].times.value
yglitchtimes=yts[youtliers].times.value
#print 'Glitches: ', [int(x) for x in xglitchtimes[1:5]]

# 4) Plot both timeseries (blue and green) and then
#    highlight what outliers were identified in
#    red timeseries overlay fragments.

times = xts.times.value
xlim = xts.span
plot = xts.plot(figsize=(12,6), color='blue', label=xchan.replace('_', r'\_'),
        linewidth=0.5)
ax = plot.gca()
ax.plot(yts, color='green', label=ychan.replace('_', r'\_'),
        linewidth=0.5)
ax.plot(xglitches, color='red',marker=".",
        label='X-glitches', linewidth=0.5)
ax.plot(yglitches.value, color='magenta',marker=".",
        label='Y-glitches', linewidth=0.5)
ax.set_ylabel('Transmitted power [unknown]')
ax.legend(loc='best')
plot.add_segments_bar(xsegs, label='Guardian-X')
plot.add_segments_bar(ysegs, label='Guardian-Y')

out_filename = os.path.join(args.outdir, '{:s}-ALSts-{:s}.png'.format(ifo,gpsstub))
plot1 = save_figure(plot, out_filename)

run_time = time.time() - start_of_run
logger.info('Runtime: {:.1f} seconds'.format(run_time))


# Function definitions


def remove_outliers(ts, N):
    outliers = find_outliers(ts, N)
    c = 1
    if outliers.any():
        print("-- Found %d outliers in %s, recursively removing"
              % (len(outliers), ts.name))
        while outliers.any():
            unit = ts.unit
            cache = outliers
            mask = numpy.ones(len(ts), dtype=bool)
            mask[outliers] = False
            spline = UnivariateSpline(ts[mask].times.value, ts[mask].value,
                                      s=0, k=3)
            ts[outliers] = (spline(ts[outliers].times.value) * unit)
            outliers = find_outliers(ts, N)
            print("Completed %d removal cycles" % c)
            if numpy.array_equal(outliers, cache):
                print("Outliers did not change, breaking recursion")
                break
            print("%d outliers remain" % len(outliers))
            c += 1



# EOF
