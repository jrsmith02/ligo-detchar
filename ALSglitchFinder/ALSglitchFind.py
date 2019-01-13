# This is first draft code to respond to the request for an 
# ALS glitch monitor in:
# https://alog.ligo-wa.caltech.edu/aLOG/index.php?callRep=45924
# Josh Smith and Joe Areeda

# 1) Loads the guardian segments for ALS XARM/YARM and ISC and 
# decides what times to include
# 2) loads raw data for H1:ALS-C_TRX_A_LF_OUT_DQ and TRY
# 3) Calculates the standard deviation of these timeseries in
# the correct guardian segments
# 4) Identifies any outliers that exceed N standard deviations,
# where N is user-specified
# 5) Plots both timeseries and highlights outliers in red

import time
start_of_run = time.time()  #nopep8

import numpy
import argparse
import os
import re
import logging
import warnings
from collections import OrderedDict
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt


from gwpy.timeseries import TimeSeries, TimeSeriesDict, StateTimeSeries
from gwpy.time import from_gps
from gwpy.plot import Plot
from matplotlib import use
use('agg')  # nopep8

from gwpy.segments import DataQualityFlag
from gwpy.segments import SegmentList
from gwpy.segments import Segment
from gwpy.time import to_gps
from matplotlib import rcParams

rcParams['text.usetex'] = False


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


def get_guardian_segs(gts, ts, iscts, wfsts):
    """Examine  a time series of Guardian states, returning a segment list
    of periods when the  green laser is on

    INPUT
    =====
    gts: TimeSeries of Guardian state
    ts: TimeSeries of data to examine.
        Used here only to resize flags appropriately

    RETURN
    ======
    segs:   SegmentList for Green laser on
    flags:  StateTimeSeries for us with input ts
    """

    f1 = -19 <= gts.value
    f2 = gts.value < 100
    flg = numpy.logical_and(f1, f2)
    flag = StateTimeSeries(flg, t0=gts.t0, dt=gts.dt,
                           name=gts.name+' green laser', channel=gts.channel)

    grn_segs = flag.to_dqflag()

    iscflag = iscts >= 12 * iscts.unit
    isc_flags =  StateTimeSeries(iscflag, t0=iscflag.t0, dt=iscflag.dt,
                                 name=iscts.name, channel=iscts.channel)
    isc_segs = isc_flags.to_dqflag()
    isc_segs.name = u'ISC Grdn > 12'

    cflg = numpy.logical_and(flg, iscflag.value)

    wfs_flg = wfsts.value > 0
    wfs_flags = StateTimeSeries(wfs_flg, t0=wfsts.t0, dt=wfsts.dt,
                                 name=wfsts.name, channel=wfsts.channel)
    wfs_segs = wfs_flags.to_dqflag()
    wfs_segs.name = u'WFS TRIG'

    cflg = numpy.logical_and(cflg, wfs_flags.value)

    combo_flag = StateTimeSeries(cflg, t0=iscflag.t0, dt=iscflag.dt,
                           name=gts.name+' green laser', channel=gts.channel)
    combo_segs = combo_flag.to_dqflag()

    combo_flag = resample_bool(combo_flag, iscflag, ts)

    return grn_segs, isc_segs, wfs_segs, combo_segs, combo_flag

def resample_bool(inflag, ints, outts):
    out_array = numpy.ndarray(len(outts), dtype=bool)
    fact = int(len(outts)/len(ints))
    for i in range(0, len(ints)):
        if inflag[i]:
            for o in range(i*fact, (i+1)*fact):
                out_array[o] = True

    outflag = StateTimeSeries(out_array, t0=outts.t0, dt=outts.dt)

    return outflag


def plotit(segs, isc_segs, wfs_segs, combo_segs, flag, ts, axis):
    """plot timeseries segments and outliers
    INPUT
    ======
    segs -  segments describing when the green laser is on for this axis
    combo_segs -  segments that combine segs  and ISC in a proper state
    flag -  a Boolean array corresponding to combo_segs
    ts -  the time series to plot and analyze during combo_seg
    axis  - 'X' or 'Y'  for labels
    """
    global ifo, gpsstub, N, args

    mean = ts[flag].mean()
    std = ts[flag].std()
    low_limit = mean - N * std

    chan = ts.channel.name
    lim = ts.span

    plot = None

    for seg in combo_segs.active:
        #  plot only the active segments as dots so we don't have lines through
        #  the unanalyzed segments
        seg = seg.protract(-1)

        if seg[1] - seg[0] > 1:
            partial_ts = ts.crop(seg[0], seg[1])

            outlier_ts = partial_ts.copy()
            omax = outlier_ts.max()

            cntr = (outlier_ts - mean).abs() <= N * std
            outlier_ts[cntr] = omax * 2
            out_count = len(partial_ts) - cntr.sum().value
            logger.info('Seg: {}, size: {:d}, outliers {:.0f}'.format(
                    seg, len(partial_ts), out_count))

            if plot == None:
                plot = partial_ts.plot(figsize=(18, 6), color='blue', marker='.',
                                        markersize=1, label=chan, linewidth=0)
                ax = plot.gca()

                plot2 = outlier_ts.plot(figsize=(18, 6), color='red', marker='.',
                                        markersize=4, label=chan, linewidth=0)
                ax2 = plot2.gca()
            else:
                ax.plot(partial_ts, color='blue', marker='.', markersize=1,
                         linewidth=0)

                ax.plot(outlier_ts,  color='red', markersize=1,
                         linewidth=0, marker='.')

                ax2.plot(outlier_ts, color='red', markersize=1,
                        linewidth=0, marker='.')

    ax2.set_ylim(0, ts.max().value * 1.1)
    ax2.set_xlim(lim[0], lim[1])
    ax2.set_title('Test {:s} and glitches'.format(chan))
    testfnane = os.path.join(args.outdir, '{:s}-Tst-{:s}-{:s}.png'.format(
            ifo, axis, gpsstub))
    plot2.add_segments_bar(combo_segs, label='Combined segs')
    plot2.add_segments_bar(wfs_segs, label='{:s}-WFS trg delay'.format(axis))

    plot2.add_segments_bar(segs, label='{:s}-grn-lsr'.format(axis))
    plot2.add_segments_bar(isc_segs, label=isc_segs.name)

    plot2.savefig(testfnane)
    logger.info('Wrote test image {:s}'.format(testfnane))

    ax.set_ylim(0, ts.max().value*1.1)
    ax.set_xlim(lim[0], lim[1])
    ax.set_title('{:s} and glitches'.format(chan))
    plot.add_segments_bar(combo_segs, label='Combined segs')
    plot.add_segments_bar(wfs_segs, label='{:s}-WFS trg delay'.format(axis))

    plot.add_segments_bar(segs, label='{:s}-grn-lsr'.format(axis))
    plot.add_segments_bar(isc_segs, label=isc_segs.name)

    out_filename = os.path.join(args.outdir, '{:s}-ALS-{:s}-{:s}.png'.format(
            ifo, axis, gpsstub))
    plot.savefig(out_filename)
    logger.info('Wrote {:s}'.format(out_filename))


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

logger.info('args: {}'.format(args))
xchan = '%s:ALS-C_TRX_A_LF_OUT_DQ' % ifo
ychan = '%s:ALS-C_TRY_A_LF_OUT_DQ' % ifo
iscgrd = '%s:GRD-ISC_LOCK_STATE_N' % ifo
xgrd = '%s:GRD-ALS_XARM_STATE_N' % ifo
ygrd = '%s:GRD-ALS_YARM_STATE_N' % ifo
xwfs = '%s:ALS-X_WFS_DOF_FM_TRIG_DELAYED' % ifo
ywfs = '%s:ALS-Y_WFS_DOF_FM_TRIG_DELAYED' % ifo

# get a list of ALL chanel names for efficient read/transfer
chan_list = [iscgrd, xgrd, ygrd, xchan, ychan, xwfs, ywfs]

# 1) Get all channels atonce for efficiency: LOCK ALS ARMS & Guardian states
tsd =  TimeSeriesDict.get(chan_list, start, end, verbose=(verbosity > 1),
                          nproc=nproc)
# move the ALs channels to convenience variables
xts = tsd[xchan]
yts = tsd[ychan]


xsegs, isc_segs, xwfs_segs, xcombo_segs, xcombo_flag = get_guardian_segs(tsd[xgrd],
                                            tsd[xchan], tsd[iscgrd], tsd[xwfs])
xsegs.name = 'X-green laser on'

ysegs, isc_segs, ywfs_segs, ycombo_segs, ycombo_flag = get_guardian_segs(tsd[ygrd],
                                            tsd[ychan], tsd[iscgrd], tsd[ywfs])
ysegs.name = 'Y-green laser on'

gpsstub = '%d-%d' % (start, end-start)

plotit(xsegs, isc_segs, xwfs_segs, xcombo_segs, xcombo_flag, xts, 'X')
plotit(ysegs, isc_segs, ywfs_segs, ycombo_segs, ycombo_flag, yts, 'Y')

run_time = time.time() - start_of_run
logger.info('Runtime: {:.1f} seconds'.format(run_time))

