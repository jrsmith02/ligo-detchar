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

import numpy 

from gwpy.timeseries import TimeSeries
from gwpy.time import from_gps
from gwpy.plot import Plot
from matplotlib import use
use('agg')  # nopep8

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
    return pngfile

# 1) Load the segment LOCK ALS ARMS

# Will do this later. 

# 2) load raw data for H1:ALS-C_TRX_A_LF_OUT_DQ and TRY
ifo = 'H1'
xchan = '%s:ALS-C_TRX_A_LF_OUT_DQ' % ifo
ychan = '%s:ALS-C_TRY_A_LF_OUT_DQ' % ifo
start=1228730598
end=1228731098
nproc=4
gpsstub = '%d-%d' % (start, end-start)
xts = TimeSeries.get(xchan, start, end,verbose=True, nproc=nproc)
yts = TimeSeries.get(ychan, start, end,verbose=True, nproc=nproc)

# 4) Plot both timeseries (blue and green) and then highlight what outliers were identified in red timeseries overlay fragments.

times = xts.times.value
xlim = xts.span
plot = Plot(figsize=(12,6))
ax = plot.gca(xscale='auto-gps', epoch=start, xlim=xlim)
ax.plot(times, xts.value, color='blue', label=xchan.replace('_', r'\_'), linewidth=0.5)
ax.plot(times, yts.value, color='green', label=ychan.replace('_', r'\_'), linewidth=0.5)
ax.legend(loc='best')
plot1 = save_figure(plot, '%s-ALSts-%s.png' % (ifo,gpsstub))

# Function definitions

def find_outliers(ts, N):
    ts = ts.value  # strip out Quantity extras
    return numpy.nonzero(abs(ts - numpy.mean(ts)) > N*numpy.std(ts))[0]


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
