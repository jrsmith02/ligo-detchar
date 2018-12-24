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

from gwpy.timeseries import TimeSeries
from gwpy.time import from_gps

# 1) Load the segment LOCK ALS ARMS

# Will do this later. 

# 2) load raw data for H1:ALS-C_TRX_A_LF_OUT_DQ and TRY
ifo = 'H1'
xchan = '{ifo}:ALS-C_TRX_A_LF_OUT_DQ'
start=1228730598
end=1228731098
xts = TimeSeries.get(xchan, start, end,verbose=True, nproc=args.nproc)

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
