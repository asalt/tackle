# rutils.py
import warnings

_rpy2 = False
try:
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    _rpy2 = True
except ModuleNotFoundError:
    warnings.warn("rpy2 needs to be installed")

if _rpy2 is True:
    grdevices = importr("grDevices")
    pandas2ri.activate()
    r_source = robjects.r["source"]
    gr_devices = {
        ".png": grdevices.png,
        ".pdf": grdevices.pdf,
        ".svg": grdevices.svg,
    }


def close_grdevice():
    grdevices.dev_off()
