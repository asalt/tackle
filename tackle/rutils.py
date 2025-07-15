# rutils.py
import warnings

from rpy2.robjects import r
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
grdevices = importr("grDevices")
#pandas2ri.activate()
r_source = robjects.r["source"]
gr_devices = {
    ".png": grdevices.png,
    ".pdf": grdevices.pdf,
    ".svg": grdevices.svg,
}



def close_grdevice():
    grdevices.dev_off()
