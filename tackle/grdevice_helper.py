from functools import partial
from rpy2.robjects.packages import importr

grdevices = importr("grDevices")
gr_devices = {
   "png": grdevices.png,
   "pdf": grdevices.pdf,
   "svg": grdevices.svg,
}

def get_grid_kws(filetype, width=7, height=7, units="in", res=300):
    gr_kws = {
       "png": dict(width=width, height=height, units="in", res=300),
       "pdf": dict(
           width=width,
           height=height,
       ),
       "svg": dict(
           width=width,
           height=height,
       ),
    }
    if filetype not in gr_kws:
        raise ValueError(f"Invalid filetype: {filetype}")
    return gr_kws[filetype]


def get_device(filetype="pdf", **kwargs):
    filetype = filetype.strip(".")
    if filetype not in gr_devices:
        raise ValueError(f"Invalid filetype: {filetype}")
    gr_device = gr_devices[filetype]
    gr_kwargs = get_grid_kws(filetype, **kwargs)
    return partial(gr_device, **gr_kwargs)




