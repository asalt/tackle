from contextlib import contextmanager
from functools import partial


@contextmanager
def grdevice(filename, width=7, height=7, units="in", res=300):
    from pathlib import Path
    from rpy2.robjects.packages import importr
    grdevices = importr("grDevices")
    filetype = Path(filename).suffix.strip(".").lower()

    device = get_device(filetype, width=width, height=height, units=units, res=res)
    device(filename)  # Open the device

    try:
        yield
    finally:
        grdevices.dev_off()

def get_gr_devices():
    from rpy2.robjects.packages import importr

    grdevices = importr("grDevices")
    gr_devices = {
       "png": grdevices.png,
       "pdf": grdevices.cairo_pdf,
       "svg": grdevices.svg,
    }
    return gr_devices

def get_grid_kws(filetype, width=7, height=7, units="in", res=300):


    from rpy2.robjects.packages import importr

    grdevices = importr("grDevices")
    gr_devices = {
       "png": grdevices.png,
       "pdf": grdevices.pdf,
       "svg": grdevices.svg,
    }

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
    gr_devices = get_gr_devices()
    if filetype not in gr_devices:
        raise ValueError(f"Invalid filetype: {filetype}")
    gr_device = gr_devices[filetype]
    gr_kwargs = get_grid_kws(filetype, **kwargs)
    return partial(gr_device, **gr_kwargs)




