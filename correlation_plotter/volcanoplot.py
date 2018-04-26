import os

import numpy as np
import pandas as pd
from .utils import get_outname

def volcanoplot(ctx, foldchange, number, scale):
    data_obj = ctx.obj['data_obj']

    group = data_obj.group #
    if group is None:
        print('Must supply a group value.')
        return

    if data_obj.col_metadata.loc[group].nunique() != 2:
        print('Error in volcanoplot, number of groups must be exactly 2.')
        return
    if data_obj.col_metadata.loc[group].value_counts().min() < 3:
        print('Each group must have at least 3 replicates.')
        return

    groups = dict()
    for grp in data_obj.col_metadata.loc[group].unique():
        samples = ((data_obj.col_metadata.loc[group] == grp)
                   .where(lambda x: x)
                   .dropna()
                   .index
        )
        groups[grp] = samples

    group0, group1 = data_obj.col_metadata.loc[group].values[[0, -1]]
    samples0, samples1 = groups[group0], groups[group1]
    # print(group0, group1)
    # print(samples0, samples1)


    values = data_obj.areas_log_shifted

    qvals   = data_obj.qvalues

    log2_fc = ((values[samples1].mean(1) - values[samples0].mean(1))
               .apply(lambda x: np.power(10, x))
               # .pipe(np.power, 10)
               .pipe(np.log2)
    )
    log2_fc.name = 'log2_Fold_Change'

    df = qvals.join(log2_fc.to_frame())
    df['GeneSymbol'] = df.index.map(lambda x: data_obj.gid_symbol.get(x, '?'))
    df['FunCats']    = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ''))
    df.index.name = 'GeneID'

    outname = get_outname('volcanoplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          batch=data_obj.batch_applied,
                          batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath)

    out = outname + '.tab'
    print("Saving", out, '...', end='', flush=True)
    df.to_csv(out, sep='\t')
    print('done', flush=True)

    try:
        from rpy2.robjects import r
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        _viaR = True
    except ModuleNotFoundError:
        _viaR = False

    if _viaR:

        pandas2ri.activate()
        r_source = robjects.r['source']
        r_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                              'R', 'volcanoplot.R')
        r_source(r_file)
        Rvolcanoplot = robjects.r['volcanoplot']

        file_fmts = ctx.obj['file_fmts']
        grdevices = importr('grDevices')
        gr_devices = {'.png': grdevices.png,
                      '.pdf': grdevices.pdf,
                      '.svg': grdevices.svg}
        gr_kws = {'.png': dict(width=5, height=5, units='in', res=300),
                  '.pdf': dict(width=5, height=5,),
                  '.svg': dict(width=5, height=5,)
        }
        for file_fmt in file_fmts:

            grdevice = gr_devices[file_fmt]
            gr_kw = gr_kws[file_fmt]
            out = outname + file_fmt
            print("Saving", out, '...', end='', flush=True)

            grdevice(file=out, **gr_kw)

            Rvolcanoplot(pandas2ri.py2ri(df.reset_index()), max_labels=number,
                         fc_cutoff=foldchange, label_cex=scale, group0=group0,
                         group1=group1)

            grdevices.dev_off()
            print('done.', flush=True)

    else:
        print('Must install rpy2')
