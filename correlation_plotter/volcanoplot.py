import os
import itertools

import numpy as np
import pandas as pd
from .utils import get_outname, parse_gid_file

def volcanoplot(ctx, foldchange, expression_data, number, only_sig=False, sig=.05,
                yaxis='pAdj', scale=1.2, highlight_geneids=None):

    data_obj = ctx.obj['data_obj']

    gids_to_highlight = None
    if highlight_geneids is not None:
        gids_to_highlight = parse_gid_file(highlight_geneids)


    if yaxis not in ('pValue', 'pAdj'):
        raise ValueError('Must choose between `pValue` and `pAdj`')

    group = data_obj.group #
    if group is None:
        print('Must supply a group value.')
        return

    if data_obj.col_metadata.loc[group].nunique() < 2:
        print('Error in volcanoplot, number of groups must be at least 2.')
        return
    # if data_obj.col_metadata.loc[group].nunique() != 2:
    #     print('Error in volcanoplot, number of groups must be exactly 2.')
    #     return
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

    values = data_obj.areas_log_shifted
    padj   = data_obj.padj
    from .utils import filter_observations
    filtering = values.copy()
    filtering.index = [filtering.index, ['area']*len(filtering)]
    values_filtered = filter_observations(filtering, 'area',
                                          data_obj.non_zeros, data_obj.nonzero_subgroup,
                                          data_obj.col_metadata)
    values = values_filtered.reset_index(level=1, drop=True)


    for group in itertools.combinations(data_obj.col_metadata.loc[group].unique(), 2):
        group0, group1 = group

        # group0, group1 = data_obj.col_metadata.loc[group].unique()[[0, -1]]
        samples0, samples1 = groups[group0], groups[group1]

        log2_fc = ((values[samples1].mean(1) - values[samples0].mean(1))
                .apply(lambda x: np.power(10, x))
                # .pipe(np.power, 10)
                .pipe(np.log2)
        )
        log2_fc.name = 'log2_Fold_Change'

        df = padj.join(log2_fc.to_frame())
        df['GeneSymbol'] = df.index.map(lambda x: data_obj.gid_symbol.get(x, '?'))
        df['FunCats']    = df.index.map(lambda x: data_obj.gid_funcat_mapping.get(x, ''))
        df.index.name = 'GeneID'
        df['highlight'] = False
        if gids_to_highlight is not None:
            df.loc[gids_to_highlight, 'highlight'] = True


        outname = get_outname('volcanoplot', name=data_obj.outpath_name, taxon=data_obj.taxon,
                              non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                              batch=data_obj.batch_applied,
                              batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                              outpath=data_obj.outpath,
                              group='{}_vs_{}'.format(group0, group1),
        )

        out = outname + '.tab'
        print("Saving", out, '...', end='', flush=True)
        export_data = df
        if only_sig:
            export_data = df.query('qValue < @sig')
        if expression_data:
            export_data = export_data.join(values)
        export_cols = [x for x in export_data.columns if x not in ('highlight', )]
        export_data[export_cols].to_csv(out, sep='\t')
        print('done', flush=True)

        try:
            from rpy2.robjects import r
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            _viaR = True
        except ModuleNotFoundError:
            _viaR = False

        if not _viaR:
            print('Must install rpy2')
            return

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

            Rvolcanoplot(pandas2ri.py2ri(df.reset_index()), max_labels=number, fc_cutoff=foldchange,
                         sig=sig, yaxis=yaxis, label_cex=scale, group0=group0, group1=group1)

            grdevices.dev_off()
            print('done.', flush=True)
