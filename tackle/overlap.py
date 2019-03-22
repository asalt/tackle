import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sb
import click

from pyupset.visualisation import DataExtractor

from .utils import get_outname, save_multiple, genefilter, filter_sra, filter_taxon, filter_observations
from .upset import make_plot as make_upset

idx = pd.IndexSlice

from itertools import combinations

def calc_combos(tot):
    return sum(len(x) for x in [list(combinations('a'*tot, i)) for i in np.arange(1, tot+1)])

def make_overlap(data_obj, group=None, file_fmts=('.png',), non_zeros=1., maxsize=15):

    outname = get_outname('overlap', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=data_obj.non_zeros, colors_only=data_obj.colors_only,
                          # batch=data_obj.batch_applied,
                          # batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath,
                          group=group
                          # **kws
    )


    overlap_dict = dict()
    if group:
        cols = data_obj.col_metadata.T.groupby(group).groups
    else:
        cols = {col: [col] for col in data_obj.data.columns}

    # ensure not more combinations than reasonable to calculate (12 experiments/groups)
    num_combos = calc_combos(min(len(cols), 20))
    if num_combos > 4095:
        if len(cols) > 20:
            s = 'more than {} combinations to calculate.'.format(num_combos)
        else:
            s = '{} combinations to calculate.'.format(num_combos)
        stdout = 'Too many experiments. With {} experiments, {}. Consider using --groups'.format(len(cols), s)
        click.echo(stdout)
        click.echo('Skipping...')
        return


    for name, col in cols.items():

        res = (filter_observations(data_obj.data[col], 'area', nonzero_value=non_zeros)
               .pipe(filter_sra)
               .loc[ idx[:, 'SRA'], :]
               .reset_index()
               .rename(columns={'level_0':'GeneID', 'level_1': 'SRA'})
        )


        # res = (data_obj.data.loc[ idx[:, 'SRA'], col ]
        #        .where(lambda x: x == 'S')
        #        .dropna()
        #        .reset_index()
        #        .rename(columns={'level_0':'GeneID', 'level_1': 'SRA'})
        # )

        overlap_dict[name] = res

    de = DataExtractor(overlap_dict, 'GeneID' )

    # ensure overlap is small enough to be reasonable for display
    maxiter, bound_min, size = 1e3, 0, np.inf
    iiter = 0
    while size > maxsize:  # no more than that! default 15
        iiter += 1
        bound_min += 1
        filtered_intersections = de.get_filtered_intersections(sort_by='size',
                                                               inters_size_bounds=(bound_min, np.inf),
                                                               inters_degree_bounds=(0, np.inf) )
        size = min(size, len(filtered_intersections[0]))
        # print(size, len(filtered_intersections[0]), min(filtered_intersections[0]))
        if iiter > maxiter:
            click.echo('Could not find small enough overlap...Skipping')
            return


    pyupset_res = make_upset(overlap_dict, unique_keys=('GeneID',), h_ratio=2, v_ratio=2.25,
                             dot_size=180, figsize=(12,10.5), bound_min=bound_min)
    fig = pyupset_res['figure']

    save_multiple(fig, outname, *file_fmts)


    # export gene membership

    membership_df = (pd.DataFrame.from_dict(
        {'|'.join(membership) : de.inters_df_dict[ membership ]['GeneID']
         for membership in filtered_intersections[1]
        },
        orient='columns'
    )
                     .melt(var_name='Sample_Set', value_name='GeneID')
                     .dropna()
    )

    membership_df['GeneID'] = membership_df['GeneID'].astype(int)

    print('Saving {}+.tsv')
    membership_df.to_csv(outname+'.tsv', index=False, sep='\t')
