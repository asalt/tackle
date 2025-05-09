import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sb
import click

# from pyupset.visualisation import DataExtractor

from .utils import get_outname, save_multiple, genefilter, filter_sra, filter_taxon, filter_observations
#from .upset import make_plot as make_upset

idx = pd.IndexSlice

from itertools import combinations

def calc_combos(tot):
    return sum(len(x) for x in [list(combinations('a'*tot, i)) for i in np.arange(1, tot+1)])

def make_overlap(data_obj, group=None, file_fmts=('.png',), non_zeros=1., maxsize=15, figsize=None):

    if figsize is None:
        figsize=(12,10.5)

    outname = get_outname('overlap', name=data_obj.outpath_name, taxon=data_obj.taxon,
                          non_zeros=non_zeros, colors_only=data_obj.colors_only,
                          # batch=data_obj.batch_applied,
                          # batch_method = 'parametric' if not data_obj.batch_nonparametric else 'nonparametric',
                          outpath=data_obj.outpath,
                          group=group
                          # **kws
    )


    overlap_dict = dict()
    if group:
        cols = data_obj.col_metadata.groupby(group).groups
    else:
        # cols = {col: [col] for col in data_obj.data.columns if col }
        cols = {col: [col] for col in data_obj.col_metadata.index}

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

        # res = (filter_observations(data_obj.data[col], 'area', nonzero_value=non_zeros)

        # res = (filter_observations(data_obj.data[['GeneID', 'Metric', *col]], 'area', nonzero_value=non_zeros)
        res = (filter_observations(data_obj.data[['GeneID', 'Metric', *col]], 'area', nonzero_value=non_zeros)
               .pipe(filter_sra)
               # .loc[ idx[:, 'SRA'], :]
               .query('Metric=="SRA"')
               # .reset_index()
               .rename(columns={'level_0':'GeneID', 'level_1': 'SRA'})
        )


        # res = (data_obj.data.loc[ idx[:, 'SRA'], col ]
        #        .where(lambda x: x == 'S')
        #        .dropna()
        #        .reset_index()
        #        .rename(columns={'level_0':'GeneID', 'level_1': 'SRA'})
        # )

        overlap_dict[name] = res.GeneID.tolist()


    from rpy2.robjects import r
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()

    grdevices = importr("grDevices")
    grid = importr("grid")
    complex_heatmap = importr('ComplexHeatmap')
    robjects.r('library(ComplexHeatmap)')
    robjects.r('library(grid)')


    robjects.r.assign('listing', robjects.ListVector(overlap_dict))
    robjects.r('cmat <- make_comb_mat(listing)')
    robjects.r.assign('maxsize', maxsize)
    robjects.r("""
      if (ncol(cmat) > maxsize){
        cmat <- cmat[, 1:maxsize]
      }
    """)
    cmat = robjects.r['cmat']



    upset_plot = robjects.r("""UpSet(
        cmat,
        top_annotation=upset_top_annotation(
            cmat,
            height=unit(2, "in"),
            annotation_name_rot=90
        ),
        right_annotation=upset_right_annotation(
            cmat,
            width=unit(1.5, "in")
        )
    )""")


    width, height = cmat.shape
    max_row_name = len(sorted(overlap_dict.keys(), key = lambda x: len(x))[-1])
    figwidth = 2 + max(width, 2)*.4 + (max_row_name*.08)
    figheight = 2 + height*.24
    print(f"widthxheight: {figwidth}x{figheight}")

    grdevices.png(outname+'.png', width=figwidth, height=figheight, units='in', res=300)
    complex_heatmap.draw(upset_plot)
    grdevices.dev_off()
    return # temporary



    # de = DataExtractor(overlap_dict, 'GeneID' )

    # # ensure overlap is small enough to be reasonable for display
    # maxiter, bound_min, size = 1e3, 0, np.inf
    # iiter = 0
    # while size > maxsize:  # no more than that! default 15
    #     iiter += 1
    #     bound_min += 1
    #     filtered_intersections = de.get_filtered_intersections(sort_by='size',
    #                                                            inters_size_bounds=(bound_min, np.inf),
    #                                                            inters_degree_bounds=(0, np.inf) )
    #     size = min(size, len(filtered_intersections[0]))
    #     # print(size, len(filtered_intersections[0]), min(filtered_intersections[0]))
    #     if iiter > maxiter:
    #         click.echo('Could not find small enough overlap...Skipping')
    #         return

    # pyupset_res = make_upset(overlap_dict, unique_keys=('GeneID',), h_ratio=2, v_ratio=2.25,
    #                          dot_size=180, figsize=figsize, bound_min=bound_min)
    # fig = pyupset_res['figure']

    # save_multiple(fig, outname, *file_fmts)


    # # export gene membership

    # all_intersections = de.get_filtered_intersections(sort_by='size',
    #                                                   inters_size_bounds=(0, np.inf),
    #                                                   inters_degree_bounds=(0, np.inf) )


    # membership_df = (pd.DataFrame.from_dict(
    #     {'|'.join(membership) : de.inters_df_dict[ membership ]['GeneID']
    #      for membership in all_intersections[1]
    #     },
    #     orient='columns'
    # )
    #                  .melt(var_name='Sample_Set', value_name='GeneID')
    #                  .dropna()
    # )

    # membership_df['GeneID'] = membership_df['GeneID'].astype(int)

    # print('Saving {}+.tsv'.format(outname))
    # membership_df.to_csv(outname+'.tsv', index=False, sep='\t')
