import os

from .containers import GeneMapper
from .utils import hgene_map, named_temp

import matplotlib.pyplot as plt
from rpy2 import robjects

def ssGSEA(expression, metadata,
           geneset='hallmark',
           norm='rank',
           combine_mode='combine.off',
           statistic='area.under.RES',
           weight=0.75,
           correl='z_score',
           nperm=1000,
           min_overlap=10,
           extended_output=True,
           globalfdr=False,
           output_score='NES',
           output_prefix='run',
           log_file='run.log',
           rank_plots=False,
           seed=1234
):


    plt.rc('font',**{'family':'sans-serif','sans-serif':["DejaVu Sans", "Arial", "Liberation Sans",
                            "Bitstream Vera Sans", "sans-serif"]})


    geneset_mapping = {'hallmark': 'h.all.v7.0.entrez.gmt',
                        'go_biological': 'c5.bp.v7.0.entrez.gmt',
                       'curated.all': 'c2.all.v7.0.entrez.gmt',
                        'curated.CGP': 'c2.cgp.v7.0.entrez.gmt',
                        'curated.CP.all': 'c2.cp.v7.0.entrez.gmt',
                        'curated.CP.BioCarta': 'c2.cp.biocarta.v7.0.entrez.gmt',
                        'curated.CP.KEGG': 'c2.cp.kegg.v7.0.entrez.gmt',
                        'curated.CP.Reactome': 'c2.cp.reactome.v7.0.entrez.gmt',
                        'curated.CP.PID': 'c2.cp.pid.v7.0.entrez.gmt',
                        'oncogenic.C6': 'c6.all.v7.0.entrez.gmt',
                        'go.All': 'c5.all.v7.0.entrez.gmt',
                        'go.Bio': 'c5.bp.v7.0.entrez.gmt',
                        'go.Cell': 'c5.cc.v7.0.entrez.gmt',
                        'go.Molecular': 'c5.mf.v7.0.entrez.gmt',
                        'motif.gene.sets': 'c3.all.v7.0.entrez.gmt'
    }

    expression = hgene_map(expression) # will map other species to human if needed

    EXCLUDE = ['recno', 'runno', 'searchno', 'label']
    metadata = metadata[[x for x in metadata if x not in EXCLUDE]]


    _genemapper = GeneMapper()

    expression['GeneSymbol'] = expression.index.astype(str).map(lambda x: _genemapper.symbol.get(x, '?'))

    # need to re-order to have GeneSymbol at front
    # order = ['GeneSymbol'] + expression.columns.tolist()
    order = ['GeneSymbol'] + metadata.index.tolist()
    expression = expression[order]



    nrow = len(expression)
    ncol = len(metadata)
    nrow_meta = 1 # only 1 row metadta info - the GeneSymbol. In the future could add more (like description)
    ncol_meta = len(metadata.columns)

    geneset_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                'GSEA', 'genesets', geneset_mapping[geneset])


    # with named_temp(suffix='.gct', mode='w') as f:
    with open(os.path.join(output_prefix, 'input_dataset.gct'), 'w') as f:


        f.write('#1.3\n')

        s = '{}\t{}\t{}\t{}\n'.format(nrow, ncol, nrow_meta, ncol_meta)
        f.write(s)

        s = 'id\tGeneSymbol\t{}\n'.format('\t'.join(metadata.index))
        f.write(s)

        for ix, row in metadata.T.iterrows():
            _s = '{}\t-666\t{}\n'.format(ix, '\t'.join(row.astype(str)))
            f.write(_s)

        # not this
        # metadata.to_csv(f.name, index=False, header=False, sep='\t', mode='a')

        # f.file, not f.name!!! f.name screws up the file for some reason..
        # expression.to_csv(f.file, header=False,  sep='\t', mode='a')
        expression.to_csv(f, header=False,  sep='\t', mode='a')

        f.close() # windows compat

        # # Deal with NAs.. NO DON'T ssGSEA.R does it for you
        # (df[['GeneID', 'GeneSymbol',*meta.index]].pipe(z_score, axis=0)
        #  .to_csv(f.name, index=False, header=False, sep='\t', mode='a')
        # )

        # load ssGSEA into R namespace
        r_source = robjects.r['source']
        rfile = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                         'ssGSEA2.0', 'src', 'ssGSEA2.0.R')
        r_source(rfile)

        ssGSEA2 = robjects.r['ssGSEA2']
        ssGSEA2(f.name,
                output_prefix=output_prefix,
                gene_set_databases=geneset_file,
                sample_norm_type=norm,
                weight=weight,
                statistic=statistic,
                output_score_type=output_score,
                nperm=nperm,
                combine_mode=combine_mode,
                min_overlap=min_overlap,
                correl_type=correl,
                global_fdr=globalfdr,
                extended_output=extended_output,
                par=True, # parallel
                spare_cores=1,
                export_signat_gct=False, # will generate gct files with expression values for each signature
                param_file=True,
                log_file=log_file,
                seed=seed
        )


        # plot results
        rfile = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                             'R', 'ssGSEA_plots.R'
        )
        r_source(rfile)
        make_plots = robjects.r['make_plots']
        make_heatmap = robjects.r['make_heatmap']
        basedir = os.path.split(output_prefix)[0]
        # make_plots(f.name, geneset_file, basedir, geneset, rank_plots=rank_plots)
        make_heatmap(f.name, geneset=geneset, basedir=basedir)




    # ssGSEA2 <- function (
    #                      input.ds,                      ## input data file in gct format, first column (Name) must contain gene symbols^M
    #                      output.prefix,                 ## prefix used for output tables^M
    #                      gene.set.databases=gsea.dbs,   ## list of genesets (in gmt format) to evaluate enrichment on^M
    #                      sample.norm.type=c("rank", "log", "log.rank", "none"),  ## sample normalization^M
    #                      weight= 0,                     ## when weight==0, all genes have the same weight; if weight>0,^M
    #                                                     ## actual values matter, and can change the resulting score^M
    #                      statistic           = c("area.under.RES", "Kolmogorov-Smirnov"), ## test statistic^M
    #                      output.score.type   = c("NES", "ES"),^M
    #                      nperm               = 1000,    ## number of random permutations for NES case^M
    #                      combine.mode        = c("combine.off", "combine.replace", "combine.add"),^M
    #                      ## "combine.off" do not combine *_UP and *_DN versions in a single score.^M
    #                      ## "combine.replace" combine *_UP and *_DN versions in a single score.^M
    #                      ## "combine.add" combine *_UP and *_DN versions in a single score and^M
    #                      ##    add it but keeping the individual *_UP and *_DN versions.^M
    #                      min.overlap         = 10,^M
    #                      correl.type         = c("rank", "z.score", "symm.rank"),  ## correlation type: "rank", "z.score", "symm.rank"^M
    #                      #fdr.pvalue          = TRUE,    ## output adjusted (FDR) p-values^M
    #                      global.fdr          = FALSE,   ## if TRUE calculate global FDR; else calculate FDR sample-by-sample^M
    #                      extended.output     = TRUE,    ## if TRUE the result GCT files will contain statistics about gene coverage etc.  ^M
    #                      par=F,^M
    #                      spare.cores=1,^M
    #                      export.signat.gct=T, ## if TRUE gct files with expression values for each signature will be generated^M
    #                      param.file=T,^M
    #                      log.file='run.log')
