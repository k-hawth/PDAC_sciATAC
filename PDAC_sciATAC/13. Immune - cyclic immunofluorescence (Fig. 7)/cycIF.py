#!/usr/bin/env python
# coding: utf-8

'''
python 3.11.14
pandas 2.3.3
numpy 2.3.5
seaborn 0.13.2
matplotlib 3.10.8
scipy 1.16.3
statsmodels 0.14.6
statannotations 0.7.1
'''

# #  cycIF results patient level


# import libraries
import importlib
import pandas as pd
import util
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings


# set path
kevin_out = './'
os.makedirs(kevin_out, exist_ok=True)
os.makedirs(f'{kevin_out}/figures', exist_ok=True)


# ## load metadata
df_atac = pd.read_csv(f'{kevin_out}/cycIF_metadata.csv',index_col=0)


# load and plot
figsize = (2.5,4)
ls_pos = [
    'results_20230326_U54-TMA_FractionPositive_byPatient_byleiden_inT cell_PDAC.csv',#just pri for atac analysis
 'results_20230326_U54-TMA_FractionPositive_byPatient_byleidencelltype5_inepithelial_PDAC.csv',
 'results_20230326_U54-TMA_FractionPositive_byPatient_byleidencelltype5_inimmune_PDAC.csv',]
test_type = 'Mann-Whitney'
linebreak = '\n'
lls_tests = [['ATAC_Subtype',['SIP', 'Classical-NRP']],
             ['ATAC_Subtype',['SIP', 'Classical-NRP','Basal Committed']],]
for s_load in ls_pos:
    print(f'loading {s_load}')
    df_myc = pd.read_csv(f'{kevin_out}/{s_load}',index_col=0)
    s_title = s_load.split('_')[-1].split('.csv')[0]
    s_end = s_load.split("_in")[-1].split('_by')[-1].split('_')[-2].split(".")[0]
    s_end2 = s_load.split("_in")[-1].split('_by')[-1].split('_')[-1].split(".")[0]
    print(s_end)
    ls_col = df_myc.columns[(df_myc.dtypes=='float64')]
    df_myc['ATAC_Subtype'] = df_myc.index.map(dict(zip(df_atac.Public_Specimen_ID,df_atac.sciATAC_Subtype)))
    for x, ls_groups in lls_tests:
        if df_myc[x].nunique() < 2:
            print(f'skipping {x}')
            continue
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            for y in ls_col:
                if y.find('DAPI') > -1 or y.find('+') > -1:
                    continue
                if len(ls_groups) == 2:
                    b_correct=False
                    fig_alpha = .05
                else:
                    b_correct=True 
                    fig_alpha = .1
                data = df_myc[(df_myc.loc[:,x].notna())]
                fig, ax, corrected, raw_pvals = util.debugged_fdr_plot(
                          data, x, y, ls_groups, test_type,figsize=figsize,b_correct=b_correct)
                if (np.array(corrected) > fig_alpha).all():
                    plt.close(fig)
                else:
                    s_end_title = f'{s_end.replace("_"," ")}'
                    ax.set_ylabel(f'Fraction {y.replace("_"," ")}+')
                    ax.set_title(f'{y.replace("_"," ")} in {s_end_title}')
                    plt.tight_layout()
                    plt.show(fig)
                    fig.savefig(f'{kevin_out}/figures/stripplot_TMA_clust_vs_{y}_{x}_{s_end}_{s_end2}.pdf',bbox_inches='tight')


# end
print('done')







