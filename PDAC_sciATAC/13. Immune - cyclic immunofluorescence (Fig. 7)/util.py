import itertools
import textwrap

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator


def debugged_fdr_plot(data, x, y, order, test_type='Mann-Whitney',figsize=(4, 4),b_correct=True,s_plot_type='box'):
    # 1. Prepare Data - Strict cleaning
    clean_df = data.dropna(subset=[x, y]).copy()
    # 2. Manual Stat Calculation (No Annotator Magic)
    # We do this FIRST so we know the p-values are real
    pairs = list(itertools.combinations(order, 2))
    raw_pvals = []
    valid_pairs = []

    for p1, p2 in pairs:
        group1 = clean_df[clean_df[x] == p1][y]
        group2 = clean_df[clean_df[x] == p2][y]

        if len(group1) > 1 and len(group2) > 1:
            # Perform the actual test
            if test_type == 't-test_ind':
                # Welch's T-test (equal_var=False) is safer for unequal N
                stat, p = ttest_ind(group1, group2, equal_var=False)
            elif test_type == 'Mann-Whitney':
                stat, p = mannwhitneyu(group1, group2)

            raw_pvals.append(p)
            valid_pairs.append((p1, p2))
    # 3. FDR Correction
    if b_correct:
        _, corrected_pvals, _, _ = multipletests(raw_pvals, method='fdr_bh')
    else:
        corrected_pvals = raw_pvals

    # 4. Plotting
    fig, ax = plt.subplots(figsize=figsize)
    if s_plot_type == 'box':
        sns.boxplot(data=clean_df, x=x, y=y, order=order, hue_order=order,hue=x,showfliers=False, ax=ax, palette='muted')
    elif s_plot_type == 'violin':
        sns.violinplot(data=clean_df, x=x, y=y, order=order, hue_order=order,hue=x,inner='quart', ax=ax, palette='muted')
    sns.stripplot(data=clean_df, x=x, y=y, order=order,hue=x,hue_order=order, ax=ax, palette='dark', alpha=0.8,s=3)

    # Add N-counts manually to be safe
    n_counts = clean_df.groupby(x,observed=False).size()
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{l}\n(N={n_counts.get(l, 0)})" for l in order], rotation=90)

    # Wraps the text every 30 characters
    wrapped_label = "\n".join(textwrap.wrap(str(y), width=30))

    plt.ylabel(wrapped_label)
    # 5. Use Annotator ONLY for drawing the lines
    # We pass our manually calculated/corrected p-values
    annotator = Annotator(ax, valid_pairs, data=clean_df, x=x, y=y, order=order)
    annotator.configure(test=None, verbose=False)
    if b_correct:
        formatted_p = [f"FDR={p:.3f}" if p > 0.001 else f"FDR={p:.2e}" for p in corrected_pvals]
    else:
        formatted_p = [f"P={p:.3f}" if p > 0.001 else f"P={p:.2e}" for p in corrected_pvals]
    annotator.set_custom_annotations(formatted_p)
    annotator.annotate()
    # Safely remove only if it exists
    legend = ax.get_legend()
    if legend:
        legend.remove()
    plt.tight_layout()
    return fig, ax, np.array(corrected_pvals), raw_pvals
