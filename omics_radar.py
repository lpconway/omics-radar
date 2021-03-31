# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 22:50:30 2021

@author: Louis
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex
import numpy as np

mpl.rcParams['font.sans-serif'] = ['Arial'] + mpl.rcParams['font.sans-serif']

def radar(fold_changes,
                     p_values,
                     labels,
                     p_threshold = 0.05,
                     fold_threshold = 3,
                     scaling_function = lambda a : a**2 / 2,
                     y_max = None,
                     cmap='ocean',
                     edgecolor='#4A4D57',
                     seed=None
                     ):
    """
    
    Parameters
    ----------
    fold_changes : array-like
        An object containing the fold changes for each protein
        (rows) in each experimental treatment (columns).
    p_values : array-like
        An object containing the p-values for each protein (rows)
        in each experimental treatment (columns).
    labels : list of strings
        The labels for each experimental treatment.
    p_threshold : numeric, optional
        p-Value significance threshold. The default is 0.05.
    fold_threshold : numeric, optional
        Fold-change significance threshold. The default is 3.
    scaling_function : function, optional
        A function relating circle size to fold-change. The default is
        lambda a : a**2 / 2.
    y_max : numeric, optional
        The maximum radius of the plot. If None, maxmimum radius is -log10 of
        the minimum non-NaN p value. The default is None.
    cmap : string, optional
        The matplotlib colormap used to color each sector. The default is
        'ocean'.
    edgecolor : string, optional
        The color of the edges of the circles. The default is '#4A4D57'.
    seed : numeric
        Random seed controlling the angular plaement of the circles

    Returns
    -------
    figure, axis

    """

    # Find angle for each sector
    n=len(labels)
    angle = np.pi*2/n
    
    # Coerce input fold and p-values into arrays
    fold_changes = np.array(fold_changes)
    p_values = np.array(p_values)
    
    # Generate list of colors for each sector
    colormap = plt.cm.get_cmap(cmap)
    colors = [rgb2hex(colormap(x)) for x in np.arange(0, 1+1/n, 1/n)]

    # Plot data
    fig = plt.figure(figsize=(4,4), dpi=300)
    axis = fig.add_subplot(111, projection='polar')

    for index in range(n):
        # Get parameters for this sector
        label = labels[index]
        theta_min = angle*index
        color = colors[index]
        fold = fold_changes[:,index]
        p = p_values[:,index]
        
        # Create a subset of 'significant' points (i.e. above the fold-change
        # threshold and below the p-value threshold)
        mask = (fold > fold_threshold) & (p < p_threshold)
        subset_n = np.sum(mask)
        subset_fold = fold[mask]
        subset_p = p[mask]
        
        # The angle of each point is random
        np.random.seed(seed=seed)
        subset_theta = np.random.random(subset_n) * angle + theta_min
        
        # Plot 
        axis.scatter(subset_theta, -np.log10(subset_p), s=scaling_function(subset_fold), alpha = 0.5, c=color, edgecolor=edgecolor)
        
        inv_mask = mask == False
        inv_n = np.sum(inv_mask)
        inv_fold = fold[inv_mask]
        inv_p = p[inv_mask]
        inv_theta = np.random.random(inv_n) * angle + theta_min
        axis.scatter(inv_theta, -np.log10(inv_p), s=scaling_function(inv_fold), alpha = 0.05, c=color)
        
        axis.set_xticks(np.arange(0, np.pi*2, angle))
        axis.tick_params(axis='x',label1On=False)
        axis.tick_params(axis='y',label1On=False)
        
    if y_max:
        axis.set_ylim(0, y_max)
        
    for y in axis.get_yticks():
        if y <= axis.get_ylim()[1]:
            axis.text(0, y, "{:2.0f}".format(y), va='top', ha='right', fontsize=6)

    axis.text(0, axis.get_ylim()[1]*1.05, "-log$_{10}(p)$", va='top', ha='left', fontsize=6)

    axis.set_yticks([-np.log10(p_threshold)], minor=True)
    axis.tick_params(which='minor', length=5, width=2, direction='in')
    axis.grid(which='minor', linewidth='1', color='black')

    top = axis.get_ylim()[1]*1.1
    for label, iota, color in zip(labels, np.arange(angle/2, np.pi*2+angle/2, angle), colors):
        axis.text(iota, top, label, va='center', ha='center', rotation=180*iota/np.pi-90, color=color)
    
    max_log_fold = int(np.log2(np.nanmax(fold)))+1
    legend_plots = []
    for log_fold in range(max_log_fold):
        legend_plot = plt.scatter([],[], s=scaling_function(2**log_fold), c=colors[0], alpha=0.5, edgecolors=edgecolor)
        legend_plots.append(legend_plot)
    
    legend_labels = [str(2**x) for x in np.arange(1, max_log_fold+1)]

    axis.legend(legend_plots, legend_labels, frameon=False, fontsize=10, handlelength=2, loc = 'right', handletextpad=1, scatterpoints = 1, bbox_to_anchor=(1.4, 0.5), title='Fold')
    
    return(fig, axis)