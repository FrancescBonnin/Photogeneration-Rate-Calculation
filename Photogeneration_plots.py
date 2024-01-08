# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 10:00:35 2019

@author: Francesc Bonnín Ripoll

Aquest és el programa definitiu pel càlcul de la photogeneracio
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import matplotlib.ticker as mtick

def G_Plot(array_z, array_G, label_plot, input_file, total_width, directory):
    
    file_name = 'G'+ input_file
    os.chdir(directory)
    sns.set()
    sns.set_style("ticks",{"xtick.direction": "in","ytick.direction": "in"})

    plt.figure()
    sns.set_palette('colorblind')

    fig = plt.scatter(array_z,array_G,label=label_plot, marker='.')
    plt.xlabel('$z$ (nm)',fontsize=20)
    plt.ylabel('$G$ (cm$^{-3}$)',fontsize=20)
    plt.legend(loc=0,ncol=1, shadow=True, fancybox=True, fontsize=15)
    plt.tick_params(labelsize=18,pad=10)
    plt.xlim(0,total_width)
    plt.savefig(file_name + '.png', format='png', dpi=1000,bbox_inches='tight')   
    plt.show()
    
    return fig
