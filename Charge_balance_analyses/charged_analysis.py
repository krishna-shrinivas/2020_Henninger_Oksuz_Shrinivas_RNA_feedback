#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 15:01:45 2019

@author: krishna
"""

import os
from Bio import SeqIO
from localcider.sequenceParameters import SequenceParameters
import argparse

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib; mpl =matplotlib

mpl.rcParams['figure.figsize'] = (12.0,9.0) # default = (6.0, 4.0)
mpl.rcParams['font.size']      = 28        # default = 10

mpl.rcParams['axes.linewidth']    = 0.75 # default = 1.0
mpl.rcParams['lines.linewidth']   = 1.5 # default = 1.0
mpl.rcParams['patch.linewidth']   = 1.0 # default = 1.0
mpl.rcParams['grid.linewidth']    = 0.5 # default = 0.5
mpl.rcParams['xtick.major.width'] = 1.0 # default = 0.5
mpl.rcParams['xtick.minor.width'] = 0.0 # default = 0.5
mpl.rcParams['ytick.major.width'] = 1.0 # default = 0.5
mpl.rcParams['ytick.minor.width'] = 0.0 # default = 0.5

def make_nice_axis(ax):
    """ Function to beautify axis"""

    ax.spines['top'].set_visible(False) # hide top axs
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 30))
    ax.spines['left'].set_position(('outward', 30))
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(pad=10)
    ax.yaxis.set_tick_params(pad=10)
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20

def return_equal_charge_concentrations(protein_file,RNA_file,protein_conc=1000,count_histidine=0,count_gfp=1):


    protein = list(SeqIO.parse(protein_file,'fasta'))
    RNAS = list(SeqIO.parse(RNA_file,'fasta'))

    pos_protein = 0;
    for p in protein:
        protein_obj = SequenceParameters(str(p.seq))
        pos_protein = protein_obj.get_countPos() -protein_obj.get_countNeg()

    GFP = SequenceParameters("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK")



    if count_histidine:
        pos_protein = pos_protein + 0.5*protein_obj.get_amino_acid_fractions()['H']*(protein_obj.get_length())
    if count_gfp:
        pos_protein = pos_protein +(GFP.get_countPos()-GFP.get_countNeg())
    counter_conc_max = {};
    for rna in RNAS:
        counter_conc_max[rna.id.split('(')[0]] = protein_conc*pos_protein/len(rna.seq);

    return(counter_conc_max)



def plot_partition(args):
    partition_data = args.i;
    df = pd.read_excel(partition_data,skiprows=1,sheet_name=args.partition)

    counter_conc_max = return_equal_charge_concentrations(args.p,args.r,protein_conc=float(args.conc))

    counter_rna_peak = list(counter_conc_max.values())[0];
    RNA = [str(x) for x  in list(df.columns.values)]
    RNA_counter = [ min(float(x),counter_rna_peak)/max(counter_rna_peak,float(x)) for x in RNA]




    fig,axes = plt.subplots(1,1,figsize=(16,12))

    axes = sns.boxplot(data=df,palette=sns.color_palette("Blues"))

    for patch in axes.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.8))

    axes = sns.swarmplot(data=df,palette=sns.color_palette("Blues"))
    axes.set_ylabel(args.partition + ' Partition',color='#377eb8')
    axes.set_xlabel('RNA concentration (nM)')
    axes.tick_params(axis='y', labelcolor='#377eb8')

    make_nice_axis(axes)

    ax2 = axes.twinx()  # instantiate a second axes that shares the same x-axis
    make_nice_axis(ax2)

    ax2.spines['top'].set_visible(False) # hide top axs
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(True)

    ax2.spines['right'].set_position(('outward', 30))
    ax2.yaxis.set_ticks_position('right')


    color = '#ff7f00'
    ax2.set_ylabel('Charge ratio', color=color)
    ax2.set_ylim((0,1))
    ax2.plot(RNA, RNA_counter, color=color,lw=3)
    ax2.tick_params(axis='y', labelcolor=color)

    file_save = False;
    if args.output is not None:
        file_save = True;
    if file_save:
        file_name =  'Output/' + args.output;
        os.makedirs(os.path.dirname(file_name),exist_ok=True);
        fig.savefig(file_name+'.svg',format='svg',dpi=600,bbox_inches='tight')
        fig.savefig(file_name+'.png',format='png',dpi=600,bbox_inches='tight')
    plt.show()


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Given 2 fasta files + partition, plot output')
  parser.add_argument('--i',help="Name of partition file", required = True);
  parser.add_argument('--p',help="Name of protein fasta", required = True);
  parser.add_argument('--r',help="Name of RNA fasta", required = True);
  parser.add_argument('--partition',help="Name of partition to plot, defaults to protein", required = False, default= 'protein');

  parser.add_argument('--conc',help="Protein concentration in nM", required = False,default=1000.0);
  parser.add_argument('--output',help="Name of output file", required = False);
  args = parser.parse_args();

  plot_partition(args)
