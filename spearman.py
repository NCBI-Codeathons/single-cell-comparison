import numpy as np
import pandas as pd
import scipy as sp
import scipy.cluster.hierarchy as hc
from scipy import stats
from os import listdir
import os.path
import glob
import string
from matplotlib import pyplot

import pdb

base = "/data/gene_count_length_files/SRP050499/"
acc1 = "SRR2013520"
acc2 = "SRR2013521"

#pdb.set_trace()

flist = glob.glob("../../data/gene_count_length_files/SRP011546/*.tsv")

#flist = fflist[0:10]

rank_mat = np.zeros(( len(flist),len(flist) ))

for ii, a in enumerate(flist):
    for jj,b in enumerate(flist):
        if (a >b ):
            #fname1 = base + acc1 + ".tsv"
            #fname2 = base + acc2 + ".tsv"
            fname1 = a
            fname2 = b
            (run1,suf) = os.path.basename(a).split(".")
            (run2,suf) = os.path.basename(b).split(".")

            df1 = pd.read_csv(fname1,sep='\t',skiprows=2,names=['g','len','count'])
            df2 = pd.read_csv(fname2,sep='\t',skiprows=2,names=['g','len','count'])

            df1['rank'] = df1['count'].rank(method="min",ascending=False)
            df2['rank'] = df2['count'].rank(method="min",ascending=False)

            #intersection of nonzero element cardinality? nb nonzero returns odd tuple, need zeroth element 
            idxp1 = pd.Series(data=df1['count'].to_numpy().nonzero()[0])
            idxp2 = pd.Series(data=df2['count'].to_numpy().nonzero()[0])

            i_both=pd.concat([idxp1,idxp2],join='inner')

            corr,p = stats.spearmanr(df1['count'],df2['count']) 
            print ( corr, p, i_both.shape, run1,run2 )
            rank_mat[ii][jj] = corr
            rank_mat[jj][ii] = corr


#hierarchical clustering of ranks
p_rank = pd.DataFrame(data=rank_mat[0:,0:])
link = hc.linkage(p_rank,method='centroid')
o1 = hc.leaves_list(link)
disp_mat = p_rank.iloc[o1,:]
disp_mat = p_rank.iloc[:, o1[::-1]]
pyplot.imshow(disp_mat)
pyplot.savefig('SRP011546.png')

print (disp_mat)
