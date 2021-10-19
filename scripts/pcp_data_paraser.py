import sys
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from math import isnan, log, exp
import scipy.linalg as la

its = "Intensity "
lfq = "LFQ intensity "
rat = "Ratio H/L "

experiment_lst = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"]

lines = open("proteinGroups.txt", 'r').readlines()
#tags
tags = {}
es = lines[0].strip().split('\t')
for i, e in enumerate(es):
  tags[e] = i

lib_fn = sys.argv[1]
map_mass_kD = {}
for record in SeqIO.parse(lib_fn, "fasta"):
    try:
        analysed_seq = ProteinAnalysis(str(record.seq))
        map_mass_kD[record.id] = int(analysed_seq.molecular_weight()/1000)
    except:
        pass

#pdf = matplotlib.backends.backend_pdf.PdfPages("outputs.pdf")

def process(name, Ls, Hs, Rs):
    mass_lst = []
    for p in name.split(';'):
        if p in map_mass_kD.keys():
            mass_lst.append( map_mass_kD[p] )
    if len(mass_lst) == 0: return
    real_mass = np.median(mass_lst)
    print name, real_mass
    #for l, h, r in zip( Ls, Hs, Rs ):
    #    print l, h, r
    ratio = 0.5
    kT = 1.0
    Nfrac = len(Ls)
    #print Nfrac, len(Hs), len(Rs)
    P_mat = np.zeros([Nfrac+1,Nfrac+1])
    LH_tr = {}
    vHs = []
    sum_tr = 0.0
    sum_H = 0.0
    for n in xrange(Nfrac):
        nf = n+1
        vH = Hs[n]
        vR = Rs[n]
        #print vH, vR
        vHs.append(vH)
        sum_H += vH
        if not isnan(vR):
            LH_tr[nf] = vR
            sum_tr += vR
    vHs = np.asarray(vHs)
    vHs = vHs / sum_H
    #fill H-LFQ prob
    for i in xrange(Nfrac):
        P_mat[i+1, 0] = 1-ratio
        for j in xrange(Nfrac):
            P_mat[i+1, j+1] = vHs[j]*ratio
    #connect with L
    Es = {}
    Ps = {}
    sum_expE = 1.0
    for k in LH_tr.keys():
        r = LH_tr[k]
        P_L = 1/(1+r)
        P_H = r/(1+r)
        E_L = -log(P_L)
        E_H = -log(P_H)
        dE = (E_H - E_L)
        Es[k] = dE
        sum_expE = sum_expE + exp(-dE/kT)
    #fill first line
    P_mat[0,0] = 1.0/sum_expE
    for k in LH_tr.keys():
        Ps[k] = exp(-Es[k]/kT)/sum_expE
        P_mat[0, k] = Ps[k]
    #limit
    P_mat = np.dot(P_mat, P_mat)
    P_mat = np.dot(P_mat, P_mat)
    P_mat = np.dot(P_mat, P_mat)
    sum_Hs = np.sum(P_mat[0,1:])
    correct_Hs = P_mat[0,1:]/sum_Hs
    for i in xrange(Nfrac):
        print "%6.4f" % correct_Hs[i]

for l in lines[1:]:
    es = l.strip().split('\t')
    pros = es[tags["Protein IDs"]]
    raw_intens_L = []
    lfq_intens_L = []
    raw_intens_H = []
    lfq_intens_H = []
    ratios = []
    sum_lfq = 0.0
    sum_its = 0.0
    for nexp in experiment_lst:
        raw_intens_L.append( float(es[tags[its+"L "+nexp]]) )
        lfq_intens_L.append( float(es[tags[lfq+"L "+nexp]]) )
        raw_intens_H.append( float(es[tags[its+"H "+nexp]]) )
        lfq_intens_H.append( float(es[tags[lfq+"H "+nexp]]) )
        ratios.append( float(es[tags[rat+nexp]]) )
        sum_lfq += lfq_intens_L[-1] + lfq_intens_H[-1]
        sum_its += raw_intens_L[-1] + raw_intens_H[-1]

    if sum_lfq > 10.0:
        process(pros, lfq_intens_L, lfq_intens_H, ratios)
    elif sum_its > 10.0:
        process(pros, raw_intens_L, raw_intens_H, ratios)

#pdf.close()
