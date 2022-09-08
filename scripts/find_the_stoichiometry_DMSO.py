# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:23:51 2022

@author: reslo
"""

import sys
import matplotlib.pyplot as plt
#from scipy.cluster.hierarchy import fcluster, median, centroid, average, weighted, dendrogram
#from scipy.spatial.distance import pdist
import numpy as np
import scipy.stats as stats
import matplotlib
import matplotlib.backends.backend_pdf

def process_my_data(my_data_list):
    a = 0
    b = 0
    c = 0
    d = 0
    dele = list()
    for my_data in my_data_list:
        if my_data["ratio"] == []:
            dele.append(my_data)
            continue
            d += 1
        max_peak = max(my_data['ratio'])
        sum_peak_area = sum(my_data['ratio'])
        mw_upper_bound = my_data['molecular_weight'] + 10000
        mw_lower_bound = my_data['molecular_weight'] - 10000
        max_index = my_data['ratio'].index(max_peak)
        predict_mw = max_index * 2500 + 3750
        if my_data['molecular_weight'] == 0 or max_peak/sum_peak_area < 0.1 or max_peak <= 0:
            dele.append(my_data)
            a += 1
        elif int(my_data['molecular_weight']) > 62000:
            dele.append(my_data)
            b += 1
        elif predict_mw >= mw_upper_bound or predict_mw <= mw_lower_bound:
            dele.append(my_data)
            c += 1
    for bad_data in dele:
        my_data_list.remove(bad_data)
path = sys.argv[1]
mw_path = sys.argv[2]
mw_file = open(mw_path,'r').readlines()
mw_dic = {}
for mw_line in mw_file:
    ipi = mw_line.strip().split()[0]
    mw = int(mw_line.strip().split()[1])
    mw_dic[ipi] = mw
text = open(path,'r').readlines()
protein = []
for line in text:
    if '|' in line:
        protein.append({})
        gene = line.strip().split()[1]
        protein[-1]['symbol'] = gene
        protein[-1]['ratio'] = []
        if gene in mw_dic.keys():
            protein[-1]['molecular_weight'] = mw_dic[gene]
        else:
            protein[-1]['molecular_weight'] = 0
    else:
        protein[-1]['ratio'] += [float(line.strip().split()[1])]
#print(protein)
process_my_data(protein)
cluster_data = []
cluster_label = []
xdata = [2.5*i + 3.75 for i in range(1,24)]
output_file = open("index_DMSO.txt",'w')
matplotlib.use('Agg')
label = sys.argv[3]
ref_file = open("ref_list_"+label+".txt",'w')
pdf = matplotlib.backends.backend_pdf.PdfPages("output_"+label+".pdf")
for p in protein:
    max_peak = max(p['ratio'])
    max_peak_index = p['ratio'].index(max_peak)
    predict_index = (p['molecular_weight'] - 3750)//2500
    stoi = []
    b = sum(p["ratio"])
    for i in range(0,23):
        p["ratio"][i] = p['ratio'][i]/b
    for i in range(0,23):
        if p['ratio'][i] < 0.1:
            p['ratio'][i] = 0.0
    b = max(p["ratio"])
    if b != 0:
        for i in range(0,23):
            p["ratio"][i] = p['ratio'][i]/b
    else:
        continue
    index = []
    a = 0
    p['ratio'].append(0.0)
    ref = 0
    for r in p["ratio"]:
        if r >= 0.1:
            ref += 1
    if ref == 1:
        print(p['symbol'],file = ref_file)
#    if p["symbol"] == "P04406":
#        print(p["ratio"])
    for r in p['ratio']:
        label = 0
        if r >= 0.1:
            ind = p['ratio'].index(r)
            if ind in range(predict_index - 2, predict_index + 3):
                if ind == 23:
                    stoi.append(1.0)
                    index.append(ind)
                else:
                    stoi.append(0.0)
                    index.append(ind)
                    if ind == 0:
                        stoi[-1] += (p["ratio"][ind])
                        ind += 1
                    while(ind < 23 and label < 5):
                        stoi[-1] += (p["ratio"][ind] + p["ratio"][ind-1])
                        if p["ratio"][ind] < 0.1:
                            label += 1
                            stoi.append(0.0)
                            index.append(ind + 1)
                        elif ind != 22:
                            if p["ratio"][ind] <= p['ratio'][ind-1] + 0.2 and p["ratio"][ind] <= p['ratio'][ind+1] + 0.2:
                                stoi.append(0.0)
                                index.append(ind + 1)
                        else:
                            label = 0
                        ind += 1
                break
    if len(stoi) != 0 and sum(stoi) != 0:
        b = sum(stoi)
        for i in range(0,len(stoi)):
            stoi[i] = stoi[i]/b
        for i in range(0,len(stoi)):
            if stoi[i] < 0.1:
                stoi[i] = 0.0
        b = sum(stoi)
        for i in range(0,len(stoi)):
            stoi[i] = stoi[i]/b
        print(p["symbol"],sep = '\t',end = '\t')
        for i in range(0,len(stoi)):
            if stoi[i] != 0.0:
                print(str(stoi[i]),sep = '\t',end = '\t')
        print('\t')
        cluster_data.append(stoi)
        cluster_label.append(p["symbol"])
        ydata =[]
        b = sum(p["ratio"])
        for ratio in p["ratio"]:
            ratio = ratio/b
            ydata.append(ratio)
        pre_mw = float(p["molecular_weight"])/1000.0
        x_exp = [-15,-5]
        x_pre = [-30,-20]
        f = plt.figure()
        print('|' + '\t' + p["symbol"],file = output_file)
        si = 2*np.sqrt(2*np.pi)
        for i in range(0,len(stoi)):
            if i != len(stoi)-1 and index[i+1] - index[i] != 1:
                lw = 6*(index[i+1] - index[i])
                mw_exp = 3.75 + 1.25*(index[i+1] + index[i])
                if stoi[i] != 0:
                    y_gaussian = np.arange(3.75 + 2.5*index[i],3.75 + 2.5*index[i+1],0.001)
                    x_gaussian = stats.norm.pdf(y_gaussian,mw_exp,(index[i+1]-index[i])/si)
                    x_g = [gau*stoi[i]*100 for gau in x_gaussian]
                    plt.plot(x_g,y_gaussian,'b')
                    plt.text(110,mw_exp,str(round(100*stoi[i],1)) + "%",verticalalignment = "center")
                    print(index[i],stoi[i],sep = '\t',file = output_file)
            else:
                lw = 6
                mw_exp = 3.75 + 2.5*index[i]
                if stoi[i] != 0:
                    y_gaussian = np.arange(3.75 + 2.5*index[i] - 2,3.75 + 2.5*index[i] + 2,0.001)
                    x_gaussian = stats.norm.pdf(y_gaussian,mw_exp,2/si)
                    x_g = [gau*stoi[i]*100 for gau in x_gaussian]
                    plt.plot(x_g,y_gaussian,'b')
                    plt.text(110,mw_exp,str(round(100*stoi[i],1)) + "%",verticalalignment = "center")
                    print(index[i],stoi[i],sep = '\t',file = output_file)
            plt.plot(x_exp,[mw_exp,mw_exp],color = (1-stoi[i],1-stoi[i],1-stoi[i]),linewidth = lw)
        plt.axis([-20,135,0,66])
        plt.xticks([0,20,40,60,80,100])
        pre_x = [pre_mw - 2, pre_mw, pre_mw + 2]
        pre_y = [0,1,0]
#        for i in range(0,len(ydata)):
#            if ydata[i] > 0.1:
#                plt.plot([1,100*ydata[i]],[xdata[i],xdata[i]],color = (0.3,0.3,1),linewidth = 6)
        plt.plot([0,0],[0,66],color = (0,0,0))
        plt.ylabel("molecular weight (kDa)")
        plt.xlabel("Density")
        plt.title(p["symbol"])
        pdf.savefig(f)
        plt.close()
pdf.close()
output_file.close()