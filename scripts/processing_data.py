import sys
import pprint

#read my data
def get_my_data (path_1, path_2):
	my_data = list()
	filehandle_1 = open(path_1,'r').readlines()
	filehandle_2 = open(path_2,'r').readlines()
	mw_dic = dict()
	for mw in filehandle_2:
		gene = mw.strip().split()[0]
		if gene not in mw_dic.keys():
			mw_dic[gene] = mw.strip().split()[1]
	for line in filehandle_1:
		line = line.strip()
		if "sp" in line:
			gene_list =[]
			for gene in line.split(';'):
				if "sp" in gene:
					gene = gene.split('|')[1]
					gene_list.append(gene)
			my_data.append({})
			my_data[-1]['gene'] = gene_list
			i = 0
			my_data[-1]["fraction"] = []
			my_data[-1]["ratio"] = []
		else:
			i += 1
			if line.split()[2] != "nan":
				my_data[-1]["fraction"].append(i)
				my_data[-1]['ratio'].append(1/float(line.split()[2]))
	data_list = []
	for data in my_data:
		for gene in data["gene"]:
			data_dic = dict()
			data_dic["gene"] = gene
			data_dic["fraction"] = data["fraction"]
			data_dic["ratio"] = data["ratio"]
			if gene in mw_dic.keys():
				data_dic["molecular_weight"] = int(mw_dic[gene])
			else:
				data_dic["molecular_weight"] = 0
			data_list.append(data_dic)
	return data_list

#process my data
def process_my_data(my_data_list):
	dele = list()
	for my_data in my_data_list:
		if my_data["ratio"] == []:
			dele.append(my_data)
			continue
		eff_length = len(my_data['fraction'])
		max_peak = max(my_data['ratio'])
		sum_peak_area = sum(my_data['ratio'])
		mw_upper_bound = my_data['molecular_weight'] + 6000
		mw_lower_bound = my_data['molecular_weight'] - 6000
		max_index = my_data['ratio'].index(max_peak)
		predict_fraction = my_data['fraction'][max_index]
		predict_mw = (predict_fraction - 1) * 2000 + 15000
		if my_data['molecular_weight'] == 0 or sum_peak_area < 0.5 or max_peak < 0.3:
			dele.append(my_data)
		elif int(my_data['molecular_weight']) > 62000  or int(my_data['molecular_weight']) < 12000:
			dele.append(my_data)
		elif predict_mw >= mw_upper_bound or predict_mw <= mw_lower_bound:
			dele.append(my_data)
	for bad_data in dele:
		my_data_list.remove(bad_data)

protein_list = []
combine_data = sys.argv[1]
mw = sys.argv[2]
protein_list = get_my_data (combine_data,mw)
#pprint.pprint (protein_list)
process_my_data (protein_list)
for protein in protein_list:
	print('| ' + protein['gene'])
	for i in range(1,24):
		if i in protein['fraction']:
			ii = protein['fraction'].index(i)
			print(str(i) + ' ' + str(protein['ratio'][ii]))
		else:
			print(str(i) + ' 0.0')

