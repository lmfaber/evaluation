import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def generate_heatmap(metrics_dataframe, filename):
	"""Generates a heatmap for a pandas dataframe, but exclues the last row. For this script the last row is the sum score of all other rows and has another colour than the heatmap.
	
	Arguments:
		metrics_dataframe {DataFrame} -- dataframe to plot
		filename {str} -- file name of output file.
	"""
	
	# Set the size of the plot to DIN A4 so values and axis are visible
	fig, ax = plt.subplots()
	# the size of A4 paper
	fig.set_size_inches(11.7, 8.27)

	# Set the mask to exclude the last row, which includes the sum score
	mask = np.zeros(metrics_dataframe.shape)
	mask[-1,:] = True
	# Heatmap without the last row which contains the sum score
	hm = sns.heatmap(metrics_dataframe,
					mask=mask, 
					cbar=True,
					annot=True,
					vmin=metrics_dataframe.values[:-1,:].ravel().min(),
					vmax=metrics_dataframe.values[:-1,:].ravel().max(),
					linewidth=0.75,
					cmap='coolwarm')

	# Set the mask to exlucde all data except of the sum score
	mask = np.ones(metrics_dataframe.shape)
	mask[-1,:] = False
	# Heatmap for the last row.
	sns.heatmap(metrics_dataframe,
					mask=mask,
					alpha=0.5,
					cbar=False,
					linewidth=0.75,
					cmap = 'Greys',
					vmax=metrics_dataframe.values[-1,:].ravel().min(),
					annot=True,
					annot_kws={'size': 8, 'color':'black', 'weight': 'bold'},
					fmt='.2f')

	plt.xticks(rotation=45, ha='right')		# Rotate the xaxis at 45 degree
	plt.tight_layout()                      # command that x and y labels are not cut off.
	plt.savefig(filename, format='svg')		# Save the figure as svg

def normalize(array, reverse=False):
	"""
	Normalizes any float/int array to the range of 0-1.
	:param array: float/int-type array
	:return: float-type array
	"""
	if reverse:
		minimum = max(array)
		maximum = min(array)
	else:
		maximum = max(array)
		minimum = min(array)
	returnArray = []
	for a in array:
		if (maximum - minimum) == 0:
			returnArray.append(0)
		else:
			returnArray.append(abs((a-minimum)/(maximum-minimum)))
	return(returnArray)

assembler = snakemake.params['assemblies'] # Used assemblers as given in the contig file

# Create empty dataframe with assemblers as columns
assemblerDict = {}
for a in assembler:
	assemblerDict[a] = []
infoDF = pd.DataFrame(data=assemblerDict)

##############################
########## rnaQUAST ##########
##############################
sensitivityFile = str(snakemake.input['rnaquast_sensitivity'])


# whitespace delimiter, but only if its 2 or more whitespaces after another
sensitivity = pd.read_csv(sensitivityFile, delimiter=r"\s{2,}", index_col=0, na_values = '*', engine='python')
sensitivity = sensitivity.dropna(axis=0)
infoDF = pd.concat([infoDF, sensitivity])

##############################
########## Transrate #########
##############################
assemblies = str(snakemake.input['transrate'])
transrateInfos = pd.read_csv(assemblies, index_col=0).transpose()
transrateInfos = transrateInfos.dropna(axis=0)
transrateInfos.columns = assembler
infoDF = pd.concat([infoDF, transrateInfos])

##############################
########## DETONATE ##########
##############################
detonateDir = str(snakemake.input['detonateDir'])
# kc_paths = snakemake.input['kc']
# rsem_eval_paths = snakemake.input['rsem_eval']

detonateInfos = pd.DataFrame()

for count, a in enumerate(assembler):
	#kc_path = kc_paths[count]
	kc_path = detonateDir + '/kc_' + a + '.txt'
	kc = pd.read_csv(kc_path, delimiter='\t', index_col=0, names=[a])

	contig_nucl_path = detonateDir + '/contig_nucl_' + a + '.txt'
	contig = pd.read_csv(contig_nucl_path, delimiter='\t', index_col=0, names=[a])
	all = kc.append(contig)

	#rsem_eval_path = rsem_eval_paths[count]
	rsem_eval_path = detonateDir + '/rsem_eval_' + a + '.score'
	rsem = pd.read_csv(rsem_eval_path, delimiter='\t', index_col=0, names=[a])
	all = all.append(rsem)


	detonateInfos = pd.concat([detonateInfos, all], axis=1)

infoDF = pd.concat([infoDF, detonateInfos])

##############################
########### BUSCO ############
##############################
busco_dir = str(snakemake.input['busco'])
buscoInfos = pd.DataFrame()
for count, a in enumerate(assembler):
	busco_results = f'{busco_dir}/run_{a}/short_summary_{a}.txt'
	with open(busco_results, 'r') as busco_reader:
		lines = busco_reader.readlines()
		complete_buscos = int(lines[9].split('\t')[1])
		complete_single_buscos = int(lines[10].split('\t')[1])
		complete_duplicated_buscos = int(lines[11].split('\t')[1])
		fragemented_buscos = int(lines[12].split('\t')[1])
		missing_buscos = int(lines[13].split('\t')[1])
	metric_names = ['Complete BUSCOs', 'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs']
	metrics = [complete_buscos, complete_single_buscos, complete_duplicated_buscos, fragemented_buscos, missing_buscos]
	col_name = [a]
	temp_busco_results = pd.DataFrame(metrics, columns=col_name, index=metric_names)
	buscoInfos = pd.concat([buscoInfos, temp_busco_results], axis=1)

infoDF = pd.concat([infoDF, buscoInfos])

##############################
########### Ex90N50 ##########
##############################
EX90N50_dir = str(snakemake.input['EX90N50'])
EX90N50_infos = pd.DataFrame()
for count, a in enumerate(assembler):
	stats_file = f'{EX90N50_dir}/{a}/ExN50.stats'
	with open(stats_file, 'r') as reader:
		for line in reader:
			if line.startswith('90'):
				EX90N50_value = int(line.split()[2])
	metric_names = ['Ex90N50']
	metrics = [EX90N50_value]
	col_name = [a]
	temp_results = pd.DataFrame(metrics, columns=col_name, index=metric_names)
	EX90N50_infos = pd.concat([EX90N50_infos, temp_results], axis=1)

infoDF = pd.concat([infoDF, EX90N50_infos])

##############################
####### Remapping rates ######
##############################

mapping_rates = str(snakemake.input['mapping_rates'])
mapping_infos = pd.read_csv(mapping_rates, sep = '\t', index_col = 0)

infoDF = pd.concat([infoDF, mapping_infos])

# Remove rows that contain NaN values
infoDF = infoDF.dropna(axis=0)
##############################
######### SAVE FILES #########
##############################

# Save unnormalized
rawFile = str(snakemake.output.raw)

# Rename some indices
infoDF.rename(index={'Score': "RSEM_EVAL",
                'kmer_compression_score': 'KC_VALUE'}, inplace=True)
infoDF.to_csv(rawFile, sep='\t')


# Save chosen values in etra file.
reducedRawFile = str(snakemake.output.minRaw)
metricsToExtract = ['Database coverage', 'Duplication ratio', '95%-assembled isoforms', 'reference_coverage','mean_orf_percent',
					'KC_VALUE', 'RSEM_EVAL', 'score', 'optimal_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 'Complete BUSCOs',
					'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs', 'Ex90N50', 'mapping_rate_per_base']
reducedRawScores = infoDF.loc[metricsToExtract]
reducedRawScores.to_csv(reducedRawFile, '\t')


highValues = ['Database coverage', '95%-assembled isoforms', 'reference_coverage', 'mean_orf_percent', 'RSEM_EVAL', 'KC_VALUE', 'score', 'optimal_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 'Complete BUSCOs', 'Complete and single-copy BUSCOs','Ex90N50', 'mapping_rate_per_base']
lowValues = ['Duplication ratio', 'Missing BUSCOs', 'Fragmented BUSCOs', 'Complete and duplicated BUSCOs']

normDF_highValues = reducedRawScores.loc[highValues].apply(normalize, axis='columns', raw=True, result_type='expand')
normDF_lowValues = reducedRawScores.loc[lowValues].apply(normalize, args=(True,), axis='columns', raw=True, result_type='expand')
normDF = normDF_highValues.append(normDF_lowValues)

colSum = pd.DataFrame(normDF.apply(sum, axis='rows')).transpose()
colSum.index = ['SUM_SCORE']
normDF = pd.concat([normDF, colSum])

# Save reduced normalized
normalizedFile = str(snakemake.output.normalized)
normDF.columns = assembler
normDF.to_csv(normalizedFile, sep='\t')
heatmapSVG = str(snakemake.output.heatmap)
generate_heatmap(normDF, filename=heatmapSVG)



# Save all normalized
# normalizedFile = str(snakemake.output.normalized)
# normDF = infoDF.apply(normalize, axis='columns', raw=True, result_type='expand')
# normDF.columns = assembler
# normDF.to_csv(normalizedFile, sep='\t')
