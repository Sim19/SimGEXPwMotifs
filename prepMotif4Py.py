#!/usr/bin/python

# ----------------------------------------------------------------------
# This file is a python script,
# that reshapes the motif file, the result from gimme motif,
# such that it can easily be used with the model.py file
# THE FILE CAN BE CALLED
# 	-VIA CMD-LINE with
#		python code/data_pred/prepMotif4Py.py -m $MOTIFFILE -o $OUTPUTFILE
#	--> saves data in $OUTPUTFILE
# 	-VIA PYTHON LIBRARY
#		import sys
#		sys.path.append('./code/data_pred')
#		import prepMotif4Py
#		prepMotif4Py.main(motif_file, cond_file)
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
#   LOAD PACKAGES
# ----------------------------------------------------------------------
# Data Science
import pandas as pd
import numpy as np

# from LINCS
from cmapPy.pandasGEXpress import parse

# OS and sys
import sys
import getopt  # to hand over arguments from command line
import os
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
#   DEFINE FUNCTIONS
# ----------------------------------------------------------------------

# Read arguments from command line
def get_args(argv):
	
	mfile = ''
	gfile = './data/l1000/genes_lincs_landmark.bed'
	cfile = './data/l1000/mcf7_landmark_level5.gctx'
	ofile = ''
	
	try:
		opts, args = getopt.getopt(argv, 'hm:g::c::o:',
								   ['mfile=', 'gfile=',
									'cfile=', 'ofile='])
	except getopt.GetoptError:
		print('test.py -m <motiffile> -g <genefile> -c <conditionfile> -o <outputfile>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('test.py -m <motiffile> -g <genefile> -c <conditionsfile> -o <outputfile>')
			sys.exit()
		elif opt in ("-m", "--mfile"):
			mfile = arg
		elif opt in ('-g', '--gfile'):
			gfile = arg
		elif opt in ('-c', '--cfile'):
			cfile = arg
		elif opt in ('-o', '--ofile'):
			ofile = arg
	print('Motif input file is %s' %str(mfile))
	print('Gene input file is %s' %gfile)
	print('Conditions file is %s' %cfile)
	print('Output file is %s' %ofile)
	
	return [mfile, gfile, cfile, ofile]


# read files based on names from command lines
def read_file(file_name, sep='\t', lineterminator='\n', header=1,
			  index_col=False):
	#file_motif = motiffile  # './data/l1000/motif_lincs.txt'
	file = pd.read_csv(file_name, sep=sep, lineterminator=lineterminator,
					   header=header, index_col=index_col)
	# motif.iloc[0]
	return file



# get motifs aligned with genes
def align_motif_w_genes(motif_file, cond_file):
	TSS_NCBI = motif_file[motif_file.columns[0]]
	N = TSS_NCBI.shape[0]
	ncbi_id = np.zeros(N)
	idx_swap = np.zeros(N)

	for idx in range(0, N):
		ncbi_id[idx] = TSS_NCBI[idx].split(' ')[1]
		idx_swap[idx] = np.where(
				cond_file.row_metadata_df.index == ncbi_id[idx]
				)[-1]

	motif_file[motif_file.columns[0]] = np.int_(ncbi_id)
	idx_swap = np.int_(idx_swap)

	motif2 = motif_file.copy()

	for idx in range(0, N):
		motif2.iloc[idx] = motif_file.iloc[
			np.where(idx_swap == idx)[-1][0]]

	motif2[motif2.columns[0]] = np.int_(motif2[motif2.columns[0]])
	motif2.rename(columns={'Unnamed: 0': 'NCBI_gene_ID'}, inplace=True)

	# set index
	motif_ind = motif2.set_index("NCBI_gene_ID", )
	motif_ind.index = motif_ind.index.map(unicode)

	return motif_ind
	
	
### MAIN ####
def main(motif_file, cond_file):
	
	# Read motif file
	motif = read_file(file_name=motif_file, header=4)
	
	# Read genes file
	# genes = read_file(file_name=gene_file, header=-1)

	# Read gene-condition matrix (from LINCS)
	conds = parse(cond_file)
	
	motifIndexed = align_motif_w_genes(motif_file=motif, cond_file=conds)
	
	return motifIndexed
	
	
# run main()
if __name__ == '__main__':
	
	# get arguments passed on from command line
	[motiffile, genefile, condfile, outputfile] = get_args(sys.argv[1:])
	
	# call the main function
	motif_indexed = main(motif_file=motiffile, cond_file=condfile)
	
	# save newly indexed motiffile to outputfile
	motif_indexed.to_csv(outputfile, header=True, sep='\t', index=True)

###################################