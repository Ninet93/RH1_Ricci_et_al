import os, sys
import re
from collections import OrderedDict
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped, generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from operator import itemgetter
from itertools import groupby
from textwrap import wrap


def Filter_CDS(MA_file, path_MA):
	os.chdir(path_MA)

	print(MA_file)

	out_f_bis = MA_file.split('.')[0] + '_length.txt'
	out_bis =  open(out_f_bis, 'w')
	out_bis.write('#ID\tlength\n')

	out_f = MA_file.split('.')[0] + '.txt'
	out = open(out_f, 'w')
	out.write('#ID\tDivisibleBy3\tNoStartCodon\tNoStopCodon\tEarlyCodonStop\tPos_EarlyCodonStop\n')


	record = SeqIO.to_dict(SeqIO.parse(MA_file, 'fasta'))

	max_seq=[]
	for ID in record.keys():
		seq_ID_wo_gaps = str(record[ID].seq).replace('-', '')

		len_seq_wo_gap = len(seq_ID_wo_gaps)

		out_bis.write('%s\t'%ID)
		out_bis.write('%d\n'%len_seq_wo_gap)

		max_seq.append(len_seq_wo_gap)

	maximum = max(max_seq)

	df_max_seq = pd.Series(max_seq)
	print('Summary statistics of CDS length:')
	print(df_max_seq.describe())

	max_calcul = float(maximum/3)

	if max_calcul.is_integer() == False:
		print('The maximal nucleotide length in MA is not divisible by 3.')
		print(maximum)


	sequences_CDS=[]
	sequences_AA=[]
	not_divisibly_by_3_YN=False
	NoStartCodon_YN=False
	NoStopCodon_YN=False
	for ID in record.keys():
		not_divisibly_by_3_YN=False
		NoStartCodon_YN=False
		NoStopCodon_YN=False

		seq_ID_wo_gaps = str(record[ID].seq).replace('-', '')

		len_seq_wo_gap = len(seq_ID_wo_gaps)

		calcul = float(len_seq_wo_gap/3)

		if calcul.is_integer() == True:

			CDS_seq = SeqRecord(Seq(str(record[ID].seq).replace('-', '')), id=ID, description='')

			AA = str(Seq(seq_ID_wo_gaps).translate())

			AA_seq = SeqRecord(Seq(AA), id=ID, description='')

		else:

			not_divisibly_by_3_YN = True



		# check if early stop codon
		#codons = [str(record[ID].seq).upper()[i:i+3] for i in range(0, len(record[ID].seq), 3)] # gaps not removed
		codons = [seq_ID_wo_gaps.upper()[i:i+3] for i in range(0, len(seq_ID_wo_gaps), 3)] # gaps removed

		stop_codons=['TAA', 'TAG', 'TGA']

		NoStartCodon_YN=False
		NoStopCodon_YN=False

		if codons[0] != 'ATG':
			NoStartCodon_YN = True
			#print('Start')
			#print(ID)
			#print(codons[0])

		if codons[len(codons)-1] not in stop_codons:
			NoStopCodon_YN = True
			#print('Stop')
			#print(ID)
			#print(codons[len(codons)-1])



		pos_codons=-2
		early_stop_codons_YN=False
		list_pos_stop_codons=[]
		for x in codons[:-1]:
			pos_codons += 3
			if x in stop_codons:
				early_stop_codons_YN = True
				list_pos_stop_codons.append(str(pos_codons))

		list_YN=[NoStartCodon_YN, NoStopCodon_YN, not_divisibly_by_3_YN, early_stop_codons_YN]

		if True in list_YN:
			out.write('%s\t'%ID)
			out.write('%s\t'%str(not_divisibly_by_3_YN))
			out.write('%s\t'%str(NoStartCodon_YN))
			out.write('%s\t'%str(NoStopCodon_YN))
			out.write('%s\t'%str(early_stop_codons_YN))


			if early_stop_codons_YN == True:
				list_pos_stop_codons_bis = ', '.join(list_pos_stop_codons)
				out.write('%s\n'%list_pos_stop_codons_bis)

			else:
				out.write('%s\n'%'')


		else:
			sequences_CDS.append(CDS_seq)
			sequences_AA.append(AA_seq)


	out.close()
	out_bis.close()

	out_CDS = MA_file.split('.')[0] + '_wo_gaps.fasta'

	SeqIO.write(sequences_CDS, out_CDS, 'fasta')


	out_AA = MA_file.split('.')[0] + '_AA_wo_gaps.fasta'

	SeqIO.write(sequences_AA, out_AA, 'fasta')


MA_file = sys.argv[1]
path_MA = sys.argv[2]
Filter_CDS(MA_file, path_MA)

