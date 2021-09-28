import os, sys
import re
from collections import OrderedDict
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import Bio.Align
from operator import itemgetter
from itertools import groupby
from textwrap import wrap

def Extract_CDS(mRNA_file, COO, path):
	os.chdir(path)


	coo = pd.read_csv(COO)
	coo.columns=['Opsin', 'Gene', 'Protein', 'Type', 'Scf_start', 'Scf_stop', 'From_mRNA_start', 'From_mRNA_stop']
	# coo starts with 1
	# get right sequences with python >> [start : stop + 1]

	coo['Scf_start'] = coo['Scf_start'].astype(int)
	coo['Scf_stop'] = coo['Scf_stop'].astype(int)
	coo['From_mRNA_start'] = coo['From_mRNA_start'].astype(int)
	coo['From_mRNA_stop'] = coo['From_mRNA_stop'].astype(int)

	list_Opsin = coo['Opsin'].unique().tolist()
	list_Gene = coo['Gene'].unique().tolist()
	list_Protein = coo['Protein'].unique().tolist()
	list_Type = coo['Type'].unique().tolist()


	os.chdir(path)
	align = AlignIO.read(mRNA_file, 'fasta')

	len_align = align.get_alignment_length() # get alignment length

	count=-1
	for record in align:
		count += 1

		if record.id in list_Gene:
			ref = record.id

			ref_count = count

			ref_CDS = coo[(coo['Gene'] == ref) & (coo['Type'] == 'CDS')]



	ref_seq = align[ref_count].seq # get ref sequence
	ref_seq_wo_gap = str(ref_seq).replace('-', '') # get ref sequence without gaps >> initial ref sequence
	len_ref_seq = len(ref_seq_wo_gap) # length initial ref sequence



	# extract MA block when alignment with RefSeq (when there is no gap)
	count=0
	for p in list(range(0, len_align)):
		if align[ref_count, p] != '-':
			count += 1

			if count == 1:
				block = align[:, p:p+1]

			else:
				block += align[:, p:p+1]

	# extract MA block when alignment with CDS of RefSeq >> allowing for translation into amino acids
	list_SeqRecord_bis=[]

	count_bis=-1
	for record in block:
		no_rev_complement=''
		count_bis += 1

		count=0
		for index, row in ref_CDS.iterrows():
			count += 1

			start = row['From_mRNA_start']
			stop = row['From_mRNA_stop']

			if count == 1:

				no_rev_complement = block[count_bis, start:stop+1].seq.reverse_complement()
				#print('1')
				#print(no_rev_complement)
				#print('-----')

			else:
				no_rev_complement += block[count_bis, start:stop+1].seq.reverse_complement()
				#print('2')
				#print(no_rev_complement)
				#print('-----')

		SeqRec_bis=SeqRecord(seq=no_rev_complement, id=block[count_bis].id, name=block[count_bis].name, description=block[count_bis].description)

		list_SeqRecord_bis.append(SeqRec_bis)

	extract = MultipleSeqAlignment(list_SeqRecord_bis)
	# extracted CDS MA with all sequences as the initial file - reverse complement of Onil scaffold - starts with ATG...


	AlignIO.write(extract, mRNA_file.split('.')[0] + '_CDS.fasta', 'fasta')



mRNA_file = sys.argv[1]
COO = sys.argv[2]
path = sys.argv[3]
Extract_CDS(mRNA_file, COO, path)

