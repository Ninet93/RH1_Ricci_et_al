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
from operator import itemgetter
from itertools import groupby
from textwrap import wrap


def Get_LT_cichlids_complete_ID_seq(FASTA, ID, path_data, path_out, OPSIN):
	os.chdir(path_data)

	record = SeqIO.to_dict(SeqIO.parse(FASTA, 'fasta'))

	sequences=[]

	for i in record.keys():
		if ID in i:
			sequences.append(record[i])

	os.chdir(path_out)
	#print(path_out)

	out_f = OPSIN + '_' + ID + '.fasta'
	SeqIO.write(sequences, out_f, 'fasta')




FASTA = sys.argv[1]
ID = sys.argv[2]
path_data = sys.argv[3]
path_out = sys.argv[4]
OPSIN = sys.argv[5]
Get_LT_cichlids_complete_ID_seq(FASTA, ID, path_data, path_out, OPSIN)

