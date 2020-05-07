import numpy as np
from uuid import uuid4 
from fast5_research import Fast5, BulkFast5, fast5

from ont_fast5_api.fast5_interface import get_fast5_file

from Bio import SeqIO

filename = 'data/fq200420_1.fasta'

records = list(SeqIO.parse(filename, "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record


record_dict = SeqIO.index(filename, "fasta")
print(record_dict[records[0].id])

filename2 = 'data/fastq_runid_191d5e13ace2d2025598d80406d55d220bb411df_0_0.fastq'

records2 = list(SeqIO.parse(filename2, "fastq"))
print(records2[0].id)  # first record
print(records2[-1].id)  # last record

record_dict = SeqIO.index(filename2, "fastq")
print(record_dict[records2[0].id])

'''
filename = 'temp123.fast5'

with get_fast5_file(filename, mode="r") as f5:
  for read in f5.get_reads():
    raw_data = read.get_raw_data()
    print(read.read_id, raw_data)
'''
'''
filename2 = 'temp123.fast5'

with Fast5(filename2) as fh:
    raw = fh.get_read(raw=True)
    summary = fh.summary()
print('Raw is {} samples long.'.format(len(raw)))
print('Summary {}.'.format(summary))
'''