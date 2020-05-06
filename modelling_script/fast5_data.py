import numpy as np
from uuid import uuid4 
from fast5_research import Fast5, BulkFast5, fast5

from ont_fast5_api.fast5_interface import get_fast5_file

filename = 'FAN19880_d8b24b7f_0.fast5'

with get_fast5_file(filename, mode="r") as f5:
  for read in f5.get_reads():
    raw_data = read.get_raw_data()
    print(read.read_id, raw_data)

'''
with Fast5(filename) as fh:
    raw = fh.get_read(raw=True)
    summary = fh.summary()
print('Raw is {} samples long.'.format(len(raw)))
print('Summary {}.'.format(summary))
'''