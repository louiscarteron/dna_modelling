from Bio import SeqIO
from statistics import mean
from Levenshtein import distance, editops
from csv import reader, DictReader
from collections import defaultdict, Counter
from munkres import Munkres, print_matrix
from itertools import islice
import json

#filename = '/mnt/c/Users/Louis/Documents/DnaModelling/42k_data/output_reads.fastq'
#filename = 'data/porechop/eins3_data_chopped.fastq'
filename = "data/42k2b/42k2b_porechop_nomiddle.fastq"

'''
Parses the input file 'filename' and filters it based on a length greater than 'length' and a quality score greater than 'score'
'''
def trim_inputs(minlength, maxlength, score):
  input_seq_iterator = SeqIO.parse(filename, "fastq")
  short_seq_iterator = (record for record in input_seq_iterator if len(record.seq) > minlength and len(record.seq) < maxlength and mean(record.letter_annotations['phred_quality']) > score)

  SeqIO.write(short_seq_iterator, "data/42k2b/trimmed/42k2b_porechop_nomiddle_q10.fastq", "fastq")

def trim_inputs_mod(length, score):
  input_seq_iterator = SeqIO.parse(filename, "fastq")

  temp = ("".join([o for i, o in enumerate(str(record.seq)) if record.letter_annotations['phred_quality'][i] > score]) for record in input_seq_iterator )
  short_temp = (record for record in temp if len(record) > length)

  #short_seq_iterator = (record for record in input_seq_iterator if len(record.seq) > length and mean(record.letter_annotations['phred_quality']) > score)
  return short_temp
  #SeqIO.write(short_temp, "data/porechop/test.fastq", "fastq")


def _read_csv(filepath):
  columns = []
  # Read csv and skip first column
  with open(filepath, "r") as f:
    columns = [row[1] for row in reader(f)]

  return columns

def _stripOligos(oligos):
  priming_5 = "CTACAACGCAGATTACAACCTCAGT"
  priming_3 = "CCATCCTTGCCAGCGTTACC"

  main_bp_length = 25

  stripped = [o[len(priming_5):][:main_bp_length] for o in oligos]

  return stripped

'''
def _read_csv(filepath):
  columns = defaultdict(list)
  # Read csv and skip first column
  with open(filepath, "r") as f:
    reader = DictReader(f)
    for row in reader:
      for (i, v) in enumerate(row):
        columns[i].append(v)

  return columns[1]
'''

def _read_input(filepath, filetype):
  input_seq_iterator = SeqIO.parse(filepath, filetype)

  trimmed_oligos = [str(record.seq) for record in input_seq_iterator]

  return trimmed_oligos

def match(input_oligos, read_oligos):
  report = []

  short = input_oligos[:100]

  for i in short:

    best_score = 1e4
    best_read = ""

    for r in read_oligos:
      score = distance(i, r)
      if score < best_score:
        best_score = score
        best_read = r

    edits = editops(i, best_read)
    counts = Counter(x[0] for x in edits)
    
    report.append({
      'input_oligo': i,
      'match': {
        'score': best_score,
        'closest_match': best_read,
        'edits': dict(counts)
      }
    })

  return report

def match2(input_oligos, read_oligos):
  report = []

  short = input_oligos[:100]

  for i in short:

    best_score = 1e4
    best_read = ""
    match25bp = []

    for r in read_oligos:
      score = distance(i, r)
      if score < best_score:
        best_score = score
        best_read = r

      if i in r:
        match25bp.append(r)

    edits = editops(i, best_read)
    counts = Counter(x[0] for x in edits)
    
    report.append({
      'input_oligo': i,
      'match': {
        'score': best_score,
        'closest_match': best_read,
        'edits': dict(counts)
      },
      '25bp': {
        'oligos': match25bp
      }
    })

  return report

class GeneratorLen(object):
  def __init__(self, gen, length):
    self.gen = gen
    self.length = length

  def __len__(self): 
    return self.length

  def __iter__(self):
    return self.gen

  def __getitem__(self, key):
    if isinstance(key, int) and key >= 0:
      return islice(self.gen, key, key + 1)
    elif isinstance(key, slice):
      return islice(self.gen, key.start, key.stop, key.step)
    else:
      raise KeyError(f'Key must be non-negative integer or slice, not {key}')

def match_try(input_oligos, read_oligos):
  distance_matrix = []

  print(len(read_oligos))

  short_input = input_oligos[:5]
  short_read = read_oligos[:20000]

  for i in short_input:
    distance_matrix.append(GeneratorLen((distance(i, r) for r in short_read), 20000))

  print("done")

  m = Munkres()
  indexes = m.compute(distance_matrix)

  print("done")

  report = []

  for row, column in indexes:
    io = short_input[row]
    ro = short_read[column]

    edits = editops(io, ro)
    counts = Counter(x[0] for x in edits)

    report.append({
      'input_oligo': io,
      'match': {
        'score': distance(io, ro),
        'closest_match': ro,
        'edits': dict(counts)
      }
    })
    
  return report


def process_report(report):
  temp = [r['match']['score'] for r in report]
  print(mean(temp))

def dump_report(report):
  with open("data/42k2b/reports/report_25bp_pc_nm_42k2b.json", "w+") as fp:
    json.dump(report, fp, indent = 2)

def main():
  #trim_inputs(90, 150, 10)
  #return
  input_oligos = _read_csv("data/3xr6.csv")
  temp = _stripOligos(input_oligos)
  read_oligos = _read_input("data/42k2b/trimmed/42k2b_porechop_nomiddle_q10.fastq", "fastq")
  #report = match(input_oligos, read_oligos)
  report = match2(temp, read_oligos)
  process_report(report)
  dump_report(report)
  

if __name__ == "__main__":
  main()
