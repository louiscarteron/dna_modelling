from Bio import SeqIO
from statistics import mean
from Levenshtein import distance, editops
from csv import reader, DictReader
from collections import defaultdict, Counter
from munkres import Munkres, print_matrix
from itertools import islice
from fuzzywuzzy import fuzz
import math
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

def match25bp(input_oligos, read_oligos):
  report = []

  reads = read_oligos.copy()

  short = input_oligos[:420]

  for i in input_oligos:

    #match25bp = [r for r in reads if i in r or isCloseSubstring(r, i, 90)]

    match25bp = [r for r in reads if i in r]

    #closeMatch = [r for r in reads if isCloseSubstring(r, i, 0.90)]

    closeMatch = []

    for m in match25bp:
      reads.remove(m)
    
    report.append({
      'input_oligo': i,
      'match': match25bp,
      'closeMatch': closeMatch
    })

  return report

def isCloseSubstring(mainstring, substring, tol):

  return fuzz.partial_ratio(mainstring, substring) > tol

  '''
  window_size = len(substring)
  search_space = [mainstring[i:i + window_size] for i in range(len(mainstring) - window_size)]

  for s in search_space:
    if distance(substring, s) < tol:
      return True

  return False
  '''

def process_25bp(report):
  forward = "CTACAACGCAGATTACAACCTCAGT"
  backwards = "CCATCCTTGCCAGCGTTACC"
  input_size = 420
  total_matches = sum([len(r['match']) for r in report])
  print(f'Total Matches: {total_matches}. Overall percentage: {total_matches * 100 / input_size}%')

  close_matches = sum([len(r['closeMatch']) for r in report])
  print(f'Close Matches: {close_matches}. Overall percentage: {close_matches * 100 / input_size}%')


def process_25bp_report(report):

  new_report = []

  for r in report:
    matches = r.get('match', [])

    if len(matches) == 0:
      continue

    in_oligo = r['input_oligo']
    splits = [align_25bp(in_oligo, m) for m in matches]

    splits_info = []

    for s in splits:
      for t in s:
        edits = editops(t, in_oligo)
        counts = Counter(x[0] for x in edits)
        splits_info.append({
          'strand': t,
          'editops': dict(counts)
        })

    new_report.append({
      'input_oligo': in_oligo,
      'splits': splits_info
    })
  
  return new_report



def align_25bp(true25bp, match):
  index = match.find(true25bp)
  bp_length = len(true25bp)
  max_len = len(match)
  search_range = math.ceil(index / bp_length)
  min_start = index - (math.floor(index / bp_length) * bp_length)

  indexes_to_cut_at = [i for i in range(min_start, max_len, bp_length)]

  lines = [match[i: i+bp_length] for i in indexes_to_cut_at]
  lines.insert(0, match[0:min_start]) #insert missing start in case it is needed

  return lines


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

def dump_report(report, filepath):
  with open(filepath, "w+") as fp:
    json.dump(report, fp, indent = 2)

def read_json(filepath):
  with open(filepath) as fp:
    data = json.load(fp)
  
  return data

def main():
  report = read_json("data/42k2b/reports/report_25bp_exact_42k2b.json")
  new_report = process_25bp_report(report)
  dump_report(new_report, "data/42k2b/reports/test.json")




  return 
  #trim_inputs(90, 150, 10)
  #return
  input_oligos = _read_csv("data/3xr6.csv")
  temp = _stripOligos(input_oligos)
  read_oligos = _read_input("data/42k2b/trimmed/42k2b_porechop_nomiddle_q10.fastq", "fastq")
  #report = match(input_oligos, read_oligos)
  report = match25bp(temp, read_oligos)
  process_25bp(report)
  dump_report(report)
  

if __name__ == "__main__":
  main()
