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

def _read_input(filepath, filetype):
  input_seq_iterator = SeqIO.parse(filepath, filetype)

  trimmed_oligos = [str(record.seq) for record in input_seq_iterator]

  return trimmed_oligos

def _round_up(i):
  return int(round(i + 5.1, -1))

def _get_length_distribution(filepath, filetype):
  input_seq_iterator = SeqIO.parse(filename, filetype)
  length_iterator = (_round_up(len(str(record.seq))) for record in input_seq_iterator)
  return Counter(length_iterator)

def _get_score_distribution(filepath, filetype):
  input_seq_iterator = SeqIO.parse(filename, filetype)
  score_iterator = (round(mean(record.letter_annotations['phred_quality'])) for record in input_seq_iterator)
  return Counter(score_iterator)
  



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

  short = input_oligos[:100]

  for i in short:

    #match25bp = [r for r in reads if i in r or isCloseSubstring(r, i, 90)]

    matched25bp = [r for r in reads if i in r]

    closeMatch = [r for r in reads if len(r) > 90 and isCloseSubstring(i, r, 90)]

    #closeMatch = []
    
    for m in matched25bp:
      reads.remove(m)
    '''
    for c in closeMatch:
      reads.remove(c)
    '''
    report.append({
      'input_oligo': i,
      'match': matched25bp,
      'closeMatch': closeMatch
    })

  return report

def isCloseSubstring(mainstring, substring, tol):

  return fuzz.partial_ratio(mainstring, substring) > tol

def process_25bp(report):
  forward = "CTACAACGCAGATTACAACCTCAGT"
  backwards = "CCATCCTTGCCAGCGTTACC"
  input_size = 42000
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
    splits_with_idx = [align_25bp(in_oligo, m) for m in matches]

    splits, split_idxes = [list(t) for t in zip(*splits_with_idx)]

    match_info = []

    for idx, val in enumerate(matches):
      curr_splits = splits[idx]
      match_index = split_idxes[idx]
      strand_ops = []

      seq_block = [100 for i in range(len(curr_splits))]

      for i, s in enumerate(curr_splits):

        #edits = editops(s, in_oligo)
        edits = editops(in_oligo, s)
        counts = Counter(x[0] for x in edits)
        strand_ops.append({
          'strand': s,
          'editops': dict(counts)
        })

        if len(s) == 25:
          seq_block[i] = sum(list(dict(counts).values()))
          if sum(list(dict(counts).values())) == 0:
            seq_block[i] = 0


      start_block = find_smallest_sum(seq_block)
      
      if abs(start_block - match_index) > 2:
        print(val)
        print(strand_ops)
      
      match_info.append({
        "sequence": val,
        "seq_block_start": start_block,
        "splits": strand_ops
      })

    new_report.append({
      'input_oligo': in_oligo,
      'matches': match_info
    })
  
  return new_report

def find_smallest_sum(arr):
  size = 3
  idx = 0
  min_sum = 1000
  for i in range(len(arr)-size + 1):
    curr_sum = sum(arr[i:i+size])
    if curr_sum < min_sum:
      min_sum = curr_sum
      idx = i
  return idx

def align_25bp(true25bp, match):
  index = match.find(true25bp)
  bp_length = len(true25bp)
  max_len = len(match)
  search_range = math.ceil(index / bp_length)
  min_start = index - (math.floor(index / bp_length) * bp_length)

  indexes_to_cut_at = [i for i in range(min_start, max_len, bp_length)]

  lines = [match[i: i+bp_length] for i in indexes_to_cut_at]
  lines.insert(0, match[0:min_start]) #insert missing start in case it is needed

  return lines, lines.index(true25bp)

def print_table_from_report(report, print_results=True):
  total_inputs = len(report)

  bp_length = len(report[0].get('input_oligo'))

  sub_count = 0
  ins_count = 0
  del_count = 0
  total_nucleotides = 0

  sub_info = {'A': {'A': 0, 'G': 0, 'C': 0, 'T': 0},
   'G': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'C': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'T': {'A': 0, 'G': 0, 'C': 0, 'T': 0}
  }

  del_info = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

  ins_info = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

  ins_after_info = {'A': {'A': 0, 'G': 0, 'C': 0, 'T': 0},
   'G': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'C': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'T': {'A': 0, 'G': 0, 'C': 0, 'T': 0},
   'start': {'A': 0, 'G': 0, 'C': 0, 'T': 0}
  }

  del_after_info = {'A': {'A': 0, 'G': 0, 'C': 0, 'T': 0},
   'G': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'C': {'A': 0, 'G': 0, 'C': 0, 'T': 0}, 
   'T': {'A': 0, 'G': 0, 'C': 0, 'T': 0},
   'start': {'A': 0, 'G': 0, 'C': 0, 'T': 0}
  }

  nucleotide_distribution = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

  read_distribution = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

  total_input_25_length = 0

  # Ugly for now
  temp = {
    0: {
    'ins': 0,
    'sub': 0,
    'del': 0
    },
    1: {
    'ins': 0,
    'sub': 0,
    'del': 0
    },
    2: {
    'ins': 0,
    'sub': 0,
    'del': 0
    },
  }

  temparr = []

  te123 = {0: 0, 1:0, 2:0}

  for r in report:
    matches = r.get('matches', [])
    input_oligo = r.get('input_oligo', "")

    for m in matches:
      splits = m.get('splits', [])
      
      seq_block_start = m.get('seq_block_start', 1) # Default to 1 in case

      seq_count = 0

      for i, s in enumerate(splits):

        if i < seq_block_start or i - 2 > seq_block_start:
          continue

        current_strand = s.get('strand', "")

        total_input_25_length += len(input_oligo)
        input_dist = Counter(input_oligo)
        
        # Populate input distributions
        for (key, val) in input_dist.items():
          nucleotide_distribution[key] += val


        #replacements = process_replacements(current_strand, input_oligo)
        replacements = process_replacements(input_oligo, current_strand)

        for rep in replacements:

          rep_type = rep.get('type')

          if rep_type == 'insert':
            inse = rep.get('modification').get('insert')
            prev = rep.get('modification').get('after')
            ins_info[inse] += 1
            ins_after_info[prev][inse] += 1

          elif rep_type == 'replace':
            og = rep.get('modification').get('original')
            nw_olg = rep.get('modification').get('replacement')
            sub_info[og][nw_olg] += 1

          else:
            dele = rep.get('modification').get('delete')
            prev = rep.get('modification').get('after')
            temparr.append(rep.get('seq3'))
            del_info[dele] += 1
            del_after_info[prev][dele] += 1

        strand_dist = Counter(current_strand)
        for (key, val) in strand_dist.items():
          read_distribution[key] += val
        
        total_nucleotides += len(current_strand)
        
        editops = s.get('editops', {})

        if editops == {}:
          te123[seq_count] += 1

        sub_count += editops.get('replace', 0)
        ins_count += editops.get('insert', 0)
        del_count += editops.get('delete', 0)
  
        temp[seq_count]['ins'] += editops.get('insert', 0)
        temp[seq_count]['del'] += editops.get('delete', 0)
        temp[seq_count]['sub'] += editops.get('replace', 0)

        seq_count += 1

        

  total_errors = sub_count + ins_count + del_count

  if print_results:

    oli_keys = ["A", "G", "C", "T"]

    print(f'Nucleotide distribution in inputs')
    for oli in oli_keys:
      print(f'  {oli}: {nucleotide_distribution[oli]} ({nucleotide_distribution[oli] * 100 / total_input_25_length:.2f}%)')

    print("")
    print(f'Nucleotide distribution in reads')
    for oli in oli_keys:
      print(f'  {oli}: {read_distribution[oli]} ({read_distribution[oli] * 100/ total_nucleotides:.2f}%)')

    print("\n")
    print(f'Total errors: {total_errors} ({total_errors * 100/total_nucleotides:.2f}%)')
    print('Breakdown as follows:')
    print(f'Substitutions: {sub_count} ({sub_count * 100 / total_errors:.2f}%)')
    for (key, val) in sub_info.items():
    
      total_for_sub = sum(val.values())
      print(f'  Original nucleotide {key} : {total_for_sub} ({total_for_sub * 100 / sub_count:.2f}%)')
      for (key, val) in val.items():
        print(f'    {key}: {val} ({val * 100/ total_for_sub:.2f}%)')


    print("\n")
    print(f'Insertions: {ins_count} ({ins_count * 100 / total_errors:.2f}%)')

    for oli in oli_keys:
      print(f'  {oli}: {ins_info[oli]} ({ins_info[oli] * 100/ ins_count:.2f}%)')
    
    for (key, val) in ins_after_info.items():
    
      total_for_ins = sum(val.values())
      print(f'  Previous nucleotide {key} : {total_for_ins} ({total_for_ins * 100 / ins_count:.2f}%)')
      for (key, val) in val.items():
        print(f'    {key}: {val} ({val * 100/ total_for_ins:.2f}%)')
    
    print("\n")
    print(f'Deletions: {del_count} ({del_count * 100 / total_errors:.2f}%)')
    for oli in oli_keys:
      print(f'  {oli}: {del_info[oli]} ({del_info[oli] * 100/ del_count:.2f}%)')

    #print(Counter(temparr))
    
    for (key, val) in del_after_info.items():
    
      total_for_del = sum(val.values())
      print(f'  Previous nucleotide {key} : {total_for_del} ({total_for_del * 100 / del_count:.2f}%)')
      for (key, val) in val.items():
        print(f'    {key}: {val} ({val * 100/ total_for_del:.2f}%)')
    
    print("\n")
    for (key, v) in temp.items():
      print(f'Errors in strand {key}')
      for (k, val) in v.items():
        print(f'{k}: {val}')

    print(te123)

  return del_info, ins_info, sub_info, nucleotide_distribution, read_distribution


def process_replacements(str1, str2):
  edits = editops(str1, str2)
  temp = []

  for e in edits:

    (op, spos, dpos) = e

    if op == 'replace':
      info = {'original': str1[spos], 'replacement': str2[dpos]}
      temp.append({
        'type': op,
        'modification': info
      })

    elif op == 'insert':
      after = str2[dpos - 1] if dpos > 0 else 'start'
      info = {'insert': str2[dpos], 'after': after}
      temp.append({
        'type': op,
        'modification': info
      })

    else:
      after = str1[spos - 1] if spos > 0 else 'start'
      info = {'delete': str1[spos], 'after': after}

      seq3 = str1[spos - 1: spos + 2] if spos > 0 else 'start'

      temp.append({
        'type': op,
        'modification': info,
        'seq3': seq3
      })

  return temp
      
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
  '''
  report = read_json("data/flowcell/report/full_flowcell1.json")
  process_25bp(report)
  new_report = process_25bp_report(report)
  #dump_report(new_report, "data/flowcell/report/full_flowcell3_errors.json")

  print_table_from_report(new_report)

  return
  '''
  

  #trim_inputs(90, 150, 10)
  #return
  input_oligos = _read_csv("data/3xr6.csv")
  temp = _stripOligos(input_oligos)
  #read_oligos = _read_input("data/42k2b/trimmed/42k2b_porechop_nomiddle_q10.fastq", "fastq")
  read_oligos = _read_input("data/flowcell/data/pc_nm_flowcell_2.fasta", "fasta")
  #report = match(input_oligos, read_oligos)
  report = match25bp(temp, read_oligos)
  process_25bp(report)
  dump_report(report, "data/flowcell/report/flowcell2_fuzzy.json")
  

if __name__ == "__main__":
  main()
