from Bio import SeqIO
from statistics import mean
from Levenshtein import distance
from csv import reader, DictReader
from collections import defaultdict
import json

filename = '/mnt/c/Users/Louis/Documents/DnaModelling/42k_data/output_reads.fastq'

'''
Parses the input file 'filename' and filters it based on a length greater than 'length' and a quality score greater than 'score'
'''
def trim_inputs(length, score):
  input_seq_iterator = SeqIO.parse(filename, "fastq")
  short_seq_iterator = (record for record in input_seq_iterator if len(record.seq) > length and mean(record.letter_annotations['phred_quality']) > score)

  SeqIO.write(short_seq_iterator, "quality_seq_output_qscore13.fastq", "fastq")


def _read_csv(filepath):
  columns = []
  # Read csv and skip first column
  with open(filepath, "r") as f:
    columns = [row[1] for row in reader(f)]

  return columns

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

def _read_input(filepath):
  input_seq_iterator = SeqIO.parse(filepath, "fastq")
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
    
    report.append({
      'input_oligo': i,
      'match': {
        'score': best_score,
        'closest_match': best_read
      }
    })

  return report


def process_report(report):
  temp = [r['match']['score'] for r in report]
  print(mean(temp))

def dump_report(report):
  with open("report.json", "w+") as fp:
    json.dump(report, fp, indent = 2)

def main():
  #trim_inputs(90, 13)
  input_oligos = _read_csv("data/3xr6.csv")
  read_oligos = _read_input("quality_seq_output_qscore13.fastq")
  report = match(input_oligos, read_oligos)
  process_report(report)
  dump_report(report)
  

if __name__ == "__main__":
  main()
