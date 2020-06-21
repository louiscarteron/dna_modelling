from random import random, randint
import math
from collections import Counter, Iterable
import data_config
import substitution
import time

debug = True

nuc2str = {
  0: 'A',
  1: 'G',
  2: 'T', 
  3: 'C'
}

str2nuc = {
  'A': 0,
  'G': 1,
  'T': 2,
  'C': 3,
}

base_oligo_length = 64

'''Temp function for testing encoding'''

def encode_binary(str):
  output = []
  prev = 0
  for i, s in enumerate(str):
    new = ''
    if s == '0':
      new = 'A' if i % 2 == 0 else 'T'
    if s == '1':
      new = 'G' if i % 2 == 0 else 'C'
    prev = s
    output.append(new)
  return ''.join(output)


def generateRandomSequence():
  rand_oligo = [[randint(0, 3) for i in range(140)] for j in range(5)]

  str_oligos = [''.join([nuc2str[s] for s in oligo]) for oligo in rand_oligo]

  with open("input.txt", "w") as f:
    for oligo in str_oligos:

      f.write("%s\n" % oligo)
    

def decode_file():
  with open("input.txt") as f:
    str_oligos = f.read().splitlines()
  

  oligos = [[str2nuc[s] for s in oligo] for oligo in str_oligos]

  return oligos 

def read_file(file_name):
  with open(file_name) as f:
    str_oligos = f.read().splitlines()
  

  oligos = [oligo for oligo in str_oligos]

  return oligos 


def encode_file(oligos):
  
  str_oligos = [''.join([nuc2str[s] for s in oligo]) for oligo in oligos]

  with open("output.txt", "w") as f:
    for oligo in str_oligos:

      f.write("%s\n" % oligo)

# TODO: Account for termination factor of 0.05% at every iteration. Start of inner loop. 
def synthesis(oligos, method):
  syn_oligos = []

  sub_rate, ins_rate, del_rate = method.values()

  total_sub = 0
  total_ins = 0
  total_del = 0
  
  # Going through each oligo to synthesis
  for i, oligo in enumerate(oligos):

    sub_counter = 0
    ins_counter = 0
    del_counter = 0

    new_oligo = []

    for nuc in oligo:

      # Early termination 
      if random() < 0.0005:
        break

      new_nuc = nuc

      if random() < sub_rate:
        new_nuc = substitution.getSubNucleotide(nuc, 'average')

        sub_counter += 1

      elif random() < ins_rate:
        new_oligo.append(0) # Adding new nucleotide. 
        ins_counter += 1

      elif random() < del_rate:
        del_counter += 1
        # About to delete the nucleotide, so we can just ignore it and not add it to the new oligo
        continue
      
      new_oligo.append(new_nuc)
    '''
    if debug:

      print(f'Sub counter for olgio {i}: {sub_counter}')
      print(f'Ins counter for olgio {i}: {ins_counter}')
      print(f'Del counter for olgio {i}: {del_counter}')
    '''
    syn_oligos.append(new_oligo)

    # Used to track metrics 
    total_sub += sub_counter
    total_ins += ins_counter
    total_del += del_counter

  if debug:

    number_of_nuc = base_oligo_length * len(oligos)

    print(f'Total sub events: {total_sub}. Proportion is: {total_sub / number_of_nuc}. Target is: {sub_rate}')
    print(f'Total ins events: {total_ins}. Proportion is: {total_ins / number_of_nuc}. Target is: {ins_rate}')
    print(f'Total del events: {total_del}. Proportion is: {total_del / number_of_nuc}. Target is: {del_rate}')

  return syn_oligos, total_sub


# Data taken from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169774
# Using data from Pacific Biosciences RSII

def pcr(oligos, cycles=10):

  pcr_oligos = []

  total_sub = 0
  total_ins = 0
  total_del = 0

  for i in range(cycles):
    
    for oligo in oligos:

      new_oligo = []

      sub_counter = 0
      ins_counter = 0
      del_counter = 0

      for nuc in oligo: 

        rand = random()

        new_nuc = nuc

        # Error happened
        if (rand < 1.8e-4):

          error_type = random()

          # Substitution error
          if (error_type < 0.973):
            new_nuc = substitution.getSubNucleotideForPCR(nuc)
            sub_counter += 1
          # Deletion error
          elif (error_type < 0.99):
            del_counter += 1
            continue
          
          # Addition error
          else:
            ins_counter += 1

        new_oligo.append(new_nuc)

      pcr_oligos.append(new_oligo)

      total_sub += sub_counter
      total_ins += ins_counter
      total_del += del_counter

  if debug:

    number_of_nuc = base_oligo_length * len(oligos)

    print(f'Total sub events: {total_sub}. Proportion is: {total_sub / number_of_nuc}')
    print(f'Total ins events: {total_ins}. Proportion is: {total_ins / number_of_nuc}')
    print(f'Total del events: {total_del}. Proportion is: {total_del / number_of_nuc}')


  return pcr_oligos




# Validate an oligo to check if valid G-C content and correct homopolymers. 
def validateOligo(oligo):
  gcCount = 0

  nucHistory = [] # Need to make max size of 4 (max of 3 homopolymers)

  for o in oligo:

    if o == 1 or o == 3:
      gcCount += 1

    if len(nucHistory) == 4:
      nucHistory.pop(0)
    
    nucHistory.append(o)

    if len(nucHistory) == 4 and nucHistory.count(nucHistory[0]) == len(nucHistory):
      return False

  return (gcCount / len(oligo) > 0.45) and (gcCount / len(oligo) < 0.55)

# Generates a legal Oligo that satisfies 45-55% G-C content and no more than 3 homopolymers.
def generateLegalOligo(length):

  isValid = False
  oligo = []

  #tempCount = 0

  while(not isValid):
    oligo = [randint(0, 3) for i in range(length)]
    isValid = validateOligo(oligo)

    #tempCount += 1
  
  #print(f'Had to iterate {tempCount} times to get a valid oligo')

  return oligo

def generateRandomOligos(length, total):
  oligos = [generateLegalOligo(base_oligo_length) for i in range(total)]

  str_oligos = [''.join([nuc2str[s] for s in oligo]) for oligo in oligos]

  with open("input2.txt", "w") as f:
    for oligo in str_oligos:

      f.write("%s\n" % oligo)

  return oligos



def monteCarloDecaySimulation(oligos, decayEvents):

  for i in range(decayEvents):

    decayed_oligo_index = randint(0, len(oligos))

    decayed_oligo = oligos[decayed_oligo_index]
    
    # TODO: Need to account for broken strands 
    if len(decayed_oligo) > 10:
      break_point = randint(1, base_oligo_length - 1)
      new_oligo = [decayed_oligo[:break_point], decayed_oligo[break_point:]]

      oligos[decayed_oligo_index] = new_oligo

    else:
      print("Should have broken an already broken strand")


# Time is in years to keep units consistent. 
# Breaks a strand in 2, meaning it can't be applified in the PCR stage. 
#https://gyazo.com/53e1f6eaa1d4e2166f1ae85279cdb7fd
def storage(oligos, time, redundancy):
  half_life = 521

  dt = 1

  # Probability of one nucleotide decaying in a time frame (1 year)
  r = math.log(2) / 521 

  k = math.exp(41.2 - 15267.6 * (1/283.15))

  final_oligos = oligos.copy()

  t = 0
  decay = 0

  while t < time:

    rnd = random()

    if rnd < k:
      # 5 is number of decay events to be simulated. Can be changed to fit different environments (hot, cold etc)
      monteCarloDecaySimulation(final_oligos, 5)

      decay += 1

    t += dt

  if debug:
    print(f'Number of decay events: {decay}')

  return final_oligos

def sequence(oligos, method):
  seq_oligos = []

  sub_rate, ins_rate, del_rate = method.values()

  total_sub = 0
  total_ins = 0
  total_del = 0
  
  # Going through each oligo to synthesis
  for i, oligo in enumerate(oligos):

    sub_counter = 0
    ins_counter = 0
    del_counter = 0

    new_oligo = []

    for nuc in oligo:

      new_nuc = nuc

      if random() < sub_rate:
        new_nuc = substitution.getSubNucleotide(nuc, 'average')

        sub_counter += 1

      elif random() < ins_rate:
        new_oligo.append(0) # Adding new nucleotide.
        ins_counter += 1

      elif random() < del_rate:
        del_counter += 1
        # About to delete the nucleotide, so we can just ignore it and not add it to the new oligo
        continue
      
      new_oligo.append(new_nuc)
    '''
    if debug:

      print(f'Sub counter for olgio {i}: {sub_counter}')
      print(f'Ins counter for olgio {i}: {ins_counter}')
      print(f'Del counter for olgio {i}: {del_counter}')
    '''
    seq_oligos.append(new_oligo)

    # Used to track metrics 
    total_sub += sub_counter
    total_ins += ins_counter
    total_del += del_counter

  if debug:

    number_of_nuc = base_oligo_length * len(oligos)

    print(f'Total sub events: {total_sub}. Proportion is: {total_sub / number_of_nuc}. Target is: {sub_rate}')
    print(f'Total ins events: {total_ins}. Proportion is: {total_ins / number_of_nuc}. Target is: {ins_rate}')
    print(f'Total del events: {total_del}. Proportion is: {total_del / number_of_nuc}. Target is: {del_rate}')

  return seq_oligos


def getOligoLengths(item):
  if (any(isinstance(i, list) for i in item)):
    return [getOligoLengths(i) for i in item]
  else:
    return len(item)

# Helper to flatten a list [a, [b, c]] to [a, b, c]
def flatten(x):
  if isinstance(x, Iterable):
      return [a for i in x for a in flatten(i)]
  else:
      return [x]

# Helper to get frequency of each oligo length after storage. 
def getDecayInformation(oligos):

  lengths = flatten(getOligoLengths(oligos))

  freq = {}

  for i in set(lengths):
    freq[i] = lengths.count(i)

  print(freq)

def test():

  freq = {}

  temp_array = [substitution.getSubNucleotide(3, 'high_pr4t') for i in range(0, 100)]

  for i in set(temp_array):
    freq[i] = temp_array.count(i)

  print(freq)

def count_errors(inputs, outputs):
  '''
  err_count = 0
  for i in range(len(inputs)):
    if inputs[i] != outputs[i]:
      err_count += 1

  return err_count
  '''
  temp = []
  for i in range(len(inputs)):
    c = Counter(inputs[i])
    c.subtract(Counter(outputs[i]))
    temp.append(c)
  
  return temp

def analyse_output(outputs):
  lens = [len(o) for o in outputs]
  print(sum(lens) / len(lens))
  counter = Counter(lens)
  print(dict(counter))

def main():

  #st = "0110111001100001011011100110111101110000011011110111001001100101"


  print("Decoding file...")
  raw_oligos = decode_file()

  print("Synthesising oligos...")
  syn_oligos, _ = synthesis(raw_oligos, data_config.average)
  
  print("Simulating storage...")
  stg_oligos = storage(syn_oligos, 521, 1)

  getDecayInformation(stg_oligos)
  #return 

  '''
  print("Applying PCR...")
  #pcr_oligos = pcr(stg_oligos)
  '''

  pcr_oligos = stg_oligos

  print("Sequencing stored oligos...")
  #seq_oligos = sequence(pcr_oligos)
  seq_oligos = sequence(stg_oligos, data_config.sequencing)

  print("Finished")

  encode_file(seq_oligos)

  err = count_errors(raw_oligos, seq_oligos)
  print(err)



if __name__ == "__main__":
  start_time = time.time()
  main()
  print(f'--- {time.time() - start_time} seconds ---')