from random import random, randint
import math
import collections

'''
high_pr = {
  'sub_rate': 0.005,
  'add_rate': ,
  'del_rate': ,
}

low_pr = {
  'sub_rate': ,
  'add_rate': ,
  'del_rate': ,
}

erlich = {
  'sub_rate': ,
  'add_rate': ,
  'del_rate': ,
}

goldman = {
  'sub_rate': ,
  'add_rate': ,
  'del_rate': ,
}

high_pr4t = {
  'sub_rate': ,
  'add_rate': ,
  'del_rate': ,
}
'''


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

base_oligo_length = 113

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


def encode_file(oligos):
  
  str_oligos = [''.join([nuc2str[s] for s in oligo]) for oligo in oligos]

  with open("output.txt", "w") as f:
    for oligo in str_oligos:

      f.write("%s\n" % oligo)


# Function to get the substituting nucleotide based on probabilities.
# TODO: Check if correct from the paper. 
def getSubNucleotide(nuc):
  rand = random()

  # TODO: Figure out how to get these mappings to work correctly. 
  # TODO: Clean up this mess of a code and make sure the numbers are correct. 

  # A2C: 0.05 -> 0.25
  # A2G: 0.1  -> 0.5
  # A2T: 0.05 -> 0.25

  # T2G: 0.05  -> 0.28
  # T2C: 0.1   -> 0.57
  # T2A: 0.025 -> 0.15

  # C2A: 0.07 -> 0.175
  # C2G: 0.03 -> 0.075
  # C2T: 0.3  -> 0.75

  # G2C: 0.02 -> 0.07
  # G2A: 0.2  -> 0.67
  # G2T: 0.08 -> 0.26

  if nuc == 0:
    if rand < 0.25:
      return 3
    elif rand < 0.75: 
      return 1
    else:
      return 2

  if nuc == 1:
    if rand < 0.07:
      return 3
    elif rand < 0.74:
      return 0
    else:
      return 2

  if nuc == 2:
    if rand < 0.28:
      return 1
    elif rand < 0.85:
      return 3
    else:
      return 0

  if nuc == 3:
    if rand < 0.175:
      return 0
    elif rand < 0.25:
      return 1
    else:
      return 2

# TODO: Account for termination factor of 0.05% at every iteration. Start of inner loop. 
def synthesis(oligos, del_rate, sub_rate, add_rate):
  syn_oligos = []
  
  # Going through each oligo to synthesis
  for i, oligo in enumerate(oligos):

    sub_counter = 0
    add_counter = 0
    del_counter = 0

    new_oligo = []
    '''
    Need to go through each nucleotide and decide what to do with them. 
    TODO: encode conditional probability for substitution and addition.
    Defaulting to A (0) for now. 
    '''
    for nuc in oligo:

      new_nuc = nuc

      if random() < sub_rate:
        new_nuc = getSubNucleotide(nuc)

        sub_counter += 1

      elif random() < add_rate:
        new_oligo.append(0) # Adding new nucleotide. TODO: need to choose which one to add from condition probability distrubution. 
        add_counter += 1

      elif random() < del_rate:
        del_counter += 1
        continue
      
      new_oligo.append(new_nuc)

    print(f'Sub counter for olgio {i}: {sub_counter}')
    print(f'Add counter for olgio {i}: {add_counter}')
    print(f'Del counter for olgio {i}: {del_counter}')
    syn_oligos.append(new_oligo)

  return syn_oligos


def pcr(oligos, cycles=10):
  pass


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

  with open("input.txt", "w") as f:
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
def storage(oligos, time, redundancy):
  half_life = 521

  dt = 1

  # Probability of one nucleotide decaying in a time frame (1 year)
  r = math.log(2) / 521 

  final_oligos = oligos.copy()

  t = 0
  decay = 0

  while t < time:

    rnd = random()

    if rnd < r:
      # 5 is number of decay events to be simulated. Can be changed to fit different environments (hot, cold etc)
      monteCarloDecaySimulation(final_oligos, 5)

      decay += 1

    t += dt

  print(decay)

  return final_oligos

def sequence(oligos):
  pass


def getOligoLengths(item):
  if (any(isinstance(i, list) for i in item)):
    return [getOligoLengths(i) for i in item]
  else:
    return len(item)

# Helper to flatten a list [a, [b, c]] to [a, b, c]
def flatten(x):
  if isinstance(x, collections.Iterable):
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

def main():

  # Used to generate 100 random oligos 
  #oligos = generateRandomOligos(base_oligo_length, 100)
  #print(len(oligos))


  print("Decoding file...")
  raw_oligos = decode_file()

  print("Synthesising oligos...")
  syn_oligos = synthesis(raw_oligos, 0.01, 0.01, 0.01)

  
  print("Simulating storage...")
  stg_oligos = storage(syn_oligos, 521, 1)

  getDecayInformation(stg_oligos)
  return 

  print("Applying PCR...")
  pcr_oligos = pcr(stg_oligos)


  print("Sequencing stored oligos...")
  fwd_primer = "ATG"
  rvs_primer = "GCA"
  seq_oligos = sequence(pcr_oligos, fwd_primer, rvs_primer)

  print("Finished")



if __name__ == "__main__":
  main()