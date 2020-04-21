from random import random, randint # soloq teammates btw 

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

def pcr(oligos, cycles=10):
  pass

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

def storage(oligos, time=0, sub_rate):
  pass

def sequence(oligos):
  pass


def main():
  #generateRandomSequence()

  print("Decoding file...")
  raw_oligos = decode_file()

  print("Synthesising oligos...")
  syn_oligos = synthesis(raw_oligos, 0.01, 0.01, 0.01)
  
  encode_file(syn_oligos)

  return

  print("Simulating storage...")
  str_oligos = storage(syn_oligos, 100)


  print("Applying PCR...")
  pcr_oligos = pcr(str_oligos)


  print("Sequencing stored oligos...")
  fwd_primer = "ATG"
  rvs_primer = "GCA"
  seq_oligos = sequence(pcr_oligos, fwd_primer, rvs_primer)

  print("Finished")



if __name__ == "__main__":
  main()