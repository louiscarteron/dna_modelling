def decode_file():
  with open("input.txt") as f:
    oligos = f.read().splitlines()
  
  return oligos


def encode_file():
  pass

def pcr(oligos, cycles=10):
  pass

def synthesis(oligos, error_rate):
  pass

def storage(oligos, time=0):
  pass

def sequence(oligos, fwd_primer, rvs_primer):
  pass


def main():
  print("Decoding file...")
  raw_oligos = decode_file()


  print("Synthesising oligos...")
  syn_oligos = synthesis(raw_oligos, 0.1)


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