from random import random

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