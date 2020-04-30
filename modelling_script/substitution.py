from random import random

'''
str2nuc = {
  'A': 0,
  'G': 1,
  'T': 2,
  'C': 3,
}
'''

high_pr = {
  0: {
    0.35: 3,
    0.9: 1,
    1: 2
  },
  1: {
    0.18: 2,
    0.25: 3,
    1: 0
  }, 
  2: {
    0.12: 1,
    0.88: 3,
    1: 0
  },
  3: {
    0.14: 0,
    0.2: 1,
    1: 2
  }
}

low_pr = {
  0: {
    0.21: 3,
    0.79: 1,
    1: 2
  }, 
  1: {
    0.26: 2,
    0.3: 3,
    1: 0
  }, 
  2: {
    0.14: 1,
    0.86: 3,
    1: 0
  },
  3: {
    0.18: 0,
    0.21: 1,
    1: 2
  }
}

'''
erlich = {
  0: {
    : 3,
    : 1,
    1: 2
  }, 
  1: {
    : 2,
    : 3,
    1: 0
  }, 
  2: {
    : 1,
    : 3,
    1: 0
  },
  3: {
    : 0,
    : 1,
    1: 2
  }
}

goldman = {
  0: {
    : 3,
    : 1,
    1: 2
  }, 
  1: {
    : 2,
    : 3,
    1: 0
  }, 
  2: {
    : 1,
    : 3,
    1: 0
  },
  3: {
    : 0,
    : 1,
    1: 2
  }
}

high_pr4t = {
  0: {
    : 3,
    : 1,
    1: 2
  }, 
  1: {
    : 2,
    : 3,
    1: 0
  }, 
  2: {
    : 1,
    : 3,
    1: 0
  },
  3: {
    : 0,
    : 1,
    1: 2
  }
}
'''

average = {
  0 : {
    0.25: 3,
    0.75: 1,
    1: 2
  },
  1 : {
    0.07: 3,
    0.74: 0,
    1: 2
  },
  2 : {
    0.28: 1,
    0.85: 3,
    1: 0
  },
  3 : {
    0.175: 0,
    0.25: 1,
    1: 2
  }
}

method_lookup = {
  'high_pr': high_pr,
  'low_pr': low_pr,
  'erlich': {},
  'goldman': {},
  'high_pr4t': {},
  'average': average
}

# Function to get the substituting nucleotide based on probabilities.
# TODO: Check if correct from the paper. 
def getSubNucleotide(nuc, method):
  rand = random()

  method_probs = method_lookup[method]

  probabilities = method_probs[nuc]

  for p in probabilities:
    if rand < p:
      return probabilities[p]

