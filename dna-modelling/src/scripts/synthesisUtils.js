export const getSubNucleotide = (nuc, method) => {
  const rand = Math.random() * 100;

  const methodProbs = methodLookup[method];
  const probabilities = methodProbs[nuc];

  for (const p in probabilities) {
    if (rand < p) {
      return probabilities[p];
    }
  }
}

export const getInsNucleotide = (nuc, method) => {
  console.log("Getting ins");
  return nuc;
}

const high_pr = {
  0: {
    35: 3,
    90: 1,
    100: 2
  },
  1: {
    18: 2,
    25: 3,
    100: 0
  }, 
  2: {
    12: 1,
    88: 3,
    100: 0
  },
  3: {
    14: 0,
    20: 1,
    100: 2
  }
}

const low_pr = {
  0: {
    21: 3,
    79: 1,
    100: 2
  }, 
  1: {
    26: 2,
    30: 3,
    100: 0
  }, 
  2: {
    14: 1,
    86: 3,
    100: 0
  },
  3: {
    18: 0,
    21: 1,
    100: 2
  }
}

const erlich = {
  0: {
    23: 3,
    59: 1,
    100: 2
  }, 
  1: {
    71: 2,
    84: 3,
    100: 0
  }, 
  2: {
    46: 1,
    75: 3,
    100: 0
  },
  3: {
    27: 0,
    50: 1,
    100: 2
  }
}

const goldman = {
  0: {
    11: 3,
    74: 1,
    100: 2
  }, 
  1: {
    17: 2,
    45: 3,
    100: 0
  }, 
  2: {
    45: 1,
    78: 3,
    100: 0
  },
  3: {
    13: 0,
    40: 1,
    100: 2
  }
}

const high_pr4t = {
  0: {
    30: 3,
    80: 1,
    100: 2
  }, 
  1: {
    8: 2,
    13: 3,
    100: 0
  }, 
  2: {
    10: 1,
    80: 3,
    100: 0
  },
  3: {
    7: 0,
    12: 1,
    100: 2
  }
}

const average = {
  0 : {
    25: 3,
    75: 1,
    100: 2
  },
  1 : {
    7: 3,
    74: 0,
    100: 2
  },
  2 : {
    28: 1,
    85: 3,
    100: 0
  },
  3 : {
    17: 0,
    25: 1,
    100: 2
  }
}

const methodLookup = {
  'high_pr': high_pr,
  'low_pr': low_pr,
  'erlich': erlich,
  'goldman': goldman,
  'high_pr4t': high_pr4t,
  'average': average
}