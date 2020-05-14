export const getSubNucleotidePCR = (nuc) => {
  const rand = Math.random();
  // A or T
  if (nuc % 2 == 0) {
    if (rand < 0.85) {
      return (nuc + 1) % 4;
    } else if (rand < 0.97) {
      return (nuc + 2) % 4;
    } else {
      return (nuc + 3) % 4;
    }
  }
  //G or C
  else {
    if (rand < 0.84) {
      return (nuc + 3) % 4;
    } else if (rand < 0.91) {
      return (nuc + 2) % 4;
    } else {
      return (nuc + 1) % 4;
    }
  }
}