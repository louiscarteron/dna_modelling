
const nuc2str = {
  0: 'A',
  1: 'G',
  2: 'T', 
  3: 'C'
}

const str2nuc = {
  'A': 0,
  'G': 1,
  'T': 2,
  'C': 3
}

export const processInput = (oligos) => {

  return oligos.split('\n');

}

export const encodeOligos = (oligos) => {

  return oligos.map(oligo => [...oligo].map(o => str2nuc[o]));

}

export const decodeOligos = (oligos) => {
  
  return oligos.map(oligo => oligo.map(o => nuc2str[o]).join(''));
  
}

export const flattenPCR = (pcrOligos) => {

  return pcrOligos.flat()
  
}
