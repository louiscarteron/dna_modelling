import { processInput, encodeOligos, decodeOligos } from "./utils";
import { getSubNucleotide, getInsNucleotide } from "./synthesisUtils";

onmessage = function(e) {
  console.log('Worker: Message received from main script, in modelling.worker file');

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);

  const synOligos = synthesis(encodedOligos);

  const decodedOligos = decodeOligos(synOligos);
  console.log(decodedOligos);
  
  postMessage("Finished in worker");
}

//TODO: consider adding report object. probs gonna do it now.
const synthesis = (oligos) => {
  postMessage("Starting Synthesis");

  const errorReport = {};

  const synOligos = [];

  // TODO: Change the following to use input objects for the parameters. 

  const insRate = 0.001;
  const subRate = 0.001;
  const delRate = 0.001;

  let totalIns = 0, totalSub = 0, totalDel = 0;

  // Needed to access loop index 
  for (const [i, oligo] of oligos.entries()) {

    //TODO: Post this to the main thread to get a progress bar working. 
    console.log(i);

    let insCounter = 0, subCounter = 0, delCounter = 0;

    const newOligo = [];

    for (const nuc of oligo) {
      
      let newNuc = nuc;

      if (Math.random() < subRate) {
        newNuc = getSubNucleotide(nuc, 'average');
        subCounter++;
      }

      else if (Math.random() < insRate) {
        newOligo.push(getInsNucleotide(nuc, 'average'));
        insCounter++;
      }
      
      else if (Math.random() < delRate) {
        delCounter++;
        continue;
      }

      newOligo.push(newNuc);
    }

    synOligos.push(newOligo);

    errorReport[`oligo${i}`] = {
      insertions: insCounter,
      substitutions: subCounter,
      deletions: delCounter
    }

    totalIns += insCounter;
    totalSub += subCounter;
    totalDel += delCounter;

  }

  errorReport['total'] = {
    insertions: totalIns,
    substitutions: totalSub,
    deletions: totalDel
  };

  postMessage("Finished Synthesis");

  return synOligos;
}