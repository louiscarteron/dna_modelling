import { processInput, encodeOligos, decodeOligos } from "./utils";
import { getSubNucleotide, getInsNucleotide } from "./synthesisUtils";

onmessage = function(e) {
  console.log('Worker: Message received from main script, in modelling.worker file');

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);

  synthesis(encodedOligos);

  const decodedOligos = decodeOligos(encodedOligos);
  console.log(decodedOligos[0]);
  
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

    console.log(i);

    let insCounter = 0, subCounter = 0, delCounter = 0;

    const newOligo = [];

    for (const nuc of oligo) {
      
      let newNuc = nuc;

      if (Math.random() < insRate) {
        newNuc = getSubNucleotide(nuc, 'average');
        insCounter++;
      }

      else if (Math.random() < subRate) {
        newOligo.push(getInsNucleotide(nuc, 'average'));
        subCounter++;
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

  console.log(errorReport);

  postMessage("Finished Synthesis");
}