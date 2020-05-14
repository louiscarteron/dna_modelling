import { processInput, encodeOligos, decodeOligos } from "./utils";
import { getSubNucleotide, getInsNucleotide } from "./synthesisUtils";
import { decaySimulation } from "./storageUtils";

onmessage = function(e) {
  console.log('Worker: Message received from main script, in modelling.worker file');

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);

  const baseOligoLength = encodedOligos[0].length;

  const synOligos = synthesis(encodedOligos);

  const storedOligos = storage(synOligos, 521, baseOligoLength);

  const decodedOligos = decodeOligos(synOligos);
  
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

//TODO: Add capabilities to account for temperature settings to decay faster
const storage = (oligos, time, baseOligoLength) => {
  //Half life of DNA in an ideal environment
  const halfLife = 521;

  const dt = 1;

  const decayRate = Math.log(2) / halfLife;

  const finalOligos = oligos.slice(); // Create copy to modify.

  let t = 0;
  let decay = 0;

  while (t < time) {
    const rand = Math.random();

    if (rand < decayRate) {
      decaySimulation(finalOligos, 5, baseOligoLength);
      decay++;

    }

    t += dt;
  }

  console.log(`Number of decay events was ${decay}`);

  return finalOligos;
}