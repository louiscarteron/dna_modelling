import { processInput, encodeOligos, decodeOligos } from "./utils";
import { getSubNucleotide, getInsNucleotide } from "./synthesisUtils";
import { decaySimulation } from "./storageUtils";
import { getSubNucleotidePCR } from "./pcrUtils";

onmessage = function(e) {
  console.log('Worker: Message received from main script, in modelling.worker file');

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);

  console.log(`Encoded length: ${encodedOligos.length}`);

  const baseOligoLength = encodedOligos[0].length;

  const synOligos = synthesis(encodedOligos);

  console.log(`Syn length: ${synOligos.length}`);

  const storedOligos = storage(synOligos, 521, baseOligoLength);

  console.log(`Storage length: ${storedOligos.length}`);

  const pcrOligos = pcr(storedOligos, 60);

  console.log(`PCR length: ${pcrOligos.length}`);

  const seqOligos = sequencing(pcrOligos);

  console.log(`Sequencing length: ${seqOligos.length}`);

  const decodedOligos = decodeOligos(seqOligos);

  console.log(`Decoded length: ${decodedOligos.length}`);
  
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
    //console.log(i);

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

//TODO: If strand is split, can not apply pcr to it. Need to apply logic to do so and report in statistics. 
//TODO: fix error with creating oligos * cycles number of final oligos. Is it really an error?
const pcr = (oligos, cycles = 10) => {

  const pcrOligos = [];

  let subCounter = 0, insCounter = 0, delCounter = 0;

  for (let i = 0; i < cycles; i++) {

    for (const oligo of oligos) {

      const newOligo = [];

      for (const nuc of oligo) {

        const rand = Math.random();

        let newNuc = nuc;

        if (rand < 1.8e-4) {

          const errorType = Math.random();
          if (errorType < 0.973) {
            newNuc = getSubNucleotidePCR(nuc);
            subCounter++;
          } else if (errorType < 0.99) {
            delCounter++;
            continue;
          } else {
            //TODO: need to do insertion;
            insCounter++;
          }

        }

        newOligo.push(newNuc);

      }

      pcrOligos.push(newOligo);

    }

  }

  console.log(`Subs: ${subCounter}`);
  console.log(`Ins: ${insCounter}`);
  console.log(`Del: ${delCounter}`);

  return pcrOligos;

}

const sequencing = (oligos) => {
  console.log("Should sequence oligos here");

  return oligos;
}