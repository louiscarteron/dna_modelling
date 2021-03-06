import { processInput, encodeOligos, decodeOligos, flattenPCR } from "./utils";
import { getSubNucleotide, getInsNucleotide } from "./synthesisUtils";
import { decaySimulation } from "./storageUtils";
import { getSubNucleotidePCR } from "./pcrUtils";
import { StatusEnum } from "./modellingStatus";
import sampleSize from "lodash/sampleSize";
import levenshtein from "js-levenshtein";
import Error from "./errors.js";

Array.prototype.extend = function (other) {
  other.forEach(function(v) {
    this.push(v)
  }, this);
}

onmessage = function(e) {

  const errors = new Error(["synthesis", "storage", "pcr", "sequencing"]);

  postMessage({status: StatusEnum.START});

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);

  const baseOligoLength = encodedOligos[0].length;

  postMessage({status: StatusEnum.SYNTHESIS});
  const synOligos = synthesis(encodedOligos, errors);

  postMessage({status: StatusEnum.STORAGE});
  const storedOligos = storage(synOligos, 1400, baseOligoLength);

  postMessage({status: StatusEnum.PCR});
  //const pcrOligos = pcr(storedOligos, 60);
  const pcrOligos = exponentialPCR(storedOligos, 20);

  const flatPcr = flattenPCR(pcrOligos);

  postMessage({status: StatusEnum.SEQUENCING});
  const seqOligos = sequencing(flatPcr);

  const decodedOligos = decodeOligos(seqOligos);
  
  postMessage({status: StatusEnum.FINISH});

  calculateLevenshteinScore(splitOligos, decodedOligos);

  errors.printReport();
}

//TODO: consider adding report object. probs gonna do it now.
const synthesis = (oligos, errors) => {

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

  errors.incrementError("synthesis", "inserts", totalIns);
  errors.incrementError("synthesis", "deletes", totalDel);
  errors.incrementError("synthesis", "substitutions", totalSub);

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
const pcrLowMemory = (oligos, cycles = 10) => {

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

//TODO: Do some smart garbage clean up to not make the size of the array blow up. 
// Probably solve the memory issue with doing some random sampling at the end to get a sub sample. 

const exponentialPCR = (oligos, cycles = 10) => {

  console.info("Using exponential PCR which is very memory intensive.");

  if (cycles > 10) {
    console.warn(`About to apply ${cycles} PCR cycles, resulting in high CPU and Memory usage.`);
    //return;
  }

  //const pcrOligos = oligos.slice();
  const pcrOligos = [];

  let subCounter = 0, insCounter = 0, delCounter = 0;

  // Generating PCR in order of input oligos. 

  const totalOligos = oligos.length;

  for (const [index, oligo] of oligos.entries()) {

    // Needed to update main thread on current progress 
    postMessage({status: StatusEnum.PCR, progress: {
      current: index,
      total: totalOligos
    }})
  
    let pcrForOligo = [oligo];

    for (let i = 0; i < cycles; i++) {

      const oligoBuffer = [];

      for (const oligo of pcrForOligo) {

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

        oligoBuffer.push(newOligo);

      }

      pcrForOligo.extend(oligoBuffer);
      
      //TODO: Figure out a good middle ground for max size. 
      const maxSize = Math.pow(2, 11);

      if (pcrForOligo.length > maxSize) {
        pcrForOligo = sampleSize(pcrForOligo, maxSize);
      }

    }

    pcrOligos[index] = pcrForOligo;

  }

  console.log(pcrOligos.length);
  console.log(pcrOligos[0].length);
  console.log(pcrOligos[0].slice(0, 10));

  console.log(`Subs: ${subCounter}`);
  console.log(`Ins: ${insCounter}`);
  console.log(`Del: ${delCounter}`);
  
  return pcrOligos; 
}

const sequencing = (oligos) => {
  console.log("Should sequence oligos here");

  return oligos;
}

const calculateLevenshteinScore = (inputOligos, outputOligos) => {

  console.log(inputOligos.length, outputOligos.length);

  const pcrDups = outputOligos.length / inputOligos.length;

  const distances = [];

  for (let i = 0; i < outputOligos.length; i++) {
    const score = levenshtein(inputOligos[~~(i/pcrDups)], outputOligos[i]);

    const temp = {
      oligo: i,
      score: score
    }
    distances.push(temp);
  }

  console.table(distances.slice(0, 10));
  
  return;
}