import { processInput, encodeOligos, decodeOligos } from "./utils";

onmessage = function(e) {
  console.log('Worker: Message received from main script, in modelling.worker file');

  const splitOligos = processInput(e.data);
  const encodedOligos = encodeOligos(splitOligos);
  const decodedOligos = decodeOligos(encodedOligos);
  console.log(decodedOligos[0]);
  
  postMessage("Finished in worker");
}

const synthesis = (oligos) => {
  postMessage("Starting Synthesis");

  const syn_oligos = [];

  // Needed to access loop index 
  for (const [i, oligo] of oligos.entries()) {

    const new_oligo = [];

  }

  postMessage("Finished Synthesis");
}