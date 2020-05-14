export const decaySimulation = (oligos, decayEvents, baseOligoLength) => {
  for(let i = 0; i < decayEvents; i++) {

    const decayedOligoIndex = Math.floor(Math.random() * oligos.length);

    const decayedOligo = oligos[decayedOligoIndex];

    if (decayedOligo.length > 10) {
      const breakPoint = Math.floor(Math.random() * baseOligoLength) + 1;

      const temp = decayedOligo.splice(breakPoint);

      oligos[decayedOligoIndex] = [decayedOligo, temp];

    } else {
      console.log("Should have broken an already broken strand");
    }
      
  }
}