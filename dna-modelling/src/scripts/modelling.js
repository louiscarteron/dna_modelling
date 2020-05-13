export default () => {


  //const utils = await import("./utils");

  onmessage = function(e) {
    console.log('Worker: Message received from main script');
    //const splitOligos = utils.processInput(e.data);
    //const encodedOligos = utils.encodeOligos(splitOligos);
    //utils.decodeOligos(encodedOligos);
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

}

/*
syn_oligos = []

  sub_rate, ins_rate, del_rate = method.values()

  total_sub = 0
  total_ins = 0
  total_del = 0
  
  # Going through each oligo to synthesis
  for i, oligo in enumerate(oligos):

    sub_counter = 0
    ins_counter = 0
    del_counter = 0

    new_oligo = []
    '''
    Need to go through each nucleotide and decide what to do with them. 
    TODO: encode conditional probability for substitution and addition.
    Defaulting to A (0) for now. 
    '''
    for nuc in oligo:

      new_nuc = nuc

      if random() < sub_rate:
        new_nuc = substitution.getSubNucleotide(nuc, 'average')

        sub_counter += 1

      elif random() < ins_rate:
        new_oligo.append(0) # Adding new nucleotide. TODO: need to choose which one to add from condition probability distrubution. 
        ins_counter += 1

      elif random() < del_rate:
        del_counter += 1
        # About to delete the nucleotide, so we can just ignore it and not add it to the new oligo
        continue
      
      new_oligo.append(new_nuc)

    if debug:

      print(f'Sub counter for olgio {i}: {sub_counter}')
      print(f'Ins counter for olgio {i}: {ins_counter}')
      print(f'Del counter for olgio {i}: {del_counter}')

    syn_oligos.append(new_oligo)

    # Used to track metrics 
    total_sub += sub_counter
    total_ins += ins_counter
    total_del += del_counter

  if debug:

    number_of_nuc = base_oligo_length * len(oligos)

    print(f'Total sub events: {total_sub}. Proportion is: {total_sub / number_of_nuc}. Target is: {sub_rate}')
    print(f'Total ins events: {total_ins}. Proportion is: {total_ins / number_of_nuc}. Target is: {ins_rate}')
    print(f'Total del events: {total_del}. Proportion is: {total_del / number_of_nuc}. Target is: {del_rate}')

  return syn_oligos
*/