import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from Levenshtein import distance, editops
from collections import Counter
import json
from process_fastq import print_table_from_report, process_25bp_report, _get_length_distribution, _get_score_distribution
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson

def read_json(filepath):
  with open(filepath) as fp:
    data = json.load(fp)
  
  return data

def stacked():

  plt.style.use('seaborn-colorblind')
  
  inserts = [4.4, 4.4, 4.5]
  deletes = [5.4, 5.2, 5.3]
  subs    = [9.6, 9.9, 10.0]

  width = 0.35

  ind = [0, 1, 2]

  temp1 = np.add(subs, deletes).tolist()

  
  p1 = plt.bar(ind, subs, width)
  p2 = plt.bar(ind, deletes, width, bottom=subs)
  p3 = plt.bar(ind, inserts, width, bottom=temp1)

  labels=['Flowcell 1', 'Flowcell 2', 'Flowcell 3']
  x = [0, 1, 2]

  plt.xticks(x, labels)
  plt.ylabel('Error Probability (%)')
  plt.title('Error Probability (%) per nucleotide for Flowcell 1, 2, 3')

  axes = plt.gca()
  axes.set_ylim([0,30])

  lines_x = [0.5, 1.5]
  for xc in lines_x:
    axes.axvline(x=xc, color="0.8")

  
  plt.legend((p1[0], p2[0], p3[0]), ('Substitutions', 'Deletes', 'Inserts'), bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

  plt.savefig("graphs/test3.png", bbox_inches="tight")

  return
  #pass

def graph_distribution(distribution, title, file_name):
  
  total_1 = sum(distribution[0].values())
  total_2 = sum(distribution[1].values())
  total_3 = sum(distribution[2].values())

  mod_data = {
    'Flowcell1': [i / total_1 for i in distribution[0].values()],
    'Flowcell2': [i / total_2 for i in distribution[1].values()],
    'Flowcell3': [i / total_3 for i in distribution[2].values()]
  }
  

  tags = distribution[0].keys()
  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)

  labels=['A', 'G', 'C', 'T']
  x = [0, 1, 2, 3]

  lines_x = [0.5, 1.5, 2.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  tick_spacing_y = 0.01

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Nucelotide Proportion')
  plt.title(title)

  plt.savefig(f'graphs/{file_name}.png', bbox_inches="tight")

def graph_insertions(insertions):
  
  total_1 = sum(insertions[0].values())
  total_2 = sum(insertions[1].values())
  total_3 = sum(insertions[2].values())

  mod_data = {
    'Flowcell1': [i / total_1 for i in insertions[0].values()],
    'Flowcell2': [i / total_2 for i in insertions[1].values()],
    'Flowcell3': [i / total_3 for i in insertions[2].values()]
  }
  

  tags = insertions[0].keys()
  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)

  labels=['A', 'G', 'C', 'T']
  x = [0, 1, 2, 3]

  lines_x = [0.5, 1.5, 2.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  tick_spacing_y = 0.01

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Conditional err prob')
  plt.title('Nucleotide Insertion for Flowcell 1, 2, 3')

  plt.savefig("graphs/insertions4.png", bbox_inches="tight")

def graph_deletions(deletions):
  
  total_1 = sum(deletions[0].values())
  total_2 = sum(deletions[1].values())
  total_3 = sum(deletions[2].values())

  mod_data = {
    'Flowcell1': [i / total_1 for i in deletions[0].values()],
    'Flowcell2': [i / total_2 for i in deletions[1].values()],
    'Flowcell3': [i / total_3 for i in deletions[2].values()]
  }
  

  tags = deletions[0].keys()
  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)

  labels=['A', 'G', 'C', 'T']
  x = [0, 1, 2, 3]

  lines_x = [0.5, 1.5, 2.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  tick_spacing_y = 0.01

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Conditional err prob')
  plt.title('Nucleotide Deletion for Flowcell 1, 2, 3')

  plt.savefig("graphs/deletions4.png", bbox_inches="tight")

def graph_subs2(substitutions):

  labels = ['A', 'G', 'C', 'T']

  total_1 = np.sum([list(a.values()) for a in substitutions[0].values()])
  total_2 = np.sum([list(a.values()) for a in substitutions[1].values()])
  total_3 = np.sum([list(a.values()) for a in substitutions[2].values()])

  '''
  for (key, val) in sub_info.items():
    
      total_for_sub = sum(val.values())
  '''

  mod_data1 = {
    'Flowcell1': [sum(val.values()) * 100/ total_1 for val in substitutions[0].values()],
    'Flowcell2': [sum(val.values()) * 100/ total_2 for val in substitutions[1].values()],
    'Flowcell3': [sum(val.values()) * 100/ total_3 for val in substitutions[2].values()]
  }

  tags = substitutions[0].keys()
  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()
  bar_plot(ax, mod_data1, total_width=.8, single_width=.9)

  x = [0, 1, 2, 3]

  lines_x = [0.5, 1.5, 2.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  #tick_spacing_y = 0.01

  #ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  #for label in ax.get_yticklabels()[1::2]:
  #  label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Proportion (%)')
  plt.title('Nucleotide substituted out for Flowcell 1, 2, 3')

  plt.savefig("graphs/subbedOut.png", bbox_inches="tight")

def graph_subs3(substitutions):

  labels = ['A', 'G', 'C', 'T']

  temp1 = [sum([list(val.values())[i] for val in substitutions[0].values()]) for i in range(4)]
  temp2 = [sum([list(val.values())[i] for val in substitutions[1].values()]) for i in range(4)]
  temp3 = [sum([list(val.values())[i] for val in substitutions[2].values()]) for i in range(4)]

  total_1 = sum(temp1)
  total_2 = sum(temp2)
  total_3 = sum(temp3)

  mod_data1 = {
    'Flowcell1': [i * 100 / total_1 for i in temp1],
    'Flowcell2': [i * 100 / total_2 for i in temp2],
    'Flowcell3': [i * 100 / total_3 for i in temp3]
  }

  tags = substitutions[0].keys()
  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()
  bar_plot(ax, mod_data1, total_width=.8, single_width=.9)

  x = [0, 1, 2, 3]

  lines_x = [0.5, 1.5, 2.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  #tick_spacing_y = 0.01

  #ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  #for label in ax.get_yticklabels()[1::2]:
  #  label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Proportion (%)')
  plt.title('Nucleotide substituted in for Flowcell 1, 2, 3')

  plt.savefig("graphs/subbedIn.png", bbox_inches="tight")

def graph_substitutions(substitutions):

  total_1 = np.sum([list(a.values()) for a in substitutions[0].values()])
  total_2 = np.sum([list(a.values()) for a in substitutions[1].values()])
  total_3 = np.sum([list(a.values()) for a in substitutions[2].values()])

  labels = ['A2G', 'A2C', 'A2T', 'G2A', 'G2C', 'G2T', 'C2A', 'C2G', 'C2T', 'T2A', 'T2G', 'T2C']
  
  mod_data = {
    'Flowcell1': [i / total_1 for a in substitutions[0].values() for i in a.values() if i != 0],
    'Flowcell2': [i / total_2 for a in substitutions[1].values() for i in a.values() if i != 0],
    'Flowcell3': [i / total_3 for a in substitutions[2].values() for i in a.values() if i != 0]
  }

  tags = substitutions[0].keys()
  plt.style.use('seaborn-colorblind')

  SMALL_SIZE = 20
  MEDIUM_SIZE = 24
  BIGGER_SIZE = 28

  plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
  plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
  plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
  plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
  plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
  plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
  plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

  fig, ax = plt.subplots(figsize=(20, 10))
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)

  x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

  lines_x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  tick_spacing_y = 0.01

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('conditionial err prob')
  plt.title('Nucleotide Substitution for Flowcell 1, 2, 3')

  plt.savefig("graphs/substitutions4.png", bbox_inches="tight")


def process_reports(reports):
  deletions = []
  insertions = []
  substitutions = []
  in_distribution = []
  read_distribution = []
  for report in reports:
    temp = process_25bp_report(report)
    dels, inss, subs, in_dist, read_dist = print_table_from_report(temp, False)
    deletions.append(dels)
    insertions.append(inss)
    substitutions.append(subs)
    in_distribution.append(in_dist)
    read_distribution.append(read_dist)

  #graph_deletions(deletions)
  #graph_insertions(insertions)
  #graph_substitutions(substitutions)
  graph_subs2(substitutions)
  graph_subs3(substitutions)

  #graph_distribution(in_distribution, "Nucleotide Distribution for Matched Inputs", "indist")
  #graph_distribution(read_distribution, "Nucleotide Distribution for Matched Reads", "redist")

  return

  '''
  total_1 = sum(insertions[0].values())
  total_2 = sum(insertions[1].values())
  total_3 = sum(insertions[2].values())

  mod_data = {
    'Flowcell1': [i * 100 / total_1 for i in insertions[0].values()],
    'Flowcell2': [i * 100 / total_2 for i in insertions[1].values()],
    'Flowcell3': [i * 100 / total_3 for i in insertions[2].values()]
  }
  '''

  total_1 = np.sum([list(a.values()) for a in substitutions[0].values()])
  total_2 = np.sum([list(a.values()) for a in substitutions[1].values()])
  total_3 = np.sum([list(a.values()) for a in substitutions[2].values()])

  labels = ['A2G', 'A2C', 'A2T', 'G2A', 'G2C', 'G2T', 'C2A', 'C2G', 'C2T', 'T2A', 'T2G', 'T2C']
  
  mod_data = {
    'Flowcell1': [i / total_1 for a in substitutions[0].values() for i in a.values() if i != 0],
    'Flowcell2': [i / total_2 for a in substitutions[1].values() for i in a.values() if i != 0],
    'Flowcell3': [i / total_3 for a in substitutions[2].values() for i in a.values() if i != 0]
  }

  colors = []

  tags = substitutions[0].keys()
  plt.style.use('seaborn-colorblind')
  fig, ax = plt.subplots(figsize=(20, 10))
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)
  #w = 0.3
  #ax.bar(tags, deletions[0].values(), width=w, color='b', align='center')
  #ax.bar(tags, deletions[1].values(), width=w, color='g', align='center')
  #ax.bar(tags, deletions[2].values(), width=w, color='r', align='center')
  #ax.autoscale(tight=True)

  #labels=['A', 'C', 'G', 'T']
  x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

  lines_x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
  for xc in lines_x:
    ax.axvline(x=xc, color="0.8")

  tick_spacing_y = 0.01
  tick_spacing_x = 4

  ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_x))

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('conditionial err prob')
  plt.title('Nucleotide Substitution for Flowcell 1, 2, 3')

  plt.savefig("graphs/substitutions.png", bbox_inches="tight")

def bar_plot(ax, data, colors=None, total_width=0.8, single_width=1, legend=True):
  """Draws a bar plot with multiple bars per data point.

  Parameters
  ----------
  ax : matplotlib.pyplot.axis
      The axis we want to draw our plot on.

  data: dictionary
      A dictionary containing the data we want to plot. Keys are the names of the
      data, the items is a list of the values.

      Example:
      data = {
          "x":[1,2,3],
          "y":[1,2,3],
          "z":[1,2,3],
      }

  colors : array-like, optional
      A list of colors which are used for the bars. If None, the colors
      will be the standard matplotlib color cyle. (default: None)

  total_width : float, optional, default: 0.8
      The width of a bar group. 0.8 means that 80% of the x-axis is covered
      by bars and 20% will be spaces between the bars.

  single_width: float, optional, default: 1
      The relative width of a single bar within a group. 1 means the bars
      will touch eachother within a group, values less than 1 will make
      these bars thinner.

  legend: bool, optional, default: True
      If this is set to true, a legend will be added to the axis.
  """

  # Check if colors where provided, otherwhise use the default color cycle
  if colors is None:
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

  # Number of bars per group
  n_bars = len(data)

  # The width of a single bar
  bar_width = total_width / n_bars

  # List containing handles for the drawn bars, used for the legend
  bars = []

  # Iterate over all data
  for i, (name, values) in enumerate(data.items()):
    # The offset in x direction of that bar
    x_offset = (i - n_bars / 2) * bar_width + bar_width / 2

    # Draw a bar for every value of that type
    for x, y in enumerate(values):
      bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=colors[i % len(colors)])

    # Add a handle to the last drawn bar, which we'll need for the legend
    bars.append(bar[0])

  # Draw legend if we need
  if legend:
    ax.legend(bars, data.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


def graph_lengths():
  
  counter = _get_length_distribution("data/flowcell/data/pc_nm_flowcell_1.fasta", "fastq")
  dict_counter = dict(counter)

  plt.style.use('seaborn-colorblind')

  SMALL_SIZE = 20
  MEDIUM_SIZE = 24
  BIGGER_SIZE = 28

  plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
  plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
  plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
  plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
  plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
  plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
  plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


  fig, ax = plt.subplots(figsize=(20, 10))

  sorted_dict = dict(sorted(dict_counter.items()))

  x_sum = np.sum(list(sorted_dict.values()))

  x = np.cumsum(list(sorted_dict.values()) / (x_sum * 0.01))[:50]

  x_labels = list(sorted_dict.keys())[:50]

  x_pos = [i for i in range(50)]

 

  ax.bar(x_pos, x)
  #plt.axis([0, 20, 0, 10000])

  tick_spacing_x = 20
  tick_spacing_y = 1

  ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_x))

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
  [l.set_visible(False) for (i,l) in enumerate(ax.yaxis.get_ticklabels()) if (i - 1) % 10 != 0]
  

  plt.xlabel("Length")
  plt.ylabel("Proportion < length (%)")
  plt.title("Cumulative lengths of Flowcell 1")

  #plt.xticks(x_pos, x)
  plt.xticks(x_pos, x_labels)
  for label in ax.get_xticklabels()[1::2]:
    label.set_visible(False)

  plt.savefig("graphs/lengths1.png", bbox_inches="tight")

def graph_scores():
  pass

def fit_function(k, lamb):
    '''poisson function, parameter lamb is the fit parameter'''
    return poisson.pmf(k, lamb)

def temp():
  s = [66, 57, 66, 64, 45, 57, 87, 66, 59, 59, 65, 82, 69, 62, 64, 55, 76, 58, 63, 62, 68, 63, 62, 68, 69, 79, 58, 63, 47, 67, 43, 53, 58, 66, 53, 52, 71, 60, 75, 53, 50, 67, 66, 65, 71, 62, 65, 60, 59, 59, 57, 68, 60, 58, 61, 56, 65, 58, 50, 53, 66, 71, 56, 68, 69, 71, 70, 62, 62, 61, 56, 59, 55, 57, 54, 54, 64, 64, 41, 61, 59, 61, 73, 67, 54, 49, 60, 58, 79, 56, 67, 67, 64, 61, 62, 60, 70, 70, 49, 68]
  count, bins, ignored = plt.hist(s)

  plt.savefig("graphs/poisson.png")

def temp2():
  #s1 = np.random.poisson(61, 100)

  t = np.arange(0, 100, 0.1)
  d = np.exp(-61.44)*np.power(61.44, t)/factorial(t)

  plt.style.use('seaborn-pastel')


  s = [66, 57, 66, 64, 45, 57, 87, 66, 59, 59, 65, 82, 69, 62, 64, 55, 76, 58, 63, 62, 68, 63, 62, 68, 69, 79, 58, 63, 47, 67, 43, 53, 58, 66, 53, 52, 71, 60, 75, 53, 50, 67, 66, 65, 71, 62, 65, 60, 59, 59, 57, 68, 60, 58, 61, 56, 65, 58, 50, 53, 66, 71, 56, 68, 69, 71, 70, 62, 62, 61, 56, 59, 55, 57, 54, 54, 64, 64, 41, 61, 59, 61, 73, 67, 54, 49, 60, 58, 79, 56, 67, 67, 64, 61, 62, 60, 70, 70, 49, 68]

  #s = [72, 61, 57, 54, 73, 68, 66, 63, 55, 49, 64, 50, 53, 55, 61, 66, 61, 66, 56, 59, 62, 53, 60, 62, 54, 61, 64, 61, 48, 52, 49, 63, 66, 73, 56, 64, 67, 63, 80, 57, 66, 54, 61, 62, 44, 59, 66, 70, 45, 64, 54, 57, 62, 66, 70, 76, 55, 62, 41, 73, 57, 60, 68, 60, 62, 61, 52, 61, 56, 67, 70, 58, 65, 65, 56, 58, 66, 60, 58, 54, 64, 67, 54, 59, 55, 49, 60, 57, 62, 51, 62, 62, 56, 67, 60, 60, 66, 63, 62, 58]

  bins = np.arange(100) - 0.5

  count, bin_edge, ignored = plt.hist(s, bins=bins, density=True, label='Substitution Events')

  bin_middles = 0.5 * (bin_edge[1:] + bin_edge[:-1])

  parameters, cov_matrix = curve_fit(fit_function, bin_middles, count)

  x_plot = np.arange(0, 100)

  plt.plot(
      x_plot,
      fit_function(x_plot, *parameters),
      linestyle='--',
      label='Fit result',
      linewidth=2.0,
      color='red'
  )

  plt.plot(t, d, '--', label='Poisson(61)', linewidth=2.0, color='black')

  plt.legend()

  plt.title('Distribution of Substitution for 100 runs over 100 sequences (128nt)')

  plt.savefig("graphs/poisson3.png", bbox_inches="tight")

def main():
  reports = [read_json(f'data/flowcell/report/full_flowcell{i}.json') for i in range(1, 4)]
  process_reports(reports)


if __name__ == "__main__":
  main()