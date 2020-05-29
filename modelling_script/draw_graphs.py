import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from Levenshtein import distance, editops
from collections import Counter
import json
from process_fastq import print_table_from_report, process_25bp_report

def read_json(filepath):
  with open(filepath) as fp:
    data = json.load(fp)
  
  return data

def stacked():
  '''
  a_scores = [deletions[0]['A'] * 100 / total_1, deletions[1]['A'] * 100 / total_2, deletions[2]['A'] * 100 / total_3]
  c_scores = [deletions[0]['C'] * 100 / total_1, deletions[1]['C'] * 100 / total_2, deletions[2]['C'] * 100 / total_3]
  g_scores = [deletions[0]['G'] * 100 / total_1, deletions[1]['G'] * 100 / total_2, deletions[2]['G'] * 100 / total_3]
  t_scores = [deletions[0]['T'] * 100 / total_1, deletions[1]['T'] * 100 / total_2, deletions[2]['T'] * 100 / total_3]

  width = 0.35

  ind = [0, 1, 2]

  temp1 = np.add(a_scores, c_scores).tolist()
  temp2 = np.add(temp1, g_scores).tolist()

  
  p1 = plt.bar(ind, a_scores, width)
  p2 = plt.bar(ind, c_scores, width, bottom=a_scores)
  p3 = plt.bar(ind, g_scores, width, bottom=temp1)
  p4 = plt.bar(ind, t_scores, width, bottom=temp2)

  labels=['Flowcell 1', 'Flowcell 2', 'Flowcell 3']
  x = [0, 1, 2]

  plt.xticks(x, labels)
  plt.ylabel('Probability (%)')
  plt.title('Nucleotide Deletion for Flowcell 1, 2, 3')

  
  plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A', 'C', 'G', 'T'), bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

  plt.savefig("test2.png", bbox_inches="tight")

  return
  '''
  pass

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

  labels=['A', 'C', 'G', 'T']
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

  plt.savefig("graphs/insertions.png", bbox_inches="tight")

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

  labels=['A', 'C', 'G', 'T']
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

  plt.savefig("graphs/deletions.png", bbox_inches="tight")

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

  plt.savefig("graphs/substitutions.png", bbox_inches="tight")


def process_reports(reports):
  deletions = []
  insertions = []
  substitutions = []
  for report in reports:
    temp = process_25bp_report(report)
    dels, inss, subs, = print_table_from_report(temp, False)
    deletions.append(dels)
    insertions.append(inss)
    substitutions.append(subs)

  graph_deletions(deletions)
  graph_insertions(insertions)
  graph_substitutions(substitutions)

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

def main():
  reports = [read_json(f'data/flowcell/report/full_flowcell{i}.json') for i in range(1, 4)]
  process_reports(reports)


if __name__ == "__main__":
  main()