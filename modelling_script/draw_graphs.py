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


  total_1 = sum(insertions[0].values())
  total_2 = sum(insertions[1].values())
  total_3 = sum(insertions[2].values())

  mod_data = {
    'Flowcell1': [i * 100 / total_1 for i in insertions[0].values()],
    'Flowcell2': [i * 100 / total_2 for i in insertions[1].values()],
    'Flowcell3': [i * 100 / total_3 for i in insertions[2].values()]
  }

  tags = insertions[0].keys()
  fig, ax = plt.subplots()
  bar_plot(ax, mod_data, total_width=.8, single_width=.9)
  #w = 0.3
  #ax.bar(tags, deletions[0].values(), width=w, color='b', align='center')
  #ax.bar(tags, deletions[1].values(), width=w, color='g', align='center')
  #ax.bar(tags, deletions[2].values(), width=w, color='r', align='center')
  #ax.autoscale(tight=True)

  labels=['A', 'C', 'G', 'T']
  x = [0, 1, 2, 3]

  tick_spacing = 1

  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
  for label in ax.get_yticklabels()[1::2]:
    label.set_visible(False)

  plt.xticks(x, labels)
  plt.ylabel('Probability (%)')
  plt.title('Nucleotide Insertion for Flowcell 1, 2, 3')

  plt.savefig("graphs/insertions.png", bbox_inches="tight")

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