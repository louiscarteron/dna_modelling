#modelling_script 

#Model
To run model, have input sequences in ```input.txt``` (line seperated)
```python3 modelling.py```
Numbers can be tuned in ```data_config.py```

#Nanopore
To run nanopore analysis, get a copy of 3xr6 data set and flowcell 1, 2 and 3 from https://www.researchsquare.com/article/rs-27205/v1. It is recommended to run Porechop on the flowcell data.

You'll need to install Levensthein, fuzzywuzzy and Bio (can be found on pip)
Run analysis in 2 parts by invoking:
```python3 process_fastq.py```
First run generate the matches in a report.
Second run generates error profile.

#Graphs
If reports have been generated, graphs can be plotted using:
```python3 draw_graphs.py```
Requires numpy, scipy to be installed (can be found on pip)