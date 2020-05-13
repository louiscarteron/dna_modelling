import csv

def process_csv(path):

  data = []

  with open(path) as f:
    csv_reader = csv.reader(f, delimiter=',')
    for r in csv_reader:
      data.append(r[1])

  return data

def remove_primers(oligos):
  

def main():
  oligos = process_csv("data/3xr6.csv")
  print(len(oligos))

if __name__ == "__main__":
  main()