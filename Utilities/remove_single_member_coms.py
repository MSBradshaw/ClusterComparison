import sys

for line in sys.stdin:
    row = line.strip().split('\t')
    if len(row) > 2:
        print(line.strip())
    else:
        quit()