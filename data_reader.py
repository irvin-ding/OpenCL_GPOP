## Import library
import csv

def data_extract(file):
    with open(file) as f:
        reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
        data = list(reader)
        f.close()

    n = len(data)

    m = 0
    for i in range(0,n):
        if not data[i]:
            m += 1

    grid = list()

    p = n/m
    p = int(p)

    for i in range(0,m):
        L = list()
        for j in range(0,p-1):
            L.append(data[i*p+j+1])
        grid.append(L)

    return(grid)

## print(data_extract())
