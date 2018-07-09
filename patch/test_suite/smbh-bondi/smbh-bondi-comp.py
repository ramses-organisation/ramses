with open('output_00002/sink_00002.csv') as f:
    for line in f:
        if line.split(' ')[1] != '#':
            print line.split(',')[12]
