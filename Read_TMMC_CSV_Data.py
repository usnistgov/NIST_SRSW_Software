#!/usr/bin/python

import csv

#-----------------------------------------------------------------
# Routine to read TMMC data from a CSV file
def Read_TMMC_CSV_Data(fileprefix,directory,tag):
    input_file = open(directory+fileprefix+tag+'.csv','r')
    reader = csv.reader(input_file)
    N = [ ]
    data = [ ]
    for row in reader:
        N.append(int(row[0]))
        data.append(float(row[1]))
    input_file.close()
    return N, data
#-----------------------------------------------------------------
