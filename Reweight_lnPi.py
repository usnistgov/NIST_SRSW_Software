#!/usr/bin/python

#Modules
import csv,math,sys,os

#-----------------------------------------------------------------
# Function to determine the number of state points in present data file
def line_count(filename):
    line_test = open(filename, mode='r')
    lines = 0
    for line in line_test:
        lines = lines + 1
    line_test.close()
    return lines
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# Routine to reweight a macrostate distribution (like a subroutine)
def reweight_lnPi(lnPi_old,N,beta,mu_old,mu_new):
    max_lnPi = lnPi_old[0]
    lnPi_star = [ ]
    for row in N:
        lnPi_temp = lnPi_old[row] + beta * float(N[row]) * (mu_new - mu_old)
        lnPi_star.append(lnPi_temp)
        if lnPi_temp > max_lnPi: max_lnPi = lnPi_temp
    sum_lnPi_star = 0.0
    for row in N:
        lnPi_star[row] = lnPi_star[row] - max_lnPi
        sum_lnPi_star = sum_lnPi_star + math.exp(lnPi_star[row])
    for row in N:
        lnPi_star[row] = lnPi_star[row] - math.log(sum_lnPi_star)
    return lnPi_star
#-----------------------------------------------------------------

#Main Program

#Unpack the data bundle
print 'Unpacking TMMC data in '+str(sys.argv[1])
os.system( 'mkdir ./tmp ; cd tmp ; tar -xzf ../'+str(sys.argv[1]) )

#Read in the file meta-data
numlines = line_count('./tmp/meta_data')
meta_data = open('./tmp/meta_data',mode='r')
line = 0
for line in range(numlines):
    tmp = meta_data.readline().rstrip() #Use the "rstrip" to get ride of \n
    data = tmp.partition(':')
    if data[0] == 'Units type':
        units_type = data[2].strip()
    elif data[0] == 'Temperature':
        temperature = float(data[2])
    elif data[0] == 'Chemical Potential':
        mu = float(data[2])
    elif data[0] == 'Volume':
        volume = float(data[2])
    elif data[0] == 'Epsilon/kB':
        eps_kB = float(data[2])
    elif data[0] == 'Sigma/Ang':
        sigma = float(data[2])
    elif data[0] == 'File Prefix':
        fileprefix = data[2].strip()
meta_data.close()

#Read in the Macrostate Distribution Data
macrostate_filename = './tmp/'+fileprefix+'lnpi.csv'
macrostate_file = open(macrostate_filename,'r')
reader = csv.reader(macrostate_file)
N = [ ]
lnPi = [ ]
sumPi = 0.0
for row in reader:
    N.append(int(row[0]))
    lnPi.append(float(row[1]))
    sumPi = sumPi + math.exp(float(row[1]))
Nmax = int(row[0])
macrostate_file.close()

#Error checking 
 #Confirm that the zero-density limit is sampled
if N[0] != 0:
    print 'error'
 #Confirm that lnPi[Nmax] is sufficiently small ( <= lnPi_max-10 )

#Variables; conversion if necessary
if units_type == 'reduced':
    beta = 1.0/temperature
else:
    print 'CODE INCOMPLETE: Invoked '+units_type+' units type'
    sys.exit()

#Pi(N) normalization
for row in N:
    lnPi[row] = lnPi[row] - math.log(sumPi)
sumPi = 1.0


#Histogram Reweighting
mu_old = mu
mu_new = mu + 0.20
lnPi_old = lnPi
lnPi_new = reweight_lnPi(lnPi_old,N,beta,mu_old,mu_new)
print N[0], lnPi_new[0]
print N[200], lnPi_new[200]

#File Cleanup
os.system( 'rm -rf ./tmp' )
