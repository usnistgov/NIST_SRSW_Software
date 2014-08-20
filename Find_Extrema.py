#!/usr/bin/python

#Modules
import csv,math,sys,os

# Function to determine the number of state points in present data file
def line_count(filename):
    line_test = open(filename, mode='r')
    lines = 0
    for line in line_test:
        lines = lines + 1

    line_test.close()
    return lines



#Main Program

#Constants
N_range = 50
max_grad2 = -4.0
max_fraction = 0.550

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
N_max = int(row[0])
N_min = N[0]
macrostate_file.close()

#Error checking 
 #Confirm that the zero-density limit is sampled
if N[0] != 0:
    print 'error'
 #Confirm that lnPi[N_max] is sufficiently small ( <= lnPi_max-10 )

#Pi(N) normalization
for row in N:
    lnPi[row] = lnPi[row] - math.log(sumPi)
sumPi = 1.0

#Log(|dlnPidN|)
dydx2 = [ ]
for row in N:
    if row == N_min:
        dydx2.append(5.0) #just make it positive and big
    elif row == N_max:
        dydx2.append(5.0) #just make it positive and big
    else:
        term = (lnPi[row+1]-lnPi[row-1])/(N[row+1]-N[row-1])
        dydx2.append(math.log(math.sqrt(term * term)))

#Identify candidate extrema
N_start = N_min+2
N_end = N_max-2
N_index = N_start
list_extrema_locations = [ ]
isearch = 0
while N_index < N_end:
    N_lo = N_index
    N_hi = N_index + N_range-1
    if N_hi > N_end: N_hi = N_end
    #Candidates in the current range
    list_extrema_locations.append(0)
    list_extrema_locations[isearch] = N_lo
    N_local = N_lo
    while N_local <= N_hi:
        if dydx2[N_local] < dydx2[list_extrema_locations[isearch]]:
            list_extrema_locations[isearch] = N_local
        N_local = N_local + 1
    #Check adjacent ranges
    if isearch > 0:
        if abs(list_extrema_locations[isearch]-list_extrema_locations[isearch-1]) < N_range:
            if dydx2[list_extrema_locations[isearch]] < dydx2[list_extrema_locations[isearch-1]]:
                list_extrema_locations[isearch-1] = list_extrema_locations[isearch]
            isearch = isearch - 1

    isearch = isearch + 1
    N_index = N_index + N_range
isearch = isearch - 1

print 'Candidate Extrema'
for index in range(isearch+1):
    print list_extrema_locations[index], dydx2[list_extrema_locations[index]]
print ''

#Finalize the list of candidate extrema
final_extrema_locations = [ ]
n_extrema = 0
print 'Finalist Extrema'
for index in range(isearch):
#    print index, list_extrema_locations[index], dydx2[list_extrema_locations[index]]
    if dydx2[list_extrema_locations[index]] < max_grad2:
        final_extrema_locations.append(list_extrema_locations[index])
        n_extrema = n_extrema + 1
        print final_extrema_locations[n_extrema-1], dydx2[final_extrema_locations[n_extrema-1]]
print ''

#Determine which extrema are maxima:
n_maxima = 0
maximum_list = [ ]
n_search = N_range + int(float(N_range)/2.0)
for index in range(n_extrema):
#    print 'Examining extremum at '+str(final_extrema_locations[index])
    max_reference = lnPi[final_extrema_locations[index]]
    max_counter = 0
    delta_n = -n_search
    while delta_n <= n_search:
        if delta_n != 0:
            N_temp = final_extrema_locations[index] + delta_n
            if N_min < N_temp < N_max:
                if max_reference > lnPi[N_temp]: max_counter = max_counter +1
        delta_n = delta_n + 1
#    print max_counter, final_extrema_locations[index],lnPi[final_extrema_locations[index]], index
    
    #Deal with noise
    max_critical = max_fraction * float(2*n_search+1)
    if max_counter >= max_critical:
        max_location_temp = final_extrema_locations[index]
        max_reference = lnPi[max_location_temp]
        delta_n = -n_search
        while delta_n <= n_search:
            if delta_n != 0:
                N_temp = final_extrema_locations[index]+delta_n
                if N_min < N_temp < N_max:
                    if max_reference < lnPi[N_temp]:
                        max_location_temp = N_temp
                        max_reference = lnPi[N_temp]

            delta_n = delta_n + 1
        n_maxima = n_maxima + 1
        maximum_list.append(max_location_temp)

#        print 'LOCAL MAXIMUM: '+str(maximum_list[n_maxima-1])
    if n_maxima > 1:
        if (maximum_list[n_maxima-1] - maximum_list[n_maxima-2]) < (N_range+1):
            if lnPi[maximum_list[n_maxima-2]] < lnPi[maximum_list[n_maxima-1]]:
                maximum_list[n_maxima-2] = maximum_list[n_maxima-1]
            n_maxima = n_maxima - 1

#Verify that there is a true minima between two maxima
if n_maxima == 0:
    print 'ERROR: No stable phase for mu = '+str(mu)
maximum_true = [ ]
maximum_true.append('T') #Start with the assumption that all are local maxima
for index in range(n_maxima-1):
    maximum_true.append('T') #have to add an entry
    N_lo = maximum_list[index]
    N_hi = maximum_list[index+1]
    ref_point = lnPi[N_lo]
    ref_location = N_lo
    N_local = N_lo
    while N_local <= N_hi:
        if lnPi[N_local] < ref_point:
            ref_point = lnPi[N_local]
            ref_location = N_local
        N_local = N_local + 1

    W_lo = lnPi[N_lo] - ref_point
    W_hi = lnPi[N_hi] - ref_point

    if W_lo <= 1.0: 
        if maximum_true[index] == 'T': maximum_true[index] == 'F'
    if W_hi <= 1.0: 
        if maximum_true[index+1] == 'T': maximum_true[index+1] == 'F'


print 'Identified the Stable Local Maxima:'
for index in range(n_maxima):
    if maximum_true[index] == 'T':
        print 'LOCAL MAXIMUM: '+str(maximum_list[index])
    


#File Cleanup
os.system( 'rm -rf ./tmp' )
