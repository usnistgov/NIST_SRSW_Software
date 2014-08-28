#!/usr/bin/python

#------------SOFTWARE DISCLAIMER AND REDISTRIBUTION CONDITIONS----------------
#   This software was developed at the National Institute of Standards and
#   Technology by employees of the Federal Government in the course of their
#   official duties. Pursuant to Title 17 Section 105 of the United States 
#   Code this software is not subject to copyright protection and is in the 
#   public domain. muVT_VLE.py is an experimental system. NIST assumes no
#   responsibility whatsoever for its use by other parties, and makes no
#   guarantees, expressed or implied, about its quality, reliability, or any
#   other characteristic. We would appreciate acknowledgement if the software
#   is used.
#
#   This software can be redistributed and/or modified freely provided that 
#   any derivative works bear some notice that they are derived from it, and 
#   any modified versions bear some notice that they have been modified.
#
#-----------------------------------------------------------------------------
#   Software Name: muVT_VLE.py
#   Brief Description: Script to identify the vapor-liquid equilibrium
#    conditions from a particle-number probability distribution
#   Authors: Daniel W. Siderius, PhD ; Vincent K. Shen, PhD
#   Contact Information: daniel.siderius@nist.gov
#
#   Computation of the equation of state is based on a method discussed by
#    Errington in J. Chem. Phys., 118(22):9915-9925, (2003).
#    The method is presented on page 9917 and in equation 13.
#
#   Version History:
#
#   Version: 1.0
#   Release Date: 27 August 2012
#   Notes: Original version
#-----------------------------------------------------------------------------

#Modules
import sys,os,shutil,numpy, csv
from scipy.signal import argrelextrema
from Reweight_lnPi import *
from Bulk_Calculate_Properties import *
from Read_TMMC_CSV_Data import *
from Renormalize_lnPi import *
from Parse_Standard_XML_Data import *

# Directories and Files
work_dir = './tmp/'
input_filename = str(sys.argv[1])
meta_data = work_dir+'metadata.xml'

# Unpack the data bundle & read the XML metadata
print 'Unpacking TMMC data in '+str(sys.argv[1])
if not os.path.exists(work_dir): os.makedirs(work_dir)
os.system( 'cd '+work_dir+' ; tar -xzf ../'+input_filename )
(temperature,lnZ,volume,fileprefix,units_type) = Parse_Standard_XML_data(meta_data)

# Useful Conversions
mu = temperature * lnZ #log(activity) -> chemical potential
beta = 1.0e0 / temperature

# Read in the Macrostate Probability Distribution Data, get N bounds, and normalize lnPi
(N,lnPi) = Read_TMMC_CSV_Data(fileprefix,work_dir,'lnpi')
N_max = N[-1]; N_min = N[0]
(lnPi,sumPi) = Renormalize(lnPi,shift=True)

# Read in the Canonical Ensemble Energy Data
(Ntemp,energy) = Read_TMMC_CSV_Data(fileprefix,work_dir,'energy')

# Set default values for search criteria
tolerance = 1.0e-14
lnPi_threshold = 10.0e0
mu_new = mu
midpoint = int( (N_max + N_min) / 2 )
# Read in optional "helper" data for locating the coexistence point, if helper file is specified
if len(sys.argv) > 2:
    print 'Reading solver-aid data from '+sys.argv[2]
    try:
        optional_dom = minidom.parse( sys.argv[2] )
        if optional_dom.getElementsByTagName('tolerance'):
            tolerance = float(optional_dom.getElementsByTagName('tolerance')[0].firstChild.data)
            print 'Using convergence tolerance = '+str(tolerance)
        if optional_dom.getElementsByTagName('lnPi_threshold'):
            lnPi_threshold = float(optional_dom.getElementsByTagName('lnPi_threshold')[0].firstChild.data)
            print 'Using delta(lnPi) threshold = '+str(lnPi_threshold)
        if optional_dom.getElementsByTagName('mu_guess'):
            mu_new = float(optional_dom.getElementsByTagName('mu_guess')[0].firstChild.data)
            print 'Using first guess of coexistence mu* = '+str(mu_new)
        if optional_dom.getElementsByTagName('midpoint'):
            midpoint = float(optional_dom.getElementsByTagName('midpoint')[0].firstChild.data)
            print 'Using '+str(midpoint)+' as first guess of midpoint'
    except:
        print sys.argv[2]+' is not a well-formed XML file. Please specify a proper XML helper file'
        sys.exit()

#Store "old" macrostate data
mu_old = mu
lnPi_old = lnPi

# Crudely Bracket the coexistence point
print 'Estimating Initial Bracket'
#  One bound based on either the default mu or something input in optional data
lnPi_new = Reweight_lnPi(lnPi_old,N,beta,mu,mu_new)
minima = argrelextrema(numpy.array(lnPi_new),numpy.less)
if len(minima[0]) == 1:
    min_N = N[minima[0]]
elif len(minima[0]) > 1:
    min_N = int(sum(minima[0])/len(minima[0]))
else:
    min_N = midpoint
(nmols, pressure, U, GPFE) = Bulk_Calculate_Properties(lnPi_new,energy,N,mu_new,volume,beta,min_N)
delta_p = [ pressure[0]-pressure[1], 0.0e0 ]
#  Identify another bound by a crude search
mu_ref = mu_new
delta_mu = mu_ref * 1.0e-4
solved = False
iterations = 0
while not solved:
    iterations = iterations + 1
    mu_new = mu_ref + delta_mu * float(iterations)
    lnPi_new = Reweight_lnPi(lnPi_old,N,beta,mu,mu_new)
    minima = argrelextrema(numpy.array(lnPi_new),numpy.less)
    if len(minima[0]) == 1:
        min_N = N[minima[0]]
    elif len(minima[0]) > 1:
        min_N = int(sum(minima[0])/len(minima[0]))
    else:
        min_N = int( (N_max + N_min) / 2 )
    (nmols, pressure, U, GPFE) = Bulk_Calculate_Properties(lnPi_new,energy,N,mu_new,volume,beta,min_N)
    delta_p[1] = pressure[0]-pressure[1]
    print min_N, mu_new, delta_p[0]*delta_p[1]
    print 
    if delta_p[0]*delta_p[1] < 0.0:
        print 'Solution Bracketed'
#        print 'mu_ref '+str(mu_ref)
#        print 'mu_new '+str(mu_new)
        print
        solved = True
    delta_mu = 1.0005e0 * delta_mu * ( -1.00e0**iterations)

# Search for the coexistence point by a refined bracketing algorithm
if delta_p[0] < delta_p[1]:
    mu_low = mu_ref
    mu_high = mu_new
else:
    mu_low = mu_new
    mu_high = mu_ref
    [ delta_p[0],delta_p[1] ] = [ delta_p[1],delta_p[0] ] #Swap values
print 'Delta p: '+str(delta_p[0])+'   '+str(delta_p[1])
print
iterations = 0
solved = False
print '{:>10}'.format('Iteration') +'{:>15}'.format('mu*') +'{:>15}'.format('rho*(1)') +'{:>15}'.format('rho*(2)') +'{:>15}'.format('p*(1)') +'{:>15}'.format('p*(2)') +'{:>15}'.format('mu_low') +'{:>15}'.format('mu_high')
while not solved:
    iterations = iterations + 1
    mu_test = (mu_low + mu_high) / 2.0e0
    lnPi_new = Reweight_lnPi(lnPi_old,N,beta,mu,mu_test)
    minima = argrelextrema(numpy.array(lnPi_new),numpy.less)
    min_N = int(sum(minima[0])/len(minima[0])) #This is a lazy way to deal with noise in lnPi
    (nmols, pressure, U, GPFE) = Bulk_Calculate_Properties(lnPi_new,energy,N,mu_new,volume,beta,min_N)
    delta_p_new = pressure[0] - pressure[1]
    if delta_p_new < 0.0:
        delta_p[0] = delta_p_new
        mu_low = mu_test
    else:
        delta_p[1] = delta_p_new
        mu_high = mu_test
    print '{:10}'.format(iterations) +'  {: 12.6e}'.format(mu_test) +'  {: 12.6e}'.format(nmols[0]/volume) +'  {: 12.6e}'.format(nmols[1]/volume) +'  {: 12.6e}'.format(pressure[0]) +'  {: 12.6e}'.format(pressure[1]) +'  {: 12.6e}'.format(mu_low) +'  {: 12.6e}'.format(mu_high)
    if abs((mu_low-mu_high)/mu_low) < tolerance:
        solved = True

#Final Checks
lnPi_new = Reweight_lnPi(lnPi_old,N,beta,mu,mu_test)
maxima = argrelextrema(numpy.array(lnPi_new),numpy.greater)
minima = argrelextrema(numpy.array(lnPi_new),numpy.less)
# Check that the saturation lnPi has a sufficiently small high-N tail
if abs(lnPi_new[maxima[0][-1]] - lnPi_new[-1]) < lnPi_threshold:
    print 'WARNING: Low Probability Tail in lnPi is not present!'
    print 'Current mu*: ', mu_test
    print 'Local Maximum at N = '+str(maxima[0][-1])+' lnPi = '+str(lnPi_new[maxima[0][-1]])
    print 'Tail at          N = '+str(N_max)+' lnPi = '+str(lnPi_new[-1])
    print 'Total Delta: '+str(abs(lnPi_new[maxima[0][-1]] - lnPi_new[-1]))
    print 'Insufficient Macrostate Probability Data to ensure minimum total delta of '+str(lnPi_threshold)
    print 'Leaving data contents in: '+work_dir
    sys.exit()
# Confirm that the free energy barrier between phases is > kT
stable = True
for entry in maxima[0]:
    if abs(lnPi_new[entry] - lnPi_new[minima[0][0]]) < 1.0e0:
        print
        print 'WARNING: Barrier between phases is smaller than kT'
        print 'Maximum located at: '+str(entry)
        print 'Minimum located at: '+str(minima[0][0])
        print ' Barrier height: '+'{:15.10e}'.format(abs(lnPi_new[entry] - lnPi_new[minima[0]]))
        stable = False
if not stable:
    print 'Phase Coexistence is not stable'
    sys.exit()

#Output the Saturation Data
# Macrostate distribution at coexistence
lnpi_output_filename = fileprefix+'sat.lnpi.csv'
lnpi_output = open(lnpi_output_filename,'wb')
writer = csv.writer(lnpi_output)
for row in N:
    writer.writerow( (row, lnPi_new[row]) )
lnpi_output.close()
print
print 'lnPi Macrostate Distribution at Coexistence output as '+str(lnpi_output_filename)
print
# Saturation conditions
saturation_output_filename = fileprefix+'sat.txt'
saturation_output = open(saturation_output_filename,'wb')
saturation_output.write('          kT*lnZ                 p(1)*                  p(2)*             rho(1)*                rho(2)*                 u(1)*             u(2)* \n')
saturation_output.write('{:20.10e}'.format(mu_test) +'  {:20.10e}'.format(pressure[0]) +'  {:20.10e}'.format(pressure[1]) +'  {:20.10e}'.format(nmols[0]/volume) +'  {:20.10e}'.format(nmols[1]/volume) +'  {:20.10e}'.format(U[0]/nmols[0]) +'  {:20.10e}'.format(U[1]/nmols[1]) +'\n')
saturation_output.close()

print '   Saturation Conditions'
print '          kT*lnZ                 p(1)*                  p(2)*             rho(1)*                rho(2)*                 u(1)*             u(2)* '
print '{:20.10e}'.format(mu_test) +'  {:20.10e}'.format(pressure[0]) +'  {:20.10e}'.format(pressure[1]) +'  {:20.10e}'.format(nmols[0]/volume) +'  {:20.10e}'.format(nmols[1]/volume) +'  {:20.10e}'.format(U[0]/nmols[0]) +'  {:20.10e}'.format(U[1]/nmols[1]) +'\n'

#File Cleanup
shutil.rmtree(work_dir)
