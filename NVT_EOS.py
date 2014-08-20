#!/usr/bin/python

#------------SOFTWARE DISCLAIMER AND REDISTRIBUTION CONDITIONS----------------
#   This software was developed at the National Institute of Standards and
#   Technology by employees of the Federal Government in the course of their
#   official duties. Pursuant to Title 17 Section 105 of the United States 
#   Code this software is not subject to copyright protection and is in the 
#   public domain. EOS.py is an experimental system. NIST assumes no
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
#   Software Name: EOS.py
#   Brief Description: Script to compute the density-pressure equation
#    of state from a particle-number probability distribution
#   Author: Daniel W. Siderius, PhD
#   Contact Information: daniel.siderius@nist.gov
#
#   Version History:
#
#   Version: 0.91
#   Release Date: 20 August 2012
#   Notes: Original version
#-----------------------------------------------------------------------------

#Modules
import csv,math,sys,os,shutil
import scipy.interpolate
import xml.dom.minidom as minidom
import CODATA_2010_Constants as CODATA

#Directories and Files
work_dir = './tmp'
input_filename = str(sys.argv[1])
meta_data = work_dir+'/meta_data.xml'

#Unpack the data bundle
print 'Unpacking TMMC data in '+str(sys.argv[1])
if not os.path.exists(work_dir): os.makedirs(work_dir)
os.system( 'cd '+work_dir+' ; tar -xzf ../'+input_filename )

#Parse the meta data
dom = minidom.parse( meta_data )
temperature = float(dom.getElementsByTagName('temperature')[0].firstChild.data)
lnZ = float(dom.getElementsByTagName('log_activity')[0].firstChild.data)
volume = float(dom.getElementsByTagName('volume')[0].firstChild.data)
fileprefix = dom.getElementsByTagName('file_prefix')[0].firstChild.data

#Convert from log(activity) to chemical potential
mu = temperature * lnZ

#Read in the Macrostate Distribution Data
macrostate_filename = work_dir+'/'+fileprefix+'lnpi.csv'
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

#Pi(N) normalization (just to be sure)
for row in N:
    lnPi[row] = lnPi[row] - math.log(sumPi)
sumPi = 1.0

#Calculation of the Helmholtz Free Energy
FofNVT = [ ]
for row in N:
    FofNVT.append( mu * float(N[row]) + temperature * (-lnPi[row] - math.log(sumPi) + lnPi[0]) )

#Calculation of the N-V-T Ensemble Chemical Potential
print 'Using SciPy spline routine for numerical derivatives'
spline_ticks = scipy.interpolate.splrep( N, FofNVT )
muofNVT = scipy.interpolate.splev( N, spline_ticks, der = 1 )
#Alternate Method for Debugging: Finite Difference Derivatives
#muofNVT = [ ]
#for row in N:
#    if row == 0:
#        muofNVT.append( (-3.0*FofNVT[0] + 4.0*FofNVT[1] - FofNVT[2])/2.0 )
#    elif row == Nmax:
#        muofNVT.append( (3.0*FofNVT[Nmax] - 4.0*FofNVT[Nmax-1] + FofNVT[Nmax-2])/2.0 )
#    else:
#        muofNVT.append( (FofNVT[row+1]-FofNVT[row-1])/2.0 )

#Calculation of the N-V-T Ensemble Pressure
pofNVT = [ ]
for row in N:
    pofNVT.append( -(FofNVT[row] - muofNVT[row]*float(row))/volume )

#Output the isotherm
isotherm_filename = fileprefix+'isotherm.csv'
isotherm_file = open(isotherm_filename,'wb')
writer = csv.writer(isotherm_file)
for row in N:
    density = float(N[row])/volume
    pressure = pofNVT[row]
    writer.writerow( (density, pressure) )
isotherm_file.close()
print 'rho-p Isotherm output as '+str(isotherm_filename)

#File Cleanup
shutil.rmtree(work_dir)
