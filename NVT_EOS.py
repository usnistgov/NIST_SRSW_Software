#!/usr/bin/python

#------------SOFTWARE DISCLAIMER AND REDISTRIBUTION CONDITIONS----------------
#   This software was developed at the National Institute of Standards and
#   Technology by employees of the Federal Government in the course of their
#   official duties. Pursuant to Title 17 Section 105 of the United States 
#   Code this software is not subject to copyright protection and is in the 
#   public domain. NVT_EOS.py is an experimental system. NIST assumes no
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
#   Software Name: NVT_EOS.py
#   Brief Description: Script to compute the density-pressure equation
#    of state from a particle-number probability distribution. This
#    method yields an isotherm that contains stable, metastable, and
#    unstable state points.
#   Author: Daniel W. Siderius, PhD
#   Contact Information: daniel.siderius@nist.gov
#
#   Computation of the equation of state is based on a method discussed by
#    Shen and Errington in J. Phys. Chem. B, 108(51):19595-19606, (2004).
#    The method is presented on page 19598 in equations 12-19.
#
#   Version History:
#
#   Version: 1.0
#   Release Date: 27 August 2012
#   Notes: Original version
#-----------------------------------------------------------------------------

#Modules
import math,sys,os,shutil, csv
import scipy.interpolate
import CODATA_2010_Constants as CODATA
from Read_TMMC_CSV_Data import *
from Renormalize_lnPi import *
from Parse_Standard_XML_Data import *

#Directories and Files
work_dir = './tmp/'
input_filename = str(sys.argv[1])
meta_data = work_dir+'metadata.xml'

# Unpack the data bundle & read the XML metadata
print 'Unpacking TMMC data in '+str(sys.argv[1])
if not os.path.exists(work_dir): os.makedirs(work_dir)
os.system( 'cd '+work_dir+' ; tar -xzf ../'+input_filename )
XML_input = Parse_Standard_XML_data(meta_data)
temperature = XML_input[0]
lnZ = XML_input[1]
volume = XML_input[2]
fileprefix = XML_input[3]
units_type = XML_input[4]
if units_type == 'Absolute':
    LJ_sigma = XML_input[5]
    LJ_epsilon_kB = XML_input[6]

#Convert from log(activity) to chemical potential
mu = temperature * lnZ

# Read in the Macrostate Probability Distribution Data, get N bounds, and normalize lnPi
(N,lnPi) = Read_TMMC_CSV_Data(fileprefix,work_dir,'lnpi')
N_max = N[-1]; N_min = N[0]
(lnPi,sumPi) = Renormalize(lnPi,shift=True) #Make sure shift=False, otherwise the EOS will be wrong

#Calculation of the Helmholtz Free Energy
FofNVT = [ ]
for row in N:
    FofNVT.append( mu * float(N[row]) + temperature * (-lnPi[row] + lnPi[0]) )

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
