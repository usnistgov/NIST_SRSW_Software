#!/usr/bin/python

import math

#-----------------------------------------------------------------
# Routine to Normalize and, if desired, shift lnPi
def Renormalize(lnPi,shift=False):
    sumPi = 0.0e0
    for entry in lnPi:
        sumPi = sumPi + math.exp(entry)
    lnPi[:] = [x - math.log(sumPi) for x in lnPi]
    if shift:
        maxval = max(lnPi)
        lnPi[:] = [x - maxval for x in lnPi]
    sumPi = sum( [math.exp(entry) for entry in lnPi] )
      
    return lnPi, sumPi
#-----------------------------------------------------------------
