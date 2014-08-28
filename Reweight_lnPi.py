#!/usr/bin/python

#Modules
import math

#-----------------------------------------------------------------
# Routine to reweight a macrostate distribution
def Reweight_lnPi(lnPi_old,N,beta,mu_old,mu_new):
    max_lnPi = lnPi_old[0]
    lnPi_star = [ ]
    for row in N:
        lnPi_temp = lnPi_old[row] + beta * float(N[row]) * (mu_new - mu_old)
        lnPi_star.append(lnPi_temp)
        if lnPi_temp > max_lnPi: max_lnPi = lnPi_temp

    lnPi_star[:] = [ entry - max_lnPi for entry in lnPi_star ]
    sum_Pi_star = sum( [math.exp(entry) for entry in lnPi_star] )
    lnPi_star[:] = [ entry - math.log(sum_Pi_star) for entry in lnPi_star ]

    return lnPi_star
#-----------------------------------------------------------------
