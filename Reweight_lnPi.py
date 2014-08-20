#!/usr/bin/python

#Modules
import csv,math,sys,os

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
