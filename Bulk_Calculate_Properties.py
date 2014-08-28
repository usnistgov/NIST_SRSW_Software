#!/usr/bin/python

#Modules
import math

#-----------------------------------------------------------------
# Routine to compute thermodynamic properties from lnPi
# Assumes existence of *at most* two phases
def Bulk_Calculate_Properties(lnPi,energy,N,mu,volume,beta,min_N):
    norm_phase = [ 0.0e0 ] * 2
    nmols = [ 0.0e0 ] * 2
    internal_energy = [ 0.0e0 ] * 2
    for Ni in N:
        Prob_i = math.exp(lnPi[Ni])
        if Ni < min_N:
            phase_ID = 0
        else:
            phase_ID = 1
        nmols[phase_ID] = nmols[phase_ID] + float(Ni) * Prob_i
        norm_phase[phase_ID] = norm_phase[phase_ID] + Prob_i
        internal_energy[phase_ID] = internal_energy[phase_ID] + energy[Ni]*Prob_i
    pressure = [ 0.0e0 ] * 2
    free_energy = [ 0.0e0 ] * 2
    for phase_ID in range(2):
        if norm_phase[phase_ID] > 0.0:
            nmols[phase_ID] = nmols[phase_ID] / norm_phase[phase_ID]
            internal_energy[phase_ID] = internal_energy[phase_ID] / norm_phase[phase_ID]
            free_energy[phase_ID] = -(math.log(norm_phase[phase_ID]) - lnPi[0]) / beta
            pressure[phase_ID] = -free_energy[phase_ID] / volume
        else:
            nmols[phase_ID] = 0.0e0
            pressure[phase_ID] = 0.0e0
            internal_energy[phase_ID] = 0.0e0
    return nmols, pressure, internal_energy, free_energy
#-----------------------------------------------------------------
