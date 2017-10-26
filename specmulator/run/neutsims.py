'''

python script to run neutrino sims 

'''
import sys

import env 
import data as Data

str_mneut = sys.argv[1]
if str_mneut == '0.0': 
    mneut = 0.0
nreal = int(sys.argv[2])
nredshift = int(sys.argv[3])

_ = Data.NeutHalos(mneut, nreal, nredshift)
