'''

test data.py

'''

import env 
import data as Data 

def runNeutParticles():
    ''' Test that the method NeutrinoHalos works properly 
    '''
    cat = Data.NeutParticles(0.0, 1, 4, clobber=True)  # mneut = 0.0eV, realization #1, z = 0 
    print type(cat)
    return None 


def runNeutHalos():
    ''' Test that the method NeutrinoHalos works properly 
    '''
    halocat = Data.NeutHalos(0.0, 1, 4, clobber=True)  # mneut = 0.0eV, realization #1, z = 0 
    print type(halocat) 
    return None 


if __name__=="__main__": 
    runNeutHalos()
