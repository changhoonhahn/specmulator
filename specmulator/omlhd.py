import numpy as np 

def omlhd(k, n):
    ''' 
    '''
    nsearch=20;             # number of different initials   
    p=15;
    tfac=0.95;
    eps=0.0000001;
    nterp=150; #//10*19*4;   // parameter for Simulated Annealing 
    bestyet=100000;

    xtry = np.empty((k, n))
    xbest =np.empty((k, n))
    x = np.empty((k, n))

    dbar = (n+1) * k * 3**-1
    c = 0.5
    p0 = 0.99
    pf = 0.01 
    t1 = (pow(dbar,(-p))*n*(n-1)) * pow(4*c*(p-1),(-1))
    t2=pow((1-c),(1-p))
    t3=pow((1+c),(1-p))
    phi0=pow((t1*(t2-t3)),(pow(p,(-1))))
    delt0=pow((pow(phi0,p)-pow((1-c)*dbar,(-p))+pow((1-c)*dbar-1,(-p))-pow((1+c)*dbar,(-p))+pow((1+c)*dbar+1,(-p))),pow(p,(-1)))

    t0=-delt0*pow((np.log(p0)),(-1))

    for loop in range(5): 
        nterp += 100 
        print('nterp=', nterp)
        print('t0=', t0)
        x = LHD(n, k)


    return None


def LHD(n, k): 
    r = np.empty(n)
    lhd = np.empty((k, n))

    for j in range(n):

    return lhd

if __name__=="__main__": 
    omlhd(3, 9)
