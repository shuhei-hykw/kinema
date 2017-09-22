# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 13:48:51 2016
@author: jyoshida-sci, MyintKyawSoe
"""

import numpy as np


if __name__ == '__main__':

    a0 = []
    a0r = []#radian
    SM = []
    
    f = open("aaaa.txt","r")#Input data is p1, p2, theta1, theta2, phi1, phi2, and their uncertainties.    
    for row in f:
        print row
        vals = row.split()
        a0.append( float(vals[0]) )
        a0r.append( np.deg2rad(float(vals[0])))
        SM.append( float(vals[1]) )
    f.close()


    A0 = np.array([a0[0], a0[1], a0[2], a0[3], a0[4], a0[5]])
    Vr0 = np.zeros((6,6))
    for i in range(6):
        Vr0[i,i] = SM[i]**2


    D = np.array([
    [
    np.sin(a0r[1])*np.cos(a0r[2]),
    a0[0]*np.cos(a0r[1])*np.cos(a0r[2]),
    -a0[0]*np.sin(a0r[1])*np.sin(a0r[2]),
    np.sin(a0r[4])*np.cos(a0r[5]),
    a0[3]*np.cos(a0r[4])*np.cos(a0r[5]),
    -a0[3]*np.sin(a0r[4])*np.sin(a0r[5])
    ],[
    np.sin(a0r[1])*np.sin(a0r[2]),
    a0[0]*np.cos(a0r[1])*np.sin(a0r[2]),
    a0[0]*np.sin(a0r[1])*np.cos(a0r[2]),
    np.sin(a0r[4])*np.sin(a0r[5]),
    a0[3]*np.cos(a0r[4])*np.sin(a0r[5]),
    a0[3]*np.sin(a0r[4])*np.cos(a0r[5])
    ],[
    np.cos(a0r[1]),
    -a0[0]*np.sin(a0r[1]),
    0.0,
    np.cos(a0r[4]),
    -a0[3]*np.sin(a0r[4]),
    0.0
    ]
    ])
            
    d =  np.array([
            [a0[0]*np.sin(a0r[1])*np.cos(a0r[2])+a0[3]*np.sin(a0r[4])*np.cos(a0r[5])],
            [a0[0]*np.sin(a0r[1])*np.sin(a0r[2])+a0[3]*np.sin(a0r[4])*np.sin(a0r[5])],
            [a0[0]*np.cos(a0r[1])+a0[3]*np.cos(a0r[4])]
            ]);

    Dt = D.transpose()
    
    Vd0 = D.dot(Vr0.dot(Dt))

    VD = np.linalg.inv(Vd0)
   
    lmbda = VD.dot(d)

    lmbda_t = lmbda.transpose()

    a = A0 - Vr0.dot(Dt.dot(lmbda))

    A = []
    Ar = []#radian
    for k in range(6):
	   A.append( a[0,k])
	   Ar.append( np.deg2rad(a[0,k]))


    Vr = Vr0 - Vr0.dot(Dt.dot(VD.dot(D.dot(Vr0))))

    chi2 = lmbda_t.dot(d)
    

    dd = np.array([
    [A[0]*np.sin(Ar[1])*np.cos(Ar[2])+A[3]*np.sin(Ar[4])*np.cos(Ar[5])],
    [A[0]*np.sin(Ar[1])*np.sin(Ar[2])+A[3]*np.sin(Ar[4])*np.sin(Ar[5])],
    [A[0]*np.cos(Ar[1])+A[3]*np.cos(Ar[4])]
    ])


    print "track1 parameters"
    print "P\t{0}\t+-{1}".format(a[0,0], np.sqrt(Vr[0,0]))
    print "theta\t{0}\t+-{1}".format(a[0,1], np.sqrt(Vr[1,1]))
    print "phi\t{0}\t+-{1}".format(a[0,2], np.sqrt(Vr[2,2]))
    print ""    
    print "track2 parameters"
    print "P\t{0}\t+-{1}".format(a[0,3],np.sqrt(Vr[3,3]))
    print "theta\t{0}\t+-{1}".format( a[0,4], np.sqrt(Vr[4,4]))
    print "phi\t{0}\t+-{1}".format(a[0,5], np.sqrt(Vr[5,5]))
    print "";
    print "the constraints equations";
    print dd
    print "chi2";
    print chi2


'''
output of original CPP code
track2 parameters
P       100     +-0.707107
theta   90      +-0.707107
phi     -1.7949e-009    +-0.707107
track3 parameters
P       100     +-0.707107
theta   90      +-0.707107
phi     180     +-0.707107
the constraints equations are
[0;
 3.527139231042815e-007;
 3.527139147759126e-007]
'''
