# -*- coding: utf-8 -*-
"""
Created on 20160806
@author: jyoshida-sci
https://www.phys.ufl.edu/~avery/fitting.html Fitting Theory I: General Least Squares Fitting Theory
"""



import numpy as np
from scipy import optimize



def chisq(m1, dispflag=False):
    B_xi = 0.0
    m_Nuclear = 11174.866 #Carbon
    m_xi = 1321.71
    err_m_xi = 0.07
    
    initial_m1 = m1	
    m2 = 3727.380
    m3 = 2808.922
    
    mom1 = 170.961
    mom2 = 82.661
    mom3 = 166.049
    
    theta1 = np.deg2rad(44.9)
    theta2 = np.deg2rad(57.7)
    theta3 = np.deg2rad(156.2)
    
    phi1 = np.deg2rad(337.5)
    phi2 = np.deg2rad(174.9)
    phi3 = np.deg2rad(143.0)
    
    
    err_mom1 = 3.506
    err_mom2 = 5.855
    err_mom3 = 0.666
    
    err_theta1 = np.deg2rad(4.72)
    err_theta2 = np.deg2rad(8.41)
    err_theta3 = np.deg2rad(1.86)
    
    err_phi1 = np.deg2rad(4.64)
    err_phi2 = np.deg2rad(7.21)
    err_phi3 = np.deg2rad(2.05)
    
    err_theta1_phi1 = 1.74231E-5
    err_theta2_phi2 = 7.87901E-4
    err_theta3_phi3 = 5.28957E-6
    
    B_xi=0 #0.13


    #a-matrix: input value
    a0 = np.array([
        [mom1],
        [theta1],
        [phi1],
        [mom2],
        [theta2],
        [phi2],
        [mom3],
        [theta3],
        [phi3],
        [m_xi]
        ])

    #constraints derivation matrix
    D11 = mom1/np.sqrt(m1**2+mom1**2)
    D21 = np.sin(theta1)*np.cos(phi1)
    D31 = np.sin(theta1)*np.sin(phi1)
    D41 = np.cos(theta1)
    
    D12 = 0.0
    D22 = mom1*np.cos(theta1)*np.cos(phi1)
    D32 = mom1*np.cos(theta1)*np.sin(phi1)
    D42 = -mom1*np.sin(theta1)
    
    D13 = 0.0
    D23 = -mom1*np.sin(theta1)*np.sin(phi1)
    D33 = mom1*np.sin(theta1)*np.cos(phi1)
    D43 = 0.0
    
    D14 = mom2/np.sqrt(m2**2+mom2**2)
    D24 = np.sin(theta2)*np.cos(phi2)
    D34 = np.sin(theta2)*np.sin(phi2)
    D44 = np.cos(theta2)
    
    D15 = 0.0
    D25 = mom2*np.cos(theta2)*np.cos(phi2)
    D35 = mom2*np.cos(theta2)*np.sin(phi2)
    D45 = -mom2*np.sin(theta2)
    
    D16 = 0.0
    D26 = -mom2*np.sin(theta2)*np.sin(phi2)
    D36 = mom2*np.sin(theta2)*np.cos(phi2)
    D46 = 0.0
    
    D17 = mom3/np.sqrt(m3**2+mom3**2)
    D27 = np.sin(theta3)*np.cos(phi3)
    D37 = np.sin(theta3)*np.sin(phi3)
    D47 = np.cos(theta3)
    
    D18 = 0.0
    D28 = mom3*np.cos(theta3)*np.cos(phi3)
    D38 = mom3*np.cos(theta3)*np.sin(phi3)
    D48 = -mom3*np.sin(theta3)
    
    D19 = 0.0
    D29 = -mom3*np.sin(theta3)*np.sin(phi3)
    D39 = mom3*np.sin(theta3)*np.cos(phi3)
    D49 = 0.0
    
    D110 = -1.0
    D210 = 0.0
    D310 = 0.0
    D410 = 0.0
    
            
    D = np.array([
    [D11,D12,D13,D14,D15,D16,D17,D18,D19,D110],
    [D21,D22,D23,D24,D25,D26,D27,D28,D29,D210],
    [D31,D32,D33,D34,D35,D36,D37,D38,D39,D310],
    [D41,D42,D43,D44,D45,D46,D47,D48,D49,D410]
     ])


    #variance-covariance matrix 
    VA = np.zeros((10, 10))
    
    VA[0,0] = err_mom1**2     
    VA[1,1] = err_theta1**2
    VA[2,2] = err_phi1**2
    VA[1,2] = err_theta1_phi1
    VA[2,1] = err_theta1_phi1
    
    VA[3,3] = err_mom2**2
    VA[4,4] = err_theta2**2
    VA[5,5] = err_phi2**2
    VA[4,5] = err_theta2_phi2
    VA[5,4] = err_theta2_phi2
    
    VA[6,6] = err_mom3**2
    VA[7,7] = err_theta3**2
    VA[8,8] = err_phi3**2
    VA[7,8] = err_theta3_phi3
    VA[8,7] = err_theta3_phi3
    
    VA[9,9] = err_m_xi**2


    #constraints
    H1 = np.sqrt(m1**2+mom1**2)+np.sqrt(m2**2+mom2**2)+np.sqrt(m3**2+mom3**2)-(m_Nuclear + m_xi - B_xi)
    H2 = mom1*np.sin(theta1)*np.cos(phi1)+mom2*np.sin(theta2)*np.cos(phi2)+mom3*np.sin(theta3)*np.cos(phi3)
    H3 = mom1*np.sin(theta1)*np.sin(phi1)+mom2*np.sin(theta2)*np.sin(phi2)+mom3*np.sin(theta3)*np.sin(phi3)
    H4 = mom1*np.cos(theta1)+mom2*np.cos(theta2)+mom3*np.cos(theta3)
    d = np.array([
    [H1],
    [H2],
    [H3],
    [H4]
    ])

    #unknown parameter  
    E = np.array([[ m1/np.sqrt(m1**2 +mom1**2)],
                  [0.0],
                  [0.0],
                  [0.0]
                    ])


    #see
    #Applied Fitting Theory I General Least Squares Theory
    #6. Solving for Unknown Parameters in the Constraint Equations
    #1. Direct solution

    #auxiliary matricies
    Dt = D.transpose()#10x4matrix
    VD = np.linalg.inv( D.dot(VA.dot(Dt)) )#4x4matrix
    Et = E.transpose()#
    VE = np.linalg.inv( Et.dot(VD.dot(E)) )
    Vz = VE
    Vlambda = VD - VD.dot(E.dot( VE.dot( Et.dot(VD))))
    
    
    
    lambda0 = VD.dot(d)
    z = -VE.dot( Et.dot(lambda0))
    lmbda = lambda0 + VD.dot(E.dot(z))

    lambdat = lmbda.transpose()
    VD_inv = np.linalg.inv(VD)
    chisquare = lambdat.dot(VD_inv.dot(lmbda))
    chisquarevalue = chisquare[0,0]



    #Applied Fitting Theory VI Formulas for Kinematic Fitting
    #Appendix III Derivation of the Vertex Constraint Solution
    #formula(6)
    #caluclated value
    V =  VA - VA.dot(Dt.dot(Vlambda.dot(D.dot(VA))))
    a = a0 - VA.dot(Dt.dot(lmbda))

    new_mom1 = a[0,0]
    new_theta1 = a[1,0]
    new_phi1 = a[2,0]
    new_mom2 = a[3,0]
    new_theta2 = a[4,0]
    new_phi2 = a[5,0]
    new_mom3 = a[6,0]
    new_theta3 = a[7,0]
    new_phi3 = a[8,0]
    new_m_xi = a[9,0]
    
    new_err_mom1 = np.sqrt(V[0,0])
    new_err_theta1 = np.sqrt(V[1,1])
    new_err_phi1 = np.sqrt(V[2,2])
    new_err_mom2 = np.sqrt(V[3,3])
    new_err_theta2 = np.sqrt(V[4,4])
    new_err_phi2 = np.sqrt(V[5,5]) 
    new_err_mom3 = np.sqrt(V[6,6])
    new_err_theta3 = np.sqrt(V[7,7])
    new_err_phi3 = np.sqrt(V[8,8])
    new_err_m_xi = np.sqrt(V[9,9])
    
    new_m1 = initial_m1 + z[0,0]
    new_err_m1 = np.sqrt(Vz[0,0])

    new_H1 = np.sqrt(new_m1**2+new_mom1**2)+np.sqrt(m2**2+new_mom2**2)+np.sqrt(m3**2+new_mom3**2)-(m_Nuclear+new_m_xi-B_xi)
    new_H2 = new_mom1*np.sin(new_theta1)*np.cos(new_phi1)+new_mom2*np.sin(new_theta2)*np.cos(new_phi2)+new_mom3*np.sin(new_theta3)*np.cos(new_phi3)
    new_H3 = new_mom1*np.sin(new_theta1)*np.sin(new_phi1)+new_mom2*np.sin(new_theta2)*np.sin(new_phi2)+new_mom3*np.sin(new_theta3)*np.sin(new_phi3)
    new_H4 = new_mom1*np.cos(new_theta1)+new_mom2*np.cos(new_theta2)+new_mom3*np.cos(new_theta3)

    delta_H1 = d[0,0]/np.sqrt(VD_inv[0,0])
    delta_H2 = d[1,0]/np.sqrt(VD_inv[1,1])
    delta_H3 = d[2,0]/np.sqrt(VD_inv[2,2])
    delta_H4 = d[3,0]/np.sqrt(VD_inv[3,3])

    new_err_BLL = np.sqrt(4*0.006**2 + new_err_m1**2)


    if dispflag == True:
        print  "new_H1      {0:.4f}   {1:.4f}".format(new_H1,delta_H1)
        print  "new_H2      {0:.4f}   {1:.4f}".format(new_H2,delta_H2)
        print  "new_H3      {0:.4f}   {1:.4f}".format(new_H3,delta_H3)
        print  "new_H4      {0:.4f}   {1:.4f}".format(new_H4,delta_H4)
        print  ""        
        print  "chi_square       {0:.4f}".format(chisquarevalue)
        print  "momentum1 = {0:.3f} +- {1:.3f}".format(new_mom1,new_err_mom1)
        print  "momentum2 = {0:.3f} +- {1:.3f}".format(new_mom2,new_err_mom2)
        print  "momentum3 = {0:.3f} +- {1:.3f}".format(new_mom3,new_err_mom3)
        print  "theta1    = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_theta1),np.rad2deg(new_err_theta1))
        print  "phi1      = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_phi1),np.rad2deg(new_err_phi1))
        print  "theta2    = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_theta2),np.rad2deg(new_err_theta2))
        print  "phi2      = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_phi2),np.rad2deg(new_err_phi2))
        print  "theta3    = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_theta3),np.rad2deg(new_err_theta3))
        print  "phi3      = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_phi3),np.rad2deg(new_err_phi3))
        print  "m_xi      = {0:.3f} +- {1:.3f}".format(new_m_xi, new_err_m_xi)
        print  "m_D       = {0:.3f}".format(m1)
        print  "BLL       = {0:.3f}".format(1115.683 * 2.0 + 3727.380 - m1)
        print  "m_D(new)  = {0:.4f} +- {1:.4f}".format(new_m1, new_err_m1)
        print  "BLL(new)  = {0:.4f} +- {1:.4f}".format(1115.683 * 2 + 3727.380 - new_m1, new_err_BLL)

    return chisquarevalue



if __name__ == '__main__':

#    x_min = optimize.brent(f)
    m1_min = optimize.brent(chisq,brack=(5950.0, 5955.0), maxiter=100 )
    print "m1_min={0}".format(m1_min)
    
    #あらためてm1_minのときの結果を表示
    chisq(m1_min, True)

#    m1 = 5950.102826
#    chisq(m1)



    chisq(5955, True)



    
'''
Mishina's output
m1_min = 5950.102826
new_H1      0.0008   -10.3660
new_H2      -0.2159   -0.8452
new_H3      0.3284   0.0261
new_H4      -0.2957   0.8892

chi_square         1.2697
momentum1 = 170.873 ± 3.362
momentum2 = 80.295 ± 5.407
momentum3 = 166.064 ± 0.665
theta1    = 49.302 ± 2.331
phi1      = 338.882 ± 3.159
theta2    = 59.448 ± 2.829
phi2      = 173.929 ± 5.454
theta3    = 156.710 ± 1.739
phi3      = 142.807 ± 2.017
m_xi      = 1321.710 ± 0.070
m_D       = 5950.103
BLL       = 8.643
m_D(new)  = 5952.0532 ± 0.1763
BLL(new)  = 6.6928 ± 0.1767



This output
m1_min = 5954.09830062
new_H1      0.0008   11.4780
new_H2      -0.2159   -0.8452
new_H3      0.3284   0.0261
new_H4      -0.2957   0.8892

chi_square       1.2697
momentum1 = 170.873 +- 3.362
momentum2 = 80.295 +- 5.407
momentum3 = 166.064 +- 0.665
theta1    = 49.302 +- 2.331
phi1      = 338.882 +- 3.159
theta2    = 59.448 +- 2.829
phi2      = 173.929 +- 5.454
theta3    = 156.710 +- 1.739
phi3      = 142.807 +- 2.017
m_xi      = 1321.710 +- 0.070
m_D       = 5954.098
BLL       = 4.648
m_D(new)  = 5952.0532 +- 0.1762
BLL(new)  = 6.6928 +- 0.1766
'''
