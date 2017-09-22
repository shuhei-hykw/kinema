# -*- coding: utf-8 -*-
"""
Created on 20160806
@author: jyoshida-sci

this code ni ha BUG ga arimasu
chisq doesn't change even if be are changed
"""

import numpy as np
from scipy import optimize



def chisq(be, dispflag=False):
    
    B_xi = be
    initial_B_xi = be

    m_C = 13040.202#;//11174.866;//14895.084;//

    m_xi = 1321.71;
    m1 = 9499.323;
    m2 = 4839.942;

    err_m_xi = 0.07;
    err_m1 = 0.22;
    err_m2 = 0.02;

    mom1 = 342.477#;//359.069;//342.478;//359.069;//342.466;//359.126;//359.184;//342.573;//358.985;//344.514;//362.939;//363.004;
    mom2 = 342.609#;//342.071;//342.613;//342.076;//342.600;//342.003;//342.003;//342.717;//342.215;//344.647;//344.131;//344.294;
    theta1 = np.deg2rad(136.974)#;//135.4;//137.785;//136.2;
    theta2 = np.deg2rad(42.951)#;//42.5;//42.131;//41.7;
    phi1 = np.deg2rad(13.381)#;//13.3;//13.122;//13.3;
    phi2 = np.deg2rad(193.386)#;//193.5;//193.103;//193.5;

    err_mom1 = 8.016#;//8.015;//7.960;//9.543;//9.347;//7.9621;//7.9711;//7.9665;//8.051;//8.056;
    err_mom2 = 1.444#;//1.444;//1.509;//3.490;//1.637;//1.4847;//1.4893;//1.3175;//1.3495;//1.351;
    err_theta1 = np.deg2rad(4.4)#;//4.34;
    err_theta2 = np.deg2rad(1.8)#;//1.77;
    err_phi1 = np.deg2rad(4.1)#;//4.10;
    err_phi2 = np.deg2rad(1.7)#;//1.58;
    
    err_theta1_phi1 = -0.000866310#;//-0.0001263;
    err_theta2_phi2 = -0.000124961#;//-0.0008770;


    #a-matrix: input value
    a0 = np.array([
        [mom1],
        [theta1],
        [phi1],
        [mom2],
        [theta2],
        [phi2],
        [m_xi],
        [m1],
        [m2]
        ])
        


    D11 = mom1/np.sqrt(m1**2+mom1**2)
    D21 = np.sin(theta1)*np.cos(phi1)
    D31 = np.sin(theta1)*np.sin(phi1)
    D41 = np.cos(theta1)
    
    D12 = 0.0
    D22 = mom1*np.cos(theta1)*np.cos(phi1)
    D32 = mom1*np.cos(theta1)*np.sin(phi1)
    D42 = -mom1*np.sin(theta1)
        
    D13 = 0;
    D23 = -mom1*np.sin(theta1)*np.sin(phi1)
    D33 = mom1*np.sin(theta1)*np.cos(phi1)
    D43 = 0;
    
    D14 = mom2/np.sqrt(m2**2+mom2**2)
    D24 = np.sin(theta2)*np.cos(phi2)
    D34 = np.sin(theta2)*np.sin(phi2)
    D44 = np.cos(theta2)
    
    D15 = 0;
    D25 = mom2*np.cos(theta2)*np.cos(phi2)
    D35 = mom2*np.cos(theta2)*np.sin(phi2)
    D45 = -mom2*np.sin(theta2)
    
    D16 = 0.0
    D26 = -mom2*np.sin(theta2)*np.sin(phi2)
    D36 = mom2*np.sin(theta2)*np.cos(phi2)
    D46 = 0.0
    
    D17 = -1;
    D27 = 0;
    D37 = 0;
    D47 = 0;
    
    D18 = m1 / np.sqrt(m1**2 + mom1**2);
    D28 = 0;
    D38 = 0;
    D48 = 0;
    
    D19 = m2 / np.sqrt( m2**2 + mom2**2);
    D29 = 0;
    D39 = 0;
    D49 = 0;

    D = np.array([
    [D11,D12,D13,D14,D15,D16,D17,D18,D19],
    [D21,D22,D23,D24,D25,D26,D27,D28,D29],
    [D31,D32,D33,D34,D35,D36,D37,D38,D39],
    [D41,D42,D43,D44,D45,D46,D47,D48,D49]
     ])


    #variance-covariance matrix 
    Va = np.zeros((9, 9))
    
    Va[0,0] = err_mom1**2     
    Va[1,1] = err_theta1**2
    Va[2,2] = err_phi1**2
    Va[1,2] = err_theta1_phi1
    Va[2,1] = err_theta1_phi1
    
    Va[3,3] = err_mom2**2
    Va[4,4] = err_theta2**2
    Va[5,5] = err_phi2**2
    Va[4,5] = err_theta2_phi2
    Va[5,4] = err_theta2_phi2
    
    Va[6,6] = err_m_xi**2
    Va[7,7] = err_m1**2
    Va[8,8] = err_m2**2
    

    #unknown parameter  
    E = np.array([
                [1.0],
                [0.0],
                [0.0],
                [0.0]
                ])

    #constraints
    H1 = np.sqrt(m1**2+mom1**2)+np.sqrt(m2**2+mom2**2)-(m_C + m_xi - B_xi)
    H2 = mom1*np.sin(theta1)*np.cos(phi1)+mom2*np.sin(theta2)*np.cos(phi2)
    H3 = mom1*np.sin(theta1)*np.sin(phi1)+mom2*np.sin(theta2)*np.sin(phi2)
    H4 = mom1*np.cos(theta1)+mom2*np.cos(theta2)
    d = np.array([
    [H1],
    [H2],
    [H3],
    [H4]
    ])


    #solve
    Dt = D.transpose()#10x4matrix
    Vd0 = D.dot(Va.dot(Dt))#4x4matrix
    Vd = np.linalg.inv(Vd0)#4x4matrix
    E_t =  E.transpose()#
    VE0 = E_t.dot(Vd.dot(E))
    VE = np.linalg.inv(VE0)
    lambda0 = Vd.dot(d)
    z = -VE.dot( E_t.dot(lambda0))
    lmbda = lambda0 + Vd.dot(E.dot(z))
    lambda_t = lmbda.transpose()
    Vd_i = np.linalg.inv(Vd)
    chi_square = lambda_t.dot(Vd_i.dot(lmbda))

    chi_square_value = chi_square[0,0]
    
    
	
    #caluclated value
    a = a0 - Va.dot(Dt.dot(lmbda))
    Vz = VE
    V_lambda = Vd - Vd.dot(E.dot( VE.dot( E_t.dot(Vd))))
    V =  Va - Va.dot(Dt.dot(V_lambda.dot(D.dot(Va))))

    new_mom1 = a[0,0]
    new_theta1 = a[1,0]
    new_phi1 = a[2,0]
    new_mom2 = a[3,0]
    new_theta2 = a[4,0]
    new_phi2 = a[5,0]
    new_m_xi = a[6,0]
    new_m1 = a[7,0]
    new_m2 = a[8,0]
    
    new_err_mom1 = np.sqrt(V[0,0])
    new_err_theta1 = np.sqrt(V[1,1])
    new_err_phi1 = np.sqrt(V[2,2])
    new_err_mom2 = np.sqrt(V[3,3])
    new_err_theta2 = np.sqrt(V[4,4])
    new_err_phi2 = np.sqrt(V[5,5]) 
    new_err_m_xi = np.sqrt(V[6,6])
    new_err_m1 = np.sqrt(V[7,7])
    new_err_m2 = np.sqrt(V[8,8])
    
    new_B_xi = initial_B_xi + z[0,0]
    new_err_B_xi = np.sqrt(Vz[0,0])

    new_H1 = np.sqrt(new_m1**2+new_mom1**2)+np.sqrt(m2**2+new_mom2**2)-(m_C + new_m_xi - new_B_xi)
    new_H2 = new_mom1*np.sin(new_theta1)*np.cos(new_phi1)+new_mom2*np.sin(new_theta2)*np.cos(new_phi2)
    new_H3 = new_mom1*np.sin(new_theta1)*np.sin(new_phi1)+new_mom2*np.sin(new_theta2)*np.sin(new_phi2)
    new_H4 = new_mom1*np.cos(new_theta1)+new_mom2*np.cos(new_theta2)

    delta_H1 = d[0,0]/np.sqrt(Vd_i[0,0])
    delta_H2 = d[1,0]/np.sqrt(Vd_i[1,1])
    delta_H3 = d[2,0]/np.sqrt(Vd_i[2,2])
    delta_H4 = d[3,0]/np.sqrt(Vd_i[3,3])

    if dispflag == True:
        print  "new_H1      {0:.4f}   {1:.4f}".format(new_H1,delta_H1)
        print  "new_H2      {0:.4f}   {1:.4f}".format(new_H2,delta_H2)
        print  "new_H3      {0:.4f}   {1:.4f}".format(new_H3,delta_H3)
        print  "new_H4      {0:.4f}   {1:.4f}".format(new_H4,delta_H4)
        print  ""        
        print  "chi_square       {0:.8f}".format(chi_square_value)
        print  "m1        = {0:.3f} +- {1:.3f}".format(new_m1,new_err_m1)
        print  "momentum1 = {0:.3f} +- {1:.3f}".format(new_mom1,new_err_mom1)
        print  "theta1    = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_theta1),np.rad2deg(new_err_theta1))
        print  "phi1      = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_phi1),np.rad2deg(new_err_phi1))
        print  "m2        = {0:.3f} +- {1:.3f}".format(new_m2,new_err_m2)
        print  "momentum2 = {0:.3f} +- {1:.3f}".format(new_mom2,new_err_mom2)
        print  "theta2    = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_theta2),np.rad2deg(new_err_theta2))
        print  "phi2      = {0:.3f} +- {1:.3f}".format(np.rad2deg(new_phi2),np.rad2deg(new_err_phi2))
        print  "m_xi      = {0:.3f} +- {1:.3f}".format(new_m_xi, new_err_m_xi)
        print  "B_xi(new) = {0:.5f} +- {1:.5f}".format(new_B_xi, new_err_B_xi);

    return chi_square_value




if __name__ == '__main__':
    
#    x_min = optimize.brent(f)
    be_min = optimize.brent(chisq, brack=(3.0, 10.0), maxiter=100)
    print "be_min={0}".format(be_min)
    chisq(be_min, True)

#    be = 7.474889
#    chisq(be, True)


    print ""
    chisq(1.0, True)
    print ""
    chisq(10.0, True)
    print ""
    chisq(50.0, True)
    print "Bug!! The same chi-square!!!"


'''
Mishina's output
be_minimum=7.474889
chi_square_minimum=0.000520

new_H1      0.0000
new_H2      -0.0002
new_H3      -0.0001
new_H4      0.0001

chi_square         0.0005
m1        = 9499.323 ± 0.000<-Mishina's bug
momentum1 = 342.605 ± 1.421
theta1    = 137.038 ± 1.658
phi1      = 13.383 ± 1.561
m2        = 4839.942 ± 0.000<-Mishina's bug
momentum2 = 342.605 ± 1.421
theta2    = 42.962 ± 1.657
phi2      = 193.383 ± 1.562
m_xi      = 1321.710 ± 0.070
B_xi(new)  = 4.35995 ± 0.27677


This output
be_min=39.5560576278 <- ????
new_H1      0.0000   91.6492
new_H2      -0.0002   0.0113
new_H3      -0.0001   0.0019
new_H4      0.0001   0.0198

chi_square       0.00051968
m1        = 9499.323 +- 0.220
momentum1 = 342.605 +- 1.421
theta1    = 137.038 +- 1.658
phi1      = 13.383 +- 1.561
m2        = 4839.942 +- 0.020
momentum2 = 342.605 +- 1.421
theta2    = 42.962 +- 1.657
phi2      = 193.383 +- 1.562
m_xi      = 1321.710 +- 0.070
B_xi(new) = 4.35995 +- 0.27677

'''