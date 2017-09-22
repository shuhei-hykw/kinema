import numpy as np
import rangeenergy
from uncertainties import ufloat
from uncertainties.umath import sin,cos,sqrt,fabs


table = [
{'S':0, 'A':1, 'Z':0, 'M':939.565, 'N':'n'},
{'S':0, 'A':1, 'Z':1, 'M':938.272, 'N':'p'},
{'S':0, 'A':2, 'Z':1, 'M':1875.613, 'N':'d'},
{'S':0, 'A':3, 'Z':1, 'M':2808.922, 'N':'t'},
{'S':0, 'A':2, 'Z':0, 'M':939.565*2.0, 'N':'2n'},
{'S':0, 'A':3, 'Z':2, 'M':2808.392, 'N':'He3'},
{'S':0, 'A':4, 'Z':2, 'M':3727.380, 'N':'He4'},
{'S':0, 'A':5, 'Z':2, 'M':4667.845, 'N':'He5'},
{'S':0, 'A':5, 'Z':3, 'M':4667.624, 'N':'Li5'},
{'S':0, 'A':6, 'Z':2, 'M':5605.541, 'N':'He6'},
{'S':0, 'A':6, 'Z':3, 'M':5601.520, 'N':'Li6'},
{'S':0, 'A':6, 'Z':4, 'M':5605.298, 'N':'Be6'},
{'S':0, 'A':7, 'Z':2, 'M':6545.550, 'N':'He7'},
{'S':0, 'A':7, 'Z':3, 'M':6533.836, 'N':'Li7'},
{'S':0, 'A':7, 'Z':4, 'M':6534.186, 'N':'Be7'},
{'S':0, 'A':8, 'Z':2, 'M':7482.542, 'N':'He8'},
{'S':0, 'A':8, 'Z':3, 'M':7471.369, 'N':'Li8'},
{'S':0, 'A':8, 'Z':4, 'M':7454.852, 'N':'Be8'},
{'S':0, 'A':8, 'Z':5, 'M':7472.321, 'N':'B8'},
{'S':0, 'A':9, 'Z':3, 'M':8406.871, 'N':'Li9'},
{'S':0, 'A':9, 'Z':4, 'M':8392.753, 'N':'Be9'},
{'S':0, 'A':9, 'Z':5, 'M':8393.310, 'N':'B9'},
{'S':0, 'A':9, 'Z':6, 'M':8409.295, 'N':'C9'},
{'S':0, 'A':10, 'Z':4, 'M':9325.507, 'N':'Be10'},
{'S':0, 'A':10, 'Z':5, 'M':9324.440, 'N':'B10'},
{'S':0, 'A':10, 'Z':6, 'M':9327.580, 'N':'C10'},
{'S':0, 'A':11, 'Z':3, 'M':10285.846, 'N':'Li11'},
{'S':0, 'A':11, 'Z':4, 'M':10264.570, 'N':'Be11'},
{'S':0, 'A':11, 'Z':5, 'M':10252.550, 'N':'B11'},
{'S':0, 'A':11, 'Z':6, 'M':10254.022, 'N':'C11'},
{'S':0, 'A':12, 'Z':4, 'M':11200.922, 'N':'Be12'},
{'S':0, 'A':12, 'Z':5, 'M':11188.746, 'N':'B12'},
{'S':0, 'A':12, 'Z':6, 'M':11174.866, 'N':'C12'},
{'S':0, 'A':12, 'Z':7, 'M':11191.693, 'N':'N12'},
{'S':0, 'A':13, 'Z':4, 'M':12142.542, 'N':'Be13'},
{'S':0, 'A':13, 'Z':5, 'M':12123.434, 'N':'B13'},
{'S':0, 'A':13, 'Z':6, 'M':12109.485, 'N':'C13'},
{'S':0, 'A':13, 'Z':7, 'M':12111.195, 'N':'N13'},
{'S':0, 'A':13, 'Z':8, 'M':12128.443, 'N':'O13'},
{'S':0, 'A':14, 'Z':5, 'M':13062.023, 'N':'B14'},
{'S':0, 'A':14, 'Z':6, 'M':13040.874, 'N':'C14'},
{'S':0, 'A':14, 'Z':7, 'M':13040.207, 'N':'N14'},
{'S':0, 'A':14, 'Z':8, 'M':13044.841, 'N':'O14'},
{'S':0, 'A':15, 'Z':6, 'M':13979.222, 'N':'C15'},
{'S':0, 'A':15, 'Z':7, 'M':13968.939, 'N':'N15'},
{'S':0, 'A':15, 'Z':8, 'M':13971.182, 'N':'O15'},
{'S':0, 'A':16, 'Z':6, 'M':14914.536, 'N':'C16'},
{'S':0, 'A':16, 'Z':7, 'M':14906.014, 'N':'N16'},
{'S':0, 'A':16, 'Z':8, 'M':14895.084, 'N':'O16'},
{'S':0, 'A':17, 'Z':8, 'M':15830.506, 'N':'O17'},
{'S':0, 'A':0, 'Z':-1, 'M':139.570, 'N':'pi-'},
{'S':0, 'A':0, 'Z':0, 'M':134.977, 'N':'pi0'},
{'S':-1, 'A':1, 'Z':0, 'M':1115.683, 'N':'L'},
{'S':-1, 'A':3, 'Z':1, 'M':2991.166, 'N':'H3L'},
{'S':-1, 'A':4, 'Z':1, 'M':3922.565, 'N':'H4L'},
{'S':-1, 'A':4, 'Z':2, 'M':3921.685, 'N':'He4L'},
{'S':-1, 'A':5, 'Z':2, 'M':4839.943, 'N':'He5L'},
{'S':-1, 'A':6, 'Z':2, 'M':5779.348, 'N':'He6L'},
{'S':-1, 'A':6, 'Z':3, 'M':5778.807, 'N':'Li6L'},
{'S':-1, 'A':7, 'Z':3, 'M':6711.623, 'N':'Li7L'},
{'S':-1, 'A':7, 'Z':4, 'M':6715.821, 'N':'Be7L'},
{'S':-1, 'A':8, 'Z':2, 'M':7654.073, 'N':'He8L'},
{'S':-1, 'A':8, 'Z':3, 'M':7642.719, 'N':'Li8L'},
{'S':-1, 'A':8, 'Z':4, 'M':7643.029, 'N':'Be8L'},
{'S':-1, 'A':9, 'Z':3, 'M':8578.552, 'N':'Li9L'},
{'S':-1, 'A':9, 'Z':4, 'M':8563.825, 'N':'Be9L'},
{'S':-1, 'A':9, 'Z':5, 'M':8579.714, 'N':'B9L'},
{'S':-1, 'A':10, 'Z':4, 'M':9499.326, 'N':'Be10L'},
{'S':-1, 'A':10, 'Z':5, 'M':9500.103, 'N':'B10L'},
{'S':-1, 'A':11, 'Z':5, 'M':10429.883, 'N':'B11L'},
{'S':-1, 'A':12, 'Z':5, 'M':11356.863, 'N':'B12L'},
{'S':-1, 'A':12, 'Z':6, 'M':11358.905, 'N':'C12L'},
{'S':-1, 'A':13, 'Z':5, 'M':12293.059, 'N':'B13L'},
{'S':-1, 'A':13, 'Z':6, 'M':12278.859, 'N':'C13L'},
{'S':-1, 'A':14, 'Z':6, 'M':13212.998, 'N':'C14L'},
{'S':-1, 'A':14, 'Z':7, 'M':13214.708, 'N':'N14L'},
{'S':-1, 'A':15, 'Z':7, 'M':14142.300, 'N':'N15L'},
{'S':-1, 'A':16, 'Z':8, 'M':15074.365, 'N':'O16L'},
{'S':-1, 'A':18, 'Z':8, 'M':16931.689, 'N':'O18L'},
{'S':-2, 'A':4, 'Z':1, 'M':4106.719, 'N':'H4LL'},
{'S':-2, 'A':5, 'Z':1, 'M':5036.208, 'N':'H5LL'},
{'S':-2, 'A':5, 'Z':2, 'M':5034.978, 'N':'He5LL'},
{'S':-2, 'A':6, 'Z':2, 'M':5952.506, 'N':'He6LL'},
{'S':-2, 'A':7, 'Z':2, 'M':6890.851, 'N':'He7LL'},
{'S':-2, 'A':7, 'Z':3, 'M':6889.990, 'N':'Li7LL'},
{'S':-2, 'A':8, 'Z':3, 'M':7821.726, 'N':'Li8LL'},
{'S':-2, 'A':8, 'Z':4, 'M':7826.344, 'N':'Be8LL'},
{'S':-2, 'A':9, 'Z':2, 'M':8762.596, 'N':'He9LL'},
{'S':-2, 'A':9, 'Z':3, 'M':8751.602, 'N':'Li9LL'},
{'S':-2, 'A':9, 'Z':4, 'M':8751.872, 'N':'Be9LL'},
{'S':-2, 'A':10, 'Z':3, 'M':9685.735, 'N':'Li10LL'},
{'S':-2, 'A':10, 'Z':4, 'M':9672.798, 'N':'Be10LL'},
{'S':-2, 'A':10, 'Z':5, 'M':9687.107, 'N':'B10LL'},
{'S':-2, 'A':11, 'Z':4, 'M':10605.899, 'N':'Be11LL'},
{'S':-2, 'A':11, 'Z':5, 'M':10606.896, 'N':'B11LL'},
{'S':-2, 'A':12, 'Z':4, 'M':11538.653, 'N':'Be12LL'},
{'S':-2, 'A':12, 'Z':5, 'M':11535.326, 'N':'B12LL'},
{'S':-2, 'A':13, 'Z':5, 'M':12461.176, 'N':'B13LL'},
{'S':-2, 'A':13, 'Z':6, 'M':12463.788, 'N':'C13LL'},
{'S':-2, 'A':14, 'Z':5, 'M':13397.372, 'N':'B14LL'},
{'S':-2, 'A':14, 'Z':6, 'M':13382.852, 'N':'C14LL'},
{'S':-2, 'A':15, 'Z':6, 'M':14316.511, 'N':'C15LL'},
{'S':-2, 'A':15, 'Z':7, 'M':14318.221, 'N':'N15LL'},
{'S':-2, 'A':16, 'Z':6, 'M':15247.900, 'N':'C16LL'},
{'S':-2, 'A':16, 'Z':7, 'M':15244.393, 'N':'N16LL'},
{'S':-2, 'A':17, 'Z':8, 'M':16177.548, 'N':'O17LL'},
{'S':-2, 'A':19, 'Z':8, 'M':18032.872, 'N':'O19LL'},
{'S':-2, 'A':13, 'Z':5, 'M':12496.576, 'N':'Xi- & C12'},
{'S':-2, 'A':15, 'Z':6, 'M':14361.917, 'N':'Xi- & N14'},
{'S':-2, 'A':17, 'Z':7, 'M':16216.794, 'N':'Xi- & O16'}
]



def ke2mom(Mass, KE):
    E = Mass+KE;
    P = np.sqrt(E*E-Mass*Mass);
    return P

'''
#def two_charged(Range1,theta1,phi1,Range2,theta2,phi2):
if __name__ == '__main__':

    reen = rangeenergy.RangeEnergy()
    
    Range1 = ufloat(100, 0.25)
    theta1 = ufloat(100, 0.25)
    phi1 = ufloat(2, 0.25)
    Range2 = ufloat(2, 0.25)
    theta2 = ufloat(2, 0.25)
    phi2 = ufloat(2, 0.25)    
    
    density = 3.815
    
    fall = open('all.txt', 'w')
    fpos = open('possible.txt', 'w')

    for N0 in table:
        if N0['Z']==0:
            continue
        for N1 in table:
            if N1['Z']==0:
                continue
            for N2 in table:
                if N2['Z']==0:
                    continue
                

                if N0['A'] != N1['A'] + N2['A']:#mass number conservation
                    continue
                if N0['Z'] != N1['Z'] + N2['Z']:#charge conservation
                    continue

                if '&' in N0['N']:
                    if N0['S'] != N1['S'] + N2['S']:#Xi & nucleus
                        continue
                else:
                    if N0['S'] != N1['S'] + N2['S'] -1 :# decay of DLH or SH
                        continue

                Mass0 = N0['M']#parent nucleus                           
                Mass1 = N1['M']
                Mass2 = N2['M']                    
                Q = (Mass0 - Mass1 - Mass2)
                if Q < 0.001:
                    continue


                ke1_c = reen.KEfromRange(Mass1, Range1.n, abs(N1['Z']), density)
                mom1_c = ke2mom(Mass1, ke1_c)
                Range1_l = Range1.n - Range1.s
                ke1_l = reen.KEfromRange(Mass1, Range1_l, abs(N1['Z']), density)
                mom1_l = ke2mom(Mass1, ke1_l)
                Range1_r = Range1.n + Range1.s
                ke1_r = reen.KEfromRange(Mass1, Range1_r, abs(N1['Z']), density)
                mom1_r = ke2mom(Mass1, ke1_r)
                err_ke1 = np.fabs(ke1_l - ke1_r)/2.0;
                err_mom1 = np.fabs(mom1_l - mom1_r)/2.0                

                ke2_c = reen.KEfromRange(Mass2, Range2.n, abs(N2['Z']), density)
                mom2_c = ke2mom(Mass2, ke2_c)
                Range2_l = Range2.n - Range2.s
                ke2_l = reen.KEfromRange(Mass2, Range2_l, abs(N2['Z']), density)
                mom2_l = ke2mom(Mass2, ke2_l)
                Range2_r = Range2.n + Range2.s
                ke2_r = reen.KEfromRange(Mass2, Range2_r, abs(N2['Z']), density)
                mom2_r = ke2mom(Mass2, ke2_r)
                err_ke2 = np.fabs(ke2_l - ke2_r)/2.0;
                err_mom2 = np.fabs(mom2_l - mom2_r)/2.0                

                #ufloat
                ke1 = ufloat(ke1_c, err_ke1)
                mom1 = ufloat(mom1_c, err_mom1)
                ke2 = ufloat(ke2_c, err_ke2)
                mom2 = ufloat(mom2_c, err_mom2)
               
               
                #uncertainties.umath sin, cos, sqrt and fabs
                p1x = mom1*sin(theta1)*cos(phi1)
                p1y = mom1*sin(theta1)*sin(phi1)
                p1z = mom1*cos(theta1)                
                p2x = mom2*sin(theta2)*cos(phi2)
                p2y = mom2*sin(theta2)*sin(phi2)
                p2z = mom2*cos(theta2)
                
                cosA = (p1x*p2x+p1y*p2y+p1z*p2z)/((sqrt(p1x*p1x+p1y*p1y+p1z*p1z))*(sqrt(p2x*p2x+p2y*p2y+p2z*p2z)))             
                total_final_energy = sqrt(Mass1**2+mom1**2) + sqrt(Mass2**2+mom2**2)
                Caluculated_Mass = sqrt(total_final_energy**2 - (mom1**2 + mom2**2 + 2*mom1*mom2*cosA))

                total_mom = sqrt( (p1x+p2x)**2 + (p1y+p2y)**2 + (p1z+p2z)**2 )

                Mass_gap = fabs(Caluculated_Mass - Mass0)
                total_ke = ke1 + ke2
                QQ = fabs(Q-total_ke)
                

                #output                                
                mystr  = "{0} -> {1} + {2}************************\n\n".format(N0['N'], N1['N'], N2['N'])                    
                mystr += "KE1 = {0:9.3f} +- {1:9.3f}\n".format( ke1.s, ke1.n)
                mystr += "p1 = {0:9.3f} +- {1:9.3f}\n".format( mom1.s, mom1.n)
                mystr += "  p1x = {0:9.3f} +- {1:9.3f}\n".format( p1x.s, p1x.n)
                mystr += "  p1y = {0:9.3f} +- {1:9.3f}\n".format( p1y.s, p1y.n)
                mystr += "  p1z = {0:9.3f} +- {1:9.3f}\n".format( p1z.s, p1z.n)
                mystr += "KE2 = {0:9.3f} +- {1:9.3f}\n".format( ke2.s, ke2.n)
                mystr += "p2 = {0:9.3f} +- {1:9.3f}\n".format( mom2.s, mom2.n)
                mystr += "  p2x = {0:9.3f} +- {1:9.3f}\n".format( p2x.s, p2x.n)
                mystr += "  p2y = {0:9.3f} +- {1:9.3f}\n".format( p2y.s, p2y.n)
                mystr += "  p2z = {0:9.3f} +- {1:9.3f}\n".format( p2z.s, p2z.n)
                mystr += "\n"
                mystr += "     Calculated_Mass = {0:9.3f} +- {1:9.3f}\n".format(Caluculated_Mass.s, Caluculated_Mass.n)
                mystr += "      Estimated_Mass = {0:9.3f}\n".format(Mass0)
                mystr += "            Mass_gap = {0:9.3f}\n".format(Mass_gap)
                mystr += "         Signigfance = {0:9.3f} sigmas\n".format(Mass_gap/Caluculated_Mass.n)
                mystr += "           Total_mom = {0:9.3f} +- {1:9.3f}\n".format(total_mom.s, total_mom.n)
                mystr += "                   Q = {0:9.3f}\n".format(Q)
                mystr += "Total_Kinetic_Energy = {0:9.3f} +- {1:9.3f}\n\n".format(total_ke.s, total_ke.n)

                fall.write(mystr)
                
                if Mass_gap <= Caluculated_Mass.s*3.0 and  total_mom <= total_mom.s*3.0 and  QQ <= total_ke.s*3.0:
                    fpos.write(mystr)


    fall.close()
    fpos.close()

#    two_charged(Range1,theta1,phi1,Range2,theta2,phi2)
'''

#cn
if __name__ == '__main__':


    reen = rangeenergy.RangeEnergy()
    
    Range1 = ufloat(100, 0.25)
    theta1 = ufloat(100, 0.25)
    phi1 = ufloat(2, 0.25)
    Range2 = ufloat(2, 0.25)
    theta2 = ufloat(2, 0.25)
    phi2 = ufloat(2, 0.25)    
    Range3 = ufloat(2, 0.25)
    theta3 = ufloat(2, 0.25)
    phi3 = ufloat(2, 0.25)    
    
    density = 3.815
    
    fall = open('all.txt', 'w')
    fpos = open('possible.txt', 'w')

    for N0 in table:
        if N0['Z']==0:
            continue
        for N1 in table:
            if N1['Z']==0:
                continue
            for N2 in table:
                if N2['Z']==0:
                    continue
                for N3 in table:
                    if N3['Z']==0:
                        continue
                    for N4 in table:
                        if N4['Z']!=0:
                            continue

                        if N0['A'] != N1['A'] + N2['A'] + N3['A'] + N4['A']:#mass number conservation
                            continue
                        if N0['Z'] != N1['Z'] + N2['Z'] + N3['Z'] + N4['Z']:#charge conservation
                            continue

                        if '&' in N0['N']:#Xi & nucleus
                            if N0['S'] != N1['S'] + N2['S'] + N3['S'] + N4['S']:
                                continue
                        else:# decay of DLH or SH
                            if N0['S'] != N1['S'] + N2['S'] + N3['S'] + N4['S'] -1 :
                                continue

                        Mass0 = N0['M']#parent nucleus                           
                        Mass1 = N1['M']#charged 1
                        Mass2 = N2['M']#charged 2
                        Mass3 = N3['M']#charged 3
                        Mass4 = N4['M']#neutral particle(s)
                        Q = (Mass0 - Mass1 - Mass2 - Mass3 - Mass4)
                        if Q < 0.001:
                            continue


                        ke1_c = reen.KEfromRange(Mass1, Range1.n, abs(N1['Z']), density)
                        mom1_c = ke2mom(Mass1, ke1_c)
                        Range1_l = Range1.n - Range1.s
                        ke1_l = reen.KEfromRange(Mass1, Range1_l, abs(N1['Z']), density)
                        mom1_l = ke2mom(Mass1, ke1_l)
                        Range1_r = Range1.n + Range1.s
                        ke1_r = reen.KEfromRange(Mass1, Range1_r, abs(N1['Z']), density)
                        mom1_r = ke2mom(Mass1, ke1_r)
                        err_ke1 = np.fabs(ke1_l - ke1_r)/2.0;
                        err_mom1 = np.fabs(mom1_l - mom1_r)/2.0                

                        ke2_c = reen.KEfromRange(Mass2, Range2.n, abs(N2['Z']), density)
                        mom2_c = ke2mom(Mass2, ke2_c)
                        Range2_l = Range2.n - Range2.s
                        ke2_l = reen.KEfromRange(Mass2, Range2_l, abs(N2['Z']), density)
                        mom2_l = ke2mom(Mass2, ke2_l)
                        Range2_r = Range2.n + Range2.s
                        ke2_r = reen.KEfromRange(Mass2, Range2_r, abs(N2['Z']), density)
                        mom2_r = ke2mom(Mass2, ke2_r)
                        err_ke2 = np.fabs(ke2_l - ke2_r)/2.0;
                        err_mom2 = np.fabs(mom2_l - mom2_r)/2.0                

                        ke3_c = reen.KEfromRange(Mass3, Range3.n, abs(N3['Z']), density)
                        mom3_c = ke2mom(Mass3, ke3_c)
                        Range3_l = Range3.n - Range3.s
                        ke3_l = reen.KEfromRange(Mass3, Range3_l, abs(N3['Z']), density)
                        mom3_l = ke2mom(Mass3, ke3_l)
                        Range3_r = Range3.n + Range3.s
                        ke3_r = reen.KEfromRange(Mass3, Range3_r, abs(N3['Z']), density)
                        mom3_r = ke2mom(Mass3, ke3_r)
                        err_ke3 = np.fabs(ke3_l - ke3_r)/2.0;
                        err_mom3 = np.fabs(mom3_l - mom3_r)/2.0    
                        
                        #ufloat
                        ke1 = ufloat(ke1_c, err_ke1)
                        mom1 = ufloat(mom1_c, err_mom1)
                        ke2 = ufloat(ke2_c, err_ke2)
                        mom2 = ufloat(mom2_c, err_mom2)
                        ke3 = ufloat(ke3_c, err_ke3)
                        mom3 = ufloat(mom3_c, err_mom3)
                       
                        #uncertainties.umath sin, cos, sqrt and fabs
                        p1x = mom1*sin(theta1)*cos(phi1)
                        p1y = mom1*sin(theta1)*sin(phi1)
                        p1z = mom1*cos(theta1)                
                        p2x = mom2*sin(theta2)*cos(phi2)
                        p2y = mom2*sin(theta2)*sin(phi2)
                        p2z = mom2*cos(theta2)            
                        p3x = mom3*sin(theta3)*cos(phi3)
                        p3y = mom3*sin(theta3)*sin(phi3)
                        p3z = mom3*cos(theta3)
        
                        p4x = -(p1x + p2x + p3x)
                        p4y = -(p1y + p2y + p3y)
                        p4z = -(p1z + p2z + p3z)                
                        mom4 = sqrt( p4x**2 + p4y**2 + p4z**2)               

                        if True:#((Z4==0)&&(A4==2)&&(S4==0)):
                            ke4 = 0.5*(sqrt(Mass4**2 + mom4**2) - Mass4)
                        else:
                            ke4 = sqrt(Mass4**2 + mom4**2)-Mass4

                        cosA = (p1x*p2x+p1y*p2y+p1z*p2z)/(mom1*mom2)
                        cosB = (p2x*p3x+p2y*p3y+p2z*p3z)/(mom2*mom3)
                        cosC = (p1x*p3x+p1y*p3y+p1z*p3z)/(mom1*mom3)
                        cosD = (p1x*p4x+p1y*p4y+p1z*p4z)/(mom1*mom4)
                        cosE = (p2x*p4x+p2y*p4y+p2z*p4z)/(mom2*mom4)
                        cosF = (p4x*p3x+p4y*p3y+p4z*p3z)/(mom4*mom3)

                        total_final_energy = sqrt(Mass1**2+mom1**2) + sqrt(Mass2**2+mom2**2) + sqrt(Mass3**2+mom3**2) + sqrt(Mass4**2+mom4**2)

                        term2 = mom1**2 + mom2**2 + mom3**2 + mom4**2 
                        + 2.0*mom1*mom2*cosA
                        + 2.0*mom2*mom3*cosB
                        + 2.0*mom1*mom3*cosC
                        + 2.0*mom1*mom4*cosD
                        + 2.0*mom2*mom4*cosE
                        + 2.0*mom3*mom4*cosF

                        Caluculated_Mass_sq = total_final_energy**2  - term2                
                        Caluculated_Mass = sqrt(Caluculated_Mass_sq)
                        
                        
                        Mass_gap = fabs(Caluculated_Mass - Mass0)
                        total_ke = ke1 + ke2 + ke3 + ke4
                        QQ = fabs(Q-total_ke)

    fall.close()
    fpos.close()
