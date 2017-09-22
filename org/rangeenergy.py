# -*- coding: utf-8 -*- 
#jyoshida-sci 2015/08/31

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


class RangeEnergy:
    #constants
    Mp = 938.272#proton mass
    LMp = np.log(Mp)#log(proton mass)
    D0 = 3.815#density of standard emulsion
    r = 0.884#parameter for E07 emulsion, default_r=1.0

    #def __init__(self):


    #from Nkzw-san's Fortran-code
    #This function returns range in standard emulsion in micron units
    #@param input Mass[MeV] of a particle
    #@param input KE[MeV] of a particle
    def RangeInStandardEmulsionNk(self, Mass, KE):
        KEM = KE/Mass
        LKEM = np.log10(KEM)
        
        if(KEM<0.0001):
            return 479.210*pow(KEM,0.675899);
        else:
            Rs  =  6.05595
            Rs +=  1.38639*LKEM
            Rs += -0.302838*LKEM**2
            Rs += -0.0602134*LKEM**3
            Rs +=  0.0359347*LKEM**4
            Rs +=  0.0195023*LKEM**5
            Rs +=  0.00348314*LKEM**6
            Rs +=  0.000185264*LKEM**7
            return 10.0**Rs
            
            
    #Mishina's fitting 2014
    #This function returns KE in standard emulsion in MeV units
    #@param input LnRange_mm
    def  ProtonKEfromRangeInStandardEmulsion_part1(self, LnRange_mm):#fitted by Mishina
        LR = LnRange_mm
        LK = -2.288460778
        LK+= +1.382747508*LR
        LK+= -0.439300692*LR**2
        LK+= +0.162697682*LR**3
        LK+= -0.037735480*LR**4
        LK+= +0.005152047*LR**5
        LK+= -0.000373872*LR**6
        LK+= +0.000010917*LR**7
        LnKE_MeV = LK
        return LnKE_MeV


    def ProtonKEfromRangeInStandardEmulsion_part2(self, LnRange_mm):#fitted by Mishina
        LR = LnRange_mm   
        LK = 12.499454326
        LK+= -12.637449190*LR
        LK+= +5.296813187*LR**2
        LK+= -1.163641812*LR**3
        LK+= +0.151898030*LR**4
        LK+= -0.011803694*LR**5
        LK+= +0.000505820*LR**6
        LK+= -0.000009219*LR**7
        LnKE_MeV = LK
        return LnKE_MeV


    def ProtonKEfromRangeInStandardEmulsion_part3(self, LnRange_mm):#fitted by Mishina
        LR = LnRange_mm
        LK = -0.52629642
        LK+= +0.31555326*LR
        LK+= +0.021856192*LR**2
        LK+= +0.0012217823*LR**3
        LK+= -0.00026892371*LR**4
        LK+= +0.00001057489*LR**5
        LnKE_MeV = LK
        return LnKE_MeV

        
        
    def RangeInStandardEmulsion(self, Mass, KE):
        if KE<=0.0:
            return 0.0
        KEM = KE/Mass
        MKEM = np.log(KE*self.Mp/Mass)
        Rs = 0

        if(KEM < 0.0001):
            Rs = 479.210*pow(KEM,0.675899)
            
        elif (MKEM < 1.930606146327):
            KEfunc = lambda LnRange_mm : self.ProtonKEfromRangeInStandardEmulsion_part1(LnRange_mm) - MKEM
            initial_guess = 3.0
            solution = fsolve(KEfunc, initial_guess)
            Rs = np.exp(solution)

        elif (MKEM < 37.156634656805):
            KEfunc = lambda LnRange_mm : self.ProtonKEfromRangeInStandardEmulsion_part2(LnRange_mm) - MKEM
            initial_guess = 6.0
            solution = fsolve(KEfunc, initial_guess)
            Rs = np.exp(solution)

        else:
            KEfunc = lambda LnRange_mm : self.ProtonKEfromRangeInStandardEmulsion_part3(LnRange_mm) - MKEM
            initial_guess = 10.0
            solution = fsolve(KEfunc, initial_guess)
            Rs = np.exp(solution)

        return Rs

        
        
    #Mishina's original function
    def FunctionRs(self, KE, Mass):
        KEM = KE/Mass
        MKEM = np.log(KE*self.Mp/Mass)

        dd = 0.00001#;//step for italation

        if(KEM < 0.0001):
            Rs = 479.210*pow(KEM,0.675899)
            
        elif (MKEM < 1.930606146327):
            d0 = 3.0000
            y0 = self.Rs_function1(d0)
            while abs(MKEM-y0)>0.00001:
                d0 = (d0 + dd) if (MKEM > y0) else (d0 - dd)
                y0 = self.Rs_function1(d0)
            Rs = np.exp(d0)

        elif (MKEM < 37.156634656805):
            d0 = 6.0000
            y0 = self.Rs_function2(d0)
            while abs(MKEM-y0)>0.00001:
                d0 = (d0 + dd) if (MKEM > y0) else (d0 - dd)
                y0 = self.Rs_function2(d0)
            Rs = np.exp(d0)

        else:
            d0 = 10.0000
            y0 = self.Rs_function3(d0)
            while abs(MKEM-y0)>0.00001:
                d0 = (d0 + dd) if (MKEM > y0) else (d0 - dd)
                y0 = self.Rs_function3(d0)
            Rs = np.exp(d0)

        return Rs



    #RsRwRatio fitted by Dr.Tovee and Dr.Gajewski
    def FunctionRsRwRatio(self, Rs):
        LRs = np.log(Rs)        
        rate = -0.107714711
        rate += -0.332543998*LRs
        rate += +0.141029694*LRs**2
        rate += -0.044679440*LRs**3
        rate += +0.008162611*LRs**4
        rate += -0.000830409*LRs**5
        rate += +0.000044038*LRs**6
        rate += -0.000000951*LRs**7
        return  np.exp(rate)



    #Cz fitted by Mishina
    def FunctionCz(self, Z, beta):
        if(Z==1):
            return 0.0
        FX = 137.0*beta/Z;

        if FX<=0.5:# //regionI: a*FX^b
            return  0.168550736771407*pow(FX,1.90707106569386)
        elif FX<=2.51:# regionII: polinominal7
            val =  0.002624371
            val += -0.081622520*FX
            val += +0.643381535*FX**2
            val += -0.903648583*FX**3
            val += +0.697505012*FX**4
            val += -0.302935572*FX**5
            val += +0.067662990*FX**6
            val += -0.006004180*FX**7
            return val
        else:# regionIII: constant
            return 0.217598079611354



    #Energy->Range calculation
    def RangeFromKE(self, Mass, KE, Z, densityEM):
        if(KE <= 0.0):
            return 0.0
 
        #range as proton in standard emulsion 
        Rs = self.RangeInStandardEmulsion(Mass, KE)#Mishina's fitting function
        #Rs = self.RangeInStandardEmulsionNk(Mass, KE)#Nakazawa-san's fitting function

        #correction for range    
        ratio = self.FunctionRsRwRatio(Rs)# Rs/Rw ratio
        F = densityEM/self.D0 + ((self.r*(self.D0-densityEM))/(self.r*self.D0-1.0))*ratio #factor for range
        Rp = Rs/F#range as proton in this emulsion
        
        #calculating Cz
        E = Mass + KE #total energy
        P = np.sqrt(E*E-Mass*Mass) #momentum norm
        beta = P/E #beta of particle
        Cz = self.FunctionCz(Z, beta)
        
        #correction factors
        CPS = 1
        CPM = 1
        CF = 1

        #Range        
        R1 = CPS * (Mass/self.Mp) / (Z*Z) * Rp
        R2 = CPM * (Mass/self.Mp) * pow(Z,2.0/3.0) * Cz #R_ext
        R = (R1+R2)/CF
        
        return R



    #This is the inverse-function of RangeFromKineticEnergy
    def KEfromRange(self, Mass, Range, Z, densityEM):
        if Range<=0.0:
            return 0.0
        
        Rfunc = lambda KE : self.RangeFromKE(Mass,KE,Z,densityEM) - Range

        initial_guess = 1.0#any positive number
        solution = fsolve(Rfunc, initial_guess)

        return solution
        
        
        

    ######
    ######for Display   

    #1   
    def DisplayRangeInStandardEmulsionNk(self,Mass):
        KEs = []
        Ranges = []
        for KE in range(100):
            KEs.append( KE )
            Ranges.append( self.RangeInStandardEmulsionNk(Mass, KE) )
            
        plt.plot(KEs, Ranges)
        plt.xlabel("KE")
        plt.ylabel("Range")
        plt.grid()
        plt.show()
        
        
    #2
    def DisplayProtonKEfromRangeInStandardEmulsion(self):
        Ranges = []  
        LRs = []  
        Rs1s = []
        Rs3s = []
        for Range in range(1,3000):
            Ranges.append(Range)
            LR = np.log(Range)
            LRs.append(LR)
            Rs1s.append( np.exp( re.ProtonKEfromRangeInStandardEmulsion_part1(LR) ) )
            Rs3s.append( np.exp( re.ProtonKEfromRangeInStandardEmulsion_part3(LR) ) )
    
        #2は定数項が+12でなので、ゼロ近傍でものすごいことになる
        Ranges2 = []
        LRs2 = []
        Rs2s = []
        for Range in range(10,3000):
            Ranges2.append(Range)
            LR = np.log(Range)
            LRs2.append(LR)
            Rs2s.append( np.exp( re.ProtonKEfromRangeInStandardEmulsion_part2(LR) ) )
    
        plt.plot(Ranges, Rs1s)
        plt.plot(Ranges2, Rs2s)
        plt.plot(Ranges, Rs3s)
        plt.xlabel("Range[mm]")
        plt.ylabel("KE")
        plt.grid()
        plt.show()
        

        
    
    #3
    def DisplayRangeInStandardEmulsion(self, Mass):
        KEs = []
        Ranges = []
        for KE in range(260):
            KEs.append( KE )
            Range = re.RangeInStandardEmulsion(Mass, KE)
            Ranges.append( Range / 10000 )#in cm units
        plt.plot(KEs, Ranges)
        plt.xlabel("KE")
        plt.ylabel("RangeInStandardEmulsion [cm]")
        plt.grid()
        plt.show()
    
    #4----
    
    #5    
    def DisplayFunctionRsRwRatio(self):    
        Rss = []
        RsRws = []
        for Rs in range(500):
            Rss.append( Rs )
            RsRw = re.FunctionRsRwRatio(Rs) 
            RsRws.append( RsRw )
            
        plt.plot(Rss, RsRws)
        plt.xlabel("Rs")
        plt.ylabel("RsRw")
        plt.grid()
        plt.show()
        
        
    
    #6    
    def DisplayFunctionCz(self):
        czs = []
        betas = []
        for i in range(50):
            beta = i*0.001
            betas.append( beta )
            czs.append( re.FunctionCz(2, beta) )
            
        plt.plot(betas, czs)
        plt.xlabel("beta")
        plt.ylabel("Cz")
        plt.grid()
        plt.show()
    
    
    #7
    def DisplayRangeFromKE(self, Mass, Z, densityEM ):
        KEs = []
        Ranges = []
        for KE in range(100):
            KEs.append( KE )
            Range = re.RangeFromKE(Mass, KE, Z, densityEM)
            Ranges.append( Range )
            
        plt.plot(KEs, Ranges)
        plt.xlabel("KE")
        plt.ylabel("range")
        plt.grid()
        plt.show()
        
    
    #8
    def DisplayKEfromRange(self, Mass, Z, densityEM):
        KEs = []
        Ranges = []
        for Range in range(500):
            Ranges.append( Range )
            KE = re.KEfromRange(Mass, Range, Z, densityEM)
            KEs.append( KE )
            
        plt.plot(Ranges, KEs)
        plt.xlabel("Range")
        plt.ylabel("KE")
        plt.grid()
        plt.show()
    
    
    
    
    #usage
    def UsageOfProtonKEfromRangeInStandardEmulsions(self):
        myr_cm = 4.1*10**-1#[cm]
        myr_mm = myr_cm*10
        myke = np.exp( re.ProtonKEfromRangeInStandardEmulsion_part1(np.log(myr_mm)) )
        print '{1:.2f} MeV  {0:.2f}cm'.format(myr_cm,myke)
    
        myr_cm = 14.4*10**-1#[cm]
        myr_mm = myr_cm*10
        myke = np.exp( re.ProtonKEfromRangeInStandardEmulsion_part1(np.log(myr_mm)) )
        print '{1:.2f} MeV  {0:.2f}cm'.format(myr_cm,myke)
        
        myr_cm = 18.7*10**-1#[cm]
        myr_mm = myr_cm*10
        myke = np.exp( re.ProtonKEfromRangeInStandardEmulsion_part1(np.log(myr_mm)) )
        print '{1:.2f} MeV  {0:.2f}cm'.format(myr_cm,myke)    
    
    
    
if __name__ == "__main__":

    re = RangeEnergy()

    
    re.DisplayRangeInStandardEmulsionNk(re.Mp)
    re.DisplayProtonKEfromRangeInStandardEmulsion()
    re.DisplayRangeInStandardEmulsion(re.Mp)
    re.DisplayFunctionRsRwRatio()
    re.DisplayFunctionCz()
    re.DisplayRangeFromKE(re.Mp, 1, re.D0)   
    re.DisplayKEfromRange(re.Mp, 1, re.D0)


    print re.KEfromRange(938.272, 0.10, 1, 3.6)
    print "0.0052755"
    print re.KEfromRange(938.272, 500.10, 1, 3.6)
    print "9.2231038"
    print re.KEfromRange(938.272, 1000.10, 1, 3.6)
    print "13.8187089"
    print re.KEfromRange(938.272, 1500.10, 1, 3.6)
    print "17.4634169"
    
    print re.RangeFromKE(938.272, 10, 1, 3.6)

