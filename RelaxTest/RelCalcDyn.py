#!/usr/bin/env python3



from RelCalc import RelCalc
from scipy.optimize import leastsq
import copy,os,pickle,numpy,sys,math

############################################################
# Code to calculate relaxation rates in a NH spin system

# To run the code, type 'RelCalcRun RelCalc_NH.py' in the terminal
# both isotropic and axially symmetric tumbling models will be analysed.
# relaxation rates for specified values will be shown

#loops over isotropic and axial diffusion
#adding all combinations of external proton and external deuterons

############################################################

#R1 = d^2 (1/10 J(wH-wN) + 3/10 J(wN) + 3/5 J(wH+wN) )
#R2 = d^2 (1/5 J(0) + 1/20 J(wH-wN) + 3/20 J(wN) + 3/10 J(wH) + 3/10 J(wH+wN) )
#sig= d^2 ( -1/10 J(wH-wN) + 3/5 J(wH+wN) )
#NOE = 1 + gammaH/gammaN sig/R1

#R1 - sig = d^2 (1/5 J(wH-wN) + 3/10 J(wN)  )
#R1 = d^2 (1/5 J(wH-wN) + 3/10 J(wN)  ) + sig

#R2 -sig/2 = d^2 (1/5 J(0) + 1/10 (wH-wN)  + 3/20 J(wN) + 3/10 J(wH)  )
#R2 = d^2 (1/5 J(0) + 1/10 (wH-wN)   3/20 J(wN) + 3/10 J(wH) ) + sig/2

#R1 = d^2 (           1/5  J(wH-wN) + 3/10 J(wN)  ) + sig
#R2 = d^2 (1/5 J(0) + 1/10 J(wH-wN) + 3/20 J(wN) + 3/10 J(wH)   ) + sig/2

#R2 -sig/2  - R1/2 = d^2 (1/5 J(0)  + 3/10 J(wH)   )  -sig
#R2 +sig/2  - R1/2 = d^2 (1/5 J(0)  + 3/10 J(wH)   )
#R2 +1/2 ( sig-R1)  = d^2 (1/5 J(0)  + 3/10 J(wH)   )  




#sig = (NOE -1 )*gammaN/gammaH *R1

#sig= d^2 ( -1/10 J(wH-wN) + 3/5 J(wH+wN) )
#NOE = 1 + gammaH/gammaN sig/R1



EXTH=False
EXTD=False


SYM=True    # Calculate commutators symbolically?
SPARSE=False  # If not calculating commutators symbolically, use sparse matrices?

CSA=True #False

runs={}

picklFile='pickle.rick'
FORCE=True

if(os.path.exists(picklFile)==False or FORCE):

    for GLOBAL in 'isotropic','axial':

        # Create correct label (tag) for analysis
        print('Calculating NH')
        tag='NH'
        if(EXTH==True):
            tag+='H'
        if(EXTD==True):
            tag+='D'
            SYM=False

        inst=RelCalc(tag,sparse=SPARSE,sym=SYM)  #initialise class.



        ROT='static'        #local motion. rot or static.
        inst.HOP = False    # No hopping motion as no local motion
        inst.MACRO=False    #set macromolecular mode
        inst.GLOBAL =GLOBAL #set global motion

        # Define angles, can use strings or floats
        #beta='beta'  #interpret dipolar symbolically for externals

        # Include desired interactions to relaxation Hamiltonian
        inst.GetDip('H2','N1','d',beta='beta',alpha=0,rot=ROT) #NH dipolar
        if(CSA):
            inst.GetCSA('N1','cA',beta='betaC',alpha=0,rot=ROT) #N CSA
            inst.GetCSA('H2','cX',beta=0,alpha=0,rot=ROT) #H CSA

        extNum=3
        if(EXTH):
            inst.GetDip('H'+str(extNum), 'H2', 'e', beta=0,alpha=0,rot=ROT)  #HHext dipolar
            inst.GetDip('H'+str(extNum), 'N1', 'f', beta=0,alpha=0,rot=ROT)  #N Hext dipolar
            extNum+=1
        if(EXTD):
            inst.GetDip('D'+str(extNum), 'H2', 'g', beta=0,alpha=0,rot=ROT)  #HDext dipolar
            inst.GetDip('D'+str(extNum), 'N1', 'h', beta=0,alpha=0,rot=ROT)  #NDext dipolar

        inst.CrossCorrHam()  #getpairs of relaxation competant rates


        ####################
        #setup parameter dictionary
        R_NH=1.02E-10                    #NH bond length
        inst.pars['d']=R_NH,'N','H'     # N-H bond length in meters
        inst.pars['beta']=0
        #inst.pars['alpha']=0
        if(CSA):
            inst.pars['betaC']=0
            inst.pars['cX']=10.0,'H'        # CSA value for H in ppm
            inst.pars['cA']=-160.0,'N'       # CSA value for N in ppm

        if(GLOBAL=='axial'):
            # Define theta and phi angles to rotate from molecular frame to diffusion tensor
            inst.Mtheta = str('pi/4')
            #inst.Mpo = str('0pi/3')
            inst.Mpo = str('3*pi/3')
            # Define ratio of tau_parallel to tau_c (tau_perpendicular) 
            inst.tshape = 10    # (tshape>1 = oblate top, tshape<1 = prolate top)
        if(EXTH or EXTD):

            Rext=3E-10  #polar distance to external proton
            beta=111.6  #polar angle to external proton
            alpha=0     #azimuthal angle to external proton
            inst.pars['beta']=beta  #polar angle for external proton.


            if(EXTH):
                nucExt='H'
                inst.pars['e']=(Rext,'beta',alpha),(R_NH,0,0),nucExt,'H'  #X extH
                inst.pars['f']=(Rext,'beta',alpha),(0,0,0),nucExt,'N'     #X extH
            if(EXTD):
                nucExt='D'
                inst.pars['g']=(Rext,'beta',alpha),(R_NH,0,0),nucExt,'H'  #X extH
                inst.pars['h']=(Rext,'beta',alpha),(0,0,0),nucExt,'N'     #X extH
        #settings for the geometry plot:
        inst.PlotScale=1.5        #overall scaling of the figure (big streches it)
        inst.PlotEllipseZFac=0.6  #Zsize for tensor ellipsoids
        inst.PlotEllipseYFac=0.4  #Ysize for tensor ellipsoids
        inst.PlotEllipseXFac=0.2  #Xsize for tensor ellipsoids
        inst.PlotEllipsePhiNo=20  #number of phis   in ellipsoid plot
        inst.PlotEllipseThetaNo=7 #number of thetas in ellipsoid plot
        inst.PlotPS=1             #size of point that labels the atoms
        inst.MakeGeoFig('figNHgeo.pdf')  #show the geometry, computed using AX dipolars.
        ######################


        #first do relaxation rates relevant for fast dynamics analysis
        BASIS=[]
        BASIS='N1p','N1z','H2z'
        inst.EvalRate(BASIS,BASIS,verb='y') # Evaluate relaxation rates

        """
        inst.SetDegen(1)            #set number of degenerate protons.
        sym='p',
        BASIS=inst.GetBasisAXn(sym) # Define N+ basis (e.g. N+a and N+b)
        inst.EvalRate(BASIS,BASIS,verb='y') # Evaluate relaxation rates


        #now make complete spin half basis.
        SHIFTBASIS=inst.MakeShiftBasis(2)
        inst.EvalRate(SHIFTBASIS,SHIFTBASIS,verb='y') # Evaluate relaxation rates
        """

        inst.LatexReportTitle()

        # Plot 15N-H relaxation rates vs rotational correlation time at a spectrometer frequency of 600 MHz and 1200MHz
        inst.SetFreq(500)     #set proton Larmor of spectrometer.
        inst.SetTauC(4E-9)
        inst.SetPars()        #update parameters.
        #print(inst.CalcRate('N1p','N1p'))
        #sys.exit(100)
            

        inst.DoTcPlot(1E-6,1E-11,30,BASIS,'figNH_1.pdf')  #from 1E-6 to 1E-11 in tc
        #inst.SetFreq(1200)     #set proton Larmor of spectrometer.
        #inst.SetPars()         #update parameters.
        #inst.DoTcPlot(1E-6,1E-11,30,BASIS,'figNH_2.pdf')     #from 1E-6 to 1E-11 in tc


        #inst.SetTauC(1E-9)
        #inst.DoBasisPlot(SHIFTBASIS,'figNH_3.pdf')   #show redfield kite for BASIS for given pars
        #inst.SetTauC(1E-7)
        #inst.DoBasisPlot(SHIFTBASIS,'figNH_4.pdf')   #show redfield kite for BASIS for given pars


        #inst.MacroTable(AXNBASIS=BASIS,SHIFTBASIS=SHIFTBASIS,CARTBASIS=BOSIS)


        #inst.MacroTable(AXNBASIS=BASIS,SHIFTBASIS=SHIFTBASIS,CARTBASIS=BOSIS,NARROW=True)


        inst.CloseLatex()  #complete Latex file

        runs[GLOBAL]=copy.deepcopy(inst)

    print('Saving library file:',picklFile)
    with open(picklFile, 'wb') as handle:
        pickle.dump(runs, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    print('Loading library file:',picklFile)
    with open(picklFile, 'rb') as handle:
        runs = pickle.load(handle)
        



        
print("RelCalc libraries:",runs.keys())
#print(runs['isotropic'].resSet)




#exercise 1. Lets try to fit the R1/R2/NOE data to both models to see what parameters we get.
#print(runs['isotropic'].CalcRate('N1p','N1p'))
#sys.exit(10)



    

print()
print("Show the format for how we save relaxation rates in RelCalc:")
print("R1, isotropic:",runs['isotropic'].resSet['N1z']['N1z'])
print("R1, axial:    ",runs['axial'].resSet['N1z']['N1z'])

print("Show how to compute rates numerically using these objects:")
print("R1, isotropic  ",runs['isotropic'].CalcRate('N1z','N1z'),"s-1")
print("R1, axial:     ",runs['axial'].CalcRate('N1z','N1z'),"s-1")

print()

print("Now lets do a simple computation to figure out effective correlation times.")


# 25th Feb 2026
#Strong risk of this turning into spaghetti code!
#Class that reads in raw data, and takes some relcalc instances
#and then fits to various relaxation models.
#each model should probably be its own class rather than choking up the namespace as I have done here.
#there is a mixture of very local (per rate), local (per residue) and global analyses performed.
#a report is created summarising the gist of what's going on.
#currently I have not added errors in the fitting parameters and things like that.
#We're working with Art Palmers 500 Ubiquitin data, R1,R2,NOE. Simple models actually do a very good job.
#But its the R1s really that don't work, so these measurements are driving the move to more complex models.
class fittyRelax():
    def __init__(self,inst,pdbfile=''):

        self.Rel=inst
        self.ReadData()
        self.pdbfile=pdbfile
        
        #self.typs=typs

        self.LATEXFILE='summary.tex'

        #setup some output files.
        self.outfile='outy.out'                   #rates, isotropic fully local fit
        self.outfileGlob='outyGlob.out'           #rates, isotropic global fit
        self.outfileRes='outyRes.out'             #parameters, isotropic residue fit
        self.outfileResDat='outyResDat.out'       #rates, isotropic residue fit
        self.outfileResAx='outyResAx.out'         #parameters, axial local fit
        self.outfileResDatAx='outyResDatAx.out'   #rates, axial local fit
        self.outfileGlobAx='outyGlobAx.out'       #rates, axial globular fit
        self.outfileGlobAxRes='outGlobAxRes.out'  #residue angles, axial globular fit.
        
        outy=open(self.outfile,'w');outy.close()  #blank output file.
        outy=open(self.outfileRes,'w');outy.close()  #blank output file.
        outy=open(self.outfileResAx,'w');outy.close()  #blank output file.

        
        self.typs='R1','R2','NOE'  #for local analysis, use all 3 data types
        self.inst=self.Rel['isotropic']  #set current relcalc instance


        if(self.pdbfile!=''):

            Xang=9.011393158549712
            Zang=-3.8982578418659197
            
            self.pdb=ReadPDB(self.pdbfile,Xang=Xang,Zang=Zang)
        #print(self.pdb.angles)


        self.InitLatex()    #initialise the report
        self.MakePlotRaw()  #make plots of raw data.


        

        self.inst=self.Rel['isotropic']  #set current relcalc instance
        self.MakeCorrPlot()

        self.MakeJplot()
        
        self.AssembleData()          #setup global data arrays from raw data, merging the types.        

        
        self.DoFitLocal()            #do a local fit, each rate gets its own tc.
        self.DoFitGlobal()           #fit with one tauC
        self.DoFitResLocal()         #let each residue have its own tauC


        self.DoFitModelFree()         #let each residue have its own tauC

        
        self.inst=self.Rel['axial']  #set current relcalc instance
        
        self.R1R2=True
        self.DoAnisoFit()  #take angles from aligned PDB file and fit to get tc and tp


        self.fitGlobModelFreeAx()      #do fit.


        #self.R1R2=False
        #self.DoAnisoFit()  #take angles from aligned PDB file and fit to get tc and tp
        
        self.inst=self.Rel['isotropic']  #set current relcalc instance


        
        #self.inst=self.Rel['axial']      #switch to axial model
        #self.DoFitLocalAx()
        #self.DoFitGlobalAx()
        
        self.CloseLatex()            #conclude the report.
        sys.exit(100)


    def MakeCorrPlot(self):
        self.outfileCorr='outyCorr.out'
        tcs=numpy.logspace(-11,-5,100)
        outy=open(self.outfileCorr,'w')
        for tc in tcs:
            self.inst.SetTauC(tc)
            R1=self.inst.CalcRate('N1z','N1z')
            R2=self.inst.CalcRate('N1p','N1p')
            sig=self.inst.CalcRate('N1z','H2z')
            NOE=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
            outy.write('%e\t%e\t%e\t%e\t%e\n' % (tc,R1,R2,NOE,sig))
        outy.close()

        gnu=open('gnu.gp','w')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Isotropic model, local timescale fit\'\n')
        gnu.write('set logscale \n')
        gnu.write('set ylabel \'R1 or R2 (s-1)\'\n')
        gnu.write('set xlabel \'tc (ns)\'\n')

        gnu.write('set border 11\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')
        gnu.write('set y2tics\n')
        gnu.write('unset logscale y2\n')
        
        gnu.write('set title \'Isotropic model, local residue timescale fit\'\n')
        gnu.write('set y2label \' NOE\'\n')
        gnu.write('set format x "10^{%L}"\n')
        gnu.write('set format y "10^{%L}"\n')


        gnu.write('set output \'figCorr.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u ($1*1E9):2 ti \'R1\'w points pt 7 ps 0.2, ' % (self.outfileCorr))
        gnu.write('\'%s\' u ($1*1E9):3 ti \'R2\'w points pt 7 ps 0.2, ' % (self.outfileCorr))
        gnu.write('\'%s\' u ($1*1E9):4 ti \'NOE\'w points pt 7 ps 0.2 axes x1y2,' % (self.outfileCorr))
        gnu.write('\'%s\' u ($1*1E9):5 ti \'sig\'w points pt 7 ps 0.2 \n' % (self.outfileCorr))
        gnu.write('\n')

        gnu.write('set ylabel \'R1/R2\'\n')
        gnu.write('set border 3\n')
        gnu.write('unset y2label\n')
        gnu.write('set output \'fig2Corr.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u ($1*1E9):($2/$3) ti \'R1/R2\'w points pt 7 ps 0.2, ' % (self.outfileCorr))
        #gnu.write('\'%s\' u ($1*1E9):3 ti \'R2\'w points pt 7 ps 0.2, ' % (self.outfileCorr))
        #gnu.write('\'%s\' u ($1*1E9):4 ti \'NOE\'w points pt 7 ps 0.2 axes x1y2\n' % (self.outfileCorr))
        gnu.write('\n')

        
        gnu.close()
        
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Computation of rates verus correlation time, Isotropic model")
        
        self.WriteReportAddplot(('figCorr','fig2Corr'),cols=2)

        
    def MakeJplot(self):

        self.outfileJ='Jfile.out'
        oms=numpy.logspace(6,12,100)
        S2=0.1
        tc1=10E-9
        tc2=1E-10
        tce=(1/tc1+1/tc2)**(-1.)
        outy=open(self.outfileJ,'w')
        for om in oms:
            J1=S2*tc1/(1+om**2*tc1**2)
            J2=(1-S2)*tce/(1+om**2*tce**2)
            J=J1+J2
            outy.write('%e\t%e\t%e\t%e\n' % (om/1E6/(2*numpy.pi),J,J1,J2))
        outy.write('\n\n')

        oms=[]

        omH=500.*1E6*2*numpy.pi
        omN=500.*1E6*2*numpy.pi/self.inst.gammaH*self.inst.gammaN

        oms.append(1E-1*2*numpy.pi*1E6)
        oms.append(omH)
        oms.append(omN)
        oms.append(omH-omN)
        oms.append(omN+omH)

        lab=[]
        lab.append('0')
        lab.append('wH')
        lab.append('wN')
        lab.append('wH-wN')
        lab.append('wH+wN')

        for cnt,om in enumerate(oms):
            J1=S2*tc1/(1+om**2*tc1**2)
            J2=(1-S2)*tce/(1+om**2*tce**2)
            J=J1+J2
            outy.write('%e\t%e\t%e\t%e\t%s\n' % (om/(2*numpy.pi)/1E6,J,J1,J2,lab[cnt]))
        outy.close()

        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        #gnu.write('set key outside\n')

        gnu.write('set title \'Anisotropic model\'\n')
        gnu.write('set xlabel \'w (MHz)\'\n')
        gnu.write('set ylabel \'J(w)\'\n')

        gnu.write('set logscale x\n')
        gnu.write('set format x "10^{%L}"\n')

        gnu.write('set output \'figJ1.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' i 0 u 1:2 ti \'total\' w li,\'\' i 0 u 1:3 ti \'global\' w li,\'\' i 0 u 1:4 ti \'local\'w li' % (self.outfileJ))
        gnu.write(',\'%s\' i 1 u 1:2 noti w points pt 7 ,\'\' i 1 u 1:($2+1E-10):5 noti w labels' % (self.outfileJ))
        gnu.write('\n')

        """
        gnu.write('set output \'figJ2.pdf\'\n')
        gnu.write('unset logscale x\n')
        gnu.write('plot ')
        gnu.write('\'%s\' i 0 u 1:2 ti \'total\' w li,\'\' i 0 u 1:3 ti \'global\' w li,\'\' i 0 u 1:4 ti \'local\'w li' % (self.outfileJ))
        gnu.write(',\'%s\' i 1 u 1:2 noti w points pt 7 ,\'\' i 1 u 1:($2+1E-10):5 noti w labels' % (self.outfileJ))
        gnu.write('\n')
        """

        gnu.close()
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Spectral density")
        
        self.WriteReportAddplot(('figJ1.pdf',),cols=2)

        outy=open(self.LATEXFILE,'a')
        outy.write('Simulated model free spectral density function with $S^2$ %.2f, global tumbling $%.2f$ ns and local tumbling $%.2f$ ns. The sampled frequencies are indicated by points and labelled.\n\n' % (S2,tc1*1E9,tc2*1E9))
        #for typ in self.typs:
        #    outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        #outy.write("$\\tau_c$ overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcAll)*1E9,numpy.std(self.tcAll)*1E9))
        outy.close()

        
    #read datasets working out files from data type lists
    def ReadData(self):
        dat={}
        for source in 'R1','R2','NOE':
            dat[source]={}
            inny=open('raw/'+source+'.data.txt')
            for line in inny.readlines():
                test=line.split()
                if(len(test)<=1):
                    continue
                res=test[0]
                val=float(test[1])
                err=float(test[2])
                dat[source][res]=val,err
        self.dats=dat
        print("Data types:",self.dats.keys())


    #Turn raw data dictionary into arrays
    def AssembleData(self):
        #now do global.
        self.residues=[]
        self.yExp=[]
        self.yErr=[]
        for typ in self.typs:
            for key,vals in self.dats[typ].items():
                self.yExp.append(vals[0]) #get raw data
                self.yErr.append(vals[1])
                if(key not in self.residues):
                    self.residues.append(key)
        self.yExp=numpy.array(self.yExp)
        self.yErr=numpy.array(self.yErr)
        self.yCalc=numpy.zeros_like(self.yExp)


    ################################

    def CalcAnisoRatio(self):
        #set globals
        self.inst.SetTauC(self.tc)
        self.inst.tshape=self.tshape
        cnt=0
        for key,vals in self.dats['R1'].items():
            #self.inst.Mtheta=self.Mtheta[key] #set local.

            theta=self.pdb.angles[float(key)][1] #get theta. radians
            phi=self.pdb.angles[float(key)][2] #get theta. radians
            
            #update theta arrays, and spherical factors from anisotropic diffusion

            #do a stripped down setpars
            self.inst.Mtheta=self.Mtheta
            self.inst.Mpo=self.Mpo-phi
            
            self.inst.Jcache={} #clear cache
            self.inst.Dcache={} #clear cache
            #I think we don't need to do this.
            #self.inst.Mc=self.inst.Mnumba(self.inst.Mtheta)  #update spherical wigners
            self.inst.pars['beta']=theta/numpy.pi*180        #in degrees;
            #self.inst.pars['alpha']=phi/numpy.pi*180
            self.inst.pars['betaC']=self.betaC[key]/numpy.pi*180 #in degrees;
            
            #which angles do we need?
            #for beta,axes in list(self.inst.angs.items()):  #go over all angles RelCalc knows are needed...
            #    for ax,sSym in list(axes.items()):          #if there are multiple symmetries (here there are none)
            #        betaVal=self.inst.pars[beta] #try and get a number out of this
            #        self.inst.CalcConstsAng(betaVal,snam=sSym,rank=2) #rank2 only  #calculate legendre polynomials

            #if we have any 'L' coefficients for static axial interactions.           
            self.inst.Lconsts={}   #now compute the spectral density functions for axial symmetry (slow)
            for Lcoeff in self.inst.Lcoeffs: #loop over StAx spectral density functions and save.
                self.inst.Lconsts[Lcoeff]=self.inst.CalcConstsAxStatSmall(Lcoeff)

                    
            #self.inst.SetPars()         #update parameters.
            R1=self.inst.CalcRate('N1z','N1z')
            R2=self.inst.CalcRate('N1p','N1p')
            sig=self.inst.CalcRate('N1z','H2z')
            
            self.r1[cnt]=R1
            self.r2[cnt]=R2
            self.n[cnt]=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)            
            
            self.yc[cnt]=R1/R2


            
            self.yfc[cnt*3]=R1
            self.yfc[cnt*3+1]=R2
            self.yfc[cnt*3+2]=self.n[cnt]

            cnt+=1
    
    def packAniso(self):
        x=[]
        x.append(self.Mtheta)
        x.append(self.Mpo)
        x.append(self.tc)
        x.append(self.tshape)
        #for res in self.residues:
        #    x.append(self.betaC[res])
        return x
    def unpackAniso(self,x):
        cnt=0
        self.Mtheta=x[cnt];cnt+=1
        self.Mpo=x[cnt];cnt+=1
        self.tc=math.fabs(x[cnt]);cnt+=1
        self.tshape=math.fabs(x[cnt]);cnt+=1
        #for res in self.residues:
        #    self.betaC[res]=x[cnt];cnt+=1
    def chiAniso(self,x):
        self.unpackAniso(x)
        self.CalcAnisoRatio()
        if(self.R1R2):
            return (self.yc-self.y)/self.e
        else:
            return (self.yfc-self.yf)/self.ef
        
    def FitAniso(self):
        x0=leastsq(self.chiAniso,self.packAniso())
        print(x0[0])
        
    ###############################
        
    #
    def DoAnisoFit(self):
        print ("Doing anisotropic analysis")

        tag='aniso'
        if(self.R1R2):
            tag+='R1R2'
        self.outfileAniso='aniso.'+tag+'.out'
            
        self.x=[]
        self.y=[]
        self.e=[]
        self.yf=[]
        self.ef=[]
        self.betaC={}
        for res in self.residues: #setup data
            R1=self.dats['R1'][res][0]
            R1err=self.dats['R1'][res][1]
            R2=self.dats['R2'][res][0]
            R2err=self.dats['R2'][res][1]
            NOE=self.dats['NOE'][res][0]
            NOEerr=self.dats['NOE'][res][1]
            
            r=self.pdb.angles[float(res)][0]
            theta=self.pdb.angles[float(res)][1]
            phi=self.pdb.angles[float(res)][2]

            self.x.append( ((3*math.cos(theta)**2-1)/2) ) #p0
            self.y.append( (R1/R2) )
            self.e.append( ( R1/R2*math.sqrt( (R1err/R1)**2 + (R2err/R2)**2) ))

            self.yf.append(R1)
            self.yf.append(R2)
            self.yf.append(NOE)
            self.ef.append(R1err)
            self.ef.append(R2err)
            self.ef.append(NOEerr)

            self.betaC[res]=theta
            

        self.inst.SetPars()         #update parameters.
        
        self.tc=4E-9
        self.tshape=1.2
        self.Mtheta=0.0
        self.Mpo=0.0
        self.x=numpy.array(self.x)
        self.y=numpy.array(self.y)
        self.e=numpy.array(self.e)
        self.yc=numpy.zeros_like(self.y)

        self.yf=numpy.array(self.yf)
        self.ef=numpy.array(self.ef)
        self.yfc=numpy.zeros_like(self.yf)


        self.r1=numpy.zeros_like(self.y)
        self.r2=numpy.zeros_like(self.y)
        self.n=numpy.zeros_like(self.y)

        
        self.FitAniso()

        print('tshape:',self.inst.tshape)
        print('tc:    ',self.inst.tc,'ns')
        print('Mtheta:',self.inst.Mtheta/numpy.pi*180,'degrees')

        while(self.Mpo>2*numpy.pi):
            self.Mpo-=2*numpy.pi
        while(self.Mpo<0):
            self.Mpo+=2*numpy.pi
        print('Mp:    ',self.Mpo/numpy.pi*180,'degrees')

        self.CalcAnisoRatio()
        print('AveChi2:',numpy.sum(self.chiAniso(self.packAniso())**2)/len(self.packAniso()))

        print("done")

        
        #self.CalcAnisoRatio()




        outy=open(self.outfileAniso,'w')
        for cnt,res in enumerate(self.residues):
            #print(res)
            #for key,vals in self.dats[typ].items():
            outy.write('%s\t' % res)
            for typ in self.typs:
                #print(self.dats[typ][res])
                outy.write('%e\t%e\t' % (self.dats[typ][res][0],self.dats[typ][res][1]))
            #print(self.pdb.angles[int(res)])

            R1=self.dats['R1'][res][0]
            R1err=self.dats['R1'][res][1]
            R2=self.dats['R2'][res][0]
            R2err=self.dats['R2'][res][1]
            r=self.pdb.angles[float(res)][0]
            theta=self.pdb.angles[float(res)][1]
            phi=self.pdb.angles[float(res)][2]

            outy.write('%e\t%e\t%e\t' % (r,theta,phi)) #8,9,10
            outy.write('%e\t%e\t%e\t%e\t' % (self.x[cnt],self.y[cnt],self.e[cnt],self.yc[cnt])) #11,12,13,14
            outy.write('%e\t%e\t%e\t' % (self.r1[cnt],self.r2[cnt],self.n[cnt])) #15,16,17
            outy.write('\n')
        outy.close()

        
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        #gnu.write('set key outside\n')

        gnu.write('set title \'Anisotropic model\'\n')
        gnu.write('set xlabel \'P0(cosbeta)\'\n')
        gnu.write('set ylabel \'R1/R2\'\n')
        gnu.write('set cblabel \'residue\'\n')
        gnu.write('set output \'fig'+tag+'.2.pdf\'\n')
        if(self.R1R2):
            gnu.write('f2(x)=m*x+c\n')
            gnu.write('fit f2(x) \'%s\' u 11:14 via m,c\n' % (self.outfileAniso))
        gnu.write('plot ')
        gnu.write('\'%s\' u 11:12:13:1 noti w err lc palette pt 7 ps 0.2,\'\' u 11:14 ti \'fit\' w points' % (self.outfileAniso))
        if(self.R1R2):
            gnu.write(',f2(x) noti ')
        gnu.write('\n')

        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'fig'+tag+'.1.pdf\'\n')
        gnu.write('set key outside\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 2:15 ti \'R1\'w points pt 7 ps 0.2,' % (self.outfileAniso))
        gnu.write('\'%s\' u 4:16 ti \'R2\'w points pt 7 ps 0.2,' % (self.outfileAniso))
        gnu.write('\'%s\' u 6:17 ti \'R2\'w points pt 7 ps 0.2,x' % (self.outfileAniso))
        gnu.write('\n')

        gnu.close()
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Axial tumbling")
        
        self.WriteReportAddplot(('fig'+tag+'.1.pdf','fig'+tag+'.2.pdf',))

        outy=open(self.LATEXFILE,'a')
        #for typ in self.typs:
        #    outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        #outy.write("$\\tau_c$ overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcAll)*1E9,numpy.std(self.tcAll)*1E9))
        

        self.CalcAnisoRatio()
        if(self.R1R2):
            outy.write('Fitting R1R2 ratio\n\n')
        else:
            outy.write('Fitting rates directly\n\n')
        outy.write('\\noindent Parameters $P$: %i\n\n'% len(self.packAniso()))
        outy.write('tshape: %.2f\n\n' % self.inst.tshape)
        outy.write('tc:     %.2f ns\n\n' % (self.inst.tc*1E9))
        outy.write('$M_\\theta$: %.2f degrees\n\n' % (self.inst.Mtheta/numpy.pi*180))
        outy.write('$M_\\phi$:     %.2f degrees\n\n' % (self.Mpo/numpy.pi*180))
        

        P=len(self.packAniso())
        N=len(self.yf)
              
        chi=self.yfc-self.yf
        chi2=numpy.sum(  (chi/self.ef)**2)

        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))

        
        outy.close()

        
        #sys.exit(100)
        #now make weighted fit
        
        
        #plot 'test.out' u ((3*cos($9)**2-1)/2):($2/$4):(($2/$4)*sqrt(($3/$2)**2+($5/$4)**2)):1 w err lc palette,f2(x)
        
    ############################

        
    #do one tauC per rate fit.
    def DoFitResLocal(self):
        #first create per residue dictionaries and datasources.
        self.tcRes=[]  #will save residue tauC values.
        for self.res in self.residues: #for each residue...
            self.typRes=[]   #what types of data do we have?
            self.yExpRes=[]  #what experimental values?
            self.yErrRes=[]  #with errors.
            for typ in self.typs:  #go over the general set of types..
                if(self.res in self.dats[typ].keys()):  #if the current type matches a residues...
                    self.yExpRes.append(self.dats[typ][self.res][0])  #save data
                    self.yErrRes.append(self.dats[typ][self.res][1])  #save err
                    self.typRes.append(typ)                           #save type.
            self.yExpRes=numpy.array(self.yExpRes)   #turn into numpy....
            self.yErrRes=numpy.array(self.yExpRes)   #turn into numpy....
            self.yCalcRes=numpy.zeros_like(self.yExpRes)  #make a computed array in numpy....
            self.fitLocRes()  #and do the fit.
            self.tcRes.append(self.tc)  #save the result.
        self.MakePlotLocRes()  #make pretty plots and add to the report.
              

    #do one tauC per rate fit.
    def DoFitModelFree(self):
        #now do global.
        self.tc=4E-9
        self.fitGlobModelFree()      #do fit.

              

        

    #do one tauC fitted per relaxation rate.
    def DoFitLocal(self):
        self.tcVals={}   #save tauC values by data type.
        self.tcAll=[]    #save tauC values
        cnt=0
        for self.t,self.typ in enumerate(self.typs): #for each data type...
            if(self.typ=='NOE'):  #if NOE...
                self.tcVals[self.typ]=[]
                for self.key,vals in self.dats[self.typ].items():
                    self.yExpLoc=self.yExp[cnt]  #get local rate
                    self.yErrLoc=self.yErr[cnt]  #get local rate
                    self.fitLocNOE()             #fit to local NOE equations...
                    self.yCalc[cnt]=copy.deepcopy(self.yCalcLoc) #set calcuated value
                    self.tcVals[self.typ].append(copy.deepcopy(self.tc)) #store values.
                    self.tcAll.append(copy.deepcopy(self.tc))            #store values
                    cnt+=1
            else:
                if(self.typ=='R1'):  #if R1
                    self.rho1='N1z'
                    self.rho2='N1z'
                elif(self.typ=='R2'): #if R2
                    self.rho1='N1p'
                    self.rho2='N1p'
                self.tcVals[self.typ]=[]
                for self.key,vals in self.dats[self.typ].items():
                    self.yExpLoc=self.yExp[cnt]  #get local rate
                    self.yErrLoc=self.yErr[cnt]  #get local rate
                    self.fitLoc()                #fit to either R1 or R2 equations...
                    self.yCalc[cnt]=copy.deepcopy(self.yCalcLoc) #set calculated value
                    self.tcVals[self.typ].append(copy.deepcopy(self.tc)) #store values.
                    self.tcAll.append(copy.deepcopy(self.tc))            #store values
                    cnt+=1
            outy=open(self.outfile,'a');outy.write('\n\n');outy.close()  #split the file per data types to make a nice plot.
            print("tauC mean: ",self.typ,numpy.average(self.tcVals[self.typ]),numpy.std(self.tcVals[self.typ]))
        print("tauC overall: ",numpy.average(self.tcAll),numpy.std(self.tcAll))
        print(numpy.sum((((self.yCalc-self.yExp))/self.yErr)**2)/len(self.yCalc)) #print overall average chi2
        self.MakePlotLoc()  #make pretty plots.


    #do one tauc for all rates.
    def DoFitGlobal(self):
        #now do global.
        self.tc=numpy.average(self.tcAll) #guess initial condition from first
        self.fitGlob()      #do fit.

        
    #do some kind of axially symmetric analysis. Needs to be refined.
    def DoFitLocalAx(self):
        print("Doing residue level axial fit")
        #first create per residue dictionaries and datasources.
        self.tcRes=[]  #will save residue tauC values.
        self.tshapeRes=[]
        self.tTheta=[]
        for self.res in self.residues: #for each residue...
            self.typRes=[]   #what types of data do we have?
            self.yExpRes=[]  #what experimental values?
            self.yErrRes=[]  #with errors.
            for typ in self.typs:  #go over the general set of types..
                if(self.res in self.dats[typ].keys()):  #if the current type matches a residues...
                    self.yExpRes.append(self.dats[typ][self.res][0])  #save data
                    self.yErrRes.append(self.dats[typ][self.res][1])  #save err
                    self.typRes.append(typ)                           #save type.
            self.yExpRes=numpy.array(self.yExpRes)   #turn into numpy....
            self.yErrRes=numpy.array(self.yExpRes)   #turn into numpy....
            self.yCalcRes=numpy.zeros_like(self.yExpRes)  #make a computed array in numpy....
            self.fitLocResAx()  #and do the fit.
            self.tcRes.append(self.tc)  #save the result.
            self.tshapeRes.append(self.tshape)
            
            self.tTheta.append(self.Mtheta)

        self.MakePlotLocResAx()  #make pretty plots and add to the report.
        #self.tc=numpy.average(self.tcAll) #guess initial condition from first
        #self.fitLocalAx()


    #do one tauc for all rates.
    def DoFitGlobalAx(self):
        print("Doing global axial fit")
        #now do global.
        self.tc=numpy.average(self.tcAll) #guess initial condition from first
        self.tshape=1.1
        self.fitGlobAx()      #do fit.

        

    
    ################################
    #LOCAL, one tauC per rate.

    def pack(self):  #for 1 tc
        x=[]
        x.append(numpy.fabs(self.tc))
        return x
    def unpack(self,x): #for 1tc
        cnt=0
        self.tc=numpy.fabs(x[cnt]);cnt+=1

    def CalcRates(self):  #for R1/R2
        self.inst.SetTauC(self.tc)
        self.yCalcLoc=self.inst.CalcRate(self.rho1,self.rho2)
    def GuessParams(self):
        self.tc=0.1E-9
    def chi(self,x):   #get chi
        self.unpack(x)
        self.CalcRates()
        return self.yCalcLoc-self.yExpLoc
    def fitLoc(self):  #fit loc
        self.GuessParams()
        x0=leastsq(self.chi,self.pack())

        outy=open(self.outfile,'a')  #write output
        outy.write('%i\t%s\t%s\t%e\t%e\t%e\t%e\n' % (self.t,self.key,self.typ,self.yExpLoc,self.yCalcLoc,self.tc,(self.yExpLoc-self.yCalcLoc)/self.yErrLoc))
        outy.close()

        
        #outy=open(self.outfileSeq,'a')
        #outy.write('%i\t%s\t%s\t%e\t%e\n' % (self.t,self.key,self.typ,self.yExp,self.yCalc))
        #outy.close()

    ################################
    
    def CalcRatesNOE(self):
        self.inst.SetTauC(self.tc)
        R1=self.inst.CalcRate('N1z','N1z')
        sig=self.inst.CalcRate('N1z','H2z')
        self.yCalcLoc=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
    def chiLocNOE(self,x):
        self.unpack(x)
        self.CalcRatesNOE()
        return self.yCalcLoc-self.yExpLoc
    def fitLocNOE(self):
        self.GuessParams()
        x0=leastsq(self.chiLocNOE,self.pack())

        outy=open(self.outfile,'a')  #write output
        outy.write('%i\t%s\t%s\t%e\t%e\t%e\t%e\n' % (self.t,self.key,self.typ,self.yExpLoc,self.yCalcLoc,self.tc,(self.yExpLoc-self.yCalcLoc)/self.yErrLoc))
        outy.close()

        

    ################################
    #LOCAL, one tauc per residue
    def CalcRatesRes(self):
        self.inst.SetTauC(self.tc)
        cnt=0
        for typ in self.typRes:
            if(typ=='NOE'):
                R1=self.inst.CalcRate('N1z','N1z')
                sig=self.inst.CalcRate('N1z','H2z')
                noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
                #for key,vals in self.dats[typ].items():
                self.yCalcRes[cnt]=noe
                cnt+=1
            else:
                if(typ=='R1'):
                    self.rho1='N1z'
                    self.rho2='N1z'
                elif(typ=='R2'):
                    self.rho1='N1p'
                    self.rho2='N1p'
                rate=self.inst.CalcRate(self.rho1,self.rho2)
                #for key,vals in self.dats[typ].items():
                self.yCalcRes[cnt]=rate
                cnt+=1
    def chiRes(self,x):
        self.unpack(x)
        self.CalcRatesRes()
        return self.yCalcRes-self.yExpRes
    def fitLocRes(self):
        self.tc=4E-9
        x0=leastsq(self.chiRes,self.pack())

        outy=open(self.outfileRes,'a') #write output
        outy.write('%s\t%e\t%e\n' % (self.res,self.tc,numpy.sum( ((self.yCalcRes-self.yExpRes)/self.yErrRes)**2.)/len(self.yCalcRes)))
        outy.close()
        
        #pour back into global ararys
        cot=0
        for cot,typRes in enumerate(self.typRes): #populate the overall rate matrix - mapping.
            cnt=0
            for typ in self.typs: #for each datapoint, type
                for key,vals in self.dats[typ].items(): #for each datapoint, datset
                    cnt+=1   #increment datapoint index
                    if(self.res==key and typRes==typ): #but if we have a match, decant.
                        self.yCalc[cnt-1]=self.yCalcRes[cot] #pour residue values into the global list


        


    ####MAKE LOCAL FIT PLOTS AND REPORT######
    def MakePlotLoc(self):
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Isotropic model, local timescale fit\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'tc (ns)\'\n')
        gnu.write('set output \'fig2.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 2:($6*1E9) ti \'%s\'w points pt 7 ps 0.2 ' % (self.outfile,i,typ))
        gnu.write('\n')
            
        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'fig1.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 4:5 ti \'%s\'w points pt 7 ps 0.2' % (self.outfile,i,typ))
        gnu.write('\n')


        #make histograms and plot.
        self.bins=30
        self.histfile='hist.out'
        outy=open(self.histfile,'w')
        for i,typ in enumerate(self.typs):
            binmax=numpy.max(self.tcVals[typ])
            binmin=numpy.min(self.tcVals[typ])
            bins=numpy.linspace(binmin,binmax,self.bins)
            hist,edges=numpy.histogram(self.tcVals[typ],bins=bins)
            print(hist.shape,edges.shape)
            edges=(edges[1:]+edges[:-1])*0.5
            for j in range(len(edges)):
                outy.write('%e\t%e\n' % (edges[j],hist[j]))
            outy.write('\n\n')
        outy.close()

        gnu.write('set xlabel \'tc (ns)\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set output \'fig3.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u ($1*1E9):2 ti \'%s\'w boxes' % (self.histfile,i,typ))
        gnu.write('\n')

        gnu.close()
        
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Isotropic model, local $\\tau_c$ fit")
        
        self.WriteReportAddplot(('fig1.pdf','fig2.pdf','fig3.pdf'))


        
        outy=open(self.LATEXFILE,'a')


        P=len(self.tcAll)
        N=len(self.yExp)
              
        chi=self.yCalc-self.yExp
        chi2=numpy.sum(  (chi/self.yErr)**2)

        outy.write('\\noindent Parameters $P$: %i\n\n'% P)
        for typ in self.typs:
            outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        outy.write("$\\tau_c$ overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcAll)*1E9,numpy.std(self.tcAll)*1E9))

        
        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))



        #outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))


        #outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        #outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-len(self.tcAll)) ) )
        #outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))

        
        outy.close()

    ####MAKE RAW PLOTS AND REPORT######
    def MakePlotRaw(self):

        os.system('paste raw/R1.data.txt raw/R2.data.txt raw/NOE.data.txt > R1R2NOE.data.txt')

        gnu=open('gnu.gp','w')
        #gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        
        gnu.write('set title \'Raw data\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'rate (s-1)\'\n')
        gnu.write('set output \'figraw.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'raw/%s.data.txt\' u 1:2:3 ti \'%s\' w err pt 1' % (typ,typ))
        gnu.write(',\'R1R2NOE.data.txt\' u 1:(($8-1)*%f*$2) ti \'sig\'  pt 1' % (self.inst.gammaN/self.inst.gammaH))

        #sig = (NOE -1 )*gammaN/gammaH *R1
        
        gnu.write('\n')


        gnu.write('set output \'fig2raw.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'R1R2NOE.data.txt\' u 1:($2/$5) ti \'R1/R2\' w err pt 1' )
        gnu.write('\n')


        gnu.close()
        
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Raw data, ubiquitin, 500 MHz")
        
        self.WriteReportAddplot(('figraw','fig2raw'),cols=2)

        """
        outy=open(self.LATEXFILE,'a')
        #for typ in self.typs:
        #    outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        outy.write("$\\tau_c$ overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcAll)*1E9,numpy.std(self.tcAll)*1E9))
        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% len(self.residues))

        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-len(self.residues)) ) )
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))

        
        outy.close()
        """
        
    ####MAKE LOCAL RESIDUE PLOTS AND REPORT######     
    def MakePlotLocRes(self):
        #pour back together.
        outy=open(self.outfileResDat,'w')
        cnt=0
        for typ in self.typs:
            for key,vals in self.dats[typ].items():
                outy.write('%s\t%e\t%e\n' % (self.yExp[cnt],self.yCalc[cnt],self.yErr[cnt]))
                cnt+=1
            outy.write('\n\n')
        outy.close()


        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 11\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        
        gnu.write('set title \'Isotropic model, local residue timescale fit\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'tc (ns)\'\n')
        gnu.write('set y2label \' Chi2/N\'\n')
        gnu.write('set format y2 "10^{%L}"\n')
        gnu.write('set output \'fig2res.pdf\'\n')
        gnu.write('set logscale y2\n')
        gnu.write('set y2tics\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 1:($2*1E9) ti \'tc\'w points pt 7 ps 0.2,\'\' u 1:3 ti \'error\'w points axes x1y2 ' % (self.outfileRes))
        gnu.write('\n')




        gnu.write('set border 3\n')
        gnu.write('unset y2label\n')

        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'fig1res.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lw 0.1' % (self.outfileResDat,i,typ))
        gnu.write('\n')


        #make histograms and plot.
        self.binsRes=30
        self.histfileRes='histRes.out'
        outy=open(self.histfileRes,'w')

        binmax=numpy.max(self.tcRes)
        binmin=numpy.min(self.tcRes)
        bins=numpy.linspace(binmin,binmax,self.binsRes)
        hist,edges=numpy.histogram(self.tcRes,bins=bins)
        edges=(edges[1:]+edges[:-1])*0.5
        for j in range(len(edges)):
            outy.write('%e\t%e\n' % (edges[j],hist[j]))
        outy.close()

        gnu.write('unset key\n')
        gnu.write('set xlabel \'tc (ns)\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set output \'fig3res.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u ($1*1E9):2 w boxes' % (self.histfileRes))
        gnu.write('\n')

        gnu.close()
        
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Isotropic model, local residue $\\tau_c$ fit")
        
        self.WriteReportAddplot(('fig1res.pdf','fig2res.pdf','fig3res.pdf'))


        
        P=len(self.residues)
        N=len(self.yExp)
              
        chi=self.yCalc-self.yExp
        chi2=numpy.sum(  (chi/self.yErr)**2)

        outy=open(self.LATEXFILE,'a')
        outy.write('\\noindent Parameters $P$: %i\n\n'% P)
        outy.write("$\\tau_c$ residue overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcRes)*1E9,numpy.std(self.tcRes)*1E9))

        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))
        outy.close()

        """
        outy=open(self.LATEXFILE,'a')
        #for typ in self.typs:
        #    outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))

        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% len(self.residues))

        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-len(self.residues)) ) )
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))
        """
        



    ####MAKE LOCAL RESIDUE PLOTS AND REPORT######     
    def MakePlotLocResAx(self):
        #pour back together.
        outy=open(self.outfileResDatAx,'w')
        cnt=0
        for typ in self.typs:
            for key,vals in self.dats[typ].items():
                outy.write('%s\t%e\t%e\n' % (self.yExp[cnt],self.yCalc[cnt],self.yErr[cnt]))
                cnt+=1
            outy.write('\n\n')
        outy.close()

        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 11\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        
        gnu.write('set title \'Isotropic model, local residue timescale fit\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'tc (ns)\'\n')
        gnu.write('set y2label \' Chi2/N\'\n')
        gnu.write('set format y2 "10^{%L}"\n')
        gnu.write('set output \'fig2resAx.pdf\'\n')
        gnu.write('set logscale y2\n')
        gnu.write('set y2tics\n')
        gnu.write('plot ')
        #gnu.write('\'%s\' u 1:($2*1E9) ti \'tc\'w points pt 7 ps 0.2,\'\' u 1:3 ti \'error\'w points axes x1y2 ' % (self.outfileResAx))
        gnu.write('\'%s\' u 1:($2*1E9) ti \'tc\'w points pt 7 ps 0.2,\'\' u 1:3 ti \'tshape\'w points,\'\' u 1:4 ti \'theta\'w points  ' % (self.outfileResAx))
        gnu.write('\n')




        gnu.write('set border 3\n')
        gnu.write('unset y2label\n')

        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'fig1resAx.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lw 0.1' % (self.outfileResDatAx,i,typ))
        gnu.write('\n')


        #make histograms and plot.
        self.binsRes=30
        self.histfileResAx='histResAx.out'
        outy=open(self.histfileResAx,'w')

        binmax=numpy.max(self.tcRes)
        binmin=numpy.min(self.tcRes)
        bins=numpy.linspace(binmin,binmax,self.binsRes)
        hist,edges=numpy.histogram(self.tcRes,bins=bins)
        edges=(edges[1:]+edges[:-1])*0.5
        for j in range(len(edges)):
            outy.write('%e\t%e\n' % (edges[j],hist[j]))
        outy.close()

        gnu.write('unset key\n')
        gnu.write('set xlabel \'tc (ns)\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set output \'fig3resAx.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u ($1*1E9):2 w boxes' % (self.histfileResAx))
        gnu.write('\n')

        gnu.close()
        
        os.system('gnuplot gnu.gp')

        self.LATEXFILE='summary.tex'

        self.AddSection("Isotropic model, local residue $\\tau_c$ fit")
        
        self.WriteReportAddplot(('fig1resAx.pdf','fig2resAx.pdf','fig3resAx.pdf'))




        outy=open(self.LATEXFILE,'a')
        #for typ in self.typs:
        #    outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        outy.write("$\\tau_c$ overall: %.3f $\\pm$ %.3f ns\n\n" %(numpy.average(self.tcAll)*1E9,numpy.std(self.tcAll)*1E9))
        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% (len(self.residues)*3,))

        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-3*len(self.residues)) ) )
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))

        
        outy.close()

        


    ###################################
    #calc rates as needed. note we're not that efficient:
    #better would be to re-order with residue, and get R1/R2 and NOE at the same time. We're double calcing R1.
    #but this is fast enough. no big issue.
    def CalcRatesGlob(self):
        self.inst.SetTauC(self.tc)  #set one tauC down.
        cnt=0
        for typ in self.typs:  #for each data type...
            if(typ=='NOE'):
                R1=self.inst.CalcRate('N1z','N1z')
                sig=self.inst.CalcRate('N1z','H2z')
                noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=noe
                    cnt+=1
            else:
                if(typ=='R1'):
                    self.rho1='N1z'
                    self.rho2='N1z'
                elif(typ=='R2'):
                    self.rho1='N1p'
                    self.rho2='N1p'
                rate=self.inst.CalcRate(self.rho1,self.rho2)
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=rate
                    cnt+=1

    def chiGlob(self,x):
        self.unpack(x)
        self.CalcRatesGlob()
        return self.yCalc-self.yExp
    def fitGlob(self):
        #self.GuessParams()
        x0=leastsq(self.chiGlob,self.pack())

        print(x0[0])

        outy=open(self.outfileGlob,'w')
        cnt=0
        for t,typ in enumerate(self.typs):
            if(typ=='R1'):
                self.rho1='N1z'
                self.rho2='N1z'
            elif(typ=='R2'):
                self.rho1='N1p'
                self.rho2='N1p'
            for key,vals in self.dats[typ].items():
                outy.write('%e\t%e\n' % (self.yExp[cnt],self.yCalc[cnt]))
                cnt+=1
            outy.write('\n\n')
        outy.close()
        self.MakePlotGlob()



        
    ###MAKE OUTPUTS FOR GLOBAL#####
    def MakePlotGlob(self):
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Isotropic model, Global timescale fit\'\n')
        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'figGlob1.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lc -1 lw 0.05' % (self.outfileGlob,i,typ))
        gnu.write('\n')
        gnu.close()
        
        os.system('gnuplot gnu.gp')


        self.AddSection("Isotropic model, global $\\tau_c$ fit")
        self.WriteReportAddplot(('figGlob1.pdf',),cols=2)



        P=1
        N=len(self.yExp)
              
        chi=self.yCalc-self.yExp
        chi2=numpy.sum(  (chi/self.yErr)**2)

        outy=open(self.LATEXFILE,'a')
        outy.write('\\noindent Parameters $P$: %i\n\n'% P)
        outy.write("$\\tau_c$ global: %.3f ns\n\n" %(self.tc*1E9))
        
        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))




        """
        print(numpy.sum((((self.yCalc-self.yExp))/self.yErr)**2)/len(self.yCalc))
        
        print(numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp)))
        


        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% 1)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-1)  ))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))
        """
        outy.close()


    ###MAKE OUTPUTS FOR GLOBAL#####
    def MakePlotGlobFree(self):
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Model free\'\n')
        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'figGlobFree1.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lc -1 lw 0.05' % (self.outfileGlobFree,i,typ))
        gnu.write('\n')

        gnu.write('set border 11\n')
        gnu.write('set tics nomirror\n')
        
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'S^2\'\n')
        gnu.write('set y2label \' taue (ns)\'\n')
        gnu.write('set y2tics\n')
        gnu.write('set logscale y2\n')
        gnu.write('set format y2 "10^{%L}"\n')
        gnu.write('set output \'figGlobFree2.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 1:2 ti \'S2\'w points pt 7 ps 0.2 lc 1' % (self.outfileGlobFreePars))
        gnu.write(',\'%s\' u 1:($3*1E9) ti \'taue\'w points pt 7 ps 0.2 lc 2 axes x1y2' % (self.outfileGlobFreePars))
        gnu.write(',\'%s\' u 1:($4*1E9) ti \'tauc\' w li lc 2 axes x1y2 ' % (self.outfileGlobFreePars))
        gnu.write('\n')


        gnu.write('set border 3\n')
        gnu.write('unset y2tics\n')
        gnu.write('set output \'figGlobFree4.pdf\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'error (s-1)\'\n')
        gnu.write('unset y2label\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 1:(abs($5-$7)+abs($8-$10)+abs($11-$13)) noti w points pt 7 ps 0.2' % (self.outfileGlobFreePars))
        gnu.write('\n')



        
        gnu.write('set border 3\n')
        gnu.write('unset y2tics\n')
        gnu.write('set output \'figGlobFree3.pdf\'\n')
        gnu.write('set xlabel \'S^2\'\n')
        gnu.write('set ylabel \'taue ns\'\n')
        gnu.write('set logscale y\n')
        gnu.write('set format y "10^{%L}"\n')
        gnu.write('set cblabel \'residue\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 2:($3*1E9):1 noti w points pt 7 ps 0.2 lc palette' % (self.outfileGlobFreePars))
        gnu.write('\n')


        
        
        gnu.close()

        
        os.system('gnuplot gnu.gp')


        self.AddSection("Model Free fit")
        self.WriteReportAddplot(('figGlobFree1.pdf','figGlobFree2.pdf'),cols=2)
        self.WriteReportAddplot(('figGlobFree3.pdf','figGlobFree4.pdf'),cols=2)



        P=len(self.packFree())
        N=len(self.yExp)
              
        chi=self.yCalc-self.yExp
        chi2=numpy.sum(  (chi/self.yErr)**2)

        outy=open(self.LATEXFILE,'a')
        outy.write('\\noindent Parameters $P$: %i\n\n'% P)
        outy.write("$\\tau_c$ global: %.3f ns\n\n" %(self.tc*1E9))
        
        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))
        outy.close()




        """
        print(numpy.sum((((self.yCalc-self.yExp))/self.yErr)**2)/len(self.yCalc))
        
        print(numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp)))
        


        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% 1)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-1)  ))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))
        """



        
    ###################################
    #calc rates as needed. note we're not that efficient:
    #better would be to re-order with residue, and get R1/R2 and NOE at the same time. We're double calcing R1.
    #but this is fast enough. no big issue.
    def CalcRatesGlobFree(self):
        self.inst.SetTauC(self.tc)  #set one tauC down.
        cnt=0
        
        R1glob=self.inst.CalcRate('N1z','N1z')
        R2glob=self.inst.CalcRate('N1p','N1p')
        sigglob=self.inst.CalcRate('N1z','H2z')

        self.rates={}
        for key,vals in self.dats['R1'].items():
            tcorr=1./(1./self.tc+1./(self.te[key]))
            self.inst.SetTauC(tcorr)  #set one tauC down.
            
            R1loc=self.inst.CalcRate('N1z','N1z')
            R2loc=self.inst.CalcRate('N1p','N1p')
            sigloc=self.inst.CalcRate('N1z','H2z')

            R2=  R2glob*self.S2[key]  + (1-self.S2[key])*R2loc
            R1=  R1glob*self.S2[key]  + (1-self.S2[key])*R1loc
            sig= sigglob*self.S2[key] + (1-self.S2[key])*sigloc

            noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
            
            self.rates[key]=R1,R2,noe

        #unpack rates 
        cnt=0
        for typ in self.typs:  #for each data type...
            if(typ=='NOE'):
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][2]
                    cnt+=1
            elif(typ=='R1'):
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][0]
                    cnt+=1
            else:
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][1]
                    cnt+=1
                        
    def packFree(self):
        x=[]
        x.append(self.tc)
        for res in self.residues:
            x.append(math.acos(self.S2[res]))
            x.append(self.te[res])
        return x
    def unpackFree(self,x):
        cnt=0
        self.tc=math.fabs(x[cnt]);cnt+=1
        for res in self.residues:
            self.S2[res]=math.fabs(math.cos(x[cnt]));cnt+=1
            self.te[res]=math.fabs(x[cnt]);cnt+=1
        
    def chiGlobFree(self,x):
        self.unpackFree(x)
        self.CalcRatesGlobFree()
        return self.yCalc-self.yExp
    def fitGlobModelFree(self):
        print("Doing model free fit...")
        #self.GuessParams()
        self.tc=4E-9
        self.S2={}
        self.te={}
        for res in self.residues:
            self.S2[res]=0.5
            self.te[res]=0.4E-9
        
        x0=leastsq(self.chiGlobFree,self.packFree())

        print(x0[0])
        self.outfileGlobFree='outy.globModelFree.out'
        
        outy=open(self.outfileGlobFree,'w')
        cnt=0
        for t,typ in enumerate(self.typs):
            for key,vals in self.dats[typ].items():
                outy.write('%e\t%e\n' % (self.yExp[cnt],self.yCalc[cnt]))
                cnt+=1
            outy.write('\n\n')
        outy.close()

        self.outfileGlobFreePars='outy.globModelFreePars.out'
        
        outy=open(self.outfileGlobFreePars,'w')
        for res in self.residues:
            outy.write('%s\t' % res)
            outy.write('%e\t' % (self.S2[res]))
            outy.write('%e\t' % (self.te[res]))
            outy.write('%e\t' % (self.tc))


            cnt=0
            for t,typ in enumerate(self.typs):
                for key,vals in self.dats[typ].items():
                    if(key==res):
                        outy.write('%e\t%e\t%e\t' % (self.yExp[cnt],self.yErr[cnt],self.yCalc[cnt]))
                    cnt+=1

            theta=self.pdb.angles[float(res)][1] #get theta. radians
            phi=self.pdb.angles[float(res)][2] #get theta. radians
            outy.write('%e\t%e\t' % (theta,phi))
                    
            outy.write('\n')
        outy.close()
        
        self.MakePlotGlobFree()




    ###################################
    #calc rates as needed. note we're not that efficient:
    #better would be to re-order with residue, and get R1/R2 and NOE at the same time. We're double calcing R1.
    #but this is fast enough. no big issue.
    def CalcRatesGlobFreeAx(self):
        """
        self.inst.SetTauC(self.tc)  #set one tauC down.
        cnt=0
        
        R1glob=self.inst.CalcRate('N1z','N1z')
        R2glob=self.inst.CalcRate('N1p','N1p')
        sigglob=self.inst.CalcRate('N1z','H2z')

        self.rates={}
        for key,vals in self.dats['R1'].items():
            tcorr=1./(1./self.tc+1./(self.te[key]))
            self.inst.SetTauC(tcorr)  #set one tauC down.
            
            R1loc=self.inst.CalcRate('N1z','N1z')
            R2loc=self.inst.CalcRate('N1p','N1p')
            sigloc=self.inst.CalcRate('N1z','H2z')

            R2=  R2glob*self.S2[key]  + (1-self.S2[key])*R2loc
            R1=  R1glob*self.S2[key]  + (1-self.S2[key])*R1loc
            sig= sigglob*self.S2[key] + (1-self.S2[key])*sigloc

            noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
            
            self.rates[key]=R1,R2,noe
        """

        #set globals
        self.inst.SetTauC(self.tc)
        self.inst.tshape=self.tshape
        self.rates={}
        for key,vals in self.dats['R1'].items():
            #self.inst.Mtheta=self.Mtheta[key] #set local.

            theta=self.pdb.angles[float(key)][1] #get theta. radians
            phi=self.pdb.angles[float(key)][2] #get theta. radians
            
            #update theta arrays, and spherical factors from anisotropic diffusion

            #do a stripped down setpars
            self.inst.Mtheta=self.Mtheta
            self.inst.Mpo=self.Mpo-phi
            
            self.inst.Jcache={} #clear cache
            self.inst.Dcache={} #clear cache
            #I think we don't need to do this.
            #self.inst.Mc=self.inst.Mnumba(self.inst.Mtheta)  #update spherical wigners
            self.inst.pars['beta']=theta/numpy.pi*180        #in degrees;
            #self.inst.pars['alpha']=phi/numpy.pi*180
            self.inst.pars['betaC']=self.betaC[key]/numpy.pi*180 #in degrees;
            
            #which angles do we need?
            #for beta,axes in list(self.inst.angs.items()):  #go over all angles RelCalc knows are needed...
            #    for ax,sSym in list(axes.items()):          #if there are multiple symmetries (here there are none)
            #        betaVal=self.inst.pars[beta] #try and get a number out of this
            #        self.inst.CalcConstsAng(betaVal,snam=sSym,rank=2) #rank2 only  #calculate legendre polynomials

            #if we have any 'L' coefficients for static axial interactions.           
            self.inst.Lconsts={}   #now compute the spectral density functions for axial symmetry (slow)
            for Lcoeff in self.inst.Lcoeffs: #loop over StAx spectral density functions and save.
                self.inst.Lconsts[Lcoeff]=self.inst.CalcConstsAxStatSmall(Lcoeff)



            #self.inst.SetPars()         #update parameters.
            R1glob=self.inst.CalcRate('N1z','N1z')
            R2glob=self.inst.CalcRate('N1p','N1p')
            sigglob=self.inst.CalcRate('N1z','H2z')

            self.inst.MODELFREE=True  #turn it on
            self.inst.tauLoc=self.te[key]

            self.inst.Jcache={} #clear Jcache
            R1loc=self.inst.CalcRate('N1z','N1z')
            R2loc=self.inst.CalcRate('N1p','N1p')
            sigloc=self.inst.CalcRate('N1z','H2z')

            self.inst.MODELFREE=False #turn it off.

            R2=  R2glob*self.S2[key]  + (1-self.S2[key])*R2loc
            R1=  R1glob*self.S2[key]  + (1-self.S2[key])*R1loc
            sig= sigglob*self.S2[key] + (1-self.S2[key])*sigloc

            noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)

            self.rates[key]=R1,R2,noe

            
        #unpack rates 
        cnt=0
        for typ in self.typs:  #for each data type...
            if(typ=='NOE'):
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][2]
                    cnt+=1
            elif(typ=='R1'):
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][0]
                    cnt+=1
            else:
                for key,vals in self.dats[typ].items():
                    self.yCalc[cnt]=self.rates[key][1]
                    cnt+=1

                    
    def packFreeAx(self):
        x=[]
        x.append(self.tc)
        x.append(self.tshape)
        #x.append(self.Mtheta)
        #x.append(self.Mpo)
        for res in self.residues:
            x.append(math.acos(self.S2[res]))
            x.append(self.te[res])
        return x
    def unpackFreeAx(self,x):
        cnt=0
        self.tc=math.fabs(x[cnt]);cnt+=1
        self.tshape=math.fabs(x[cnt]);cnt+=1
        #self.Mtheta=math.fabs(x[cnt]);cnt+=1
        #self.Mpo=math.fabs(x[cnt]);cnt+=1
        
        for res in self.residues:
            self.S2[res]=math.fabs(math.cos(x[cnt]));cnt+=1
            self.te[res]=math.fabs(x[cnt]);cnt+=1
        
    def chiGlobFreeAx(self,x):
        self.unpackFreeAx(x)
        self.CalcRatesGlobFreeAx()
        return self.yCalc-self.yExp
    def fitGlobModelFreeAx(self):
        print("Doing model free fit axial...")

        print('tauc:   ',self.tc)
        print('tshape: ',self.tshape)
        print('Mtheta: ',self.Mtheta)
        print('Mpo:    ',self.Mpo)
        #self.GuessParams()
        #self.tc=4E-9
        #self.tshape=1.18
        #self.Mtheta=0.0
        #self.Mpo=0.0

        #should e already initialised
        #self.S2={}
        #self.te={}
        #self.betaC={}
        #for res in self.residues:
        #    self.S2[res]=0.5
        #    self.te[res]=0.4E-9

            #self.betaC[res]=theta

            
        x0=leastsq(self.chiGlobFreeAx,self.packFreeAx())

        print(x0[0])
        self.outfileGlobFreeAx='outy.globModelFreeAx.out'
        
        outy=open(self.outfileGlobFreeAx,'w')
        cnt=0
        for t,typ in enumerate(self.typs):
            for key,vals in self.dats[typ].items():
                outy.write('%e\t%e\n' % (self.yExp[cnt],self.yCalc[cnt]))
                cnt+=1
            outy.write('\n\n')
        outy.close()

        self.outfileGlobFreeParsAx='outy.globModelFreeParsAx.out'
        
        outy=open(self.outfileGlobFreeParsAx,'w')
        for res in self.residues:
            outy.write('%s\t' % res)
            outy.write('%e\t' % (self.S2[res]))
            outy.write('%e\t' % (self.te[res]))
            outy.write('%e\t' % (self.tc))
            outy.write('%e\t' % (self.tshape))


            cnt=0
            for t,typ in enumerate(self.typs):
                for key,vals in self.dats[typ].items():
                    if(key==res):
                        outy.write('%e\t%e\t%e\t' % (self.yExp[cnt],self.yErr[cnt],self.yCalc[cnt]))
                    cnt+=1

            theta=self.pdb.angles[float(res)][1] #get theta. radians
            phi=self.pdb.angles[float(res)][2] #get theta. radians
            outy.write('%e\t%e\t' % (theta,phi))
                    
            outy.write('\n')
        outy.close()
        
        self.MakePlotGlobFreeAx()

    ###MAKE OUTPUTS FOR GLOBAL#####
    def MakePlotGlobFreeAx(self):
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Model free\'\n')
        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'figGlobFreeAx1.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lc -1 lw 0.05' % (self.outfileGlobFreeAx,i,typ))
        gnu.write('\n')

        gnu.write('set border 11\n')
        gnu.write('set tics nomirror\n')
        
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'S^2\'\n')
        gnu.write('set y2label \' taue\'\n')
        gnu.write('set y2tics\n')
        gnu.write('set logscale y2\n')
        gnu.write('set format y2 "10^{%L}"\n')
        gnu.write('set output \'figGlobFreeAx2.pdf\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 1:2 ti \'S2\'w points pt 7 ps 0.2 lc 1' % (self.outfileGlobFreeParsAx))
        gnu.write(',\'%s\' u 1:3 ti \'taue\'w points pt 7 ps 0.2 lc 2 axes x1y2' % (self.outfileGlobFreeParsAx))
        gnu.write(',\'%s\' u 1:4 ti \'tauc\' w li lc 2 axes x1y2 ' % (self.outfileGlobFreeParsAx))
        gnu.write('\n')


        gnu.write('set border 3\n')
        gnu.write('unset y2tics\n')
        gnu.write('set output \'figGlobFreeAx4.pdf\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'error (s-1)\'\n')
        gnu.write('unset y2label\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 1:(abs($6-$8)+abs($9-$11)+abs($12-$14)) noti w points pt 7 ps 0.2' % (self.outfileGlobFreeParsAx))
        gnu.write('\n')


        
        gnu.write('set border 3\n')
        gnu.write('unset y2tics\n')
        gnu.write('set output \'figGlobFreeAx3.pdf\'\n')
        gnu.write('set xlabel \'S^2\'\n')
        gnu.write('set ylabel \'taue ns\'\n')
        gnu.write('set logscale y\n')
        gnu.write('set format y "10^{%L}"\n')
        gnu.write('set cblabel \'residue\'\n')
        gnu.write('plot ')
        gnu.write('\'%s\' u 2:($3*1E9):1 noti w points pt 7 ps 0.2 lc palette' % (self.outfileGlobFreeParsAx))
        gnu.write('\n')

        
        gnu.close()

        
        os.system('gnuplot gnu.gp')


        self.AddSection("Model Free Axial fit")
        self.WriteReportAddplot(('figGlobFreeAx1.pdf','figGlobFreeAx2.pdf',),cols=2)
        self.WriteReportAddplot(('figGlobFreeAx3.pdf','figGlobFreeAx4.pdf'),cols=2)


        P=len(self.packFreeAx())
        N=len(self.yExp)
              
        chi=self.yCalc-self.yExp
        chi2=numpy.sum(  (chi/self.yErr)**2)

        outy=open(self.LATEXFILE,'a')
        outy.write('\\noindent Parameters $P$: %i\n\n'% P)
        outy.write("$\\tau_c$ global: %.3f ns\n\n" %(self.tc*1E9))
        outy.write("shape global: %.3f\n\n" %(self.tshape))
        
        outy.write('\\noindent Datapoints $N$: %i\n\n'% N)
        outy.write('$\\chi^2/N$: %.3f\n\n' % (chi2/N))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (chi2/(N-P )))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sqrt(numpy.sum(chi**2.)/N)))
        outy.close()


        

    ###MAKE OUTPUTS FOR GLOBAL#####
    def MakePlotGlobAx(self):

        outy=open(self.outfileGlobAxRes,'w')
        for i in range(len(self.residues)):
            ang=self.Mtheta[self.residues[i]]
            while(ang>numpy.pi):
                ang-=numpy.pi
            while(ang<0):
                ang+=numpy.pi
            outy.write('%e\t%e\n' % (self.residues[i],ang))
        outy.close()
            
        gnu=open('gnu.gp','w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key outside\n')

        gnu.write('set title \'Isotropic model, Global timescale fit\'\n')
        gnu.write('set xlabel \'Rate(exp)\'\n')
        gnu.write('set ylabel \'Rate(calc)\'\n')
        gnu.write('set output \'figGlob1Ax.pdf\'\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w points pt 7 ps 0.2,x noti lc -1 lw 0.05' % (self.outfileGlobAx,i,typ))
        gnu.write('\n')
        gnu.close()
        
        os.system('gnuplot gnu.gp')


        

        self.AddSection("Axial model, global $\\tau_c$ fit")
        self.WriteReportAddplot(('figGlob1Ax.pdf',),cols=2)

        print(numpy.sum((((self.yCalc-self.yExp))/self.yErr)**2)/len(self.yCalc))
        
        print(numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp)))
        
        outy=open(self.LATEXFILE,'a')
        outy.write("$\\tau_c$ global: %.3f ns\n\n" %(self.tc*1E9))
        outy.write("$\\tau_\\mathrm{shape}$ global: %.3f \n\n" %(self.tshape))
        #print out angles.
        outy.write('Datapoints $N$: %i\n\n'% len(self.yExp))
        outy.write('Parameters $P$: %i\n\n'% len(self.packGlobAx()))
        outy.write('$\\chi^2/N$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/len(self.yExp)  ))
        outy.write('$\\chi^2/(N-P)$: %.3f\n\n' % (numpy.sum(((self.yCalc-self.yExp)/self.yErr)**2.)/(len(self.yExp)-len(self.packGlobAx()))  ))
        outy.write('RMSD:   %.4f s$^{-1}$\n\n' % (numpy.sum(((self.yCalc-self.yExp))**2.)/(len(self.yExp))))
        outy.close()

        
    ##########################################
    #this is where we are with axial:

    def CalcRatesGlobAx(self):
        #set globals
        self.inst.SetTauC(self.tc)
        self.inst.tshape=self.tshape
        cnt=0
        for typ in self.typs:
            if(typ=='NOE'):
                for key,vals in self.dats[typ].items():
                    self.inst.Mtheta=self.Mtheta[key] #set local.
                    self.inst.SetPars()         #update parameters.
                    R1=self.inst.CalcRate("N1z","N1z")
                    sig=self.inst.CalcRate("N1z","H2z")
                    self.yCalc[cnt]=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
                    cnt+=1
                
            else:
                if(typ=='R1'):
                    rho1='N1z'
                    rho2='N1z'
                if(typ=='R2'):
                    rho1='N1p'
                    rho2='N1p'
                for key,vals in self.dats[typ].items():
                    self.inst.Mtheta=self.Mtheta[key] #set local.
                    self.inst.SetPars()         #update parameters.
                    self.yCalc[cnt]=self.inst.CalcRate(rho1,rho2)
                    cnt+=1
        
    def GuessParamsGlobAx(self):
        self.tc=1E-9
        self.tshape=1.1
        self.Mtheta={}
        for res in self.residues:
            self.Mtheta[res]=numpy.random.rand()*2*numpy.pi

        
    def packGlobAx(self):
        x=[]
        x.append(self.tc)
        x.append(self.tshape)
        for res in self.residues:
            x.append(self.Mtheta[res])
        return x
    
    def unpackGlobAx(self,x):
        cnt=0
        self.tc=x[cnt];cnt+=1
        self.tshape=x[cnt];cnt+=1
        for res in self.residues:
            self.Mtheta[res]=x[cnt];cnt+=1
            
    def chiGlobAx(self,x):
        self.unpackGlobAx(x)
        self.CalcRatesGlobAx()
        print(numpy.sum((self.yCalc-self.yExp)**2)/len(self.yCalc))
        return self.yCalc-self.yExp
    def fitGlobAx(self):

        #for global axially symmetric model.
        # we need tc and shape factor.
        # lets put Mtheta to zero.
        # then for each residue, they have a theta and a phi to move around.
        # we don't have enough data to do this.

            
        self.GuessParamsGlobAx()
        x0=leastsq(self.chiGlobAx,self.packGlobAx())


        print(x0[0])
        print("Average chi2:",numpy.sum((self.chi(x0[0])/self.yErr)**2)/len(self.yExp))

        outy=open(self.outfileGlobAx,'w')
        for i in range(len(self.yExp)):
            outy.write('%e\t%e\n' % (self.yExp[i],self.yCalc[i]))
        outy.close()
        
        #print(x0[0])

    #####################


    ################################
    #LOCAL, one tauc per residue
    def CalcRatesResAx(self):
        self.inst.SetTauC(self.tc)
        self.inst.tshape=self.tshape
        self.inst.Mtheta=self.Mtheta
        self.inst.SetPars()         #update parameters.

        cnt=0
        for typ in self.typRes:
            if(typ=='NOE'):
                R1=self.inst.CalcRate('N1z','N1z')
                sig=self.inst.CalcRate('N1z','H2z')
                noe=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
                #for key,vals in self.dats[typ].items():
                self.yCalcRes[cnt]=noe
                cnt+=1
            else:
                if(typ=='R1'):
                    self.rho1='N1z'
                    self.rho2='N1z'
                elif(typ=='R2'):
                    self.rho1='N1p'
                    self.rho2='N1p'
                rate=self.inst.CalcRate(self.rho1,self.rho2)
                #for key,vals in self.dats[typ].items():
                self.yCalcRes[cnt]=rate
                cnt+=1

    def packAx(self):
        x=[]
        x.append(self.tc)
        x.append(self.tshape)
        x.append(self.Mtheta)
        return x
    def unpackAx(self,x):
        cnt=0
        self.tc=x[cnt];cnt+=1
        self.tshape=x[cnt];cnt+=1
        self.Mtheta=x[cnt];cnt+=1
                
    def chiResAx(self,x):
        self.unpackAx(x)
        self.CalcRatesResAx()
        return self.yCalcRes-self.yExpRes
    def fitLocResAx(self):
        self.tc=4E-9
        self.tshape=1.1
        self.Mtheta=0.2
        x0=leastsq(self.chiResAx,self.packAx())

        #fold Mtheta back.
        while(self.Mtheta>2*numpy.pi):
            self.Mtheta-=2*numpy.pi
        while(self.Mtheta<0):
            self.Mtheta+=2*numpy.pi

        
        outy=open(self.outfileResAx,'a') #write output
        outy.write('%s\t%e\t%e\t%e\t%e\n' % (self.res,self.tc,self.tshape,self.Mtheta,numpy.sum( ((self.yCalcRes-self.yExpRes)/self.yErrRes)**2.)/len(self.yCalcRes)))
        outy.close()
        
        #pour back into global ararys
        cot=0
        for cot,typRes in enumerate(self.typRes): #populate the overall rate matrix - mapping.
            cnt=0
            for typ in self.typs: #for each datapoint, type
                for key,vals in self.dats[typ].items(): #for each datapoint, datset
                    cnt+=1   #increment datapoint index
                    if(self.res==key and typRes==typ): #but if we have a match, decant.
                        self.yCalc[cnt-1]=self.yCalcRes[cot] #pour residue values into the global list



    



        

    ##############################################
    #HELPER FUNCTIONS FOR REPORT GENERATION
    def AddSection(self,sec):
        outy=open(self.LATEXFILE,'a')
        outy.write('\\section{%s}\n\n' % sec)
        outy.close()

    
    ### Initialise a latex report
    def InitLatex(self):
        if(os.path.exists(self.LATEXFILE)):
            os.system('rm '+self.LATEXFILE)
        outy=open(self.LATEXFILE,'w')
        #outy.write('\\documentclass[showkeys,aps,prb,prepreint,amssymb, amsmath,nobibnotes]{revtex4}\n')
        outy.write('\\documentclass[showkeys,aps,prb,prepreint,amssymb, amsmath,nobibnotes]{article}\n')
        outy.write('\\usepackage{bm,setspace}\n')
        outy.write('\\usepackage{graphicx,graphics,booktabs}\n')
        outy.write('\\usepackage[a4paper,margin=2cm]{geometry}\n')
        #outy.write('\\usepackage{accents}\n')
        outy.write('\\newlength{\\dhatheight}\n')
        #outy.write('\\newcommand{\doublehat}[1]{%\settoheight{\dhatheight}{\ensuremath{\hat{#1}}}%\addtolength{\dhatheight}{-0.35ex}%\hat{\vphantom{\rule{1pt}{\dhatheight}}%\smash{\hat{#1}}}}\n')
        outy.write('\\begin{document}\n')


        outy.close()

    ### called from dotcplot
    def WriteReportAddplot(self,outfigs,cols=3):
        outy=open(self.LATEXFILE,'a')

        if(cols==''):
            width=1./len(outfigs)*0.99
        else:
            width=1./cols*0.99
        
        #outy.write('\\center\n')
        for i in range(len(outfigs)):
            outy.write('\\includegraphics[width=%f\\textwidth]{%s}\n' % (width, outfigs[i]))
        outy.write('\\end{figure} \n')

        
    ### Complete a latex report.
    def CloseLatex(self,tig=''):
        print("Compiling latex file")
        outy=open(self.LATEXFILE,'a')
        outy.write('\\end{document}\n')
        outy.close()
        #os.system('pdflatex '+self.LATEXFILE)
        os.system('pdflatex -interaction=nonstopmode '+self.LATEXFILE+' > latex.log')
                

        os.system('rm latex.aux latex.log')

    
            
        
"""
class fittyAnIso():
    def __init__(self,inst,dats,typs,FitTheta=False):
        self.dats=dats
        self.typs=typs
        self.inst=inst
        self.FitTheta=FitTheta
        self.residues=[]
        self.yExp=[]
        self.yErr=[]
        for typ in self.typs:
            for key,vals in self.dats[typ].items():
                self.yExp.append(vals[0]) #get raw data
                self.yErr.append(vals[1])
                if(key not in self.residues):
                    self.residues.append(key)
        self.yExp=numpy.array(self.yExp)
        self.yErr=numpy.array(self.yErr)
        self.yCalc=numpy.zeros_like(self.yExp)
        self.fit()

    def CalcRates(self):
        #set globals
        self.inst.SetTauC(self.tc)
        self.inst.tshape=self.tshape
        cnt=0
        for typ in self.typs:
            if(typ=='NOE'):
                for key,vals in self.dats[typ].items():
                    if(self.FitTheta):
                        self.inst.Mtheta=self.Mtheta[key] #set local.
                        self.inst.SetPars()         #update parameters.
                    R1=self.inst.CalcRate("N1z","N1z")
                    sig=self.inst.CalcRate("N1z","H1z")
                    self.yCalc[cnt]=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
                    cnt+=1
                
            else:
                if(typ=='R1'):
                    rho1='N1z'
                    rho2='N1z'
                if(typ=='R2'):
                    rho1='N1p'
                    rho2='N1p'
                for key,vals in self.dats[typ].items():
                    if(self.FitTheta):
                        self.inst.Mtheta=self.Mtheta[key] #set local.
                        self.inst.SetPars()         #update parameters.
                    self.yCalc[cnt]=self.inst.CalcRate(rho1,rho2)
                    cnt+=1
        
    def GuessParams(self):
        self.tc=1E-9
        self.tshape=1.1
        if(self.FitTheta):
            self.Mtheta={}
            for res in self.residues:
                self.Mtheta[res]=numpy.random.rand()*2*numpy.pi

        
    def pack(self):
        x=[]
        x.append(self.tc)
        x.append(self.tshape)
        if(self.FitTheta):
            for res in self.residues:
                x.append(self.Mtheta[res])
        return x
    def unpack(self,x):
        cnt=0

        self.tc=x[cnt];cnt+=1
        self.tshape=x[cnt];cnt+=1
        if(self.FitTheta):
            for res in self.residues:
                self.Mtheta[res]=x[cnt];cnt+=1
            
    def chi(self,x):
        self.unpack(x)
        self.CalcRates()
        return self.yCalc-self.yExp
    def fit(self):
        self.GuessParams()
        x0=leastsq(self.chi,self.pack())


        print(x0[0])
        print("Average chi2:",numpy.sum((self.chi(x0[0])/self.yErr)**2)/len(self.yExp))

        outy=open('outy.out','w')
        for i in range(len(self.yExp)):
            outy.write('%e\t%e\n' % (self.yExp[i],self.yCalc[i]))
        outy.close()
        
        #print(x0[0])
"""

class ReadPDB():
    def __init__(self,infile,Xang=0,Zang=0):
        #add protons
        #obabel -ipdb 1UBQ.pdb -opdb -O output_h.pdb -h
        #centre of inertia reference frame
        #./pdb_inertia_write.py 1UBQ_h.pdb 1UBQ_align.pdb
        self.Xang=Xang/180.*numpy.pi
        self.Zang=Zang/180.*numpy.pi
        self.infile=infile
        self.extract_pdb()     #get N and H atoms from PDB file
        self.GetOrientation()  #calculate angles of NHs to Z axis
        

    def GetOrientation(self):
        print("Getting angles...")

        #Do an XZ rotation.
        print(self.Xang,self.Zang)
        RX=self.rodrigues_Mat([1,0,0],self.Xang)
        RZ=self.rodrigues_Mat([0,0,1],self.Zang)


        R=RZ@RX


        self.angles={}
        for res,atoms in self.atoms.items():
            if(len(atoms)!=2): #make sure we have two entries
                print('skipping',res,atoms[list(atoms.keys())[0]]['residue_name'])
                continue
            if('N' not in atoms):
                print('shit, no N')
                print(atoms)
            if('H' not in atoms):
                print('shit, no H')
                print(atoms)
            nx=atoms['N']['x']
            ny=atoms['N']['y']
            nz=atoms['N']['z']
            ex=atoms['H']['x']
            ey=atoms['H']['y']
            ez=atoms['H']['z']
            
            xx=nx-ex
            yy=ny-ey
            zz=nz-ez

            #xx,yy,zz=R@[xx,yy,zz ]
            [xx,yy,zz]=numpy.dot(R,[xx,yy,zz])
            
            r,t,p=self.xyz_to_spherical(xx,yy,zz)
            #print(r,t,p)
            #for atom,info in atoms.items():
            #    if(info['residue_name']=='PRO'):
            #        continue
            self.angles[res]=r,t,p
                
                #print(res,atom,info)
                #print(len(atoms))
        print("Residues with an NH orientation: ",len(self.angles))

    def rodrigues_Mat(self,axis, theta):
        """
        Rotates a point using Rodrigues' rotation formula.
        point: np.array([x, y, z])
        axis: np.array([x, y, z]) (rotation axis)
        theta: float (angle in radians)
        """
        # 1. Normalize the rotation axis
        axis = axis / numpy.linalg.norm(axis)
        
        # 2. Skew-symmetric cross-product matrix K of unit axis
        K = numpy.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
    ])
        
        # 3. Rodrigues' formula: R = I + sin(theta)K + (1-cos(theta))K^2
        I = numpy.eye(3)
        R = I + numpy.sin(theta) * K + (1 - numpy.cos(theta)) * numpy.dot(K, K)
        
        return R
        # 4. Apply rotation
        #return np.dot(R, point)
    
    # Example Usage:
    #point = np.array([1.0, 0.0, 0.0]) # Point on X-axis
    #axis = np.array([0.0, 0.0, 1.0])  # Rotate around Z-axis
    #theta = np.radians(90)            # 90 degrees
   
    #rotated_point = rodrigues_rotate(point, axis, theta)
        
    def xyz_to_spherical(self,x, y, z):
       r = math.sqrt(x**2 + y**2 + z**2)

       if r == 0:
           theta = 0.0
           phi = 0.0
       else:
           theta = math.acos(z / r)     # polar angle
           phi = math.atan2(y, x)       # azimuthal angle

       return r, theta, phi

                
        
    def extract_pdb(self):
        print("Reading ",self.infile)
        self.atoms={}       
        with open(self.infile, 'r') as file:
            for line in file:
                # ATOM and HETATM lines contain atomic coordinates
                if line.startswith(("ATOM",)):
                    #print(line)
                    try:
                        atom_name = line[12:16].strip()      # Atom name
                        if(atom_name=='N' or atom_name=='H'):

                            residue_number = int(line[22:26].strip())  # Residue sequence number
                            residue_name = line[17:20].strip() 
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())

                            if(residue_number not in self.atoms.keys()):
                                self.atoms[residue_number]={}

                            entry={
                                    "atom_name": atom_name,
                                    "residue_number": residue_number,
                                    "residue_name": residue_name,
                                    "x": x,
                                    "y": y,
                                    "z": z
                                    }

                            #print(entry)
                            if(atom_name=='N'):
                                self.atoms[residue_number][atom_name]=entry
                            elif(atom_name=='H'):
                                nx=self.atoms[residue_number]['N']['x']
                                ny=self.atoms[residue_number]['N']['y']
                                nz=self.atoms[residue_number]['N']['z']

                                ex=entry['x']
                                ey=entry['y']
                                ez=entry['z']
                                rr=((nx-ex)**2+(ny-ey)**2+(nz-ez)**2)**0.5
                                if(rr<1.5):
                                    entry['r']=rr
                                    self.atoms[residue_number][atom_name]=entry
                                #if(residue_number==52):
                                #    print()
                                #   print(self.atoms[residue_number]['N'])
                                #    print(entry)
                                #    print(rr)
                                    
                    except ValueError:
                        continue  # Skip malformed lines
        print("Number of residues found:",len(self.atoms.keys()))
        #sys.exit(10)

        

fit=fittyRelax(runs,pdbfile='1UBQ_h_align.pdb') #merge relcalc results with a module for doing various types of fitting.
#fit=fittyRelax(runs,pdbfile='1UBQ_h.pdb') #merge relcalc results with a module for doing various types of fitting.




