#!/usr/bin/env python3

import os,sys,pandas as pd,numpy as np

from Bio.SeqUtils.ProtParam import ProteinAnalysis

#take an assignment row with chemical shifts
#parse it and save class variables
#so this is a single atom assignment entry
class ass():
    def __init__(self,test):
        self.resn=test[4]           #residue name
        self.resi=float(test[5])    #residue number
        self.delta=float(test[10])  #chemical shift
        self.atom=test[20]          #literal atom (H, C, N etc)
        self.at=test[7]             #atom type (CA, CB2 etc)
        pass



#major database entry: takes a folder and makes sense of it.
#relaxation data parsed out by lucas in xlsx files convered to
#numpy arrays and saved.
#the raw BMRB file itself is parsed, but only partially.
class bmrbEntry():
    #initialise entry with an index
    def __init__(self,bmrbid):
        print()
        print("Parsing bmrbID:",bmrbid)
        self.bmrbid=bmrbid  #bmrb index
        self.name=''        #default protein name: will try to scrape this from BMRB file later
        self.bmrbAss=False  #assume no bmrb assignment
        self.bmrbFile=False  #assume there is no bmrb file
        self.T=False        #assume we have no temperature
        self.seq=''
        self.bmrb={}        #initialise bmrb database entry
        self.bmrb['seq']=False
        self.bmrb['PDB']=False
        self.bmrb['T']=False
        self.bmrb['mw']=False
        
    #test to see if entry is a folder, that contains files,
    #parse contents.
    def Setup(self):
        try:  #fails if the file is just a file not a folder (hack! but mostly works. I've had to put some folders in as formal exceptions. Can be done better.)
            self.files=os.listdir(self.bmrbid)
        except:
            print('Not a real entry') 
            return False  #this is just a file. abort.
             
        self.CountExcelFiles() #work out how mnay xlsx files are in the folder. Figure out Lucas' main excel file.
        print(self.xlfile)     #excel files.
        print("Is there a Lucas order parameter file?",self.orderParameters)  #not parsing this yet, but this is Lucas' analysis in extended model free.

        print('Excel data:',self.excelData)  #the main relaxation file edited and curated by Lucas, with different fields on different sheets.
        self.ParseData()                     #parse main relaxation data file and crate data with numpy arrays

        self.bmrbfile=self.excelData.replace('.xlsx','.str')  #should be able to deduce the raw file from the relaxation data file...
        if(os.path.exists(self.bmrbfile)):   #did that work?
            self.bmrbFile=True               #yes! we have one.
            print('bmrbfile found:',self.bmrbfile)
            self.ParseBMRB()                 #parse it.
        else:
            print("cannot find the bmrb file. Looking harder and guessing....")
            possible=[]
            for file in self.files:
                if(file[-4:]=='.str'):
                    possible.append(file)
            if(len(possible)==1):
                print("The file is probably:",possible[0])
                self.bmrbFile=True
                self.bmrbfile=possible[0]
                #print(self.bmrbFile)
                #print(self.bmrbid)
                #sys.exit(100)
                    
            print("no bmrb file found:",self.bmrbid)
        return True  #by gosh we've done it. We have an entry and it appears not complete shite


    def TableParse(self,line):
        #print("Analysing:")
        #print(line)
        line=line.split('\n')[0]
        test=[]
        e=''
        w=1
        for c in line:
            if((c=='\'' or c=='\"') and w==1): #disable whitespace
                delim=c #set the delimiter
                w=0
                continue
            if(w==0 and c==delim): #turns whitespace back on.
                w=1
                continue

            if(w==0):#whitespace is off. append toe current entry
                e+=c
                continue

            #whitespace is on
            if(c==' '):
                if(len(e)>0):
                    test.append(e)
                e=''
                continue

            e+=c
                
        if(len(e)>0):
            test.append(e)
        #print("achieved")
        #print(test)
        return test
        
    #written a half arsed parsed to look at the raw BMRB file to get things needed.
    #not completely sure of its structure. There are sections marked by save_XXX ended by save_
    #have focused currently on getting the assignment out.
    #there are some useful entires like spectromter type we might end up going fishing for.
    #this will only incresae with time as we get a better fishing expedition.
    #the bool bmrbAss is set to true if we've found an assignment list.
    #currently if there's more than one for some reason, this could trip things.
    #also i've hard coded the residue entires: a better job would be to read the IDs as they are entered
    #in the file, and use those to make sense of the index. This is belt and braces.
    def ParseBMRB(self):

        self.bmrbAssIndx=''
        #self.PDBs=[]



        self.loop=0 #counter used to tell if we're in a table loop.
        self.header=0
        self.comment=0
        self.indxSave=''
        seqReady=0  #ready to read in a sequence?  
        
        inny=open(self.bmrbfile)
        indx=''
        for line in inny.readlines():
            test=line.split()
            if(len(test)==1):  #the entires are bookended in sections by save_. Retain the index of the section.
                if('save_' in test[0]):
                    entry=test[0].split('save_')  #we are on a bookend.
                    if(len(entry[1])>0):  #we are a bookend start!
                        #print('NewEntry:',entry[1])
                        indx=entry[1]     #save the index
                        self.bmrb[indx]={}  #blank our dictionary. should check if there is already an entry here.
                    else:
                        #print('EndEntry') #we are the end of a bookend. We don't care about this.
                        pass

            if('PDB' in line):
                #should be able to fish PDBs out from the files.
                #a straight search for PDB finds PDB ascention numbers.
                #but it looks irregular. some thinking and reading needed.
                #lukas has some PDBs downloaded already, so can link to those
                pass
            
                #print(indx,line)
            if('assigned_chem_shift_list' in indx):  #this seems to be the indx for an assignment. It lets you have more than one, which could get annoying. Should check.
                if(len(test)<14):
                    continue
                #print(indx,test)
                resi=int(test[5])                   #get the residue number (a helpful index for the assignment)
                if('ass' not in self.bmrb[indx]):   #is this residue in our assignment?
                    self.bmrb[indx]['ass']={}       #blank an entry if not.
                    self.bmrbAss=True               #we definitly have a populated assignment list for this entry 
                    self.bmrbAssIndx=indx           #this is the section index for the assignemnt in this BMRB file.
                self.bmrb[indx]['ass'][resi]=ass(test)  #add an assignment entry to the assignment linked to this residue index.

            if(len(test)==0):
                continue

            if(len(test)>1 and test[0]=='_Assembly.Name'):
                #first. hack to get the overall name
                if(len(test)==2):     #if just two entires on white-space splitting, we have a single string as our name
                    self.name=test[1]
                    continue
                self.name=line.split('\'')[1]  #I'm seeing entires that are delimited by quotation marks and so the names are strings. Go fish this out.
                continue


                
            #added the last one out of desparation for bmr4364. maybe search first for which index has entity assembly? bmr4365 also.
            if(indx=='assembly' or 'system_' in indx or 'stromelysin_PNU-99533'==indx or 'assembly' in indx or '_system' in indx or 'stromelysin' in indx):
                #now, if there are any loops, like the entity and DB lists, read them in, index dictionaries by BMRB entries.
                self.ReadLoopTables(test,line)
                              
            if('conditions' in indx or 'condition_1' in indx or 'condition_one' in indx or 'cond_set_1' in indx or 'Ex-cond_1' in indx or 'Condition_1' in indx): #could be a problem if there is more than one entry.
                self.ReadLoopTables(test,line)

            if('_Entity_assembly' in self.bmrb.keys()): #if we have read in already an entity, go looking for entries
                #if(self.bmrbid=='bmr4364'):
                #    print(self.bmrb['_Entity_assembly'])
                #    sys.exit(100)
                          
                
                print(self.bmrb['_Entity_assembly']['Entity_label'])
                #unique entity labels
                unis=np.unique(self.bmrb['_Entity_assembly']['Entity_label'])
                for u in unis: #for all named entities...
                    entity=u[1:] #get the entity take (dump the dollar sign, for fuck sack, people).
                    if(indx==entity):  #if the section index maps to a named identity, lets save the sequence
                        print('hello!')

                        if('Polymer_seq_one_letter_code' in line):
                            seqReady=1
                        
                        print(indx,entity)
                        if(test[0]==';' and self.loop==0 and seqReady==1):
                            self.loop=1
                            seq=''
                            continue
                        if(test[0]==';' and self.loop==1 and seqReady==1):
                            self.loop=0
                            if(self.bmrb['seq']==False):
                                self.bmrb['seq']={}
                            self.bmrb['seq'][indx]=seq
                            seq=''
                            seqReady=0
                            continue
                        if(self.loop==1):
                            seq+=test[0]
                            


        #Matthias buck: karplus JACS : NH bod lengths 1998 maybe.
        #Hbonding :
        #Bruschweiler order parameter from structure.



        print(self.bmrb.keys())               #these are all the section indices found in the file (i'm not populating most of these yet
        #print(self.bmrb['assembly.entity'])
        #print(self.bmrb['seq'])
        #print()
        #dig out the sequence and molecular weight.
        #shoudl check if there are more than one sequences here.
        if(self.bmrb['seq']!=False):
            seqPoss=[]
            seqLen=[]
            kill=[]
            for key,vals in self.bmrb['seq'].items():
                skip=0  #flag sequences that look like amino acid sequences
                if('{' in vals or '}' in vals or ',' in vals):
                    skip=1
                    kill.append(key)
                if(skip==0):
                    seqPoss.append(vals)
                    seqLen.append(len(vals))
            
            #if(self.bmrbid=='bmr51413'):
            #    print(seqPoss)
            #    #sys.exit(100)
            argy=np.argsort(seqLen)  #take longest sequence
            seq=seqPoss[argy[-1]]
            #print(self.bmrbid,key,seq)
            #if(self.bmrbid=='bmr50001'):
            #    print(seqLen)
            #    print(argy)
            #    print(self.bmrb['seq'])
            #    sys.exit(100)
            analysed_seq = ProteinAnalysis(seq)
            mw = analysed_seq.molecular_weight()
            print(f"Molecular weight: {mw:.2f} Da")
            self.bmrb['mw']=mw

            for key in kill:  #get rid of seqs that don't look like amino acid sequences
                del self.bmrb['seq'][key]
            
            #sys.exit(100)
        #print(self.bmrb['sample_condition'])

        #get linked PDBs in assembly.db
        #for i,typ in enumerate(self.bmrb['assembly.db'):
        print(self.bmrb.keys())

        if('_Assembly_db_link' in self.bmrb.keys()):
            for i,database in enumerate(self.bmrb['_Assembly_db_link']['Database_code']):
                print(database)
                if(database=='PDB'):
                    self.bmrb['PDB']=self.bmrb['_Assembly_db_link']['Accession_code'][i].upper()


                

        #try to set temperature
        if('_Sample_condition_variable' in self.bmrb.keys()):
            key='_Sample_condition_variable'
            #print(self.bmrbid,'looking at entry',key)
            #print(key)
            #print(self.bmrb['_Sample_condition_variable'])
            for i,typ in enumerate(self.bmrb[key]['Type']):
                if(typ=='temperature'):
                    self.T=self.bmrb[key]['Val'][i]
                    self.Tunit=self.bmrb[key]['Val_units'][i]


        if(self.T!=False):
            print("Sample temperature:",self.T,self.Tunit)

        #if(self.bmrbid=='bmr4364'):
        #    print(self.bmrb.keys())
        #    print(self.bmrb['_Sample_condition_variable'])
        #    sys.exit(100)
            
        print('Assignment?',self.bmrbAss)     #did we find an assignment?
        if(self.bmrbAss):                     #if we found an assignment, how many residues are covered?
            #print(self.bmrb[self.bmrbAssIndx]['ass'].keys())
            print('Number of assigned residues:',len(self.bmrb[self.bmrbAssIndx]['ass'].keys()))

    #manually naming this shit.
    #loopMap={}
    #loopMap['_Entity_assembly.']='assembly.entity'
    #loopMap['_Assembly_db_link.']='assembly.db'
    #loopMap['_Sample_condition_variable.']='sample_condition'
    #loopMap['_Assembly_bio_function.']='bio_function'
    #loopMap['_Bond.']='bond'
    #read a loop table.
    def ReadLoopTables(self,test,line):
        print("line:",line,self.indxSave,self.loop)
        if(test[0]=='loop_'):
            self.loop=1
            self.header=1
            self.comment=0
            self.listy={} #reset the dictionary, and save it.
            self.tabLine=[] #initialise table line
            return
        if(test[0]=='stop_'):
            self.loop=0
            self.header=0
            self.bmrb[self.indxSave]=self.listy
            self.listy={} #blank the dictionary having saved it
            if(len(self.tabLine)>0):
                print('shit: still have an unparsed line')
                print(self.bmrbid)
                print(self.tabLine)

                print('bmrbid:',self.bmrbid)
                print(self.listy.keys())
                print('keys:',len(self.listy.keys()))
                print('test:',len(test))
                print(test)
                print(len(self.tabLine))
                print('exiting.')
                sys.exit(100)
            return
        if(self.loop!=1):#skip if we're not in a loop section
            return
        print(len(test),self.header,test)
        if(len(test)==1 and self.header==1 and test[0]!=';'): #read in the header
            #saving an entity assembly
            #print(self.bmrbid,test)
            self.indxSave=test[0].split('.')[0] #get the bit before the dot as main index
            tast=test[0].split('.')[1]          #get the bit after the dot as column title
            self.listy[tast]=[]   
            #if(self.bmrbid=='bmr4365'):
            #    print('header:',self.listy.keys())

            #sys.exit(100)
            """
            for key,vals in self.loopMap.items():
                tast=test[0].split(key)
                if(len(tast)>1):
                    self.indxSave=
                    self.indxSave=vals
                    self.listy[tast[1]]=[]
                    return
            """
            return
        if(self.header==1): #if we get here, first time, we've completed the header (no breaks in header allowed)
            self.header=0
        
        #parse the line. it's whitespace delimited, but some text strings are delimited by apostrophes. a pain.
        #oh fuck. the entry can be multi line (6474). it permits comments. WHY!!!!

        if(test[0]==';'): #we are a comment.
            if(self.comment==0): #turn on save comment mode.
                self.comment=1   
                self.commentEntry=''
                return
            if(self.comment==1): #turn off save comment mode.
                self.tabLine.append(self.commentEntry)  #add the comment to the line.
                if(len(test)>1):
                    for i in range(len(test)-1): #add trailing entries after comment 
                        self.tabLine.append(test[i+1])
                self.comment=0   #turn off save comment mode.

                self.TrySaveTableLine()

                return
            
        if(self.comment==1):
            self.commentEntry+=line.split('\n')[0]  #add comment to the comment string
            return
            
        testy=self.TableParse(line)  #use new parser to joint unscramble whitespace + apostrophe delimits. 
        for t in testy:
            self.tabLine.append(t)  #add entry to the line.

        self.TrySaveTableLine()
            
        #print(line)
        #print(testy)
        """
        if(len(testy)!=len(self.listy.keys())):
            print('Oh no, struggling to parse this.')

            print('bmrbid:',self.bmrbid)
            print(self.listy.keys())
            print('keys:',len(self.listy.keys()))
            print('test:',len(test))
            print(line)

            if('Accession_code' in self.listy.keys()):#sometimes not enough entires are there
                pass
            else:
                sys.exit(100)
        """

    def TrySaveTableLine(self):
        print("Adding tabline:",len(self.tabLine),len(self.listy.keys()),self.tabLine)
        if(len(self.tabLine)==len(self.listy.keys())): #line is complete. Add.
            keys=list(self.listy.keys())
            for j in range(len(self.tabLine)):
                self.listy[keys[j]].append(self.tabLine[j])
            self.tabLine=[]
            

    #how many files are excel files?
    #the shortest xlsx file seems to be the one with the relaxation data...
    def CountExcelFiles(self):
        self.xlfile=[]
        self.xlfiles=0
        for file in self.files:
            if(file[-4:]=='xlsx'):
                self.xlfiles+=1
                self.xlfile.append(os.path.join(self.bmrbid,file))
        print('Excel files:',self.xlfiles)

        self.orderParameters=False
        self.orderParameterFile=''
        for file in self.xlfile:
            if('orderParameters' in file):
                self.orderParameters=True
                self.orderParameterFile=file

        lens=[]  #sort by length of file name
        for file in self.xlfile:
            lens.append(len(file))
        lens=np.array(lens)
        #print(lens)
        argy=np.argsort(lens)  #

        #taking file with shortest name as raw data edited by Lucas.
        self.excelData=self.xlfile[argy[0]] 

        
    #having identified lucas' excel file, go fishing to get the contents.
    #save in our own dictionary, indexed by entries, with numpy arrays for easy
    #access later on.
    def ParseData(self):
        xl=pd.ExcelFile(self.excelData)
        #print(xl.sheet_names)
        self.fields=[]
        for s in xl.sheet_names:
            self.fields.append(float(s))
        #print(self.fields)

        self.data={}
        for sheet in xl.sheet_names:
            data=xl.parse(sheet)
            field=float(sheet)

            dat={}
            #get this shit out of the fucking panda
            #dat dictionary: indexed by column headers
            #values will be numpy arrays
            for row in data.iterrows():
                for key, val in row[1].items():
                    if(key not in dat.keys()):
                        dat[key]=[]
                    dat[key].append(val)
            #turn entries into numpy arrays
            for key,vals in dat.items():
                dat[key]=np.array(vals)
            self.data[field]=dat

        self.entries=[]
        for field,dat in self.data.items():
            #print(field,list(dat.keys()))
            vs=[]
            for k,vals in dat.items():
                vs.append(len(vals))
            uni=np.unique(vs)
            if(len(uni)!=1):    #if we have more than one unique list lengths, then we probably have a problem. There is some missing data.
                print('shit')
                print('unequal list lengths')
                print(unis)
                sys.exit(100)   #hard abort if this is true: check the entry. why don't we have complete data for this residue? some kind of foobar.
            self.entries.append(int(uni[0]))
        #print(self.entries)
        for i,field in enumerate(self.data.keys()):    #this is how many fields we've found
            print('field:',field,'MHz entries:',self.entries[i])                

    #make a field string
    def MakeFieldStr(self):
        fieldStr=''
        for i,field in enumerate(self.fields):
            if(i!=0):
                fieldStr+=','
            fieldStr+='%.1f' % field
        return fieldStr
            
#a simple reporting class to make latex files
#start it, then add sections, titles and figures until you're ready to close it.
class reporter():
    def __init__(self):
        self.LATEXFILE='summary.tex'
        pass

    ### Initialise a latex report
    def InitLatex(self):
        if(os.path.exists(self.LATEXFILE)):
            os.system('rm '+self.LATEXFILE)
        outy=open(self.LATEXFILE,'w')
        outy.write('\\documentclass[showkeys,aps,prb,prepreint,amssymb, amsmath,nobibnotes]{article}\n')
        outy.write('\\usepackage[a4paper, margin=2cm]{geometry}\n')
        outy.write('\\usepackage{bm,setspace}\n')
        outy.write('\\usepackage{graphicx,graphics,booktabs}\n')
        outy.write('\\newlength{\\dhatheight}\n')
        #outy.write('\\newcommand{\doublehat}[1]{%\settoheight{\dhatheight}{\ensuremath{\hat{#1}}}%\addtolength{\dhatheight}{-0.35ex}%\hat{\vphantom{\rule{1pt}{\dhatheight}}%\smash{\hat{#1}}}}\n')
        outy.write('\\begin{document}\n')
        outy.write('\\section{BMRB report}\n');
        outy.close()

    ### Complete a latex report.
    def CloseLatex(self,tig=''):

        print("Compiling latex file")
        outy=open(self.LATEXFILE,'a')
        outy.write('\\end{document}\n')
        outy.close()
        #os.system('pdflatex '+self.LATEXFILE)
        os.system('pdflatex -interaction=nonstopmode '+self.LATEXFILE+' > latex.log')
        #if(self.GLOBAL== 'isotropic'):
        #    outpdf='pdf/'+self.baseTag+'.'+'iso'+'.pdf'
        #else:
        #    outpdf='pdf/'+self.baseTag+'.'+'ax'+'.pdf'
        
        #if(tig!=''):
        #    outpdf=outpdf.replace('.pdf',tig+'.pdf')
        #os.system('mv latex.pdf '+outpdf)
        #print("Created ",outpdf)
        os.system('rm latex.aux latex.log')

    ### called from dotcplot
    def AddFigs(self,outfigs,trim=(3,0,3,0)):
        outy=open(self.LATEXFILE,'a')
        #outy.write('\\begin{figure}[h]\n')
        for i in range(len(outfigs)):
            if(os.path.exists(outfigs[i])): #trim : left, bottom, right, top
                outy.write('\\noindent \\includegraphics[trim= %fcm %fcm %fcm %fcm ,width=%f\\textwidth]{%s' % (trim[0],trim[1],trim[2],trim[3],1./(len(outfigs)*0.99),outfigs[i]) +'}\n')
            else:
                pass
        #outy.write('\\caption[]{Relaxation rates versus correlation time $\\tau_c$ for spin system %s. macro=%s. %s }\n' % (self.baseTag,self.MACRO,self.plotTxt))
        #outy.write('\\end{figure} \n')

    def AddText(self,text):
        outy=open(self.LATEXFILE,'a')
        outy.write(text)
        outy.close()

    def ClearPage(self):
        outy=open(self.LATEXFILE,'a')
        outy.write('\\clearpage\n\n')
        outy.close()        

    def MakeBMRBsummary(self,entry):
        outy=open(self.LATEXFILE,'a')

        outy.write('\\noindent \\subsection{%s: %s}\n\n' % (entry.bmrbid,entry.name))
        
        #outy.write('\\begin{minipage}[t]{0.49\\textwidth}\n')
        outy.write('\\begin{tabular}{ll}\n')


        print(entry.bmrbid)
        print(entry.fields)

        
        outy.write('\\textbf{Fields:} &  %s MHz\\\\\n' % (entry.MakeFieldStr()))
        outy.write('\\textbf{DataPts:} &  %i\\\\\n' % (entry.FieldFits['N']))

        if(entry.T!=False):
            outy.write(' \\textbf{T:} & %s %s\\\\\n' % (entry.T,entry.Tunit))
        if(entry.bmrb['mw']!=False):
            outy.write('\\textbf{Mw:} & %.2f kDa\\\\\n' % (entry.bmrb['mw']/1000))



        if(entry.bmrb['seq']!=False):
            outy.write('\\textbf{Residues: } &  ')
            for i,(key,vals) in enumerate(entry.bmrb['seq'].items()):
                if(i!=0):
                    outy.write(',')
                outy.write(' %i' % (len(vals)))
            outy.write('\\\\\n')

        
            
        outy.write('\\end{tabular}\n\n')
        outy.close()

        
    def MakePlotLocReport(self,entry):

        self.typeMan='$R_2$  (s$^{-1}$)','NOE','$R_1$ (s$^{-1}$)','$\\sigma_{zz}$  (s$^{-1}$)','$\\frac{R_1}{R_2}$'
        

        outy=open(self.LATEXFILE,'a')

        outy.write('\\begin{small}\n')
        
        colStr='l'
        for field in entry.fields:        
            colStr+='|rlrl'
        outy.write('\\noindent \\begin{tabular}{%s}\n' % colStr)

        #outy.write('\\textbf{Field:} &  %i MHz\\\\\n' % (field))

        for field in entry.fields:
            #outy.write('& %.2f & & &  ' % field)
            outy.write('& \\multicolumn{4}{c}{%s MHz}' % field)
        outy.write('\\\\\n')

        for field in entry.fields:
            outy.write('&  \\multicolumn{2}{c}{rate} & \\multicolumn{2}{c}{$\\tau_c$ (ns)}   ')

        outy.write('\\\\\n')
        outy.write('\\hline')
        
        for i,typ in enumerate(self.typeMan):
            outy.write("%s  " % (typ))

            #N=len(fitty.yExp)
            #P=N
            #chi=self.yCalc-self.yExp
            #chi2=numpy.sum(  (chi/self.yErr)**2)
            
            for field in entry.fields:
                #figs=entry.FieldFits[field]['figs']
                GausFitsRaw=entry.FieldFits[field]['raw']
                GausFits=entry.FieldFits[field]['tcs']

                outy.write("&  %.3f & $\\pm$ %.3f   " % (GausFitsRaw[i][1],GausFitsRaw[i][2]))

                if(i<3):# we have fitted these.
                    outy.write("&  %.2f & $\\pm$ %.2f " % (GausFits[i][1]*1E9,GausFits[i][2]*1E9))
                else:
                    outy.write('& & ')
            outy.write('\\\\\n')
        outy.write('\\end{tabular}\n')

        outy.write('\\end{small}\n')
        outy.write('\n\n')


        #outy.write('\\noindent Datapoints/Parameters $N$/$P$: (%i/%i)\n\n'% (N,P))
        #for i,typ in enumerate(self.typs):
        #    #outy.write("$\\tau_c$ mean from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,numpy.average(self.tcVals[typ])*1E9,numpy.std(self.tcVals[typ])*1E9))
        #    outy.write("$\\tau_c$ from %s : %.3f $\\pm$ %.3f ns\n\n" %(typ,self.GausFits[i][1]*1E9,self.GausFits[i][2]*1E9))
        #outy.write('\\end{tabular}\n\n')
        outy.close()


        for field in entry.fields:
            figs=entry.FieldFits[field]['figs']
            self.AddFigs(figs,trim=(3,0,3,0))


        
############################3


    

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
#
#I've hacked this to work with the BMRB. This means stripping it back to the basics and interfacing with the BMRB data.
#As of 20th May 2026 i'm just running the straight transform from isotropic tumbling to take all rates and HetNOEs to a single correlation time,
#very crude, but gives a view of the protein.
class fittyRelax():
    def __init__(self,inst,pdbfile=''):
        self.Rel=inst  #save the RelCalc instance

        self.LATEXFILE='summary.tex'
        
        self.typs='R2','NOE','R1'  #for local analysis, use all 3 data types
        self.inst=self.Rel['isotropic']  #set current relcalc instance


        #self.InitLatex()    #initialise the report

        self.inst=self.Rel['isotropic']  #set current relcalc instance
        #self.MakeCorrPlot()
        #self.MakeJplot()

        """        
        #self.DoFitGlobal()           #fit with one tauC
        #self.DoFitResLocal()         #let each residue have its own tauC
        #self.DoFitModelFree()         #let each residue have its own tauC
        """
        #self.CloseLatex()            #conclude the report.


    #make a correlation plot, showing rates versus tauC
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

    #plot J (spectarl density) for a current field, tauC and order parameter.
    def MakeJplot(self):

        self.outfileJ='Jfile.out'
        oms=numpy.logspace(6,12,100)
        S2=0.1
        tc1=10E-9   #set global timescale
        tc2=1E-10   #set local timescale
        tce=(1/tc1+1/tc2)**(-1.)  #inverse timescale for model free
        outy=open(self.outfileJ,'w')
        for om in oms:
            J1=S2*tc1/(1+om**2*tc1**2)
            J2=(1-S2)*tce/(1+om**2*tce**2)
            J=J1+J2
            outy.write('%e\t%e\t%e\t%e\n' % (om/1E6/(2*numpy.pi),J,J1,J2))
        outy.write('\n\n')

        oms=[]

        omH=500.*1E6*2*numpy.pi     #proton frequency at given frq (rad s-1)
        omN=500.*1E6*2*numpy.pi/self.inst.gammaH*self.inst.gammaN  #nitrogen frequency at given frq (rad s-1)

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

        
    #set spectrometer frequency of relcalc instance
    def SetFrq(self,sfrq):
        self.inst.SetFreq(sfrq)
        self.inst.SetPars()

    #set paths and output files for current bmrb entry
    def SetBMRB(self,bmrb):
        self.outPath=bmrb
        self.outfile=os.path.join(self.outPath,'LocalTc.'+str(self.inst.sfrq)+'.out')                  #rates, isotropic fully local fit
        self.rawfile=os.path.join(self.outPath,'Raw.'+str(self.inst.sfrq)+'.out')                  #rates, isotropic fully local fit
        self.rawhistfile=os.path.join(self.outPath,'Rawhist.'+str(self.inst.sfrq)+'.out')                  #rates, isotropic fully local fit        

        self.figs=[]  #we're going to make 3 figures
        for i in range(5):
            self.figs.append(self.outfile.replace('.out','.fig'+str(i)+'.pdf'))
        self.GausFile=self.outfile.replace('.out','.gaus')
        self.GausFileRaw=self.outfile.replace('.out','.raw.gaus')

        
        #self.AddSection("%s" % bmrb,sub=True)
        
        
    #data data entry from lucas, turn it into one of our data objects.
    #note: need to sort this by field also, currently fields are being kept separate.
    def SetData(self,data):
        dat={}
        for source in 'R1','R2','NOE':
            dat[source]={}

            for i,res in enumerate(data['nmrResidueCode']):
                if(source=='R1'):
                    val=data['R1'][i]
                    err=data['R1_err'][i]
                if(source=='R2'):
                    val=data['R2'][i]
                    err=data['R2_err'][i]
                if(source=='NOE'):
                    val=data['HetNoe'][i]
                    err=data['HetNoe_err'][i]
                dat[source][res]=val,err

        self.dats=dat


    #Turn raw data dictionary into arrays
    #sorted by type (R1/R2/HetNOE)
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


    """
        
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
    """

    def ReadGausFile(self,filey):
        arr=[]
        inny=open(filey)
        for line in inny.readlines():
            test=line.split()
            l=[]
            for t in test:
                l.append(float(t))
            arr.append(l)
        inny.close()
        return arr
    
    def CheckIfAlreadyDone(self):
        for fig in self.figs: #do the figures exist?
            if(os.path.exists(fig)==False):
                print('missing:',fig)
                return 1
        #print('go:',go, os.path.exists(self.GausFile))
        if(os.path.exists(self.GausFile) and os.path.exists(self.GausFileRaw)): #if all figure files are there...
            self.GausFits=self.ReadGausFile(self.GausFile)
            self.GausFitsRaw=self.ReadGausFile(self.GausFileRaw)
            return 0

        if(os.path.exists(self.GausFile)==False ):
            print('missing',self.GausFile)
        if(os.path.exists(self.GausFileRaw)==False ):
            print('missing',self.GausFileRaw)
                  
        return 1
        
    #do one tauC fitted per relaxation rate.
    def DoFitLocal(self,FORCE=False):

        #first: are all the files and figures present? if yes, just assemble
        #if no, do complete calc.
        if(FORCE): #if forcing, then do all the things.
            go=1
        else:
            go=self.CheckIfAlreadyDone()

        if(go==1): #only do this is go=1
            self.RunLocalFits() #do the local fits
            self.MakePlotRaw()  #make pretty plots of raw data
            self.MakePlotLoc()  #make pretty plots.
        return go

    #execute local fits.
    def RunLocalFits(self):
        outy=open(self.outfile,'w');outy.close()  #blank output file.

        self.tcGuess=0.1E-9  #initialise

        self.tcVals={}   #save tauC values by data type.
        self.tcAll=[]    #save tauC values
        cnt=0
        for self.t,self.typ in enumerate(self.typs): #for each data type...
            if(self.typ=='NOE'):  #if NOE...
                self.tcVals[self.typ]=[]
                for self.key,vals in self.dats[self.typ].items(): #go over each datapoint one at a time
                    self.yExpLoc=self.yExp[cnt]  #get local rate
                    self.yErrLoc=self.yErr[cnt]  #get local rate

                    self.tcGuess=0.1E-9
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
                for jj,(self.key,vals) in enumerate(self.dats[self.typ].items()):  #go over each datapoint one at a time
                    self.yExpLoc=self.yExp[cnt]  #get local rate
                    self.yErrLoc=self.yErr[cnt]  #get local rate

                    if(self.typ=='R1'):
                        self.tcGuess=self.tcVals['R2'][jj]
                    else:
                        self.tcGuess=0.1E-9
                    self.fitLoc()                #fit to either R1 or R2 equations...
                    self.yCalc[cnt]=copy.deepcopy(self.yCalcLoc) #set calculated value
                    self.tcVals[self.typ].append(copy.deepcopy(self.tc)) #store values.
                    self.tcAll.append(copy.deepcopy(self.tc))            #store values
                    cnt+=1
            outy=open(self.outfile,'a');outy.write('\n\n');outy.close()  #split the file per data types to make a nice plot.
            print("tauC mean: ",self.typ,numpy.average(self.tcVals[self.typ]),numpy.std(self.tcVals[self.typ]))
        print("tauC overall: ",numpy.average(self.tcAll),numpy.std(self.tcAll))
        print(numpy.sum((((self.yCalc-self.yExp))/self.yErr)**2)/len(self.yCalc)) #print overall average chi2

        
    
    ################################
    #LOCAL, one tauC per rate.

    def pack(self):  #for 1 tc
        x=[]
        x.append(numpy.fabs(self.tc))
        return x
    def unpack(self,x): #for 1tc
        cnt=0
        self.tc=numpy.fabs(x[cnt]);cnt+=1
    def CalcRates(self):  #for R1 or R2
        self.inst.SetTauC(self.tc)
        self.yCalcLoc=self.inst.CalcRate(self.rho1,self.rho2)
    def GuessParams(self): # aclever function to get first guess for tauC
        self.tc=self.tcGuess
    def chi(self,x):   #get chi for leastsq
        self.unpack(x)
        self.CalcRates()
        return self.yCalcLoc-self.yExpLoc
    def fitLoc(self):  #fit loc
        self.GuessParams()
        x0=leastsq(self.chi,self.pack())

        outy=open(self.outfile,'a')  #write output for plotting.
        outy.write('%i\t%s\t%s\t%e\t%e\t%e\t%e\t%e\n' % (self.t,self.key,self.typ,self.yExpLoc,self.yCalcLoc,self.tc,(self.yExpLoc-self.yCalcLoc)/self.yErrLoc,self.yErrLoc))
        outy.close()

    ################################

    #calculate a steady-state NOE (needs R1 and NOE)
    def CalcRatesNOE(self):
        self.inst.SetTauC(self.tc)
        R1=self.inst.CalcRate('N1z','N1z')
        sig=self.inst.CalcRate('N1z','H2z')
        self.yCalcLoc=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)
    #get chi for an NOE
    def chiLocNOE(self,x):
        self.unpack(x)
        self.CalcRatesNOE()
        return self.yCalcLoc-self.yExpLoc
    def fitLocNOE(self): #fit a steady state NOE
        self.GuessParams()
        x0=leastsq(self.chiLocNOE,self.pack())
        
        outy=open(self.outfile,'a')  #write output
        outy.write('%i\t%s\t%s\t%e\t%e\t%e\t%e\t%e\n' % (self.t,self.key,self.typ,self.yExpLoc,self.yCalcLoc,self.tc,(self.yExpLoc-self.yCalcLoc)/self.yErrLoc,self.yErrLoc))
        outy.close()

    ##################
    #Fit a Gaussian
        
    def PackGaus(self):
        x=[]
        x.append(self.Agaus)
        x.append(self.X0gaus)
        x.append(self.Siggaus)
        #print(x)
        return x
    def UnPackGaus(self,x):
        self.Agaus=x[0]
        self.X0gaus=x[1]
        self.Siggaus=x[2]
    def CalcGaus(self):
        self.Ycgaus=self.Agaus*numpy.exp(-(self.X0gaus-self.Xgaus)**2/(2*self.Siggaus**2))
    def ChiGaus(self,x):
        self.UnPackGaus(x)
        self.CalcGaus()
        return self.Ycgaus-self.Ygaus
    def FitGaus(self,edges,hist):
        self.Xgaus=edges
        self.Ygaus=hist
        argy=numpy.argsort(hist) #sort low to high
        self.Agaus=numpy.max(hist)
        self.X0gaus=self.Xgaus[argy[-1]]
        self.Siggaus=self.X0gaus*0.1
        x0=leastsq(self.ChiGaus,x0=self.PackGaus())
        self.GausFits.append(self.PackGaus()) #save values


    def SetupGnu(self,gnuFile):  #setup a gnuplot file
        gnu=open(gnuFile,'w')
        gnu.write('set size square\n')
        gnu.write('set border 3\n')
        gnu.write('set tics nomirror\n')
        gnu.write('set term pdf\n')
        gnu.write('set key top left\n')
        return gnu
    def GnuFinalise(self,gnu):  #complete gnuplot file
        gnu.close()
        os.system('gnuplot gnu.gp')

    def MakeLogHist(self,vals):
        histMin=np.min(vals)
        if(histMin<0):
            histMin=0.01
        histMax=np.max(vals)
        histbins=np.logspace(np.log10(histMin),np.log10(histMax),30)
        hist,edges=np.histogram(vals,bins=histbins)
        edges=(edges[1:]+edges[:-1])*0.5
        return hist,edges

    
    ####MAKE LOCAL FIT PLOTS AND REPORT######
    def MakePlotRaw(self):  #make pretty plots.
        gnu=self.SetupGnu('gnu.gp')

        outy=open(self.rawfile,'w')
        for i,res in enumerate(self.residues):
            outy.write('%i\t' % res)
            for typ in self.typs:
                yd=self.dats[typ][res][0]
                ye=self.dats[typ][res][1]
                outy.write('%e\t%e\t' % (yd,ye))
            outy.write('%e\t' % ( (self.dats['NOE'][res][0]-1)*(self.inst.gammaN/self.inst.gammaH)*self.dats['R1'][res][0]))
            outy.write('\n')
        outy.close()


        outy=open(self.rawhistfile,'w')
        self.GausFits=[]
        for typ in self.typs:
            histRaw=[]
            for i,res in enumerate(self.residues):
                histRaw.append(self.dats[typ][res][0])
            hist,edges=self.MakeLogHist(histRaw)
            self.FitGaus(edges,hist)
            for i in range(len(edges)):
                outy.write('%e\t%e\t%e\n' % (edges[i],hist[i],self.Ycgaus[i]))
            outy.write('\n\n')
        histRaw=[]
        for i,res in enumerate(self.residues):
            histRaw.append(  (self.dats['NOE'][res][0]-1)*(self.inst.gammaN/self.inst.gammaH)*self.dats['R1'][res][0] )
        hist,edges=self.MakeLogHist(histRaw)
        self.FitGaus(edges,hist)
        for i in range(len(edges)):
            outy.write('%e\t%e\t%e\n' % (edges[i],hist[i],self.Ycgaus[i]))
        outy.write('\n\n')

        histRaw=[]
        for i,res in enumerate(self.residues):
            histRaw.append(  (self.dats['R1'][res][0]/self.dats['R2'][res][0]) )
        Hist,edges=self.MakeLogHist(histRaw)
        self.FitGaus(edges,hist)
        for i in range(len(edges)):
            outy.write('%e\t%e\t%e\n' % (edges[i],hist[i],self.Ycgaus[i]))
        outy.write('\n\n')
        outy.close()

        self.GausFitsRaw=self.GausFits

        #dump file
        outy=open(self.GausFileRaw,'w')
        for i in range(len(self.GausFitsRaw)):
            for j in range(3):
                outy.write('%e\t' % (self.GausFitsRaw[i][j]))
            outy.write('\n')

        
        #NOE=1+(self.inst.gammaH/self.inst.gammaN)*(sig/R1)

        #####plot tauc versus residue######
        gnu.write('set title \'Raw data %.1f MHz\'\n' % self.inst.sfrq)
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'rate (s-1)\'\n')
        gnu.write('set pointintervalbox 0\n')
        gnu.write('set output \'figraw.pdf\'\n')
        gnu.write('set output \'%s\'\n' % self.figs[1])
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' u 1:%i:%i ti \'%s\'w err pt 7 ps 0.2 ' % (self.rawfile,i*2+2,i*2+3,typ))
        gnu.write(',\'%s\' u 1:%i ti \'%s\'w points pt 7 ps 0.2 ' % (self.rawfile,8,'sig'))
        gnu.write('\n')

        """
        gnu.write('set output \'%s\'\n' % self.figs[1])
        gnu.write('plot ')
        for i,typ in enumerate(('NOE','R1')):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' u 1:%i:%i ti \'%s\'w err pt 7 ps 0.2 ' % (self.rawfile,i*2+4,i*2+5,typ))
        gnu.write(',\'%s\' u 1:%i ti \'%s\'w points pt 7 ps 0.2 ' % (self.rawfile,8,'sig'))
        gnu.write(',\'%s\' u 1:($6/$2) ti \'%s\'w points pt 7 ps 0.2 ' % (self.rawfile,'R1/R2'))


        gnu.write('\n')
        """

        """
        gnu.write('set output \'%s\'\n' % self.figs[1])
        gnu.write('plot ')
        #for i,typ in enumerate(('NOE','R1')):
        #    if(i!=0):
        #        gnu.write(',')
        #    gnu.write('\'%s\' u 1:%i:%i ti \'%s\'w err pt 7 ps 0.2 ' % (self.rawfile,i*2+4,i*2+5,typ))
        gnu.write('\'%s\' u 1:%i ti \'%s\'w points pt 7 ps 0.2 ' % (self.rawfile,8,'sig'))
        gnu.write(',\'%s\' u 1:($6/$2) ti \'%s\'w points pt 7 ps 0.2 ' % (self.rawfile,'R1/R2'))
        gnu.write('\n')
        """
        

        gnu.write('set output \'%s\'\n' % self.figs[0])
        gnu.write('set logscale x\n')
        gnu.write('set xlabel \'value\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set xrange[0.01:*]\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 1:2 ti \'%s\'w boxes lc %i,\'\' i %i u 1:3 w li noti lc %i' % (self.rawhistfile,i,typ,i+1,i,i+1))
        gnu.write(',\'%s\' i %i u 1:2 ti \'%s\'w boxes lc %i,\'\' i %i u 1:3 w li noti lc %i ' % (self.rawhistfile,3,'sig',4,3,4))
        gnu.write(',\'%s\' i %i u 1:2 ti \'%s\'w boxes lc %i,\'\' i %i u 1:3 w li noti lc %i ' % (self.rawhistfile,4,'R1/R2',5,4,5))
        gnu.write('\n')
        

        self.GnuFinalise(gnu)



    def MakePlotLoc(self):
        fig0=2
        
        gnu=self.SetupGnu('gnu.gp')
        #####plot tauc versus residue######
        gnu.write('set title \'Isotropic model, local timescale fit\'\n')
        gnu.write('set xlabel \'residue\'\n')
        gnu.write('set ylabel \'tc (ns)\'\n')
        gnu.write('unset key\n')
        gnu.write('set logscale y\n')
        gnu.write('set output \'%s\'\n' % self.figs[fig0+1])
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 2:($6*1E9) ti \'%s\'w points pt 7 ps 0.2 ' % (self.outfile,i,typ))
        gnu.write('\n')

        ######plot calculated rate versus expeirmental rate######
        gnu.write('set xlabel \'Rate(calc)\'\n')
        gnu.write('set ylabel \'Rate(exp)\'\n')
        gnu.write('set key top left\n')
        gnu.write('unset logscale y\n')
        gnu.write('set pointintervalbox 0\n')
        gnu.write('set output \'%s\'\n' % self.figs[fig0])
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u 4:5:8 ti \'%s\'w err pt 7 ps 0.2' % (self.outfile,i,typ))
        gnu.write('\n')

        #make histograms of taucs, fit to Gaussians and and plot.
        self.bins=30
        self.histfile='hist.out'
        outy=open(self.histfile,'w')
        tcBins=numpy.logspace(-10,-7,50)
        self.GausFits=[]
        for i,typ in enumerate(self.typs):
            binmax=numpy.max(self.tcVals[typ])
            binmin=numpy.min(self.tcVals[typ])
            #bins=numpy.linspace(binmin,binmax,self.bins)
            hist,edges=numpy.histogram(self.tcVals[typ],bins=tcBins)
            edges=(edges[1:]+edges[:-1])*0.5
            #quickly fit a gaussian and plot that also.
            self.FitGaus(edges,hist)
            for j in range(len(edges)):
                outy.write('%e\t%e\t%e\n' % (edges[j],hist[j],self.Ycgaus[j]))
            outy.write('\n\n')
        outy.close()
        gnu.write('set xlabel \'tc (ns)\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set output \'%s\'\n' % self.figs[fig0+2])
        gnu.write('set logscale x\n')
        gnu.write('unset logscale y\n')
        gnu.write('plot ')
        for i,typ in enumerate(self.typs):
            if(i!=0):
                gnu.write(',')
            gnu.write('\'%s\' i %i u ($1*1E9):2 ti \'%s\'w boxes lc %i,\'\' i %i u ($1*1E9):3 noti w li lc %i' % (self.histfile,i,typ,i+1,i,i+1))
        gnu.write('\n')


        self.GnuFinalise(gnu)

        outy=open(self.GausFile,'w')
        for i in range(len(self.GausFits)):
            for j in range(3):
                outy.write('%e\t' % (self.GausFits[i][j]))
            outy.write('\n')



    ##############################################
    #HELPER FUNCTIONS FOR REPORT GENERATION
    def AddSection(self,sec,sub=False):
        outy=open(self.LATEXFILE,'a')
        if(sub):
            outy.write('\\subsection{%s}\n\n' % sec)
        else:
            outy.write('\\section{%s}\n\n' % sec)
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
            outy.write('\\noindent\\includegraphics[width=%f\\textwidth]{%s}\n' % (width, outfigs[i]))
        #outy.write('\\end{figure} \n')

        


############################

        
def MakeOverallHistograms():
    R2list=[]
    R1list=[]
    NOElist=[]
    R1min=0
    R1max=5
    R2min=0
    R2max=100
    NOEmin=0
    NOEmax=1.05
    bins=30
    R1bins=np.linspace(R1min,R1max,bins)
    R2bins=np.linspace(R2min,R2max,bins)
    NOEbins=np.linspace(NOEmin,NOEmax,bins)

    for b,entry in bmrb.items():
        for field,data in entry.data.items():
            #print(b,data.keys())
            #print(data['R1'])
            #print(data['R2'])
            #print(data['HetNoe'])
            histR1,edgesR1=np.histogram(data['R1'],bins=R1bins)
            histR2,edgesR2=np.histogram(data['R2'],bins=R2bins)
            histNOE,edgesNOE=np.histogram(data['HetNoe'],bins=NOEbins)
            #print(edgesR1.shape,histR1.shape)

            R1list.append(histR1)
            R2list.append(histR2)
            NOElist.append(histNOE)
    R1list=np.array(R1list)
    R2list=np.array(R2list)
    NOElist=np.array(NOElist)

    #print(edgesR1,histR1)

    R1mid=(edgesR1[1:]+edgesR1[:-1])*0.5
    R2mid=(edgesR2[1:]+edgesR2[:-1])*0.5
    NOEmid=(edgesNOE[1:]+edgesNOE[:-1])*0.5


    outy=open('outy.out','w')
    for i in range(bins-1):
        outy.write('%e\t%e\t' % (R1mid[i],np.sum(R1list[:,i])))
        outy.write('%e\t%e\t' % (R2mid[i],np.sum(R2list[:,i])))
        outy.write('%e\t%e\t' % (NOEmid[i],np.sum(NOElist[:,i])))
        outy.write('\n')
    outy.write('\n\n')
    for j in range(len(R1list)):
        for i in range(bins-1):
            outy.write('%e\t%e\t' % (R1mid[i],(R1list[j,i])))
            outy.write('%e\t%e\t' % (R2mid[i],(R2list[j,i])))
            outy.write('%e\t%e\t' % (NOEmid[i],(NOElist[j,i])))
            outy.write('\n')
        outy.write('\n\n')
    outy.close()

    gnu=open('gnu.gp','w')
    gnu.write('set term pdf\n')
    gnu.write('set size square\n')
    gnu.write('set border 3\n')
    gnu.write('set tics nomirror\n')
    gnu.write('unset key\n')
    gnu.write('set ylabel \'count\'\n')

    gnu.write('set output \'%s\'\n' % ('R1fig.pdf'))
    gnu.write('set title \'R1s: %i\'\n' % np.sum(R1list))
    gnu.write('set xlabel \'R1(s-1)\'\n')
    gnu.write('plot ')
    gnu.write('\'%s\' i 0 u 1:2 w li lw 2 lc 0' % ('outy.out'))
    for j in range(len(R1list)):
        gnu.write(',\'%s\' i %i u 1:2 w li' % ('outy.out',j+1))
    gnu.write('\n')

    gnu.write('set output \'%s\'\n' % ('R2fig.pdf'))
    gnu.write('set title \'R2s: %i\'\n' % np.sum(R2list))
    gnu.write('set xlabel \'R2(s-1)\'\n')
    gnu.write('plot ')
    gnu.write('\'%s\' i 0 u 3:4 w li lw 2 lc 0' % ('outy.out'))
    for j in range(len(R1list)):
        gnu.write(',\'%s\' i %i u 3:4 w li' % ('outy.out',j+1))
    gnu.write('\n')

    gnu.write('set output \'%s\'\n' % ('NOEfig.pdf'))
    gnu.write('set title \'NOEs: %i\'\n' % np.sum(NOElist))
    gnu.write('set xlabel \'HetNOE\'\n')
    gnu.write('plot ')
    gnu.write('\'%s\' i 0 u 5:6 w li lw 2 lc 0' % ('outy.out'))
    for j in range(len(R1list)):
        gnu.write(',\'%s\' i %i u 5:6 w li' % ('outy.out',j+1))
    gnu.write('\n')

    os.system('gnuplot gnu.gp')

def WriteHeader():
    print()
    print('Total bmrb entries:',len(bmrb.keys()))
    cnt=0
    fieldHist={}
    fields={}
    for b,entry in bmrb.items():

        #work out how many have assignments
        if(entry.bmrbAss):
            cnt+=1

        #work out how many fields we have
        if(len(entry.fields) not in fieldHist.keys()):
            fieldHist[len(entry.fields)]=0
        fieldHist[len(entry.fields)]+=1

        for field in entry.fields:
            v=round(int(field),-1)
            if(v not in fields):
                fields[v]=0
            fields[v]+=1


    print('Number with assignments:',cnt)
    print("Field histogram: number of entries with XXX number of fields:")
    print(fieldHist)

    report.AddText('\\noindent \\textbf{Total bmrb entries:} %i \n\n' % (len(bmrb.keys())))
    report.AddText('\\noindent \\textbf{Number with assignments:} %i \n\n' % (cnt))
    report.AddText('\\noindent \\textbf{Number of fields in dataset:} \n\n')
    report.AddText('\\noindent \\begin{tabular}{cc}\n\n')
    report.AddText('Fields & Count \\\\\n')
    keys=list(fieldHist.keys())
    keys=np.sort(keys)
    for key in keys:
        report.AddText('%i & %i \\\\\n' % (key,fieldHist[key]))
    report.AddText('\\end{tabular}\n\n\n')
    report.AddText('\\noindent \\textbf{Fields used:} \n\n')
    report.AddText('\\noindent \\begin{tabular}{cc}\n\n')
    report.AddText('Field & Count \\\\\n')
    keys=list(fields.keys())
    keys=np.sort(keys)
    for key in keys:
        report.AddText('%i & %i \\\\\n' % (key,fields[key]))
    report.AddText('\\end{tabular}\n\n\n')
    report.AddText('\\noindent \\textbf{Histogram of values, R1, R2, HetNOE} \n\n')
    report.AddFigs( ('R2fig.pdf','NOEfig.pdf','R1fig.pdf'))



def WriteBigTable(argy=[]):    
    colStr='lll' #bmrb, name,fields
    colStr+='llll' #T residues Mw pdb

    report.AddText('\\begin{tiny}')
    report.AddText('\\noindent \\begin{tabular}{%s}\n\n' % colStr)
    report.AddText('\\textbf{BMRB} & \\textbf{Name} & \\textbf{Fields}')
    report.AddText('& \\textbf{T} & \\textbf{Residues} & \\textbf{$M_w$ (kDa)} & \\textbf{PDB}')
    report.AddText('\\\\\n')
    report.AddText('\\hline\n')

    keys=list(bmrb.keys())
    if(len(argy)!=0):
        keys=np.array(keys)[argy]

    for b in keys:
        entry=bmrb[b]
        
        
        fieldStr='(%i) ' % len(entry.fields)
        for i,f in enumerate(entry.fields):
            if(i!=0):
                fieldStr+=','
            fieldStr+='%i' % int(f)
            name=entry.name.replace('N-terminal domain','NTD').replace('of Tissue Inhibitor of','').replace('Metalloproteinases-','MP').replace('MatrixMetalloProteinase-','MMP').replace('Bromodomain ','BD').replace('(BD1)','').replace('(BD2)','').replace('kinase associated protein phosphatase','KAPP').replace('Human Neutrophil Gelatinase-Associated Lipocalin','')
        report.AddText('%s & \\verb|%s| & %s' % (b,name,fieldStr))
        #try:
        #    entry.bmrb
        #except:
        #    report.AddText(' no bmrb& & & ')
        #    report.AddText('\\\\\n')
        #    continue

        if(entry.T!=False):
            report.AddText('& %s' % (entry.T))
        else:
            report.AddText('& ')

        if(entry.bmrb['seq']!=False):
            #print(len(entry.bmrb['seq']))  ##CHECK IF MORE THAN ONE SEQ!
            report.AddText('& ')


            for i,(key,vals) in enumerate(entry.bmrb['seq'].items()):
                if(i!=0):
                    report.AddText(',')
                report.AddText(' %i' % (len(vals)))
        else:
            report.AddText('& ')


        if(entry.bmrb['mw']!=False):
            report.AddText('& %.2f' % (entry.bmrb['mw']/1000.))
        else:
            report.AddText('& ')



        if(entry.bmrb['PDB']!=False):
            report.AddText('& %s' % (entry.bmrb['PDB']))
        else:
            report.AddText('& ')


        ####add in tauC values


        report.AddText('\\\\\n')

    report.AddText('\\end{tabular}\n\n\n')
    report.AddText('\\end{tiny}')

#do single rate correlation time analysis.
#if files exist already, just copy up the values
def DoSingleCorrelationAnalysis(FORCE=False,stop=-1):
    fitty=fittyRelax(runs) #setup instance of data fitter
    cnt=0
    go=0
    for b,entry in bmrb.items():
        entry.FieldFits={}
        entry.FieldFits['N']=0

        for field,data in entry.data.items():

            print(b,data.keys())

            fitty.SetBMRB(b)
            fitty.SetFrq(field)
            fitty.SetData(data)
            fitty.AssembleData()          #setup global data arrays from raw data, merging the types.        
            goCurrent=fitty.DoFitLocal(FORCE=FORCE)            #do a local fit, each rate gets its own tc.

            entry.FieldFits[field]={}        
            entry.FieldFits[field]['figs']=fitty.figs
            entry.FieldFits[field]['raw']=fitty.GausFitsRaw
            entry.FieldFits[field]['tcs']=fitty.GausFits
            entry.FieldFits['N']+=len(fitty.yExp)
            if(goCurrent==1): #we have done a fit! so downstream we need to edit things.
                go=1

        #break
        cnt+=1
        if(cnt==stop):
            break
    return go
#write single correlation analysis results to report
def WriteSingleCorrelationAnalysis(stop=-1,argy=[]):
    cnt=0
    report.AddText("\\section{Individual fits}\n\n")
    keys=list(bmrb.keys())
    if(len(argy)!=0):
        keys=np.array(keys)[argy]

    for b in keys:
        entry=bmrb[b]
        report.MakeBMRBsummary(entry)            
        report.MakePlotLocReport(entry)
        #break
        cnt+=1
        if(cnt==stop):
            break
    
    
####################    
#### REAL START ####
####################
#setup relcalc instance

from RelCalc import RelCalc
from scipy.optimize import leastsq
import copy,os,pickle,numpy,sys,math

############################################################
# Code to calculate relaxation rates in a NH spin system
############################################################

EXTH=False
EXTD=False

SYM=True    # Calculate commutators symbolically?
SPARSE=False  # If not calculating commutators symbolically, use sparse matrices?
CSA=True #False

runs={}

picklFile='pickle.rick'
FORCE=False

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
        
        
###########################################################
#off we go!

#some settings
stop=-1      #if want the program to stop early
FORCE=False  #if we want to re-do the analysis from scratch


bmrb={}  #initialise the conents of what will be the master 
contents=os.listdir('./')
for folder in contents:
    if(folder=='pdf' or folder=='__pycache__' or folder=='tmp'): 
        continue
    entry=bmrbEntry(folder)
    success=entry.Setup()
    if(success):
        bmrb[folder]=entry

MakeOverallHistograms() #lets make a histogram of experimental measurements

report=reporter()       #start up a report.
report.InitLatex()

WriteHeader()           #write overall BMRB statistics for the analysis

go=DoSingleCorrelationAnalysis(FORCE=FORCE,stop=stop) #either do single tc fits, or copy values from files


if(go==1):
    #make an output file for correlations and plotting
    outy=open('test.out','w')
    for b,entry in bmrb.items():
        for field,data in entry.data.items():
            print(b,data.keys())
            try:
                if(entry.T!=False and entry.bmrb['mw']!=False):
                    outy.write('%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n' % (entry.T,entry.bmrb['mw'],entry.FieldFits[field]['tcs'][0][1],entry.FieldFits[field]['tcs'][1][1],entry.FieldFits[field]['tcs'][2][1],entry.FieldFits[field]['raw'][0][1],entry.FieldFits[field]['raw'][1][1],entry.FieldFits[field]['raw'][2][1]))
            except:
                pass
    outy.close()

tcFigs=[]
tcFigs.append('tc1.pdf')
tcFigs.append('tc2.pdf')
tcFigs.append('tc3.pdf')

gnu=open('gnu.gp','w')
gnu.write('set term pdf\n')
gnu.write('set border 3\n')
gnu.write('set size square\n')
gnu.write('set tics nomirror\n')
gnu.write('set xlabel \'tc(fit) (ns)\'\n')
gnu.write('set ylabel \'tc(estimated) (ns)\'\n')
gnu.write('unset key\n')

gnu.write('set output \'%s\'\n' % tcFigs[0])
gnu.write('set title \'R2\'\n')
gnu.write('plot \'test.out\' u ($3*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 1,x\n')

gnu.write('set output \'%s\'\n' % tcFigs[1])
gnu.write('set title \'NOE\'\n')
gnu.write('plot \'test.out\' u ($4*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 2,x\n')

gnu.write('set output \'%s\'\n' % tcFigs[2])
gnu.write('set title \'R1\'\n')
gnu.write('plot \'test.out\' u ($5*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 3,x\n')
gnu.write('\n')
gnu.close()
os.system('gnuplot gnu.gp')

#report.AddFigs(tcFigs,trim=(3,0,3,0))
report.AddFigs(tcFigs)

report.AddText('$\\tau_c$ estimated from $\\tau_c=\\frac{\eta V}{k_b T}$. We can substitute for molecular weight and obtain $\\tau_c=\\frac{\eta}{\\rho}\\frac{M_w}{k_b T}$ where $M_w$ is in Da and $\\rho$ is the density. The line takes the density to be water, and the viscosity to be 0.1 cP. The estimated value is compared to the fitted values from the three measurements.')

#18477 700Mhz has had datatypes switched

"""
density=mass/vol

tc= nu V/ kT
tc= (nu/density) mass/ kT
"""

report.ClearPage()      #gap.


def SortBMRB(sortArg='mw'):
    arr=[]
    if(sortArg=='mw'):
        for b,entry in bmrb.items():
            arr.append(entry.bmrb['mw'])
    else:
        print('Give a sensible suggestion for sorting the bmrb entries.')
        sys.exit(100)
    return(numpy.argsort(arr))
            

argy=SortBMRB()
WriteBigTable(argy=argy)         #write overview table

report.ClearPage()      #gap.



WriteSingleCorrelationAnalysis(stop=stop,argy=argy)  #write single tc fits, one entry per protein

report.CloseLatex()

