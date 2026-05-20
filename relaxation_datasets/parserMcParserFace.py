#!/usr/bin/env python3

import os,sys,pandas as pd,numpy as np

#take an assignment row with chemical shifts
#parse it and save class variables
class ass():
    def __init__(self,test):
        self.resn=test[4]
        self.resi=float(test[5])
        self.delta=float(test[10])
        self.atom=test[20]
        self.at=test[7]

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
        self.bmrbid=bmrbid

        self.bmrbAss=False  #assume no bmrb assignment
        
    #test to see if entry is a folder, that contains files,
    #parse contents.
    def Setup(self):
        try:  #fails if the file is just a file not a folder (hack! but works...)
            self.files=os.listdir(self.bmrbid)
        except:
            print('Not a real entry') 
            return False  #this is just a file. abort.
             
        self.CountExcelFiles()#work out how mnay xlsx files
        print(self.xlfile)    #write the main data excel file with raw relaxation data
        print("Is there a Lucas order parameter file?",self.orderParameters)  #not parsing this yet, but this is Lucas' analysis

        print('Excel data:',self.excelData)  #the main relaxation file
        self.ParseData()                     #parse main relaxation data file and crate data with numpy arrays


        self.bmrbFile=False  #assume there is no bmrb file
        self.bmrbfile=self.excelData.replace('.xlsx','.str')  #should be able to deduce the raw file from the relaxation data file...
        if(os.path.exists(self.bmrbfile)):   #did that work?
            self.bmrbFile=True               #yes!
            print('bmrbfile found:',self.bmrbfile)
            self.ParseBMRB()                 #parse it.
        
        return True


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
        self.bmrb={}

        self.bmrbAssIndx=''
        self.PDBs=[]
        inny=open(self.bmrbfile)
        indx=''
        for line in inny.readlines():
            test=line.split()
            if(len(test)==1):
                if('save_' in test[0]):
                    entry=test[0].split('save_')
                    if(len(entry[1])>0):
                        #print('NewEntry:',entry[1])
                        indx=entry[1]
                        self.bmrb[indx]={}
                    else:
                        #print('EndEntry')
                        pass

            if('PDB' in line):
                #should be able to fish PDBs out from the files.
                #a straight search for PDB finds PDB ascention numbers.
                #but it looks irregular. some thinking and reading needed.
                #lukas has some PDBs downloaded already, so can link to those
                pass
            
                #print(indx,line)
            if('assigned_chem_shift_list' in indx):
                if(len(test)<14):
                    continue
                #print(indx,test)
                resi=int(test[5])
                if('ass' not in self.bmrb[indx]):
                    self.bmrb[indx]['ass']={}
                    self.bmrbAss=True
                    self.bmrbAssIndx=indx
                self.bmrb[indx]['ass'][resi]=ass(test)

                
        print(self.bmrb.keys())
        print('Assignment?',self.bmrbAss)
        if(self.bmrbAss):
            #print(self.bmrb[self.bmrbAssIndx]['ass'].keys())
            print('Number of assigned residues:',len(self.bmrb[self.bmrbAssIndx]['ass'].keys()))


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

        #taking file with shortest name as raw data.
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
            if(len(uni)!=1):
                print('shit')
                print('unequal list lengths')
                print(unis)
                sys.exit(100)
            self.entries.append(int(uni[0]))
        #print(self.entries)
        for i,field in enumerate(self.data.keys()):
            print('field:',field,'MHz entries:',self.entries[i])                
        #sys.exit(100)
            
       

#################################
#off we go!
bmrb={}
contents=os.listdir('./')
for folder in contents:
    entry=bmrbEntry(folder)
    success=entry.Setup()
    if(success):
        bmrb[folder]=entry


print()
print('Total bmrb entries:',len(bmrb.keys()))
cnt=0
fieldHist={}
for b,entry in bmrb.items():

    #work out how many have assignments
    if(entry.bmrbAss):
        cnt+=1

    #work out how many fields we have
    if(len(entry.fields) not in fieldHist.keys()):
        fieldHist[len(entry.fields)]=0
    fieldHist[len(entry.fields)]+=1
    
print('Number with assignments:',cnt)
print("Field histogram: number of entries with XXX number of fields:")
print(fieldHist)


