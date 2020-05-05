from ROOT import *
import sys, os, time

t0 = time.time() 

arg1 = sys.argv[1]  

if arg1 == 'Data': 
        input_dir = '../2lep/Data/'
elif arg1 == 'Signal': 
        input_dir = '../2lep/MC/BSM_Signal_Samples/'
elif arg1 == 'Background': 
        input_dir = '../2lep/MC/SM_Backgrounds/'


myChain = TChain('mini') 


for filename in os.listdir(input_dir):
        if not '.root' in filename: continue 
        print (filename)  
        myChain.Add(input_dir+filename) 


entries = myChain.GetEntries() 

print "-------------------------------------------"
print "Number of events to process: %d" %entries
print "-------------------------------------------"

if arg1 == 'Data': 
        myChain.Process("MySelector.C", "Data")
else: 
        myChain.Process("MySelector.C", "MC") 

t = int( time.time()-t0 )/60  

print "Time spent: %d min" %t 
