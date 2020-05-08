#!/usr/bin/env python3

import sys #use interactive with user
import pandas as pd #use pandas dataframe

def get_dnatxt(filename): #run this chunk with the (filename) given by user
    with open(filename, "r") as f: #open the file as entered by the user
      dnatxt = f.read() #reading text in the file store it to dnatxt
    return(dnatxt) #return the value of this function as "dnatxt"

def get_dnalst(dnatxt): #run this chunk with (dnatxt) from previous function
    """
    Take a string of text and removes unnecessary symbols
      Args:
        dnatxt: string
    Returns:
      string
    """
    dnalst = dnatxt #set dnalst as dnatxt for filtering unnecessary mark that will disturb data processing
    unlsted = [".",",",";",":",'"'," ","'","@","#","*"] #unlsted contains symbols that need to be removed
    for mark in unlsted: #loop through punctuation list
      Tseq = dnalst.replace(mark,"") #replace each punctuation mark in the dnatext and save the new string in dnalist
    return(Tseq) #return the value of this function as "Tseq"

###########################################################################################
###main function

def kMers(Tseq, k):
    seq = []
    obs = []
    pos = []
    for x in range(k): ###loop from 1 to k
        count = x+1 #because python base from 0, so add 1 to start from 1
        kFreq ={} ###setting up dictionary
        seqc =len(Tseq)-count+1
        for i in range(seqc):#loop to know observation and possibility for each sequence
            kmer = Tseq[i:i + count]
            if kmer in kFreq:
                kFreq[kmer] +=1
            else:
                kFreq[kmer] = 1
            obsv = len(kFreq)
            posb = seqc
        seq.append(kFreq)
        obs.append(obsv)
        pos.append(posb)
    obs.append(sum(obs))
    pos.append(sum(pos))
    p='done'###add new row to match the total row in the dataframe
    seq.append(p)
    df_dna = pd.DataFrame({
        'k': range(1,k+2),###add first column, add+2 because need a new row to match total from obs and poss
        'Observasi':obs,###add second column
        'Possible':pos,###add third column
        'Sequence' :seq,###add fourth column
    })
    
    df_dna.loc[k,'k'] = 'Total' #change the last row of 'k' column to "Total"  
    return(df_dna) #return the value of this function as dataframe "df_dna"

def ling_compx(df_dna,k):
    a = df_dna.loc[k].at['Observasi'] #set "a" as the value of total observed count, taken from the last row of column 'observed kmers'
    b = df_dna.loc[k].at['Possible']#set "b" as the value of total observed count, taken from the last row of column 'observed kmers'
    ling_compx = int(a)/int(b) #Calculate the linguistic complexity as a ratio of obs/pos
    return(ling_compx)
  
    ########################################
if __name__ == '__main__': #if this script run, run the following
    myfile = sys.argv[1] #load myfile in the command line
    k = int(sys.argv[2]) #load the k (sequence number) in the command line
    try:
      dnatxt = get_dnatxt(myfile) #get list of DNA from user's DNA file input
    except IOError: #if error in reading the file, print the following messages
      print("Could not read file: ",sys.argv[1],", write the correct filename follow by the sequence") 
      print("for example: python3[space]pythonassign.py[space]filename.txt[space]7")
      
    
    Tseq= get_dnalst(dnatxt)
    get_kMers= kMers(Tseq, k)#get kMers using get_kMers function
    linguistic= ling_compx(get_kMers,k)#get linguistic complexity using linguistic function
    
     #the following is to print the table and linguistic complexity on the screen
    print(get_kMers)
    print("Linguistic Complexity is ",linguistic)
    
      #the following print the table in file, named kmers_table.txt
    with open('kmers_table.txt', 'w') as f:
      print(get_kMers, file=f)
