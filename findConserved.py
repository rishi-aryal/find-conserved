#!/usr/bin/env python

##findConserved.py
##Rishi Aryal, 2021('reseearyal[at]gmail[dot]com') 

from __future__ import print_function
import sys
import getopt
import re
import textwrap




def tab_generator(fastafile):
    '''generates name and sequence for each fasta sequence for further processing'''
    seq = []
    title=None
    with open (fastafile,'r') as f:
        for line in f:
           if line.startswith('>'):
              if seq:
                 yield title, ''.join(seq)
              title, seq = line.strip().split(' ')[0][1:], []
           else:
              seq.append(line.strip())
        yield   title, ''.join(seq)   




def kb(seq, ksize):
    '''generates list of kmers of given size for each sequence'''
    kmer=[]
    nmer=len(seq)-ksize+1
    for i in range(nmer):
        kmer.append(seq[i:i+ksize])
    return kmer




def main(argv):
    helpdoc='''
    Finds conserved regions on multiple sequnces in a fasta file

    Usage: findConserved.py [-h,-f <fastafile>,-k<kmer_size>, -c <coverage>, -q,-l]

         -h/--help:      this useful help
         -f/--fastafile: Fasta file containing input sequences
         -l/--length:    minimum lenght of conserved sequence to be produced. 
                         Should at least be 2*kemr_size; less than that will be 
                         automatically adjusted to 2*kemr_size.[default=100] 
         -c/--coverage:  mimimum percentage of input sequence the 
                         conserved region should be present [default=100]
         -p/--position:  if location of the conserved region on each 
                         genome is desired. The postitions will be added on the name of 
                         respective sequence. [default=False]
         -k/--kmer_size: kmer size to be used for genome comparision; 
                         smaller value tend ot produce more false positives; 
                         larger value require more disc space but reduces 
                         false positives. More useful with the '--quick' option. [default=21]
         -q/--quick:     produce quick and dirty output with mapping the 
                         conserved region in each genome. Runs quick with less computational 
                         resources but may produce some false positive. good for quick output 
                         on big genomes. [default=False]

     '''

    fastafile=''
    ksize=21
    length=100
    position=False
    percent=1
    quick=False

    if len (sys.argv[1:])==0:
        print ('Usage: findConserved.py [-h,-f <fastafile>,-k<kmer_size>, -c <coverage>, -q,-l]')
        sys.exit(2)
   
    try:
        opts,args=getopt.getopt(argv,"hk:f:l:pqc:",["help","kmer_size=","fastafile=", "lenght=", "position",\
        "quick", "coverage"])
    except getopt.GetoptError:
       print ('Usage: findConserved.py [-h,-f <fastafile>,-k<kmer_size>, -c <coverage>, -q,-l]')
       sys.exit(2)

    for opt, arg in opts:
       if opt in ("-h", "--help"):
          print (helpdoc)
          sys.exit()
       elif opt in ("-f", "--fastafile"):
          fastafile = arg
       elif opt in ("-k", "--kmer_size"):
          ksize = int(arg)
       elif opt in ("-l", "--length"):
          length=int(arg)
       elif opt in ("-p", "--position"):
          position=True
       elif opt in ("-c", "--coverage"):
          percent=int(arg)/100
       elif opt in ("-q", "--quick"):
          quick=True

    seq_dict={}
    seq_list=[]
    k_dict={}
    for title, seq in tab_generator(fastafile):
        seq_dict[title]=seq
        seq_list.append(title)
        k_dict[title]=kb(seq,ksize)
    
    ## get common kmers in all sequences
    allset=[]
        
    for key, value in k_dict.items():
        if allset:
            allset=list(set.intersection(set(allset),set(value)))
        else:
            allset=value
    
    
    ## find the range of contiguous kmers and locate the position on the first sequence
    k_list1=k_dict[seq_list[0]]
    k_seq1=seq_dict[seq_list[0]]

    k_ind1=[]
    for i,v in enumerate(k_list1): 
        if v in allset:
            k_ind1.append(i)
            
                    
    k_range1=sum((list(t) for t in zip(k_ind1, k_ind1[1:]) if t[0]+1 != t[1]), [])
    k_iter1=iter(k_ind1[0:1]+k_range1+k_ind1[-1:])
    k_range_l=[[n, next(k_iter1)] for n in k_iter1]

   
    if quick:
         ## quick and dirty output, may produce some false positive, if the genome to compare is highly repetitive 
         ## uses less computational time, useful for big genomes
        seq_num=1
        for value in k_range_l:
            if value[1]-value[0]>length:
                print('>conserved_locus-'+str(seq_num)+'\n'+('\n'.join(textwrap.wrap(k_seq1[value[0]:value[1]+ksize],70)).strip())+'\n')
                seq_num+=1

    else:
         ## use regular expression to find the candidate sequence on each genome
        seq_num=1
        con_seq=[]
        for i in k_range_l:
            if i[1] + ksize - i[0] >= length:
                con_seq.append(k_seq1[i[0]:i[1]+ksize]) # make a list of conserved sequences on the first fasta sequence
           
        ## search the conserved sequence on all the input sequences on fasta file using regular expression. 
        
        for seq in con_seq:
            pos_list=[]
            for genome in seq_list: # for each genome provided in the fasta file, genome here indicated genome name to search in the name:seq dictionalry           
                m=re.compile(seq).search(seq_dict[genome])
                if not m:
                    continue        
                else:
                    pos=str(m.span()[0])+"-"+ str(m.span()[1])
                    pos_list.append(genome+':'+pos)
            if len(pos_list)/len(seq_list)>=percent:
                if position:
                    print('>conserved_locus-'+str(seq_num)+' '+' '.join(pos_list))
                else:
                    print('>conserved_locus-'+str(seq_num))
                print ('\n'.join(textwrap.wrap(seq, 70)).strip() +'\n')
    
                seq_num+=1





      

if __name__ == "__main__":
    main(sys.argv[1:])            
