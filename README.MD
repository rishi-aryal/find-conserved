This program finds the conserved regions on multiple sequnces in a fasta file

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

example run:
./findConserved.py -f example-data/example-input-sequences.fasta > example-data/example-output-conserved-100bp.fasta