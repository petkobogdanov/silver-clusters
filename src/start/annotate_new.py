#!/usr/bin/python

import sys
from merciscript import *

suff=sys.argv[1]
sequences=sys.argv[2]
#custom_motifs = create_pairs(['A', 'C', 'G', 'T'], 2)
   
runmerci('brightness' + suff + '.fa','darkness' + suff + '.fa')
#motifs1 = parseresult()
motifs1 = parse_output_file("+")

runmerci('darkness' + suff + '.fa','brightness' + suff + '.fa')
#motifs2 = parseresult()
motifs2 = parse_output_file("-")

motifs3 = merge_motifs(motifs1, motifs2)

print "num motifs " + str(len(motifs3))

cntb = read_sequence_data(sequences)
feature_vector_generator(motifs3,-1) 
save_feature_vectors('merci_combinedfeatures'+suff+'-newtest.csv', 'merci_combinedfeatures'+suff+'-newtest-occ.csv')     



