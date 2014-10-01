#!/usr/bin/python

import sys
from merciscript import *

if len(sys.argv) == 1:
  suff=""
else:
  suff=sys.argv[1]
#custom_motifs = create_pairs(['A', 'C', 'G', 'T'], 2)
   
runmerci('brightness' + suff + '.fa','darkness' + suff + '.fa')
#motifs1 = parseresult()
motifs1 = parse_output_file("+")

runmerci('darkness' + suff + '.fa','brightness' + suff + '.fa')
#motifs2 = parseresult()
motifs2 = parse_output_file("-")

motifs3 = merge_motifs(motifs1, motifs2)

print "num motifs " + str(len(motifs3))

cntb = read_sequence_data('brightness'+suff+'.csv')
read_sequence_data('darkness'+suff+'.csv')
feature_vector_generator(motifs3, cntb) 
save_feature_vectors('merci_combinedfeatures'+suff+'.csv', 'merci_combinedfeatures'+suff+'-occ.csv') 


#delete_current_data()
#read_sequence_data('brightness.csv')
#feature_vector_generator(motifs1,cntb) 
#save_feature_vectors('merci_brightfeatures.csv') 

#delete_current_data()
#read_sequence_data('darkness.csv')
#feature_vector_generator(motifs2,0) 
#save_feature_vectors('merci_darkfeatures.csv')       



