#author yasemin sahin aydin
#!/usr/bin/python

from merciscript import *


runmerci('darkness.fa', 'brightness.fa')
motifs1  = parse_occurence_file()

runmerci('darkness.fa', 'brightness.fa')
motifs2  = parse_occurence_file()

motifs3 = merge_motifs(motifs1, motifs2)

write_occurences_to_file(motifs3)
    
read_sequence_data('brightness.csv')
read_sequence_data('darkness.csv')
feature_vector_generator(motifs3) 
save_feature_vectors('combinedfeatures.csv') 



