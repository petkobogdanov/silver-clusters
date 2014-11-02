import env, os
import csv
import shlex
import math
from subprocess import Popen, PIPE
from random import shuffle
import re

global occur, data, output, gl

data = []
occur = {}
output = None
gap_length = 0
S_LEN=10
S_REG=3
posfp_bright = {}
posfp_dark = {}
pos_size = {}
pos_sticky = {}
posit = {}
add_position_features=True


def initGlobals():
    global data
    data = []
    global occur
    occur = {}
    global output
    output = None

# converts csv file into fasta format
def csvToFasta(infile, outfile):
    print "-- Extracting data from: " + infile
   
    infile = os.path.abspath(os.path.join(env._env['DATA_PATH'] , infile) )
    outfile = os.path.abspath(os.path.join(env._env['DATA_PATH'] , outfile) )
   
    outf = open(outfile, 'w')
    
    with open(infile, 'r') as csvfile:
        seqreader = csv.reader(csvfile, delimiter=',', quotechar="\"")
        seqreader.next()
        for row in seqreader:
            outf.write('>' + row[0] + '\n' + row[1] + '\n')
    
    outf.close()
    
    print "-- Fatsa file is created: " + outfile
    
# run merci to collect motif data
#def runmerci(positivefile, negativefile, topK='ALL', length=5, fp=10, fn=10, g=1, gl=1):
def runmerci(positivefile, negativefile, topK='ALL', length=10, fp=10, fn=10, g=1, gl=1):

    global output  
    global gap_length
    gap_length = gl
    from sys import platform as _platform
    slash = "/"
    if _platform == "win32":
        slash = "\\"
    cmd = "perl " + env._env['MERCI_EXE_PATH'] + " "  
    cmd += "-p " +env._env['TMP_FILES_PATH'] + slash +"{} ".format(positivefile)
    cmd += "-n " +env._env['TMP_FILES_PATH'] + slash +"{} ".format(negativefile)
    cmd += "-k {} -l {} ".format(topK, length)
    cmd += "-fp {} -fn {} ".format(fp, fn)
    cmd += "-g {} -gl {} ".format(g, gl)
    cmd += "-c " + env._env['MERCI_CLASSIFICATION_PATH'] + " " 
    cmd += "-o " + env._env['MERCI_OUTPUT_FILE']

    if env.VERBOSE:
        print "-- Running command: "    
        print cmd
    args = cmd
    if _platform != "win32":
        args = shlex.split(cmd)
    # print "\r\n\r\n",args
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    output = out
    return out
   
class Motif:
    def __init__(self, name='', regex='', width=0, sites=0, llr=0, evalue=0, suffix = ""):
        self.name = name
        self.regex = regex
        self.width = width
        self.sites = sites
        self.llr = llr
        self.evalue = evalue
        self.suffix = ""
        self.num_of_positive_occurence = 0
        self.num_of_negative_occurence = 0
    def __str__(self):
        return self.name + " : " + self.regex + " :w " + str(self.width) + " :s " + str(self.sites) + " :e " + str(self.evalue) + " :llr " + str(self.llr) 

    def to_string(self):
        return self.name + " --> " + self.regex + " - positive occ: " + self.num_of_positive_occurence + " negative occ: " + self.num_of_negative_occurence

def match_motif_regex(string):
    return re.match(r"\s*MOTIF:\s+([AGTC\[\] p]+)\s*",string,re.IGNORECASE)
    
def match_num_of_positive_occurence(string):
    return re.match(r"(\d+) positive sequences contain the motif\s*", string, re.IGNORECASE)    

def match_num_of_negative_occurence(string):
    return re.match(r"(\d+) negative sequences contain the motif\s*", string, re.IGNORECASE)    

def match_motif_regex_header(string):
    return re.match(r"\s*(Motifs:)\s*",string,re.IGNORECASE)

def match_top_motif_regex(string):
    return re.match(r"\s*([AGTC\[\] p]+)\s*",string,re.IGNORECASE)
 
def match_gl_param_regex(string):
    return re.match(r"\s+-gl\s+(?P<param_gl>\d+)\s*", string, re.IGNORECASE)
# extract motif data from merci results    
def parseresult(result=None):
    global output
    if result is None:
        if output is None:
            raise Exception("You should run merci to get output data, then parse result")
        else:
            result = output
            
    lines = result.strip().split('\n') 
    motifs = []
    motif = None
    for line in lines:
        motif_regex_match = match_motif_regex(line)
        if not motif_regex_match:
            continue
        motif = Motif()
        motif.regex = motif_regex_match.group(1).replace(" ", "")
        motifs.append(motif)
        motif = None
    return motifs
    
# read output data from merci
def parse_occurence_file():
    global gap_length
    resultfile =  env._env['MERCI_OCCURENCE_OUTPUT_FILE']
    motifs = []
    motif = None
    state = 0;
 
    with open(resultfile, 'r') as f:
        for line in f:
            if state == 0:
                motif_regex_match = match_motif_regex(line)
                if not motif_regex_match:
                    continue
                motif = Motif()
                motif.regex = ''
                motif.name = motif_regex_match.group(1).replace(" ", "")
                gaps = motif_regex_match.group(1).split()

                for ch in gaps:
                    if ch == 'gap':
                        motif.regex += "([AGTC]){0," + str(gap_length) + "}"
                    else:
                        motif.regex += ch
                
                motifs.append(motif)
                state = 1
            if state == 1:
                num_of_positive_match = match_num_of_positive_occurence(line)
                if not num_of_positive_match:
                    continue
                motif.num_of_positive_occurence = num_of_positive_match.group(1)
                state = 2
            elif state == 2:
                num_of_negative_match = match_num_of_negative_occurence(line)
                if not num_of_negative_match:
                    continue
                motif.num_of_negative_occurence = num_of_negative_match.group(1)
                state = 0
    return motifs 
    
# read output data from merci
def parse_output_file(suffix):
    if env.VERBOSE:
        print "DATA_PATH " + env._env['DATA_PATH']
    resultfile =  env._env['MERCI_OUTPUT_FILE']
    motifs = []
    motif = None
    state = 0;
    gl = 0;
    with open(resultfile, 'r') as f:
        for line in f:
            if state == 0:
                param_match = match_gl_param_regex(line)
                if not param_match:
                    continue
                gl = param_match.group('param_gl')
                state = 1
            if state == 1:
                motif_match = match_motif_regex_header(line)
                if not motif_match:
                    continue
                state = 2
            elif state == 2:
                motif_regex_match = match_top_motif_regex(line)
                if not motif_regex_match:
                    continue
                motif = Motif()
                motif.regex = ''
                motif.name = motif_regex_match.group(1).replace(" ", "")
                gaps = motif_regex_match.group(1).split()

                for ch in gaps:
                    if ch == 'gap':
                        motif.regex += "([AGTC]){0," + str(gl) + "}"
                    else:
                        motif.regex += ch
                motif.suffix = suffix 
                motifs.append(motif)
                motif = None
    return motifs  
    
def merge_motifs(motif1, motif2):
    return motif1 + motif2

def print_motifs(motifs):
    for m in motifs:
        print m.to_string()
        #print m.__str__()    
    
def match_motif_with_seq(motif_regex, sequence):
    motif_regex = ".*" + motif_regex + ".*"
    pattern = re.compile(motif_regex, re.IGNORECASE)
    return re.match(pattern, sequence)

def check_positions_in_seq_file(regex, sequences, isPos, mname):
    m = regex
    p = re.compile(m, re.IGNORECASE)
    # count the occurances
    for row in sequences[1:]:
        ps = ""
        for mch in p.finditer(row[0]):
            if ps == "": 
                ps = str(mch.start())
            else: 
                ps = ps + "_" + str(mch.start()) 
            middle = (mch.start() + mch.end())/2.0
            #print str(mch.start()) + "\t" + str(mch.end()) + "\t" + str(middle) + "\t" + str(int(math.floor(S_REG*middle/S_LEN)))
            idx = int(math.floor(S_REG*middle/S_LEN))
            if isPos: posfp_bright[row[0] + "_" +str(idx)] += 1
            else: posfp_dark[row[0] + "_" + str(idx)] +=1
        key = row[0] + "_" + mname
        if isPos: key = key + "+"
        else: key = key + "-"
        if ps != "": 
            posit[key] = ps

def check_motif_in_seq_file(motif_regex, sequences):

    for row in sequences[1:]:
        # print "row",row
        regex_match = match_motif_with_seq(motif_regex,row[0])
        if regex_match:
            row.append("1")
        else:
            row.append("0")
            
    return sequences 

# read sequence file into memory to generate feature vectors
# if you do not delete current data, it appends each file to the end 
def read_sequence_data(filename):
    cnt=0
    if env.VERBOSE:
        print "-- Reading sequence data from: " + filename
    filename = os.path.abspath(os.path.join(env._env['TMP_FILES_PATH'] , filename) )
    with open(filename, 'r') as csvfile:
        seqreader = csv.reader(csvfile, delimiter=',', quotechar="\"")
        header = seqreader.next()
        if not data:
            data.append(header)
        for row in seqreader:
            # print row    
            data.append(row)
            cnt = cnt + 1
    # print data
    return cnt

# print current data you have        
def print_current_data():
    for d in data:
        print ",".join(d)


def delete_current_data():
    global data
    data = []

# deletes feature vectors from data
def reset_current_data():
    for row in data:
        del row[2:]
           

def shuffle_current_data():
    global data
    header = data[0]
    rows = data[1:]
    shuffle(rows)
    data = []
    data.append(header)
    data = data + rows

# get some part of data    
def slice_current_data(begin=None, end=None):
    global data
    if begin is None:
        if end is None:
            return data
        else:
            data = data[:end]
    elif end is None:
        data = data[begin:]
    else:
        data = data[begin:end]

    
def feature_vector_generator(motifs,numbright, class1='bright', class2='dark'):
    if env.VERBOSE:
        print "-- Generating features"
    # print data
    found = {} 
    # add motif occurence features
    # print "data: ", data
    for m in motifs:
        if not m.name in found:
            # data[0] is feature headers
            data[0].append(m.name.replace("gap","_") + m.suffix) 
            check_motif_in_seq_file(m.regex,data)
            found[m.name]=m.name
    # print data
    # add position features
    found = {} 
    # init positions to 0
    # print "data: ",data
    for row in data[1:]:
        for i in range(S_REG):
            posfp_bright[row[0] + "_" + str(i)] = 0 
            posfp_dark[row[0] + "_" + str(i)] = 0
        # print posfp_bright
    for m in motifs:
        if not m.name in found:
           check_positions_in_seq_file(m.regex,data,m.suffix=="+",m.name)
           found[m.name]=m.name
    # print "posfp_bright: ", posfp_bright
    # print "posfp_dark: ", posfp_dark       
    if add_position_features:
        for i in range(len(data)):
            for pos in range(S_REG):
                if i==0: data[i].append(class1 + "Pos" + str(pos))
                else: data[i].append(str(posfp_bright[data[i][0]+ "_" + str(pos)])) 
            for pos in range(S_REG):
                if i==0: data[i].append(class2 + "Pos" + str(pos))
                else: data[i].append(str(posfp_dark[data[i][0]+ "_" + str(pos)]))

    # print data
    # add base size and stickiness features
    for row in data[1:]:
        for i in range(S_REG):
            pos_size[row[0] + "_" + str(i)] = 0 
            pos_sticky[row[0] + "_" + str(i)] = 0
    # print data
    for pos in range(S_REG):
        data[0].append("BaseSize" + str(pos))
    for pos in range(S_REG):
        data[0].append("BaseSticky" + str(pos))
    for i in range(len(data)):
        if i>0:
            compute_base_position_features(data[i][0]) 
            for pos in range(S_REG):
                data[i].append(str(pos_size[data[i][0]+ "_" + str(pos)])) 
            for pos in range(S_REG):
                data[i].append(str(pos_sticky[data[i][0]+ "_" + str(pos)])) 
    # print data
    # add order features
    for i in range(len(data)):
        if i > 0:
            rep = ""
            for m in motifs:
                if data[i][0] + "_" + m.name + m.suffix in posit:
                    if rep =="": 
                        rep = m.name + m.suffix + ":" + posit[data[i][0] + "_" + m.name + m.suffix]
                    else: 
                        rep = rep + ";" + m.name + m.suffix + ":" + posit[data[i][0] + "_" + m.name + m.suffix] 
            occur[data[i][0]] = rep
    for i in range(len(data)):
        if i==0: 
            data[i].append("class")
        elif numbright<0:
            if i==1: data[i].append(class1)
            else: data[i].append(class2)
        elif i<=numbright: 
            data[i].append(class1)
        else: 
            data[i].append(class2)
    # print "len of data[0] is " + str(len(data[0]))
    # print data

def compute_base_position_features(sequence):
    chunk, rem = divmod(len(sequence), S_REG)
    for pos in range(S_REG):
        for i in range(pos*chunk,(pos+1)*chunk+rem):
            sz=0
            st=0
            if(sequence[i]=='A'):
               sz=1.0
               st=0.5
            if(sequence[i]=='C'):
               sz=0.0
               st=1.0
            if(sequence[i]=='T'):
               sz=0.0
               st=0.0
            if(sequence[i]=='G'):
               sz=1.0
               st=1.0
            pos_size[sequence + "_" + str(pos)] += sz*1.0/(chunk+rem)
            pos_sticky[sequence + "_" + str(pos)] += st*1.0/(chunk+rem)
 
def save_feature_vectors(filename = "features.csv", filename_occur = "occur.csv"):

    # print data
    filename = os.path.abspath(os.path.join(env._env['TRAINING_DATA_PATH'] , filename) )
    if env.VERBOSE:
        print "saving feature vectors to " + filename
    outf = open(filename, 'w')
    filename_occur = os.path.abspath(os.path.join(env._env['TRAINING_DATA_PATH'] , filename_occur) )
    # outfocc = open(filename_occur, 'w')
    for row in data:
        # print row
        # exit(0)
        outf.write( ",".join(row) + '\n')     
        # if row[1]=="Sequence": outfocc.write(row[0] + "," + row[1] + "," + row[4] + ",SeqConfig" + '\n')
        # outfocc.write(row[0] + "," + row[1] + "," + row[4] + "," + occur[row[1]] + '\n')    
    outf.close()
    # outfocc.close()
    
def write_occurences_to_file(motifs, filename="occurences.csv"):
    filename = os.path.abspath(os.path.join(env._env['MERCI_DATA_PATH'] , filename) )
    outf = open(filename, 'w')
    for m in motifs:
        outf.write(m.name + ',' + str(m.num_of_positive_occurence) + ',' + str(m.num_of_negative_occurence) + '\n')    
    outf.close()


            
    


        

    
