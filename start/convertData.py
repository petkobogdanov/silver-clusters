import csv
import os
import env

def makeClassFromThreshold(high, low, fileName):
    fullName = os.path.join(env._env['INPUT_DATA'], fileName)
    outFile = os.path.join(env._env['INPUT_DATA'], fileName+".classes")
    # print outFile
    with open(fullName, "r") as inf, open(outFile, "w") as outf:
        csvReader = csv.reader(inf)
        for line in csvReader:
            # print line[0], line[1], float(line[1])
            if float(line[1]) < low:
                clas = 0
            elif float(line[1]) > high: 
                clas = 1
            else:
                clas = -1
            # print clas
            if clas is not -1:
                outf.write(line[0]+","+str(clas)+"\n")
    
def generateFiles(fileName):
    # print "filName",fileName
    "Convert input CSV file to fasta for MERCI"
    fullName = os.path.join(env._env['INPUT_DATA'], fileName)
    # print fileName
    baseName = os.path.splitext(fileName)[0]
    if env.USE_THRESHOLD:
        baseName = os.path.splitext(baseName)[0]
    # print baseName
    outPosfasta = os.path.join(env._env['TMP_FILES_PATH'], baseName+".POS.fa")
    outNegfasta = os.path.join(env._env['TMP_FILES_PATH'], baseName+".NEG.fa")
    outPosCSV = os.path.join(env._env['TMP_FILES_PATH'], baseName+".POS.csv")
    outNegCSV = os.path.join(env._env['TMP_FILES_PATH'], baseName+".NEG.csv")
    try:
        os.remove(outPosfasta)
        os.remove(outNegfasta)
    except OSError:
        pass
    with open(fullName, "r") as inCSV, open(outPosfasta, "w") as outPosfa, open(outNegfasta, "w") as outNegfa, \
         open(outPosCSV, "w") as outPosCsv, open(outNegCSV, "w") as outNegCsv:
        csvReader = csv.reader(inCSV)
        outPosCsv.write("Sequence\n")
        outNegCsv.write("Sequence\n")
        for line in csvReader:
            try:
                if int(line[1]) == 1:
                    outPosfa.write(">\n"+line[0]+"\n")
                    outPosCsv.write(line[0]+"\n")
                if int(line[1]) == 0:
                    outNegfa.write(">\n"+line[0]+"\n")
                    outNegCsv.write(line[0]+"\n")
            except ValueError:
                continue
    
