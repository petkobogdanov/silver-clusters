import random, sys, re, math, os, env, csv

# Globals
BASE_STICKY = {"A" : 0.5, "T" : 0.0, "C" : 1.0, "T" : 0.0}
RAND = random.seed()

def getNumLines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def ceilingKey(d, key):
    ceilK = []
    if key in d:
        return key
    ceilK = [k for k in d if k >= key]
    if ceilK:
        return min(ceilK)
    return 0.0      # TODO is returning 0 here aproblem??
    # return min(k for k in d if k >= key)            

def floorKey(d, key):
    floorK = []
    if key in d:
        return key
    floorK = [k for k in d if k <= key]
    if floorK:    
        return max(floorK)
    return 0.0
    
# Generator 
class Generator:
    def __init__(self, fn):
        self.fileName = fn
        self.n = getNumLines(self.fileName)
        self.motifs = []
        self.seqs = [""]*self.n
        self.distBnum = []
        self.distDnum = []
        self.distSticky = []
        self.distMotif = []
        self.seqValHash = {}

    def getRegions(self, ind, length, numreg):
        chunk = length / numreg
        rem = length % numreg
        res = []
        for pos in range(numreg):
            if ind >= pos*chunk and ind < (pos+1)*chunk + rem:
                res.append(pos)
        return res

    def computeAvgStickyness(self, strnew):
        cnt = 0
        sticky = 0.0
        for i in range(len(strnew)):
            if strnew[i] == "_":
                cnt = cnt +1
                sticky += BASE_STICKY[strnew[i:i+1]]
        if cnt > 0:
            sticky /= cnt
        return sticky
        
    def fillGapsBasedOnStickyness(self, seq):
        try:
            ind = seq.index('_')
        except ValueError:
            return seq
        # print "ind: ",ind
        regs = self.getRegions(ind, len(seq), self.nreg)
        # print "regs: ", regs
        base2score ={}
        totScore = 0.0
        for base in BASE_STICKY:
            base2score[base] = 0.0
        for reg in regs:
            s = reg * (len(seq)/self.nreg)
            e = (reg+1)*(len(seq)/self.nreg) + len(seq)%self.nreg
            str = seq[s:e]
            for base in base2score:
                replace = ind - s
                assert(str[replace] == "_")
                strnew = str[0:replace] + base + str[replace+1:len(str)]
                avgSticky = self.computeAvgStickyness(strnew)
                sc = 0
                if floorKey(self.distSticky[reg], avgSticky) is None:
                    sc = firstKey(self.distSticky[reg])
                else:
                    if ceilingKey(self.distSticky[reg], avgSticky) is None:
                        sc = self.distSticky[reg].keys()[-1]
                    else:
                        if floorKey(self.distSticky[reg], avgSticky - 0.000000001) is None:
                            sub = 0
                        else:
                            sub = floorKey(self.distSticky[reg], avgSticky - 000000001)
                        sc = ceilingKey(self.distSticky[reg], avgSticky) - sub
                base2score[base] = base2score[base] + sc
                totScore += sc
        prob2base = {}
        acc = 0.0
        for base in base2score:
            prob2base[base2score[base]/totScore + acc] = base
            acc += base2score[base] / totScore
        base = prob2base[ceilingKey(prob2base, random.random())]
        # print "base: ",base    
        seq = seq[0:ind] + base + seq[ind+1:len(seq)]
        return seq
    
    def generate(self, numNew):
        cnt = 0
        atMostNGaps = 1
        generated = []
        while cnt < numNew:
            seq = "__________"
            # print "calling placeMotifs", seq, atMostNGaps
            seq = self.placeMotifs(seq, atMostNGaps)
            # print "before :",seq
            seq = self.fillGapsBasedOnStickyness(seq)
            # print "after :", seq, "\n"
            generated.append(seq)
            cnt = cnt + 1
        return generated

    def getPossiblePositions(self, m, seq, region):
        resArrayList = []
        half = len(m.regex) / 2.0
        # print seq, m.regex
        # print "seq before (getpossiblepositions)", seq
        for st in range(len(seq) - len(m.regex) +1):
            if math.floor((st + half)*self.nreg*1.0/len(seq)) == region and self.fits(m, seq, st):
                resArrayList.append(st)
        # print "seq after (getpossiblepositions)", seq
        return resArrayList

    def fits(self, m, seq, position):
        for i in range(len(m.regex)):
            # print "seq before (fits)", seq
            if m.regex[i] is not "_" and seq[i+position] is not "_" and m.regex[i] is not seq[i+position]:
                return False
            # print "seq after (fits)", seq
        return True

    def insert(self, m, seq, pos):
        # print "in insert, seq: ", seq    
        assert(self.fits(m, seq, pos))
        res = ""
        for i in range(len(seq)):
            if i < pos or i >= pos + len(m.regex):
                res += seq[i]
            elif seq[i] == "_":
                res += m.regex[i-pos]
            elif m.regex[i-pos] == "_":
                res += seq[i]
            else:
                assert(seq[i] == m.regex[i-pos])
                res += seq[i]
        return res
        
    def placeMotifs(self, seq, atMostNGaps):
        pos = self.getShuffledIndices()
        numGaps = len(seq)
        # sample motifs in 1 then 2 and 3 positions
        while(numGaps > atMostNGaps):
            for j in range(self.nreg):
                i = pos[j]
                m = self.distMotif[i][ceilingKey(self.distMotif[i], random.random())]
                # print "calling getPossiblePositions", m, seq, i
                possible_pos = self.getPossiblePositions(m, seq, i)
                if len(possible_pos) > 0:
                    position = possible_pos[random.randint(0, len(possible_pos)-1)]
                    seq = self.insert(m, seq, position)
            numGaps = 0
            for i in range(len(seq)):
                if seq[i] == "_":
                    numGaps = numGaps + 1
        return seq 
                    
    def getShuffledIndices(self):
        pos = [0] * self.nreg
        for i in range(self.nreg):
            pos[i] = i
        for i in range(10):
            a = random.randint(0, self.nreg-1)
            b = random.randint(0, self.nreg-1)
            sub = pos[a]
            pos[a] = pos[b]
            pos[b] = sub
        return pos
        
    def makeArrays(self, nreg):
        self.nreg = nreg    
        self.bcnt = [[0 for x in range(self.nreg)] for x in range(self.n)]
        self.dcnt = [[0 for x in range(self.nreg)] for x in range(self.n)]
        self.bsize = [[0 for x in range(self.nreg)] for x in range(self.n)]
        self.bsticky = [[0 for x in range(self.nreg)] for x in range(self.n)]
        self.classes = [""] * self.n
        self.intensity = [0] * self.n 

    def readIntensityFile(self, fn):
        # Read intensity file; Values could also be binary (0/1)
        with open(fn, "r") as csvFile:
            csvreader = csv.reader(csvFile)
            # build Hash of Sequence, Value
            for line in csvreader:
                # line[0] is sequence, line[1] is intensity value
                self.seqValHash[line[0]] = line[1]
        # print self.seqValHash
                
    def readTrainingFile(self):
        # get relevant info. from training file
        with open(self.fileName, "r") as trainingFile:
            line = trainingFile.readline()
            sline = line.split(",")
            colFeat = 1
            colSeq = 0
            # colII = 4
            colBPos = colDPos = colBsize = colBSticky = colClass = -1
            for i in range(0, len(sline)):
                # print sline[i]
                if sline[i].lower() == "brightpos0": colBPos = i
                if sline[i].lower() == "darkpos0": colDPos = i
                if sline[i].lower() == "basesize0": colBsize = i
                if sline[i].lower() == "basesticky0": colBSticky = i
                if sline[i].lower() == "class": colClass = i
            # print colBPos, colDPos, colBsize, colBSticky, colClass
            nreg = colDPos - colBPos
            self.makeArrays(nreg)
            for i in range(colFeat, colBPos):
                self.motifs.append(Motif(sline[i], nreg))
                if "_" in sline[i]:
                    self.motifs.append(Motif(sline[i].replace("_", ""), nreg))
            nr = 0
            # print bcnt
            for line in trainingFile:
                # print sline
                sline = line.split(",")
                self.seqs[nr] = sline[colSeq]
                self.classes[nr] = sline[colClass][:-1] # truncate newline
                # print "seq",self.seqs[nr]
                # exit(0)
                # self.intensity[nr] = float(sline[colII])
                self.intensity[nr] = float(self.seqValHash[self.seqs[nr]])
                for i in range(nreg):
                    self.bcnt[nr][i] = int(sline[colBPos+i])
                    self.dcnt[nr][i] = int(sline[colDPos+i])
                    self.bsize[nr][i] = float(sline[colBsize+i])
                    self.bsticky[nr][i] = float(sline[colBSticky+i])
                nr = nr + 1
    def getTotalIntensity(self, clas):
        sumII = 0.0
        # print self.classes
        for i in range(self.n):
            if self.classes[i] == clas:
                sumII = sumII + self.intensity[i]
        # print sumII
        return sumII
    def computeDist(self, clas, sumII, pos, vals):
        m = {}
        for i in range(self.n):
            if self.classes[i] == clas:
                if vals[i][pos] not in m:
                    m[vals[i][pos]] = 0.0
                m[vals[i][pos]] = m[vals[i][pos]] + self.intensity[i] / sumII
        prob2num = {}
        acc = 0.0
        for num in m:
            prob2num[m[num]+acc] = num
            acc += m[num]
        # print prob2num, "\n"
        return prob2num

    def matchInReg(self, m, seq, reg, numreg):
        assert "_" not in seq
        regex = m.regex.replace("_", "[ACTG]{1}")
        pattern = re.compile(regex)
        for match in pattern.finditer(seq):
            # print match.start(), match.group(), seq, regex, reg, numreg
            middle = (match.start() + len(match.group())) / 2.0
            if math.floor(middle*numreg*1.0/len(seq)) == reg:
                # print "yes"
                return True
        # print "no"
        return False
            
    def computeDistMotif(self, clas):
        res = []
        totalII = [0.0] * self.nreg
        for reg in range(self.nreg):
            for m in self.motifs:
                for i in range(self.n):
                    if self.classes[i] == clas:
                        if self.matchInReg(m, self.seqs[i], reg, self.nreg):
                            m.prob[reg] += self.intensity[i]
                            totalII[reg] += self.intensity[i]
            if totalII[reg] > 0:
                for m in self.motifs:
                    m.prob[reg] /= totalII[reg]
            prob2mot = {}
            acc = 0.0
            for m in self.motifs:
                prob2mot[m.prob[reg]+acc] = m
                acc += m.prob[reg]
            res.append(prob2mot)
        return res            
            
    def computeDists(self,clas):
        # compute total integrated intensity in bright
        sumII = self.getTotalIntensity(clas)
        for pos in range(self.nreg):
            self.distBnum.append(self.computeDist(clas, sumII, pos, self.bcnt))
            self.distDnum.append(self.computeDist(clas, sumII, pos, self.dcnt))
            self.distSticky.append(self.computeDist(clas, sumII, pos, self.bsticky))
        self.distMotif = self.computeDistMotif(clas)

    def editDist(self, seq1, seq2):
        d = 0
        # print seq1, seq2
        # assert(len(seq1) == len(seq2))
        if len(seq1) is not len(seq2):
            return
        for i in range(len(seq1)):
            if seq1[i] is not "_" and seq2[i] is not "_" and seq1[i] is not seq2[i]:
                d = d + 1
        return d
        
    def isClosestKnownWithin(self, minD, maxD, s):
        closestD = len(s)
        d = 0
        assert(len(s) == len(self.seqs[0]))
        for i in range(len(self.seqs)):
            d = self.editDist(self.seqs[i], s)
            if d < closestD:
                closestD = d
        return minD <= closestD  and maxD >= closestD
        
    def writeSequences(self, fileName, sequences, MinD, MaxD):
        ident = 0
        with open(fileName, "w") as outf:
            nseq = set()
            for s in sequences:
                # print "sequence s :",s
                nseq.add(s)
            outf.write("Name,Sequence\n")
            # print nseq    
            for s in nseq:
                # if self.isClosestKnownWithin(MinD, MaxD, s):
                outf.write("10basegen" + str(ident) + "," + s + "\n")
                ident = ident + 1
            print "wrote", ident, "new sequences"
        
    def printALL(self):
        print "printing ALL: \n", \
            "Motifs: ", str(self.motifs), "\n\n", \
            "Seqs: ", str(self.seqs), "\n\n", \
            "bcnt: ", str(self.bcnt), "\n\n", \
            "dcnt: ", str(self.dcnt), "\n\n"
        for m in self.motifs:
            print m.regex
# Motif
class Motif:
    def __init__(self, nm, nreg):    
        self.name = nm
        self.isPos = nm.endswith("+")
        self.regex = nm[0:len(nm)-1]
        self.prob = [0]*nreg
    def printAll(self):
        print "printing ALL: " + self.name, self.isPos, self.regex, self.prob

def doAll(trainingFn, intensityFn, numNew, Tmin, Tmax, useThresholdp, threshList=None):
    if threshList:
        print "threshlist",threshList
    print "intensityFn",intensityFn
    # exit(0)
    outFileName = os.path.join(env._env['NEW_FEATURES_PATH'], trainingFn)
    outFileName = os.path.splitext(outFileName)[0] + "-new.csv"
    print "outFileName",outFileName
    # exit(0)
    infileName = os.path.join(env._env['TRAINING_DATA_PATH'], trainingFn)
    print "infileName",infileName
    # exit(0)
    generator = Generator(infileName) # instantiate object; get num lines in file
    # build in memory hash for Sequence,Value; 
    generator.readIntensityFile(intensityFn) 
    # exit(0)
    generator.readTrainingFile()
    # exit(0)
    generator.computeDists("bright")
    # exit(0)
    newseqs = generator.generate(numNew)
    # print outFileName
    generator.writeSequences(outFileName, newseqs, Tmin, Tmax)



    
if __name__ == "__main__":
    fileName = sys.argv[1]
    numNew = int(sys.argv[2])
    Tmin = int(sys.argv[3])
    Tmax = int(sys.argv[4])    
    generator = Generator(fileName)
    generator.readSequences()
    generator.computeDists("bright")
    newseqs = generator.generate(numNew)
    generator.writeSequences(os.path.splitext(fileName)[0] + "-new.csv", newseqs, Tmin, Tmax)
