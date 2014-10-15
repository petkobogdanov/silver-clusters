# This file is not part of the original libSVM package
from svmutil import *
import csv
import sys

def readData(fileName):
    "Read CSV file into memory for libSVM "
    x = []
    y = []
    with open(fileName, "r") as csvFile:
        reader = csv.reader(csvFile)
        reader.next()
        for line in reader:
            isPositive = True if line[-1] == "bright" else False
            line = map(float, line[1:-2])
            x.append(line)
            if isPositive is True:
                y.append(1)
            else:
                y.append(0)
    # print x
    # print y
    # print "length of x is: " + str(len(x))
    # print "length of y is: " + str(len(y))
    return x, y

if __name__ == "__main__": 
    # print sys.argv[1]
    featureFile = sys.argv[1]
    x, y = readData(featureFile)
    m = svm_train(y, x, '-v 10')
