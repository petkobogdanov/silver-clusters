import csv, os

# throw-away script to convert our raw_data file to the
# csv file which is our format
if __name__ == "__main__":
    infile = "integrated-intensity.csv"
    outfile = "integrated-intensity2.csv"
    try:
        os.remove(outfile)
    except OSError:
        pass
    print "converting data file to new format (simple csv)"
    with open(infile, "rb") as inf, open(outfile, "w") as outf:
        csvReader = csv.reader(inf)
        # outf.write(",".join(csvReader.next()[1:6]) + "\n")
        csvReader.next()
        for line in csvReader:
            outf.write(line[2]+","+line[5]+"\n")
                
                    
