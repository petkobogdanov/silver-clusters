import convertData
import merciscript
import util, generator
import sys, env, os, argparse

ROOTDIR= env._env['MAIN']
DATA= env._env['DATA_PATH']
INPUTDATA = env._env['INPUT_DATA']
# SUFFS = ["-ii"]   #, "-pi", "-gfp"] 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-th', '--thresholdhigh', type=float, help="If above this threshold, sequence in positive class")
    parser.add_argument('-tl', '--thresholdlow', type=float, help="If below this threshold, sequence in negative class")
    parser.add_argument('-g', '--generate', help="Generate new features", action="store_true")
    parser.add_argument('-v', "--verbose", help="Output debug info to stdout", dest='verbose', action="store_true")
    parser.add_argument('filename', help="Name of file to run on")
    args = parser.parse_args()
    if args.thresholdhigh and args.thresholdlow:
        env.USE_THRESHOLD = True
        # print "Using threshold: ", "H",args.thresholdhigh, "L", args.thresholdlow
        convertData.makeClassFromThreshold(args.thresholdhigh, args.thresholdlow, args.filename)
    else:
        env.USE_THRESHOLD = False
    if args.generate:
        env.GENERATE = True    
    if args.verbose:
        env.VERBOSE = True
        env.print_env()
    sys.path.append(env._env['LIBSVM_PYTHON_PATH'])
    import libsvm_wrapper, svmutil
    if args.verbose:
        util.printSeparator(True)
    # convert input files to fasta format for MERCI
    # files = [f for f in os.listdir(INPUTDATA) if os.path.isfile(os.path.join(INPUTDATA, f))]
    if env.USE_THRESHOLD:
        files = [args.filename+".classes"]
    else:
        files = [args.filename]
    numFiles = len(files)
    # print files
    # exit(0)
    if args.verbose:
        print "Running on " + str(numFiles) + " files:\n" + str(files)
    if not os.path.exists(env._env['TMP_FILES_PATH']):
        os.makedirs(env._env['TMP_FILES_PATH'])
    if not os.path.exists(env._env['TRAINING_DATA_PATH']):
        os.makedirs(env._env['TRAINING_DATA_PATH'])
    for f in files:
        # generate two csvs and two fastas, one of each for pos. neg.
        convertData.generateFiles(f) 
    # Run MERCI; save feature vectors
    for file in files:
        basename = os.path.splitext(file)[0]
        if env.USE_THRESHOLD:
            basename = os.path.splitext(basename)[0]
        posFastaFile = basename+".POS.fa"
        negFastaFile = basename+".NEG.fa"
        posCSVFile = basename+".POS.csv"
        negCSVFile = basename+".NEG.csv"
        merciscript.initGlobals()
        if args.verbose:
            util.printSeparator()
        # puts output in merciresult file
        merciscript.runmerci(posFastaFile,negFastaFile)
        # store results in list in memory
        motifs1 = merciscript.parse_output_file("+")
        # puts output in merciresult file
        merciscript.runmerci(negFastaFile,posFastaFile)
        # store results in memory
        motifs2 = merciscript.parse_output_file("-")
        motifs3 = merciscript.merge_motifs(motifs1, motifs2)
        # exit(0)
        if args.verbose:
            print "num motifs " + str(len(motifs3))
        # stores data in global variable 'data' in merciscript
        cntb = merciscript.read_sequence_data(posCSVFile) # returns count
        merciscript.read_sequence_data(negCSVFile)
        merciscript.feature_vector_generator(motifs3, cntb)
        # exit(0)
        merciscript.save_feature_vectors(basename+".combinedfeatures.csv", basename+".combinedfeatures-occ.csv")
        if args.verbose:
            util.printSeparator(True)
    # Run LibSVM; save accuracy report
    TRAINPATH=env._env['TRAINING_DATA_PATH']
    training_files = [os.path.join(TRAINPATH, files[0])]
    for tfile in training_files:
        # Run LibSVM on training data.
        basename_tfile = os.path.splitext(tfile)[0]
        if env.USE_THRESHOLD:   
            basename_tfile = os.path.splitext(basename)[0]
        basename_tfile = basename + ".combinedfeatures.csv"
        if args.verbose:
            util.printSeparator()
            print "Training from feature file:", basename_tfile
        fulltfile = os.path.join(env._env['TRAINING_DATA_PATH'], basename_tfile)
        x, y = libsvm_wrapper.readData(fulltfile)
        m = svmutil.svm_train(y, x, "-v 5")
        if args.verbose:
            util.printSeparator(True)
        # Generate new features if -g is set 
        if args.generate:
            print "Generating NEW features"
            if env.USE_THRESHOLD:
                intensityFn = os.path.join(env._env['INPUT_DATA'], basename+".csv")
                # tfileFn = os.path.join(env._env['TRAINING_DATA_PATH'], basename)
                generator.doAll(basename_tfile, intensityFn, 1000, 2, 3, True, threshList=[args.thresholdhigh, args.thresholdlow])
            else:
                intensityFn = os.path.join(env._env['INPUT_DATA'], basename+".csv")
                generator.doAll(basename_tfile, intensityFn, 1000, 2, 3, False)
    util.deleteFiles(env._env['TMP_FILES_PATH'])
