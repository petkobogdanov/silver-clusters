import convertData
import merciscript
import util, generator
import sys, env, os, argparse

if __name__ == "__main__":
    # Execution Starts Here

    # Import libSVM libary
    sys.path.append(env._env['LIBSVM_PYTHON_PATH'])
    import libsvm_wrapper, svmutil
    from svm import *

    # Add Command Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-th', '--thresholdhigh', type=float, help="If above this threshold, sequence in positive class")
    parser.add_argument('-tl', '--thresholdlow', type=float, help="If below this threshold, sequence in negative class")
    parser.add_argument('-g', '--generate', help="Generate new features", action="store_true")
    parser.add_argument('-v', "--verbose", help="Output debug info to stdout", dest='verbose', action="store_true")
    parser.add_argument('filename', help="Name of file to run on (located input_data directory)")
    parser.add_argument("-c", "--crossValidation", help="Perform cross-validation testing then stop", action="store_true")
    parser.add_argument("-Mfp", "--MERCIfp", type=int, help="Argument -fp to MERCI; the minimal frequency (absolute number) \
    for the positive sequence (defaults to 10 if not specified)")
    parser.add_argument("-Mfn", "--MERCIfn", type=int, help="Argument -fn to MERCI; the maximal frequency (absolute number) \
    for the negative sequences (defaults to 10 if not specified)")
    parser.add_argument("-Mg", "--MERCIg", type=int, help="Argument -g to MERCI; maximal number of gaps \
    (defaults to 1 if not specified)")
    parser.add_argument("-Mgl", "--MERCIgl", type=int, help="Argument -gl to MERCI; maximal gap length (defaults to 1 if not spcified)")
    args = parser.parse_args()

    # Process command line arguments
    if args.thresholdhigh and args.thresholdlow:
        env.USE_THRESHOLD = True
        # Create classes file (1's and 0's) from file with real numbers using threshold
        convertData.makeClassFromThreshold(args.thresholdhigh, args.thresholdlow, args.filename)
    else:
        env.USE_THRESHOLD = False
    if args.generate:
        env.GENERATE = True    
    if args.verbose:
        env.VERBOSE = True
        env.print_env()
    if args.MERCIfp:
        MERCIfp = args.MERCIfp        
    else:
        MERCIfp = 10
    if args.MERCIfn:
        MERCIfn = args.MERCIfn
    else:
        MERCIfn = 10
    if args.MERCIg:
        MERCIg = args.MERCIg
    else:
        MERCIg = 1
    if args.MERCIgl:
        MERCIgl = args.MERCIgl
    else:
        MERCIgl = 1
            
    util.verbose_print(util.separatorNL)    
    # Name of classes file will be different when using threshold
    if env.USE_THRESHOLD:
        classesFileName = args.filename+".classes"
    else:
        classesFileName = args.filename

    util.verbose_print(util.separator)
    util.verbose_print("Running on " + "file: " + args.filename)    
    util.verbose_print("Name of classes file: " + classesFileName)    
    # Get Basename of input file
    basename = util.getBaseName(args.filename)
    util.verbose_print(util.separatorNL)

    # Create data directories if they don't already exist
    util.mk_dirs()

    # Generate positive file and negative file (2 total) from classes file
    convertData.generateFiles(classesFileName) 
    
    # Names of input files to MERCI/libSVM
    posFastaFile = basename+".POS.fa"
    negFastaFile = basename+".NEG.fa"
    posCSVFile = basename+".POS.csv"
    negCSVFile = basename+".NEG.csv"
    
    merciscript.initGlobals()   

    # Run MERCI on positive file; puts output in merciresult file
    merciscript.runmerci(posFastaFile,negFastaFile, fp=MERCIfp, fn=MERCIfn, g=MERCIg, gl=MERCIgl)
    
    # Parse results from MERCI; store data in memory
    motifs1 = merciscript.parse_output_file("+")
    
    # Run MERCI, this time on negative file; puts output in merciresult file
    merciscript.runmerci(negFastaFile,posFastaFile, fp=MERCIfp, fn=MERCIfn, g=MERCIg, gl=MERCIgl)
    
    # Parse results from MERCI; store data in memory
    motifs2 = merciscript.parse_output_file("-")

    # Merge results from positive and negative runs; keep data in memory
    motifs3 = merciscript.merge_motifs(motifs1, motifs2)

    util.verbose_print("Number of motifs: " + str(len(motifs3)))

    # Generate features
    cntb = merciscript.read_sequence_data(posCSVFile) # returns count
    merciscript.read_sequence_data(negCSVFile)
    merciscript.feature_vector_generator(motifs3, cntb)

    # Save features for input to libSVM
    merciscript.save_feature_vectors(basename+".combinedfeatures.csv", basename+".combinedfeatures-occ.csv")

    util.verbose_print(util.separatorNL)

    # Get name of training file
    TRAINPATH=env._env['TRAINING_DATA_PATH']
    training_file = os.path.join(TRAINPATH, basename)
    training_file = training_file + ".combinedfeatures.csv"
    util.verbose_print("training_file " + training_file)

    x, y = libsvm_wrapper.readData(training_file)

    # If -c (or --crossValidation), perform cross validation and stop
    if args.crossValidation:
        svmutil.svm_train(y, x, "-v 5")
        exit(0)

    # Train and save model
    m = svmutil.svm_train(y, x)
    modelName = os.path.join(env._env['TRAINED_SVM_PATH'] , basename+".model")
    svmutil.svm_save_model(modelName, m)
    util.verbose_print(util.separatorNL)

    # Generate new features if -g is set 
    if args.generate:
        util.verbose_print("Generating NEW features")
        if env.USE_THRESHOLD:
            intensityFn = os.path.join(env._env['INPUT_DATA'], args.filename)
            generator.doAll(basename + ".combinedfeatures.csv", intensityFn, 1000, 2, 3, True, threshList=[args.thresholdhigh, args.thresholdlow])
        else:
            # Not sure if generation with classes only works "correctly"
            print "WARNING: GENERATING ON CLASSES ONLY"
            intensityFn = os.path.join(env._env['INPUT_DATA'], basename+".csv")
            generator.doAll(basename_tfile, intensityFn, 1000, 2, 3, False)

    # Cleanup
    util.deleteFiles(env._env['TMP_FILES_PATH'])
