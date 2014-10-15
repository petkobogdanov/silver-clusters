import os, inspect

global _env

USE_THRESHOLD = False
VERBOSE = False
GENERATE = False



_env = {}

_env['SCRIPT_DIRECTORY'] = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # directory this file is in

#_env['MAIN'] = os.path.abspath(os.path.join(os.getcwd(), os.path.join(os.pardir, os.pardir))) 

#_env['SRC_PATH'] = os.path.abspath(os.path.join(_env['MAIN'], "src"))

_env['SRC_PATH'] = os.path.abspath(os.path.join(_env['SCRIPT_DIRECTORY'], os.pardir))

_env['DATA_PATH'] = os.path.abspath( os.path.join(os.path.join(_env['SCRIPT_DIRECTORY'], os.pardir), "data"))

# _env['MEME_BIN_PATH'] = os.path.abspath(
#     os.path.join(_env['MEME_PATH'] , 'bin') )
# _env['MEME_EXE_PATH'] = os.path.abspath(
#     os.path.join(_env['MEME_BIN_PATH'] , 'meme') )
    
_env['NEW_FEATURES_PATH'] = os.path.abspath(
    os.path.join(_env['DATA_PATH'], "new_features")
)

_env['TRAINED_SVM_PATH'] = os.path.abspath(
    os.path.join(_env['DATA_PATH'], "trained_svm")
)

_env['MERCI_PATH'] = os.path.abspath(os.path.join(_env['SRC_PATH'], "merci"))

_env['TMP_FILES_PATH'] = os.path.abspath(
    os.path.join(_env['DATA_PATH'] , 'tmp_files'))

_env['INPUT_DATA'] = os.path.abspath(
    os.path.join(_env['DATA_PATH'], 'input_data'))

# _env['RAW_DATA'] = os.path.abspath(
    # os.path.join(_env['DATA_PATH'], "raw_data"))

_env['TRAINING_DATA_PATH'] = os.path.abspath(
    os.path.join(_env['DATA_PATH'] , 'training_data'))

_env['MERCI_EXE_PATH'] = os.path.abspath(
    os.path.join(_env['MERCI_PATH'] , 'MERCI.pl') )
    
_env['MERCI_CLASSIFICATION_PATH'] = os.path.abspath(
    os.path.join(_env['MERCI_PATH'] , 'classification') )
    
_env['MERCI_OUTPUT_FILE'] = os.path.abspath(
    os.path.join(_env['TMP_FILES_PATH'] , 'merciresult') )
    
_env['MERCI_OCCURENCE_OUTPUT_FILE'] = os.path.abspath(
    os.path.join(_env['TMP_FILES_PATH'] , 'merciresult.occurrences') )


LIBSVM_VERSION_DIR = "libsvm-3.18"
_env['LIBSVM_PYTHON_PATH'] = os.path.abspath(
    os.path.join(os.path.join(os.path.join(_env['SRC_PATH'], "libSVM"), LIBSVM_VERSION_DIR), "python")
)

envList = {"SCRIPT_DIRECTORY": _env['SCRIPT_DIRECTORY'], "DATA_PATH": _env['DATA_PATH'], "SRC_PATH": _env['SRC_PATH'],"MERCI_PATH" : _env['MERCI_PATH'],
           "TMP_FILES_PATH": _env['TMP_FILES_PATH'], "INPUT_DATA" : _env['INPUT_DATA'],
           "TRAINING_DATA_PATH": _env['TRAINING_DATA_PATH'], "MERCI_EXE_PATH": _env['MERCI_EXE_PATH'],
           "MERCI_CLASSIFICATION_PATH": _env['MERCI_CLASSIFICATION_PATH'], "MERCI_OUTPUT_FILE":_env['MERCI_OUTPUT_FILE'],
           "MERCI_OCCURENCE_OUTPUT_FILE" : _env['MERCI_OCCURENCE_OUTPUT_FILE'], "LIBSVM_PYTHON_PATH":_env['LIBSVM_PYTHON_PATH'],
           "NEW_FEATURES_PATH": _env['NEW_FEATURES_PATH']
}

def print_env():
    print "Environment:"
    for k, v in envList.iteritems():
        print "{0:30} {1:10}".format(k, v)
