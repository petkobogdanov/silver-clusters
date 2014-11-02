import os, env

separator = "------------------------------------------------------------------"
separatorNL = "------------------------------------------------------------------\n"

def runFtn(fn, args):
    print fn.__name__
    fn(args)        

# Remove all files in directory 'path'
# This is used to clean up temporary files
def deleteFiles(path):
    for f in os.listdir(path):
        fp = os.path.join(path, f)
        try:
            if os.path.isfile(fp):
                os.unlink(fp)
        except Exception, e:
            print e    

# Prints s when the -v (or --verbose) command line option is specified
def verbose_print(s):
    if env.VERBOSE is True:
        print s
    else:
        return

# Create data directories if they don't already exist
def mk_dirs():
    if not os.path.exists(env._env['TMP_FILES_PATH']):
        os.makedirs(env._env['TMP_FILES_PATH'])
    if not os.path.exists(env._env['TRAINING_DATA_PATH']):
        os.makedirs(env._env['TRAINING_DATA_PATH'])
    if not os.path.exists(env._env['NEW_FEATURES_PATH']):
        os.makedirs(env._env['NEW_FEATURES_PATH'])
    if not os.path.exists(env._env['TRAINED_SVM_PATH']):
        os.makedirs(env._env['TRAINED_SVM_PATH'])        


def getBaseName(fn):
    return os.path.splitext(fn)[0]        
