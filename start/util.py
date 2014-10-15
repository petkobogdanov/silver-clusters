import os

def printSeparator(nl=None):
    print "------------------------------------------------------------------"
    if nl is not None:
        print ""

def runFtn(fn, args):
    print fn.__name__
    fn(args)        
        
def deleteFiles(path):
    for f in os.listdir(path):
        fp = os.path.join(path, f)
        try:
            if os.path.isfile(fp):
                os.unlink(fp)
        except Exception, e:
            print e    
            
