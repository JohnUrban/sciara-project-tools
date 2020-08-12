import gzip
import cPickle as pickle

def load_pickle(inpickle):
    if inpickle.endswith('.gz'):
        with gzip.open(inpickle,'rb') as pkfile:
            f = pickle.load(pkfile)
    else:
        with open(inpickle,'rb') as pkfile:
            f = pickle.load(pkfile)
    return f
