import contextlib
import itertools
import sys
import numpy as np
from triarray import TriMatrix

def match_arg(x, lst):
    l = [el for el in lst if x in el]
    if len(l) == 0:
        sys.exit("x should be one of %s." % lst)
    else:
        return x

def set_doc(docum):
    def set_documentation(func):
        func.__doc__ = docum
        return func
    return set_documentation

@contextlib.contextmanager
def local_seed(seed):
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)

def convert_to_list(*args):
    l = args
    n = len(l)
  
    # Case "No input samples"
    if n == 0:
        return None
  
    # Case of distance matrix
    if isinstance(l[0], TriMatrix):
        if n == 1:
            return l
        coherent_inputs = True
        for i in range(1, n):
            if not isinstance(l[i], integer) or len(l[i]) != 1:
                coherent_inputs <- FALSE
                break
        if not coherent_inputs:
            sys.exit("When the first input is a distance matrix, all subsequent inputs should be integers specifying sample sizes.")
        return l
  
    if isinstance(l[0], np.ndarray):
        if l[0].ndim == 1:
            # Case of univariate data
            if n > 1:
                coherent_inputs = True
                for i in range(1, n):
                    if not isinstance(l[i], np.ndarray):
                        coherent_inputs = False
                        break
                    if l[i].ndim != 1:
                        coherent_inputs = False
                        break
                if not coherent_inputs:
                    sys.exit("When the first input is univariate data, all subsequent inputs should be univariate data as well.")
            return [[x[i,:] for i in range(x.shape[0])] for x in l]
        if l[0].ndim == 2:
            # Case of multivariate data
            if n > 1:
                coherent_inputs = True
            for i in range(1, n):
                if not isinstance(l[i], np.ndarray):
                    coherent_inputs = False
                    break
                if l[i].ndim != 2:
                    coherent_inputs = False
                    break
                if l[i].shape[1] != l[0].shape[1]:
                    coherent_inputs = False
                    break
                if not coherent_inputs:
                    sys.exit("When the first input is multivariate data, all subsequent inputs should be multivariate data as well.")  
            return [[x[i,:] for i in range(x.shape[0])] for x in l]
  
    coherent_inputs = True
    for i in range(n):
      if not isinstance(l[i], list):
          coherent_inputs = False
          break
    if not coherent_inputs:
        sys.exit("When the first input is list data, all subsequent inputs should be list data as well.")
  
    return l

def combn(n, k):
    x = list(itertools.combinations(range(n), k))
    return np.array(x).transpose()
