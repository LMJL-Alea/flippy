import contextlib
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
    n = len(args)
  
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
            sys.exit("When the first input is a distance matrix, all subsequent 
                inputs should be integers specifying sample sizes.")
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
                    sys.exit("When the first input is univariate data, all 
                        subsequent inputs should be univariate data as well.")
            return purrr::map(l, purrr::array_tree, margin = 1)
  
    # Case of multivariate data
    if (is.matrix(l[[1]])) {
      if (n > 1) {
        coherent_inputs <- TRUE
        for (i in 2:n) {
          if (!is.matrix(l[[i]]) || (ncol(l[[i]]) != ncol(l[[1]]))) {
            coherent_inputs <- FALSE
            break
          }
        }
        stopifnot(coherent_inputs)
      }
      return(purrr::map(l, purrr::array_tree, margin = 1))
    }
  
    coherent_inputs <- TRUE
    for (i in 1:n) {
      if (!is.list(l[[i]])) {
        coherent_inputs <- FALSE
        break
      }
    }
    stopifnot(coherent_inputs)
  
    return l
