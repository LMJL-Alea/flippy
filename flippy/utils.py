import sys

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
