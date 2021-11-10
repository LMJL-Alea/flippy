import numpy as np
import pandas as pd

class InferredParameter(object):
    """
    """
    # Public
    def __init__(self):
        self.estimate = None
        self.bounds = [None, None]
        self.inclusive = [True, True]
        self.label = None
    
    def __repr__(self):
        return "InferredParameter()"
    
    def __str__(self):
        return """
        {label} parameter
        
        Estimate: {estimate}
        Bounds  : {lbr}{lb}, {ub}{ubr}
        """.format(
            label = "Unnamed" if self.label is None else self.label, 
            estimate = "?" if self.estimate is None else self.estimate, 
            lbr = "[" if self.inclusive[0] else "(", 
            lb = "?" if self.bounds[0] is None else self.bounds[0], 
            ub = "?" if self.bounds[1] is None else self.bounds[1], 
            ubr = "]" if self.inclusive[1] else ")"
        )
    
    def set_estimate(self, val):
        if self.__validate_estimate(val):
            self.estimate = val
    
    def set_bounds(self, val):
        if self.__validate_bounds(val):
            self.bounds = val
    
    def set_inclusive(self, val):
        if self.__validate_inclusive(val):
            self.inclusive = val
    
    def set_label(self, val):
        if self.__validate_label(val):
            self.label = val
    
    # Private
    def __check_bounds(self, estimate, bounds):
        if not estimate is None:
            if not bounds[0] is None:
                out_of_bounds = False
                if self.inclusive[0] and estimate < bounds[0]:
                    out_of_bounds = True
                if not self.inclusive[0] and estimate <= bounds[0]:
                    out_of_bounds = True
                if out_of_bounds:
                    print("Estimate cannot be smaller than lower bound.")
                    return False
            if not bounds[1] is None:
                out_of_bounds = False
                if self.inclusive[1] and estimate > bounds[1]:
                    out_of_bounds = True
                if not self.inclusive[1] and estimate >= bounds[1]:
                    out_of_bounds = True
                if out_of_bounds:
                    print("Estimate cannot be larger than upper bound.")
                    return False
        return True
    
    def __validate_estimate(self, val):
        if not isinstance(val, (int, float)) and not val is None:
            print("The `estimate` parameter should be a real scalar or None.")
            return False
        return self.__check_bounds(val, self.bounds)
    
    def __validate_bounds(self, val):
        if not isinstance(val, list):
            print("The `bounds` parameter should be a list.")
            return False
        if len(val) != 2:
            print("The `bounds` parameter should be of length 2.")
            return False
        if not all([e is None or isinstance(e, (int, float)) for e in val]):
            print("The `bounds` parameter should contain only real values or None.")
            return False
        if all([isinstance(e, (int, float)) for e in val]):
            if val[0] > val[1]:
                print("The `bounds` parameter should be an interval. As such, the first value should be smaller than the second value.")
                return False
        return self.__check_bounds(self.estimate, val)
    
    def __validate_inclusive(self, val):
        if not isinstance(val, list):
            print("The `inclusive` parameter should be a list.")
            return False
        if len(val) != 2:
            print("The `inclusive` parameter should be of length 2.")
            return False
        if not all([isinstance(e, bool) for e in val]):
            print("The `inclusive` parameter should contain only booleans.")
            return False
        return True
    
    def __validate_label(self, val):
        if not isinstance(val, str) and not val is None:
            print("The `label` parameter should be a string or None.")
            return False
        return True

class InferredParameterSetIterator(object):
    '''
    Iterator class
    '''
    def __init__(self, param_set):
        self._param_set = param_set
        self._index = 0
    
    def __next__(self):
        '''Returns the next value from team object's lists'''
        if (self._index < len(self._param_set)):
            old_index = self._index
            self._index += 1
            return self._param_set[old_index]
        raise StopIteration

class InferredParameterSet(object):
    """
    """
    def __init__(self):
        self.parameters = []
    
    def __repr__(self):
        return "InferredParameterSet()"
    
    def __str__(self):
        l = ["""
        {label} parameter
        
        Estimate: {estimate}
        Bounds  : {lbr}{lb}, {ub}{ubr}
        """.format(
            label = "Unnamed" if param.label is None else param.label, 
            estimate = "?" if param.estimate is None else param.estimate, 
            lbr = "[" if param.inclusive[0] else "(", 
            lb = "?" if param.bounds[0] is None else param.bounds[0], 
            ub = "?" if param.bounds[1] is None else param.bounds[1], 
            ubr = "]" if param.inclusive[1] else ")"
        ) for param in self.parameters]
        return "\n\n".join(l)
    
    def __len__(self):
        return len(self.parameters)
    
    def __getitem__(self, item):
        return self.parameters[item]
    
    def __iter__(self):
        ''' Returns the Iterator object '''
        return InferredParameterSetIterator(self)
    
    def add_parameter(self, param):
        if not isinstance(param, InferredParameter):
            print("The object to be added is not an InferredParameter().")
            return
        self.parameters += [param]
    
    def create_grid(self, n = 21):
        if n % 2 == 0:
            n += 1
        half_n = int((n - 1) / 2)
        grids = [
            list(np.linspace(param.bounds[0], param.estimate, half_n + 1))[:-1] + \
            [param.estimate] + \
            list(np.linspace(param.estimate, param.bounds[1], half_n + 1))[1:] \
            for param in self.parameters
        ]
        M = np.array(np.meshgrid(*grids)).reshape(len(grids), n**len(grids)).T
        return pd.DataFrame(M, columns = [param.label for param in self.parameters])
