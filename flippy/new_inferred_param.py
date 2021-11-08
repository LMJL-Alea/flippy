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
