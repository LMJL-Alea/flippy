from .new_inferred_param import InferredParameter, InferredParameterSet
from .two_sample_test import two_sample_test
from .utils import convert_to_list

#' R6 Class representing a plausibility function
#'
#' @description A plausibility function is...
#'
#' @export
class PlausibilityFunction(object):
    """
    """
    
    def __init__(self, *args, null_spec, stat_functions, stat_assignments, seed = None):
        """
        Create a new plausibility function object.
        
        Parameters:
            null_spec (function): For one-sample problems, it should 
                transform the `x` sample (provided as first argument) using 
                the parameters (as second argument) to make its distribution 
                centered symmetric. For two-sample problems, it should 
                transform the `y` sample (provided as first argument) using 
                the parameters (as second argument) to make it exchangeable 
                with the `x` sample under a null hypothesis.
            stat_functions (list of functions): It specifies the set of test 
                statistics that should be used.
            stat_assignments (dictionary): A dictionary with named entry 
                whose values are integer vectors specifying which test 
                statistic should be associated with each parameter. The 
                length of this dictionary should match the number of 
                parameters under investigation and is thus used to set it. 
                The key of each element of the dictionary should be named 
                after the parameter it identifies.
            *args: Vectors, matrices or lists providing the observed samples.
            seed (double): Seed to be used. Defaults to `None` in which case 
                `seed = 1234` is used and the user is informed of this setting.
        
        Value:
            A new `PlausibilityFunction` object.
        """
        if not callable(null_spec):
            print("The `null_spec` argument should be of type `function`.")
            return
        self.__null_spec = null_spec
  
        self.__stat_functions = stat_functions
  
        if not isinstance(stat_assignments, dict):
            print("The `stat_assignements` argument should be of type `dict`.")
            return
        self.stat_assignments = stat_assignments
        self.nparams = len(stat_assignments)
  
        self.set_data(*args)
        
        self.parameters = InferredParameterSet()
        for name in self.stat_assignments.keys():
            param = InferredParameter()
            param.set_label(name)
            self.parameters.add_parameter(param)
        self.point_estimate = [param.estimate for param in self.parameters]
  
        if seed is None:
            print("Setting the seed for sampling permutations is mandatory for obtaining a continuous p-value function. Using `seed = 1234`.")
            seed = 1234
        self.__seed = seed
        
        self.nperms = 1000
        self.nperms_max = None
        self.alternative = "two_tail"
        self.__alternative_choices = ["two_tail", "left_tail", "right_tail"]
        self.aggregator = "tippett"
        self.__aggregator_choices = ["fisher", "tippett"]
        self.formula = "exact"
        self.__formula_choices = ["estimate", "exact", "upper_bound"]
    
    def set_data(self, *args):
        self.__data = convert_to_list(*args)
        self.__nsamples = len(self.__data)
        if self.__nsamples > 2:
            sys.exit("The PlausibilityFunction class currently only support one- and two-sample problems.")
    
    def set_nperms(self, val):
        self.nperms = val
    
    def set_nperms_max(self, val):
        self.nperms_max = val

    def set_alternative(self, val):
        if not val in self.__alternative_choices:
            sys.exit("The `alternative` argument should be one of" + ", ".join(self.__alternative_choices) + ".")
        self.alternative = val

    def set_aggregator(self, val):
        if not val in self.__aggregator_choices:
            sys.exit("The `aggregator` argument should be one of" + ", ".join(self.__aggregator_choices) + ".")
        self.aggregator = val

    def set_formula(self, val):
        if not val in self.__formula_choices:
            sys.exit("The `formula` argument should be one of" + ", ".join(self.__formula_choices) + ".")
        self.formula = val

    #' @description Computes an indicator of the plausibility of specific values
    #'   for the parameters of interest in the form of a p-value of an
    #'   hypothesis test against these values.
    #'
    #' @param parameters A vector whose length should match the `nparams` field
    #'   providing specific values of the parameters of interest for assessment
    #'   of their plausibility in the form of a p-value of the corresponding
    #'   hypothesis test.
    #' @param keep_null_distribution A boolean specifying whether the empirical
    #'   permutation null distribution should be returned as well. Defaults to
    #'   `FALSE`.
    #' @param keep_permutations A boolean specifying whether the list of sampled
    #'   permutations used to compute the empirical permutation null
    #'   distribution should be returned as well. Defaults to `FALSE`.
    #' @param ... Extra parameters specific to some statistics.
    #'
    #' @examples
    #' x <- rnorm(10)
    #' y <- rnorm(10, mean = 2)
    #' null_spec <- function(y, parameters) {purrr::map(y, ~ .x - parameters[1])}
    #' stat_functions <- list(stat_t)
    #' stat_assignments <- list(mean = 1)
    #' pf <- PlausibilityFunction$new(
    #'   null_spec = null_spec,
    #'   stat_functions = stat_functions,
    #'   stat_assignments = stat_assignments,
    #'   x, y
    #' )
    #' pf$set_nperms(50)
    #' pf$get_value(2)
    def get_value(self, 
                  parameters, 
                  keep_null_distribution = False, 
                  keep_permutations = False, 
                  **kwargs):
        if len(parameters) != self.nparams:
            sys.exit("The plausibility function has been defined to infer {nparams} parameters and you are trying to evaluate it for a vector of parameters of length {lenparam}.".format(nparams = self.nparams, lenparam = len(parameters)))
        if self.__nsamples == 1:
            x = self.__null_spec(self.__data[0], parameters)
            test_result = one_sample_test(
                x = x,
                stats = self.__stat_functions,
                B = self.nperms,
                M = self.nperms_max,
                alternative = self.alternative,
                formula = self.formula,
                combine_with = self.aggregator,
                seed = self.__seed, 
                **kwargs
            )
        else:
            y = self.__null_spec(self.__data[1], parameters)
            test_result = two_sample_test(
                x = self.__data[0],
                y = y,
                stats = self.__stat_functions,
                B = self.nperms,
                M = self.nperms_max,
                alternative = self.alternative,
                formula = self.formula,
                combine_with = self.aggregator,
                seed = self.__seed, 
                **kwargs
            )
        
        if keep_null_distribution and keep_permutations:
          return test_result
        elif keep_null_distribution:
          del test_result["permutations"]
          return test_result
        elif keep_permutations:
          del test_result["null_distribution"]
          return test_result
        else:
          return test_result["pvalue"]
