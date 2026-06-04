import numpy as np
import scipy.stats as stats
from os import remove
from os.path import exists
from bidict import bidict
import argparse


HASH_TO_PARAMETERS_FILE = "./hash_to_parameters.tsv"
CONFIG_FILE = "parameter_search_config.tsv"

'''
The Parameter class holds information for one giraffe parameter
It knows the range of potential values and can sample a value using different sampling functions

name is the name of the parameter getting called in giraffe (--name value)
value_range is a tuple of [range_start, range_end)
'''
class Parameter:
    def __init__(self, name, datatype, min_val, max_val, default, sampling_strategy, mean):
        self.name = name
        self.datatype = datatype.lower()
        self.min_val = min_val
        self.max_val = max_val
        self.default = default
        self.sampling_strategy = sampling_strategy.lower()
        if mean=="none":
            self.mean = default
        else:
            self.mean = mean
    
    def uniform(self):
        """
        Pick a sample from a uniform distribution bounded by min_val and max_val.
        """
        return np.random.randint(self.min_val, self.max_val)
    
    def log(self):
        """
        Pick a sample from a logarithmic distribution bounded by log(min_val) and log(max_val).
        """
        log_sample = np.random.uniform(np.log(self.min_val), np.log(self.max_val))
        return np.exp(log_sample)
    
    def lognormal(self):
        """
        Pick a sample from a lognormal distribution with a mean of self.mean and a standard dev of 1.
        Resample if outside of min_val and max_val.
        """
        if self.mean < self.min_val or self.mean > self.max_val:
            raise ValueError(f'Mean value of lognormal distribution ({self.mean}) cannot be outside of min/max bounds.')
        while True:
            ln_sample = np.random.lognormal(np.log(self.mean), 1) # change stddev here!
            if self.min_val <= ln_sample <= self.max_val:
                return ln_sample

    def truncated_normal(self):
        """
        Pick a sample from a truncated normal distribution (positive values only) with a mean of self.mean.
        Standard deviation is defined by a transformed version of min_val and max_val.
        Distribution is bounded by  a transformed version of min_val and max_val.
        """
        if self.mean < self.min_val or self.mean > self.max_val:
            raise ValueError(f'Mean value of truncated normal distribution ({self.mean}) cannot be outside of min/max bounds.')
        mu = self.mean
        sigma = (self.max_val - self.min_val) / 6 # adjust stddev here!
        # truncations are in number of stdv from median, below converts numeric min and max into the correct format
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html#scipy.stats.truncnorm
        a, b = (self.min_val - mu) / sigma, (self.max_val - mu) / sigma  
        truncated_normal = stats.truncnorm(a, b, loc=mu, scale=sigma)
        tn_sample = truncated_normal.rvs(size=1)[0]
        return tn_sample
    
    def sample(self):
        """
        Return a value randomly sampled from the distribution of choice in the correct datatype.
        """
        if self.min_val == self.max_val:
            return self.min_val

        value = -1
        if self.sampling_strategy == "uniform":
            value = self.uniform()
        elif self.sampling_strategy == "log":
            value = self.log()
        elif self.sampling_strategy == "lognormal":
            value = self.lognormal()
        elif self.sampling_strategy == "truncated_normal":
            value = self.truncated_normal()
        else:
            raise RuntimeError(f"Requested sampling strategy {self.sampling_strategy} not available.")
        
        if self.datatype == "int":
            value = int(value)
        elif self.datatype == "float":
            decimal_places = max(2, int(1/(self.max_val - self.min_val)))
            value = round(value, decimal_places)
        else:
            raise RuntimeError(f"Requested datatype {self.datatype} not supported.")

        return value
    
    def __repr__(self):
        return self.name + ":\n\ttype:" + self.datatype + "\n\trange: " + str(self.min_val) + "-" + str(self.max_val) + "\n\tdefault value: " + str(self.default) + "\n\tsampling strategy: " + self.sampling_strategy
    def __str__(self):
        return self.name + ":\n\ttype:" + self.datatype + "\n\trange: " + str(self.min_val) + "-" + str(self.max_val) + "\n\tdefault value: " + str(self.default) + "\n\tsampling strategy: " + self.sampling_strategy


'''
ParameterSearch is used to store information about a set of parameters.
Define the parameters to be searched in the parameter config file (default is CONFIG_FILE)
This must be a tsv with values:
#name   type    min_val max_val default sampling_strategy   mean(OPTIONAL)
Where name is the name of the flag that giraffe uses
type is the data type (int or float)
min and max val are the range of values that the parameter can take
default is the default value from giraffe. This is used to unify old runs missing parameters
sampling_strategy is how we sample the values from the range ("uniform", "log")
you can also sample from the truncated normal or lognormal distributions ("lognormal", "truncated_normal")
if desired, you can add a 7th column to the tsv titled mean for use in normal distributions

Randomly sample the parameter space with sample_parameter_space(), giving it the number of sets to return.
  This will write the sampled parameters to hash_to_parameters_file
  This gets exposed to the command line with add_random_parameters.py
Load previously generated parameter sets from hash_to_parameters_file with load_parameters_from_file().
  This can also be used to run a specific parameter set. Manually write to the hash_to_parameters_file and 
  use "." as a placeholder for the hash value. It will get filled in automatically.
  This automatically gets run every time ParameterSearch is initialized.
get_hashes_and_parameter_strings() is a generator for returning a tuple of hash value and parameter string for running
  in giraffe. It returns everything stored in the ParameterSearch from sample_parameter_space() and 
  load_parameters_from_file()
get_hashes() is a generator just for the hashes
'''
class ParameterSearch:
    def __init__(self, append, config=CONFIG_FILE, hash_to_parameters_file = HASH_TO_PARAMETERS_FILE):

        self.hash_to_parameters_file = hash_to_parameters_file

        #This defines all parameters that can be searched and their potential values
        #It gets loaded from the config file
        self.parameters = []

        f = open(config)
        for line in f:
            if line[0] != "#":
                l = line.split()
                if len(l) == 7:
                    if l[6] != "none":
                        mean_value = int(l[6]) if l[1] == "int" else float(l[6])
                    else:
                        mean_value = "none"
                else:
                    mean_value = "none"

                self.parameters.append(Parameter(l[0], l[1], 
                                                 int(l[2]) if l[1] == "int" else float(l[2]), 
                                                 int(l[3]) if l[1] == "int" else float(l[3]), 
                                                 int(l[4]) if l[1] == "int" else float(l[4]), 
                                                 l[5],
                                                 mean_value),
                                                 )
        f.close()

        #This maps a hash string to the set of parameters it represents, as a list of parameter values,
        # one for each parameter in self.parameters
        #Note that this is a two way dictionary. This is in case we later add parameters that didn't exist
        #before, we still want to be able to use the old hash value and fill in the default parameter values,
        #instead of re-hashing and losing the old results
        self.hash_to_parameters = bidict()

        # check if the file for it exists already
        file_exists = exists(self.hash_to_parameters_file)

        # if we want to add more parameters, load them
        if append and file_exists:
            self.load_parameters_from_file()
        else:
            # if we want to start fresh, remove the old one
            if not append and file_exists:
                remove(self.hash_to_parameters_file)
            # if we want to append but there is no file, warn
            if append and not file_exists:
                print(f"No previous hash_to_parameters file to append to. Creating new file.")
            
            # finally, make new file
            f = open(self.hash_to_parameters_file, "w")
            f.write("#hash\t" + '\t'.join(param.name for param in self.parameters))
            f.close()
        
    def save_parameters_to_file(self):
        """
        Write whatever is in the hash_to_parameters bidict to the output file.
        Bidict automatically removes repeated entries.
        """
        f = open(self.hash_to_parameters_file, "w")
        f.write("#hash\t" + '\t'.join(param.name for param in self.parameters))
        for k, v in self.hash_to_parameters.items():
            f.write('\n' + k + '\t' + '\t'.join(str(x) for x in v))
        f.close()

    def load_parameters_from_file(self):
        """
        Called when --append=True.
        Load parameter sets from the existing hash_to_parameters file,
        and create a new file with the new parameters appended.
        """
        #TODO make it so that you can also remove parameters from the config but still append correctly.

        f = open(self.hash_to_parameters_file, "r")

        f.readline() # skip header
        for line in f:
            l = line.split()
            new_params = tuple([(int(l[i+1]) if self.parameters[i].datatype == "int" else float(l[i+1])) if i < len(l)-1 else self.parameters[i].default for i in range(len(self.parameters))])
            #If there wasn't a hash value for the parameter set, then make one and rewrite everything
            hash_val = l[0]
            if hash_val == ".":
                hash_val = self.parameter_tuple_to_hash(new_params) 
            self.hash_to_parameters[hash_val] = new_params
        
        f.close()
        remove(self.hash_to_parameters_file)
        
        self.save_parameters_to_file()

    
    def parameter_tuple_to_hash(self, parameter_tuple):
        """
        Given a tuple representing a set of parameters, return the hash as a string
        TODO: Idk about this...
        """
        if parameter_tuple in self.hash_to_parameters.inverse:
            return self.hash_to_parameters.inverse[parameter_tuple]
        else:
            return str(abs(hash(parameter_tuple)))[:20];
        
    def parameter_tuple_to_parameter_string(self, parameter_tuple):
        """
        Given a tuple representing a set of parameters, return a string of options to be run in giraffe.
        """
        assert(len(parameter_tuple) == len(self.parameters))
        param_string = ""
        for i in range(len(parameter_tuple)):
            param_string+="--" + self.parameters[i].name
            param_string+=" " + str(parameter_tuple[i])
            if ( i != len(parameter_tuple)-1):
                param_string+=" "
        return param_string

    def hash_to_parameter_string(self, hash_val):
        return self.parameter_tuple_to_parameter_string(self.hash_to_parameters[hash_val])


    def sample_parameter_space(self, count, benchmark_default, benchmark_mean, static):
        """
        Sample the parameter space and save the new parameters to HASH_TO_PARAMETERS.
        """

        if benchmark_default:
            #add benchmark of default values
            benchmark_tuple = tuple(param.default for param in self.parameters)
            hash_val = self.parameter_tuple_to_hash(benchmark_tuple)
            self.hash_to_parameters[hash_val] = benchmark_tuple

            #f.write("\n" + hash_val + "\t" + '\t'.join([str(x) for x in benchmark_tuple]))

        if benchmark_mean:
            #add benchmark of mean values
            mean_tuple = tuple(param.mean for param in self.parameters)
            hash_val = self.parameter_tuple_to_hash(mean_tuple)
            self.hash_to_parameters[hash_val] = mean_tuple
            #f.write("\n" + hash_val + "\t" + '\t'.join([str(x) for x in mean_tuple]))

        if static:
            params = [
                {'chain-score-threshold': 234, 'min-chain-score-per-base': 0.24},
                {'chain-score-threshold': 888, 'min-chain-score-per-base': 0.29653052746724856},
                {'chain-score-threshold': 655, 'min-chain-score-per-base': 0.1562669394283844},
                {'chain-score-threshold': 816, 'min-chain-score-per-base': 0.3496272108061814},
                {'chain-score-threshold': 131, 'min-chain-score-per-base': 0.6171943193443219}
            ]
            for param in params:
                parameter_tuple = tuple(val for key, val in param.items())
                hash_val = self.parameter_tuple_to_hash(parameter_tuple)
                self.hash_to_parameters[hash_val] = parameter_tuple
        else:
            for i in range(count):
                parameter_tuple = tuple([param.sample() for param in self.parameters])
                hash_val = self.parameter_tuple_to_hash(parameter_tuple)
                self.hash_to_parameters[hash_val] = parameter_tuple

        self.save_parameters_to_file()
        
    
    def get_hashes(self):
        hashes = []
        for hash_val, parameter_tuple in self.hash_to_parameters.items():
            hashes.append(hash_val)
        return hashes

def main():
    parser = argparse.ArgumentParser(description="Add randomly sampled parameters to the file of parameters to search")
    parser.add_argument('--config-file', default=CONFIG_FILE, help="Config file for which parameters to sample and how") 
    parser.add_argument('--output-file', default=HASH_TO_PARAMETERS_FILE, help="File holding the parameter sets to search and their identifying hash value")
    parser.add_argument('--append', action="store_true", help="Whether or not to append new parameter sets to hash_to_parameters.tsv instead of rewriting it")
    parser.add_argument('--count', type=int, default=1000, help="How many parameters sets to sample [1000]")
    parser.add_argument('--benchmark-default', action="store_true", help="Whether or not to additonally run a benchmark of default parameters")
    parser.add_argument('--benchmark-mean', action="store_true", help="Whether or not to additionally run a benchmark of the mean parameters")
    parser.add_argument('--static', action="store_true", help="Whether to simply load in some parameter sets from the code")

    args = parser.parse_args()

    param_search = ParameterSearch(args.append, args.config_file, args.output_file)
    param_search.sample_parameter_space(args.count, args.benchmark_default, args.benchmark_mean, args.static)

if __name__ == "__main__":
    main()
