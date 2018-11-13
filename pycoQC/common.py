# -*- coding: utf-8 -*-

# Standard library imports
from os import access, R_OK
from inspect import signature, isfunction, ismethod

# Third party imports
import pandas as pd

#~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
class pycoQCError (Exception):
    """ Basic exception class for pycoQC package """
    pass

class pycoQCWarning (Warning):
    """ Basic Warning class for pycoQC package """
    pass

##~~~~~~~ FUNCTIONS ~~~~~~~#

def is_readable_file (fp, **kwargs):
    """Verify the readability of a file or list of file"""
    return access(fp, R_OK)

def sequencing_summary_file_sample (infile, outfile=None, n_seq=10000, **kwargs):
    """
    Sample a number read lines in infile and write the output_over_time in output_file
    If the file contains several runids the function will sample proportionally to the
    * infile: STR
    Path to a sequencing_summary input file
    * outfile: STR (default None)
        Path to a sequencing_summary output file. If not given, will return a dataframe instead
    * n_seq: STR (default 10000)
        Overall number of sequence lines to sample
    """
    df = pd.read_csv(infile, sep ="\t")
    df.dropna (inplace=True)
    total = len(df)
    print ("{} sequences".format(total))

    l = []
    for runid, runid_df in df.groupby("run_id", sort=False):
        n_to_sample = int (round (len (runid_df)/total*n_seq, 0))
        if n_to_sample == 0:
            n_to_sample=1
        print ("{} = {} seq, to sample = {}".format (runid, len(runid_df), n_to_sample))
        sdf = runid_df.sample(n_to_sample)
        sdf.sort_values("start_time", inplace=True)
        l.append(sdf)

    df = pd.concat(l)
    df.reset_index(inplace=True, drop=True)
    if outfile:
        df.to_csv(outfile, index=False, sep="\t")
    else:
        return df

def print_help (function):
    """Print a nice looking help string based on the name of a declared function"""
    if isfunction(function) or ismethod(function):
        print ("{} {}\n{}".format(function.__name__, signature(function), function.__doc__))
    else:
        print("{} is not a function".format(function))
