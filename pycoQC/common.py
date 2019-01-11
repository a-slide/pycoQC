# -*- coding: utf-8 -*-

# Standard library imports
from os import access, R_OK, listdir, path
from inspect import signature, isfunction, ismethod
from glob import iglob
import sys

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
        if outfile.endswith("gz"):
            df.to_csv(outfile, index=False, sep="\t", compression="gzip")
        else:
            df.to_csv(outfile, index=False, sep="\t", compression=None)
    else:
        return df

def print_help (function):
    """Print a nice looking help string based on the name of a declared function"""
    if isfunction(function) or ismethod(function):
        print ("{} {}\n{}".format(function.__name__, signature(function), function.__doc__))
    else:
        print("{} is not a function".format(function))

def recursive_file_gen (dir, ext, **kwargs):
    """
    create a generator listing all files with a particular extension in a folder arborescence
    The recursivity is broken when at least 1 file with a particular extenssion is found.
    """
    # In the case where the folder is a file
    if path.isdir(dir):

        # If matching files in the folder
        file_found=False
        for fn in iglob (path.join(dir, "*."+ext)):
            yield fn
            file_found=True

        # If no matching file go deeper until a leaf containing fast5 is found
        if not file_found:
            for item in listdir(dir):
                for fn in recursive_file_gen (path.join(dir, item), ext):
                    yield fn


def stderr_print (*args):
    """reproduce print with stderr.write
    """
    sys.stderr.write(" ".join(str(a) for a in args))
    sys.stderr.flush()

def counter_to_str (c):
    """Transform a counter dict to a tabulated str"""
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m
