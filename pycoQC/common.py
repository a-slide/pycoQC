# -*- coding: utf-8 -*-

# Standard library imports
from os import access, R_OK, listdir, path
import inspect
from glob import iglob, glob
import sys
import logging
from collections import *

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

def is_readable_file (fn):
    """Verify the readability of a file or list of file"""
    return path.isfile (fn) and access (fn, R_OK)

def check_arg (arg_name, arg_val, required_type, allow_none=True, min=None, max=None, choices=[]):
    """Check argument values and type"""
    if allow_none and arg_val == None:
        return arg_val

    if not isinstance(arg_val, required_type):
        try:
            arg_val = required_type(arg_val)
        except:
            raise Exception ("Argument `{}` value `{}` is not in correct type: `{}` and cannot be coerced".format(arg_name, arg_val, required_type.__name__))

    if required_type in [float, int]:
        if min and arg_val < min:
            raise Exception ("Argument `{}` value `{}` is too low. Minimal value: {}".format(arg_name, arg_val, min))
        if max and arg_val > max:
            raise Exception ("Argument `{}` value `{}` is too high. Maximal value: {}".format(arg_name, arg_val, max))

    if choices and arg_val not in choices:
        raise Exception ("Argument `{}` value `{}` is not in the list of possible choices. Choices: {}".format(arg_name, arg_val, ", ".join(choices)))

    return arg_val

def sequencing_summary_file_sample (infile, outfile=None, n_seq=10000):
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

def dict_to_str (c, prefix="\t", suffix="\n"):
    """ Transform a dict to a tabulated str """
    m = ""
    if type(c) == Counter:
        for i, j in c.most_common():
            m += "{}{}: {:,}{}}".format(prefix, i, j, suffix)
    else:
        for i, j in c.items():
            m += "{}{}: {}{}".format(prefix, i, j, suffix)
    return m

def recursive_file_gen (dir, ext):
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

def get_logger (name=None, verbose=False, quiet=False):
    """Set logger to appropriate log level"""

    logging.basicConfig(format='%(message)s')
    logger = logging.getLogger(name)

    # Define overall verbose level
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    return logger

def doc_func (func):
    """Parse the function description string"""

    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        l = l.strip()
        if l:
            if l.startswith("*"):
                break
            else:
                docstr_list.append(l)

    return " ".join(docstr_list)

def make_arg_dict (func):
    """Parse the arguments default value, type and doc"""

    # Init method for classes
    if inspect.isclass(func):
        func = func.__init__

    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if not p.name in ["self","cls"]: # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if not name in ["kwargs","args"]: # Include but skip default required and type
                    # Get Annotation
                    if p.annotation != inspect._empty:
                        d[name]["type"] = p.annotation
                    # Get default value if available
                    if p.default == inspect._empty:
                        d[name]["required"] = True
                    else:
                        d[name]["default"] = p.default

        # Parse the docstring in a dict
        docstr_dict = OrderedDict()
        lab=None
        for l in inspect.getdoc(func).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    docstr_dict[lab] = []
                elif lab:
                    docstr_dict[lab].append(l)

        # Concatenate and copy doc in main dict
        for name in d.keys():
            if name in docstr_dict:
                d[name]["help"] = " ".join(docstr_dict[name])
        return d


def arg_opt (func, arg, **kwargs):
    """Get options corresponding to argumant name and deal with special cases"""
    arg_dict = make_arg_dict(func)[arg]

    if "default" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " (default: %(default)s)"

    if "type" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " [%(type)s]"

    # Special case for boolean args
    if arg_dict["type"] == bool:
        if arg_dict["default"] == False:
            arg_dict["action"] = 'store_true'
            del arg_dict["type"]
        elif arg_dict["default"] == True:
            arg_dict["action"] = 'store_false'
            del arg_dict["type"]

    # Special case for lists args
    elif arg_dict["type"] == list:
        arg_dict["nargs"]='*'

    return arg_dict

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    The docstring synthax should follow the same synthax as the one used for this function
    * f
        Function or method to display the help message for
    """
    # Private import as this is only needed if using jupyter
    from IPython.core.display import display, Markdown

    f_doc = doc_func(f)
    arg_doc = make_arg_dict(f)

    # Signature and function documentation
    s = "**{}** ({})\n\n{}\n\n---\n\n".format(f.__name__, ", ".join(arg_doc.keys()), f_doc)

    # Args doc
    for arg_name, arg_val in arg_doc.items():
        # Arg signature section
        s+= "* **{}**".format(arg_name)
        if "default" in arg_val:
            if arg_val["default"] == "":
                  arg_val["default"] = "\"\""
            s+= " (default: {})".format(arg_val["default"])
        if "required" in arg_val:
            s+= " (required)"
        if "type" in arg_val:
            if type(list) == type:
                s+= " [{}]".format(arg_val["type"].__name__)
            else:
                s+= " [{}]".format(arg_val["type"])
        s+="\n\n"
        # Arg doc section
        if "help" in arg_val:
            s+= "{}\n\n".format(arg_val["help"])

    # Display in Jupyter
    display (Markdown(s))

def head (fp, n=10, sep="\t", comment=None):
    """
    Emulate linux head cmd. Handle gziped files and bam files
    * fp
        Path to the file to be parse.
    * n
        Number of lines to print starting from the begining of the file (Default 10)
    """
    line_list = []

    # Get lines
    try:
        with open(fp) as fh:
            line_num = 0
            while (line_num < n):
                l= next(fh).strip()
                if comment and l.startswith(comment):
                    continue
                if sep:
                    line_list.append (l.split(sep))
                else:
                    line_list.append (l)
                line_num+=1

    except StopIteration:
        pass

    # Add padding if sep given
    if sep:
        try:
            # Find longest elem per col
            col_len_list = [0 for _ in range (len(line_list[0]))]
            for ls in line_list:
                for i in range (len(ls)):
                    len_col = len(ls[i])
                    if len_col > col_len_list[i]:
                        col_len_list[i] = len_col

            # Add padding
            line_list_tab = []
            for ls in line_list:
                s = ""
                for i in range (len(ls)):
                    len_col = col_len_list[i]
                    len_cur_col = len(ls[i])
                    s += ls[i][0:len_col] + " "*(len_col-len_cur_col)+" "
                line_list_tab.append(s)
            line_list = line_list_tab

        # Fall back to non tabulated display
        except IndexError:
            return head (fp=fp, n=n, sep=None)

    for l in line_list:
        print (l)
    print()
