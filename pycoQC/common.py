# -*- coding: utf-8 -*-

# Standard library imports
from os import access, R_OK, listdir, path
import inspect
from glob import iglob
import sys
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

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    Undocumented options are not displayed.
    The docstring synthax should follow the markdown formated convention below
    * f
        Function or method to display the help message for
    """
    # For some reason signature is not always importable. In these cases the build-in help is called instead
    try:
        from IPython.core.display import display, Markdown, HTML
    except (NameError, ImportError) as E:
        NanocomporeWarning ("jupyter notebook is required to use this function. Please verify your dependencies")
        return

    if inspect.isclass(f):
        f = f.__init__

    if inspect.isfunction(f) or inspect.ismethod(f):

        # Parse arguments default values and annotations
        sig_dict = OrderedDict()
        for name, p in inspect.signature(f).parameters.items():
            sig_dict[p.name] = []
            # Get Annotation
            if p.annotation != inspect._empty:
                sig_dict[p.name].append(": {}".format(p.annotation))
            # Get default value if available
            if p.default == inspect._empty:
                sig_dict[p.name].append("(required)")
            else:
                sig_dict[p.name].append("(default = {})".format(p.default))

        # Parse the docstring
        doc_dict = OrderedDict()
        descr = []
        lab=None
        for l in inspect.getdoc(f).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    doc_dict[lab] = []
                elif lab:
                    doc_dict[lab].append(l)
                else:
                    descr.append(l)

        # Reformat collected information in Markdown synthax
        s = "---\n\n**{}.{}**\n\n{}\n\n---\n\n".format(f.__module__, f.__name__, " ".join(descr))
        for k, v in doc_dict.items():
            s+="* **{}** *{}*\n\n{}\n\n".format(k, " ".join(sig_dict[k]), " ".join(v))

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
