# -*- coding: utf-8 -*-

"""       
  ___              ___   ___ 
 | _ \_  _ __ ___ / _ \ / __|
 |  _/ || / _/ _ \ (_) | (__ 
 |_|  \_, \__\___/\__\_\\___|
      |__/      
                                __   __     ___ 
 /\  _| _. _ _   |   _ _  _ _    _) /  \ /|   / 
/--\(_|| |(-| )  |__(-(_)(-|    /__ \__/  |  /  
                      _/                        
"""

# Standard library imports
from IPython.core.display import display, HTML, Markdown
from inspect import signature, isfunction
from pkg_resources import Requirement, resource_filename
from os import access, R_OK

##~~~~~~~ FUNCTIONS ~~~~~~~#

def print(*args, **kwargs):
    """
    Format a string in HTML and print the output. Equivalent of print, but highly customizable. Many options can be passed to the function.
    * args
        One or several objects that can be cast in str
    ** kwargs
        Formatting options to tweak the html rendering
        Boolean options : bold, italic, highlight, underlined, striked, subscripted, superscripted
        String options: font, color, size, align, background_color
    """

    # Join the different elements together and cast in string
    s =  " ".join([str(i) for i in args])

    # Replace new lines and tab by their html equivalent
    s = s.replace("\n", "<br>").replace("\t", "&emsp;")

    # For boolean options
    if "bold" in kwargs and kwargs["bold"]: s = "<b>{}</b>".format(s)
    if "italic" in kwargs and kwargs["italic"]: s = "<i>{}</i>".format(s)
    if "highlight" in kwargs and kwargs["highlight"]: s = "<mark>{}</mark>".format(s)
    if "underlined" in kwargs and kwargs["underlined"]: s = "<ins>{}</ins>".format(s)
    if "striked" in kwargs and kwargs["striked"]: s = "<del>{}</del>".format(s)
    if "subscripted" in kwargs and kwargs["subscripted"]: s = "<sub>{}</sub>".format(s)
    if "superscripted" in kwargs and kwargs["superscripted"]: s = "<sup>{}</sup>".format(s)

    # for style options
    style=""
    if "font" in kwargs and kwargs["font"]: style+= "font-family:{};".format(kwargs["font"])
    if "color" in kwargs and kwargs["color"]: style+= "color:{};".format(kwargs["color"])
    if "size" in kwargs and kwargs["size"]: style+= "font-size:{}%;".format(kwargs["size"])
    if "align" in kwargs and kwargs["align"]: style+= "text-align:{};".format(kwargs["align"])
    if "background_color" in kwargs and kwargs["background_color"]: style+= "background-color:{};".format(kwargs["background_color"])

    # Format final string
    if style: s = "<p style=\"{}\">{}</p>".format(style,s)
    else: s = "<p>{}</p>".format(s)

    display(HTML(s))

def help(function, full=False):
    """
    Print a nice looking help string based on the name of a declared function. By default print the function definition and description 
    * full
        If True, the help string will included a description of all arguments
    """
    if isfunction(function):
        name = function.__name__.strip()
        sig = str(signature(function)).strip()
        display(HTML ("<b>{}</b> {}".format(name, sig)))
        
        for line in function.__doc__.split("\n"):
            line = line.strip()
            if not full and line.startswith("*"):
                break
            display(Markdown(line.strip()))
    else:
        jprint("{} is not a function".format(function))

def get_sample_file (package, path):
    """
    Verify the existence and return a file from the package data
    * package
        Name of the package
    * path
        Relative path to the file in the package. Usually package_name/data/file_name 
    """
    sample_file = resource_filename(Requirement.parse(package), path)
    if not access(sample_file, R_OK):
        print ("Sample file {} from package {} cannot be read".format(path, package))
        return
    else:
        return sample_file
