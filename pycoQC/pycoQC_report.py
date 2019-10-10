#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import json
from pkg_resources import resource_filename
import datetime
import os

# Third party imports
import plotly.offline as py
import jinja2

# Local imports
from pycoQC.common import *
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot
from pycoQC import __version__ as package_version
from pycoQC import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class pycoQC_report ():

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INIT METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __init__ (self,
        parser:pycoQC_parse,
        plotter:pycoQC_plot,
        verbose:bool=False,
        quiet:bool=False):
        """
        * parser
            A pycoQC_parse object
        * plotter
            A pycoQC_plot object
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """
        # Set logging level
        self.logger = get_logger (name=__name__, verbose=verbose, quiet=quiet)

        # Check that parser is a valid instance of pycoQC_parse
        if not isinstance(parser, pycoQC_parse):
            raise pycoQCError ("{} is not a valid pycoQC_parse object".format(parser)) ##################################### Use parser to get file list for report
        self.parser = parser

        # Check that plotter is a valid instance of pycoQC_plot
        if not isinstance(plotter, pycoQC_plot):
            raise pycoQCError ("{} is not a valid pycoQC_plot object".format(plotter))
        self.plotter = plotter

    def __repr__(self):
        return "[{}]\n".format(self.__class__.__name__)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def html_report( self,
        outfile:str,
        config_file:str="",
        template_file:str="",
        report_title:str="PycoQC report"):
        """"""
        self.logger.info("Generating HTML report")

        # Parse configuration file
        self.logger.info("\tParsing html config file")
        config_dict = self._get_config(config_file)
        self.logger.debug(config_dict)

        # Loop over configuration file and run the pycoQC functions defined
        plots = list()
        titles = list()
        for method_name, method_args in config_dict.items ():
            try:
                self.logger.info("\tRunning method {}".format(method_name))
                self.logger.debug ("\t{} ({})".format(method_name, method_args))

                # Store plot title for HTML title and remove from data passed to plotly
                plot_title = method_args["plot_title"]
                method_args["plot_title"]=""

                # Get method and generate plot
                method = getattr(self.plotter, method_name)
                fig = method(**method_args)
                plot = py.plot(
                    fig,
                    output_type='div',
                    include_plotlyjs=False,
                    image_width='',
                    image_height='',
                    show_link=False,
                    auto_open=False)

                plots.append(plot)
                titles.append(plot_title)

            except AttributeError as E:
                self.logger.info("\t\t{}".format("{} is not a valid plotting method".format(method_name)))
                self.logger.info("\t\t{}".format(E))

            except pycoQCError as E:
                self.logger.info("\t\t{}".format(E))

        # Load HTML template for Jinja
        self.logger.info("\tLoading HTML template")
        template = self._get_jinja_template(template_file)

        # Set a subtitle for the HTML report
        report_subtitle="Generated on {} with {} {}".format( datetime.datetime.now().strftime("%d/%m/%y"), package_name, package_version)

        # Define source files list
        src_files = ""
        for files_list, name in (
            (self.parser.summary_files_list,"summary"),
            (self.parser.barcode_files_list,"barcode"),
            (self.parser.bam_file_list,"bam")):

            if files_list:
                src_files += "<h4>Source {} files</h4><ul>".format(name)
                for f in files_list:
                    f = os.path.abspath(f)
                    src_files += "<li>{}</li>".format(f)
                src_files += "</ul>"

        # Render plots
        self.logger.info("\tRendering plots in d3js")
        rendering = template.render(
            plots=plots,
            titles=titles,
            plotlyjs=py.get_plotlyjs(),
            report_title=report_title,
            report_subtitle=report_subtitle,
            src_files=src_files)

        # Write to HTML file
        self.logger.info("\tWriting to HTML file")
        with open(outfile, "w") as fp:
            fp.write(rendering)

    def json_report(self,
        outfile:str):
        """"""
        self.logger.info("Generating JSON report")
        self.logger.info("\tRunning summary_stats_dict method")
        res_dict = self.plotter.summary_stats_dict ()

        self.logger.info("\tWriting to JSON file")
        with open (outfile, "w") as fp:
            json.dump(res_dict, fp, indent=2)

    #~~~~~~~~~~~~~~PRIVATE FUNCTION~~~~~~~~~~~~~~#

    def _get_config(self, config_file=None):
        """"""
        # First, try to read provided configuration file if given
        if config_file:
            self.logger.debug ("\tTry to read provided config file")
            try:
                with open(config_file, 'r') as cf:
                    return json.load(cf)
            except (FileNotFoundError, IOError, json.JSONDecodeError):
                self.logger.debug ("\t\tConfiguration file not found, non-readable or invalid")

        # Last use the default harcoded config_dict
        self.logger.debug ("\tRead default configuration file")
        config_file = resource_filename("pycoQC", "templates/pycoQC_config.json")
        with open(config_file, 'r') as cf:
            return json.load(cf)

    def _get_jinja_template(self, template_file=None):
        """"""
        # First, try to read provided configuration file if given
        if template_file:
            self.logger.debug("\tTry to load provided jinja template file")
            try:
                with open(template_file) as fp:
                    template = jinja2.Template(fp.read())
                    return template
            except (FileNotFoundError, IOError, jinja2.exceptions.TemplateNotFound, jinja2.exceptions.TemplateSyntaxError):
                self.logger.debug ("\t\tFile not found, non-readable or invalid")

        # Last use the default harcoded config_dict
        self.logger.debug ("\tRead default jinja template")
        env = jinja2.Environment (loader=jinja2.PackageLoader('pycoQC', 'templates'), autoescape=False)
        template = env.get_template('spectre.html.j2')
        return template
