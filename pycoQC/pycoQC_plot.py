# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import *
import warnings
import datetime

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# Local lib import
from pycoQC.common import *
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC import __name__ as package_name
from pycoQC import __version__ as package_version

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL SETTINGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set seed for deterministic random sampling
SEED = 42
np.random.RandomState(seed=SEED)

# Silence futurewarnings
warnings.filterwarnings("ignore", category=FutureWarning)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class pycoQC_plot ():

    def __init__ (self,
        parser:pycoQC_parse,
        min_pass_qual:int=7,
        min_pass_len:int=0,
        sample:int=100000,
        verbose:bool=False,
        quiet:bool=False):
        """
        * parser
            A pycoQC_parse object
        * min_pass_qual
            Minimum quality to consider a read as 'pass'
        * min_pass_len
            Minimum read length to consider a read as 'pass'
        * sample
            If not None a n number of reads will be randomly selected instead of the entire dataset for plotting function (deterministic sampling)
        """

        # Set logging level
        self.logger = get_logger (name=__name__, verbose=verbose, quiet=quiet)
        self.logger.warning ("Loading plotting interface")

        # Save args to self values
        self.min_pass_qual = min_pass_qual
        self.sample = sample

        # Check that parser is a valid instance of pycoQC_parse
        if not isinstance(parser, pycoQC_parse):
            raise pycoQCError ("{} is not a valid pycoQC_parse object".format(parser))
        self.parser = parser

        # Extract values from parser object
        self.all_df = parser.reads_df
        if self.has_alignment:
            self.ref_len_dict = parser.ref_len_dict
            self.alignments_df = parser.alignments_df
        self.logger.info ("\tFound {:,} total reads".format(len(self.all_df)))

        # Save df wiews and compute scaling factors
        if sample and len(self.all_df)>sample:
            self.all_sample_df = self.all_df.sample(n=sample, random_state=SEED)
            self.all_scaling_factor = len(self.all_df)/sample
        else:
            self.all_sample_df = self.all_df
            self.all_scaling_factor = 1

        self.pass_df = self.all_df.query ("mean_qscore>={} and read_len>={}".format(min_pass_qual, min_pass_len))
        if sample and len(self.pass_df)>sample:
            self.pass_sample_df = self.pass_df.sample(n=sample, random_state=SEED)
            self.pass_scaling_factor = len(self.pass_df)/sample
        else:
            self.pass_sample_df = self.pass_df
            self.pass_scaling_factor = 1
        self.logger.info ("\tFound {:,} pass reads (qual >= {} and length >= {})".format(len(self.pass_df), min_pass_qual, min_pass_len))

    def __str__(self):
        m = ""
        m+= "\tBarcode: {}\n".format(self.has_barcodes)
        m+= "\tAlignment: {}\n".format(self.has_alignment)
        m+= "\tPromethion: {}\n".format(self.is_promethion)
        m+= "\tAll reads: {:,}\n".format(len(self.all_df))
        m+= "\tAll bases: {:,}\n".format(int(self.all_df["read_len"].sum()))
        m+= "\tAll median read length: {:,}\n".format(np.median(self.all_df["read_len"]))
        m+= "\tPass reads: {:,}\n".format(len(self.pass_df))
        m+= "\tPass bases: {:,}\n".format(int(self.pass_df["read_len"].sum()))
        m+= "\tPass median read length: {:,}\n".format(np.median(self.pass_df["read_len"]))
        return m

    def __repr__(self):
        return "[{}]\n".format(self.__class__.__name__)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    @property
    def has_barcodes (self):
        return "barcode" in self.all_df
    @property
    def has_alignment (self):
        return "ref_id" in self.all_df
    @property
    def has_identity_freq (self):
        return "identity_freq" in self.all_df
    @property
    def is_promethion (self):
        return self.all_df["channel"].max() > 512
    @property
    def total_ref_len (self):
        if self.has_alignment:
            return np.sum(list(self.ref_len_dict.values()))

    #~~~~~~~SUMMARY_STATS_DICT METHOD AND HELPER~~~~~~~#
    def summary_stats_dict (self):
        """
        Return a dictionnary containing exhaustive information about the run.
        """
        self.logger.info ("\tCompute overall summary statistics")
        d = OrderedDict ()

        d["pycoqc"] = OrderedDict ()
        d["pycoqc"]["version"] = package_version
        d["pycoqc"]["date"] = datetime.datetime.now().strftime("%d/%m/%y")

        for df, lab in ((self.all_df, "All Reads"), (self.pass_df, "Pass Reads")):
            d[lab] = self._compute_stats(df)
        return d

    def _compute_stats (self, df):
        d = OrderedDict ()

        # run information
        d["run"] = OrderedDict ()
        d["run"]["run_duration"] = float(np.ptp(df["start_time"])/3600)
        d["run"]["active_channels"] = int(df["channel"].nunique())
        d["run"]["runid_number"] = int(df["run_id"].nunique())
        if self.has_barcodes:
            d["run"]["barcodes_number"] = int(df["barcode"].nunique())

        # General run information = OrderedDict ()
        d["basecall"] = OrderedDict ()
        d["basecall"]["reads_number"] = len(df)
        d["basecall"]["bases_number"] = int(df["read_len"].sum())
        d["basecall"]["N50"] = self._compute_N50(df["read_len"])
        d["basecall"]["len_percentiles"] = self._compute_percentiles (df["read_len"])
        d["basecall"]["qual_score_percentiles"] = self._compute_percentiles (df["mean_qscore"])

        if self.has_alignment:
            d["alignment"] = OrderedDict ()
            d["alignment"]["mean_coverage"] = df["align_len"].sum()/self.total_ref_len
            d["alignment"]["reads_number"] = len(df["align_len"].dropna())
            d["alignment"]["bases_number"] = int(df["align_len"].sum())
            d["alignment"]["N50"] = self._compute_N50(df["align_len"])
            d["alignment"]["len_percentiles"] = self._compute_percentiles (df["align_len"])
            if self.has_identity_freq:
                d["alignment"]["identity_freq_percentiles"] = self._compute_percentiles (df["identity_freq"])
                d["alignment"]["insertion_rate"] = df["insertion"].sum()/d["alignment"]["bases_number"]
                d["alignment"]["deletion_rate"] = df["deletion"].sum()/d["alignment"]["bases_number"]
                d["alignment"]["mismatch_rate"] = df["mismatch"].sum()/d["alignment"]["bases_number"]

        return d

    #~~~~~~~SUMMARY METHOD AND HELPER~~~~~~~#
    def summary (self,
        groupby:str = None,
        width:int = None,
        height:int = None,
        plot_title:str="Run summary"):
        """
        Plot an interactive summary table per runid
        * groupby
            Value of field to group the data in the table
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\t\tComputing plot")

        # Prepare all data
        lab1, dd1 = self.__summary_data (df_level="all", groupby=groupby)
        lab2, dd2 = self.__summary_data (df_level="pass", groupby=groupby)

        # Plot initial data
        data = [go.Table(
            header = dd1["header"][0],
            cells = dd1["cells"][0])]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.05, y=1, xanchor='right', yanchor='top', buttons = [
                dict (label=lab1, method='restyle', args=[dd1]),
                dict (label=lab2, method='restyle', args=[dd2])])]

        # Autodefine height depending on the numbers of run_ids
        if not height:
            height=300+(30*self.all_df[groupby].nunique()) if groupby else 300

        # tweak plot layout
        layout = go.Layout (
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"})

        return go.Figure (data=data, layout=layout)

    def barcode_summary (self,
        width:int = None,
        height:int = None,
        plot_title:str="Run summary by barcode"):
        """
        Plot an interactive summary table per barcode (if available)
        """
        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError ("No barcode information available")

        return self.summary(groupby="barcode", width=width, height=height, plot_title=plot_title)

    def run_id_summary (self,
        width:int = None,
        height:int = None,
        plot_title:str="Run summary by Run ID"):
        """
        Plot an interactive summary table per run_id
        """
        if self.all_df["run_id"].nunique() == 1:
            raise pycoQCError ("There is only one run_id")

        return self.summary(groupby="run_id", width=width, height=height, plot_title=plot_title)

    def __summary_data (self, df_level, groupby=None):
        """
        Private function preparing data for summary
        """
        self.logger.debug ("\t\tPreparing data for {} reads".format(df_level))

        # Get data
        df = self.pass_df if df_level == "pass" else self.all_df

        # Group by barcode if required. Otherwise fall back to runid
        cells = []
        if groupby:
            header=[groupby.capitalize(), "Reads", "Bases", "Med Read Length", "N50 Length", "Med Read Quality", "Active Channels", "Run Duration (h)"]
            for id, sdf in df.groupby (groupby):
                cells.append (self.__df_to_cell(sdf, id))
        else:
            header=["Reads", "Bases", "Med Read Length", "N50 Length", "Med Read Quality", "Active Channels", "Run Duration (h)"]
            cells.append (self.__df_to_cell(df))

        # Transpose list of list
        cells = [*zip(*cells)]

        data_dict = dict (
            header = [{"values":header, "align":"center", "fill":{"color":"grey"}, "font":{"size":14, "color":"white"}, "height":40}],
            cells  = [{"values":cells, "align":"center", "fill":{"color":"whitesmoke"}, "font":{"size":12}, "height":30}])

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict)

    def __df_to_cell (self, df, id=None):
        """Extract information from sub-dataframes and return a list of values"""
        l = []
        if id:
            l.append (id)
        l.append ("{:,}".format(len(df)))
        l.append ("{:,}".format(df["read_len"].sum()))
        l.append ("{:,.2f}".format(df["read_len"].median()))
        l.append ("{:,.2f}".format(self._compute_N50(df["read_len"])))
        l.append ("{:,.2f}".format(df["mean_qscore"].median()))
        l.append ("{:,}".format(df["channel"].nunique()))
        l.append ("{:,.2f}".format(np.ptp(df["start_time"])/3600))
        return l

    #~~~~~~~1D DISTRIBUTION METHODS AND HELPER~~~~~~~#
    def read_len_1D (self,
        color:str="lightsteelblue",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Basecalled reads length"):
        """
        Plot a distribution of read length (log scale)
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        fig = self.__1D_density_plot (
            field_name = "read_len",
            plot_title = plot_title,
            x_lab = "Basecalled length",
            color = color,
            x_scale = "log",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def read_qual_1D (self,
        color:str="salmon",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Basecalled reads PHRED quality"):
        """
        Plot a distribution of quality scores
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        fig = self.__1D_density_plot (
            field_name = "mean_qscore",
            plot_title = plot_title,
            x_lab = "Read quality scores",
            color = color,
            x_scale = "linear",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def align_len_1D (self,
        color:str="mediumseagreen",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Aligned reads length"):
        """
        Plot a distribution of read length (log scale)
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_alignment:
            raise pycoQCError ("No Alignment information available")

        fig = self.__1D_density_plot (
            field_name = "align_len",
            plot_title = plot_title,
            x_lab = "Alignment length",
            color = color,
            x_scale = "log",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def identity_freq_1D (self,
        color:str="sandybrown",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Aligned reads identity"):
        """
        Plot a distribution of alignments identity
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_identity_freq:
            raise pycoQCError ("No identity frequency information available")

        fig = self.__1D_density_plot (
            field_name = "identity_freq",
            plot_title = plot_title,
            x_lab = "Identity frequency",
            color = color,
            x_scale = "linear",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def __1D_density_plot (self, field_name, plot_title, x_lab, color, x_scale, nbins, smooth_sigma, width, height):
        """Private function generating density plots for all 1D distribution functions"""
        self.logger.info ("\t\tComputing plot")

        # Prepare all data
        lab1, dd1, ld1 = self.__1D_density_data ("all", field_name, x_scale, nbins, smooth_sigma)
        lab2, dd2, ld2 = self.__1D_density_data ("pass" ,field_name, x_scale, nbins, smooth_sigma)

        # Plot initial data
        line_style = {'color':'gray','width':1,'dash': 'dot'}
        data = [
            go.Scatter (x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], fill='tozeroy', fillcolor=color, mode='none', showlegend=True),
            go.Scatter (x=dd1["x"][1], y=dd1["y"][1], name=dd1["name"][1], text=dd1["text"][1], mode="lines+text", hoverinfo="skip", textposition='top center', line= line_style),
            go.Scatter (x=dd1["x"][2], y=dd1["y"][2], name=dd1["name"][2], text=dd1["text"][2], mode="lines+text", hoverinfo="skip", textposition='top center', line= line_style),
            go.Scatter (x=dd1["x"][3], y=dd1["y"][3], name=dd1["name"][3], text=dd1["text"][3], mode="lines+text", hoverinfo="skip", textposition='top center', line= line_style),
            go.Scatter (x=dd1["x"][4], y=dd1["y"][4], name=dd1["name"][4], text=dd1["text"][4], mode="lines+text", hoverinfo="skip", textposition='top center', line= line_style),
            go.Scatter (x=dd1["x"][5], y=dd1["y"][5], name=dd1["name"][5], text=dd1["text"][5], mode="lines+text", hoverinfo="skip", textposition='top center', line= line_style)]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.2, y=0, xanchor='left', yanchor='bottom', buttons = [
                dict (label=lab1, method='update', args=[dd1, ld1]),
                dict (label=lab2, method='update', args=[dd2, ld2])])]

        # tweak plot layout
        layout = go.Layout (
            hovermode = "closest",
            plot_bgcolor="whitesmoke",
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"},
            xaxis = {"title":x_lab, "type":x_scale, "zeroline":False, "showline":True},
            yaxis = {"title":"Read density", "zeroline":False, "showline":True, "fixedrange":True, "range":ld1["yaxis.range"]})

        return go.Figure (data=data, layout=layout)

    def __1D_density_data (self, df_level, field_name, x_scale, nbins, smooth_sigma):
        """Private function preparing data for reads_1D"""

        self.logger.debug ("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        # Get data
        df = self.pass_sample_df if df_level=="pass" else self.all_sample_df
        data = df[field_name].dropna().values

        # Count each categories in log or linear space
        min = np.nanmin(data)
        max = np.nanmax(data)
        if x_scale == "log":
            count_y, bins = np.histogram (a=data, bins=np.logspace (np.log10(min), np.log10(max)+0.1, nbins))
        elif x_scale == "linear":
            count_y, bins = np.histogram (a=data, bins= np.linspace (min, max, nbins))

        # Remove last bin from labels
        count_x = bins[1:]

        # Smooth results with a savgol filter
        if smooth_sigma:
            count_y = gaussian_filter1d (count_y, sigma=smooth_sigma)

        # Get percentiles percentiles
        stat = np.percentile (data, [10,25,50,75,90])
        y_max = count_y.max()

        data_dict = dict (
            x = [count_x, [stat[0],stat[0]], [stat[1],stat[1]], [stat[2],stat[2]], [stat[3],stat[3]], [stat[4],stat[4]]],
            y = [count_y, [0,y_max], [0,y_max], [0,y_max], [0,y_max], [0,y_max]],
            name = ["Density", "10%", "25%", "Median", "75%", "90%"],
            text = ["",
                ["", "10%<br>{:,.2f}".format(stat[0])],
                ["", "25%<br>{:,.2f}".format(stat[1])],
                ["", "Median<br>{:,.2f}".format(stat[2])],
                ["", "75%<br>{:,.2f}".format(stat[3])],
                ["", "90%<br>{:,.2f}".format(stat[4])]])

        # Make layout dict = Off set for labels on top
        layout_dict = {"yaxis.range": [0, y_max+y_max/6]}

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict, layout_dict)

    #~~~~~~~2D DISTRIBUTION METHOD AND HELPER~~~~~~~#
    def read_len_read_qual_2D (self,
        colorscale = [[0.0,'rgba(255,255,255,0)'], [0.1,'rgba(255,150,0,0)'], [0.25,'rgb(255,100,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(70,0,0)']],
        x_nbins:int=200,
        y_nbins:int=100,
        smooth_sigma:float=2,
        width:int=None,
        height:int=600,
        plot_title:str="Basecalled reads length vs reads PHRED quality"):
        """
        Plot a 2D distribution of quality scores vs length of the reads
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in (x axis)
        * y_nbins
            Number of bins to divide the read quality values in (y axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        fig = self.__2D_density_plot (
            x_field_name = "read_len",
            y_field_name = "mean_qscore",
            x_lab = "Basecalled length",
            y_lab = "PHRED quality scores",
            x_scale = "log",
            y_scale = "linear",
            x_nbins = x_nbins,
            y_nbins = y_nbins,
            colorscale = colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title = plot_title)
        return fig

    def read_len_align_len_2D (self,
        colorscale = [[0.0,'rgba(255,255,255,0)'], [0.1,'rgba(255,150,0,0)'], [0.25,'rgb(255,100,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(70,0,0)']],
        x_nbins:int=200,
        y_nbins:int=100,
        smooth_sigma:float=1,
        width:int=None,
        height:int=600,
        plot_title:str="Basecalled reads length vs alignments length"):
        """
        Plot a 2D distribution of length of the reads vs length of the alignments
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in (x axis)
        * y_nbins
            Number of bins to divide the read quality values in (y axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_alignment:
            raise pycoQCError ("No Alignment information available")

        fig = self.__2D_density_plot (
            x_field_name = "read_len",
            y_field_name = "align_len",
            x_lab = "Basecalled length",
            y_lab = "Alignment length",
            x_scale = "log",
            y_scale = "log",
            x_nbins = x_nbins,
            y_nbins = y_nbins,
            colorscale = colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title = plot_title)
        return fig

    def align_len_identity_freq_2D (self,
        colorscale = [[0.0,'rgba(255,255,255,0)'], [0.1,'rgba(255,150,0,0)'], [0.25,'rgb(255,100,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(70,0,0)']],
        x_nbins:int=200,
        y_nbins:int=100,
        smooth_sigma:float=2,
        width:int=None,
        height:int=600,
        plot_title:str="Aligned reads length vs alignments identity"):
        """
        Plot a 2D distribution of alignments length vs alignments identity
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in (x axis)
        * y_nbins
            Number of bins to divide the read quality values in (y axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_identity_freq:
            raise pycoQCError ("No identity frequency information available")

        fig = self.__2D_density_plot (
            x_field_name = "align_len",
            y_field_name = "identity_freq",
            x_lab = "Alignment length",
            y_lab = "Identity frequency",
            x_scale = "log",
            y_scale = "linear",
            x_nbins = x_nbins,
            y_nbins = y_nbins,
            colorscale = colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title = plot_title)
        return fig

    def read_qual_identity_freq_2D (self,
        colorscale = [[0.0,'rgba(255,255,255,0)'], [0.1,'rgba(255,150,0,0)'], [0.25,'rgb(255,100,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(70,0,0)']],
        x_nbins:int=200,
        y_nbins:int=100,
        smooth_sigma:float=1,
        width:int=None,
        height:int=600,
        plot_title:str="Reads PHRED quality vs alignments identity"):
        """
        Plot a 2D distribution of read quality vs alignments identity
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in (x axis)
        * y_nbins
            Number of bins to divide the read quality values in (y axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_identity_freq:
            raise pycoQCError ("No identity frequency information available")

        fig = self.__2D_density_plot (
            x_field_name = "mean_qscore",
            y_field_name = "identity_freq",
            x_lab = "PHRED quality",
            y_lab = "Identity frequency",
            x_scale = "linear",
            y_scale = "linear",
            x_nbins = x_nbins,
            y_nbins = y_nbins,
            colorscale = colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title = plot_title)
        return fig

    def __2D_density_plot (self, x_field_name, y_field_name, x_lab, y_lab, x_scale, y_scale, x_nbins, y_nbins, colorscale, smooth_sigma, width, height, plot_title):
        """Private function generating density plots for all 2D distribution functions"""
        self.logger.info ("\t\tComputing plot")

        # Prepare all data
        lab1, dd1 = self.__2D_density_data ("all", x_field_name, y_field_name, x_nbins, y_nbins, x_scale, y_scale, smooth_sigma)
        lab2, dd2 = self.__2D_density_data ("pass", x_field_name, y_field_name, x_nbins, y_nbins, x_scale, y_scale, smooth_sigma)

        # Plot initial data
        data = [
            go.Contour (x=dd1["x"][0], y=dd1["y"][0], z=dd1["z"][0], contours=dd1["contours"][0],
                name="Density", hoverinfo="name+x+y", colorscale=colorscale, showlegend=True, connectgaps=True, line={"width":0}),
            go.Scatter (x=dd1["x"][1], y=dd1["y"][1],
                mode='markers', name='Median', hoverinfo="name+x+y", marker={"size":12,"color":'black', "symbol":"x"})]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.2, y=0, xanchor='left', yanchor='bottom', buttons = [
                dict (label=lab1, method='restyle', args=[dd1]),
                dict (label=lab2, method='restyle', args=[dd2])])]

        # tweak plot layout
        layout = go.Layout (
            hovermode = "closest",
            plot_bgcolor="whitesmoke",
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"},
            xaxis = {"title":x_lab, "showgrid":True, "zeroline":False, "showline":True, "type":x_scale},
            yaxis = {"title":y_lab, "showgrid":True, "zeroline":False, "showline":True, "type":y_scale})

        return go.Figure (data=data, layout=layout)

    def __2D_density_data (self, df_level, x_field_name, y_field_name, x_nbins, y_nbins, x_scale, y_scale, smooth_sigma):
        """ Private function preparing data for 2D_density_plot """

        self.logger.debug ("\t\tPreparing data for {} reads".format(df_level))

        # Extract data field from df
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        df = df[[x_field_name, y_field_name]].dropna()

        # Prepare data for x
        x_data = df[x_field_name].values
        x_min, x_med, x_max = np.percentile (x_data, (0, 50, 100))
        if x_scale == "log":
            x_bins = np.logspace (start=np.log10((x_min)), stop=np.log10(x_max)+0.1, num=x_nbins, base=10)
        else:
            x_bins = np.linspace (start=x_min, stop=x_max, num=x_nbins)

        # Prepare data for y
        y_data = df[y_field_name].values
        y_min, y_med, y_max = np.percentile (y_data, (0, 50, 100))
        if y_scale == "log":
            y_bins = np.logspace (start=np.log10((y_min)), stop=np.log10(y_max)+0.1, num=y_nbins, base=10)
        else:
            y_bins = np.linspace (start=y_min, stop=y_max, num=y_nbins)

        # Compute 2D histogram
        z, y, x = np.histogram2d (x=y_data, y=x_data, bins=[y_bins, x_bins])
        if smooth_sigma:
            z = gaussian_filter(z, sigma=smooth_sigma)
        z_min, z_max = np.percentile (z, (0, 100))

        # Extract label and values
        data_dict = dict (
            x = [x, [x_med]], y = [y, [y_med]], z = [z, None],
            contours = [dict(start=z_min, end=z_max, size=(z_max-z_min)/15),None])

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict)

    #~~~~~~~OUTPUT_OVER_TIME METHODS AND HELPER~~~~~~~#
    def output_over_time (self,
        cumulative_color:str="rgb(204,226,255)",
        interval_color:str="rgb(102,168,255)",
        time_bins:int=500,
        width:int=None,
        height:int=500,
        plot_title:str="Output over experiment time"):
        """
        Plot a yield over time
        * cumulative_color
            Color of cumulative yield area (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * interval_color
            Color of interval yield line (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\t\tComputing plot")

        # Prepare all data
        lab1, dd1, ld1 = args=self.__output_over_time_data (df_level="all", count_level="reads", time_bins=time_bins)
        lab2, dd2, ld2 = args=self.__output_over_time_data (df_level="pass", count_level="reads", time_bins=time_bins)
        lab3, dd3, ld3 = args=self.__output_over_time_data (df_level="all", count_level="bases", time_bins=time_bins)
        lab4, dd4, ld4 = args=self.__output_over_time_data (df_level="pass", count_level="bases", time_bins=time_bins)

        # Plot initial data
        line_style = {'color':'gray','width':1,'dash':'dot'}
        data = [
            go.Scatter (x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], fill='tozeroy', fillcolor=cumulative_color, mode='none'),
            go.Scatter (x=dd1["x"][1], y=dd1["y"][1], name=dd1["name"][1], mode='lines', line={'color':interval_color,'width':2}),
            go.Scatter (x=dd1["x"][2], y=dd1["y"][2], name=dd1["name"][2], text=dd1["text"][2], mode="lines+text", hoverinfo="skip", textposition='top center', line=line_style),
            go.Scatter (x=dd1["x"][3], y=dd1["y"][3], name=dd1["name"][3], text=dd1["text"][3], mode="lines+text", hoverinfo="skip", textposition='top center', line=line_style),
            go.Scatter (x=dd1["x"][4], y=dd1["y"][4], name=dd1["name"][4], text=dd1["text"][4], mode="lines+text", hoverinfo="skip", textposition='top center', line=line_style),
            go.Scatter (x=dd1["x"][5], y=dd1["y"][5], name=dd1["name"][5], text=dd1["text"][5], mode="lines+text", hoverinfo="skip", textposition='top center', line=line_style),
            go.Scatter (x=dd1["x"][6], y=dd1["y"][6], name=dd1["name"][6], text=dd1["text"][6], mode="lines+text", hoverinfo="skip", textposition='top center', line=line_style)]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.06, y=0, xanchor='right', yanchor='bottom', buttons = [
                dict (label=lab1,  method='update', args=[dd1, ld1]),
                dict (label=lab2, method='update', args=[dd2, ld2]),
                dict (label=lab3,  method='update', args=[dd3, ld3]),
                dict (label=lab4, method='update', args=[dd4, ld4])])]

        # tweak plot layout
        layout = go.Layout (
            plot_bgcolor="whitesmoke",
            width = width,
            height = height,
            updatemenus = updatemenus,
            legend = {"x":-0.05, "y":1,"xanchor":'right',"yanchor":'top'},
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"},
            xaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":True},
            yaxis = {"title":"Count", "zeroline":False, "showline":True, "fixedrange":True, "range":ld1["yaxis.range"]})

        return go.Figure (data=data, layout=layout)

    def __output_over_time_data (self, df_level, count_level, time_bins=500):
        """Private function preparing data for output_over_time"""
        self.logger.debug ("\t\tPreparing data for {} {}".format(df_level, count_level))

        # Get data and scaling factor
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        sf = self.pass_scaling_factor if df_level == "pass" else self.all_scaling_factor

        # Bin data in categories
        t = (df["start_time"]/3600).values
        x = np.linspace (t.min(), t.max(), num=time_bins)
        t = np.digitize (t, bins=x, right=True)

        # Count reads or bases per categories
        if count_level == "reads":
            y = np.bincount(t)
        elif count_level == "bases":
            y = np.bincount(t, weights=df["read_len"].values)

        # Scale counts in case of downsampling
        y = y*sf

        # Transform to cummulative distribution
        y_cum = np.cumsum(y)
        y_cum_max = y_cum[-1]

        # Smooth and rescale interval trace
        y = gaussian_filter1d (y, sigma=1)
        y = y*y_cum_max/y.max()

        # Find percentages of data generated
        lab_text = []
        lab_name = []
        lab_x = []
        for lab in (50, 75, 90, 99, 100):
            val = y_cum_max*lab/100
            idx = (np.abs(y_cum-val)).argmin()
            lab_text.append(["", '{}%<br>{}h<br>{:,} {}'.format(lab, round(x[idx],2), int(y_cum[idx]), count_level)])
            lab_x.append ([x[idx], x[idx]])
            lab_name.append ("{}%".format(lab))

        # make data dict
        data_dict = dict(
            x = [x, x]+lab_x,
            y = [y_cum, y, [0,y_cum_max], [0,y_cum_max], [0,y_cum_max], [0,y_cum_max], [0,y_cum_max]],
            name = ["Cumulative", "Interval"]+lab_name,
            text = ["", ""]+lab_text)

        # Make layout dict = offset for labels on top
        layout_dict = {"yaxis.range": [0, y_cum_max+y_cum_max/6]}

        label = "{} {}".format(df_level.capitalize(), count_level.capitalize())
        return (label, data_dict, layout_dict)

    #~~~~~~~QUAL_OVER_TIME METHODS AND HELPER~~~~~~~#
    def read_len_over_time (self,
        median_color:str="rgb(102,168,255)",
        quartile_color:str="rgb(153,197,255)",
        extreme_color:str="rgba(153,197,255,0.5)",
        smooth_sigma:float=1,
        time_bins:int=500,
        width:int=None,
        height:int=500,
        plot_title:str="Read length over experiment time"):
        """
        Plot a read length over time
        * median_color
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        fig = self.__over_time_plot (
            field_name = "read_len",
            plot_title = plot_title,
            y_lab = "Alignment length",
            y_scale = "log",
            median_color = median_color,
            quartile_color = quartile_color,
            extreme_color = extreme_color,
            smooth_sigma = smooth_sigma,
            time_bins = time_bins,
            width = width,
            height = height)
        return fig

        return go.Figure (data=data, layout=layout)

    def read_qual_over_time (self,
        median_color:str="rgb(250,128,114)",
        quartile_color:str="rgb(250,170,160)",
        extreme_color:str="rgba(250,170,160,0.5)",
        smooth_sigma:float=1,
        time_bins:int=500,
        width:int=None,
        height:int=500,
        plot_title:str="Read quality over experiment time"):
        """
        Plot a mean quality over time
        * median_color
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        fig = self.__over_time_plot (
            field_name = "mean_qscore",
            plot_title = plot_title,
            y_lab = "Mean read PHRED quality",
            y_scale = "linear",
            median_color = median_color,
            quartile_color = quartile_color,
            extreme_color = extreme_color,
            smooth_sigma = smooth_sigma,
            time_bins = time_bins,
            width = width,
            height = height)
        return fig

    def align_len_over_time (self,
        median_color:str="rgb(102,168,255)",
        quartile_color:str="rgb(153,197,255)",
        extreme_color:str="rgba(153,197,255,0.5)",
        smooth_sigma:float=1,
        time_bins:int=500,
        width:int=None,
        height:int=500,
        plot_title:str="Aligned reads length over experiment time"):
        """
        Plot a aligned reads length over time
        * median_color
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_alignment:
            raise pycoQCError ("No Alignment information available")

        fig = self.__over_time_plot (
            field_name = "align_len",
            plot_title = plot_title,
            y_lab = "Aligned reads length",
            y_scale = "log",
            median_color = median_color,
            quartile_color = quartile_color,
            extreme_color = extreme_color,
            smooth_sigma = smooth_sigma,
            time_bins = time_bins,
            width = width,
            height = height)
        return fig

    def identity_freq_over_time (self,
        median_color:str="rgb(250,128,114)",
        quartile_color:str="rgb(250,170,160)",
        extreme_color:str="rgba(250,170,160,0.5)",
        smooth_sigma:float=1,
        time_bins:int=500,
        width:int=None,
        height:int=500,
        plot_title:str="Aligned reads identity over experiment time"):
        """
        Plot the alignment identity scores over time
        * median_color
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        # Verify that alignemnt information are available
        if not self.has_identity_freq:
            raise pycoQCError ("No identity frequency information available")

        fig = self.__over_time_plot (
            field_name = "identity_freq",
            plot_title = plot_title,
            y_lab = "Identity frequency",
            y_scale = "linear",
            median_color = median_color,
            quartile_color = quartile_color,
            extreme_color = extreme_color,
            smooth_sigma = smooth_sigma,
            time_bins = time_bins,
            width = width,
            height = height)
        return fig

    def __over_time_plot (self,
        field_name,
        plot_title,
        y_lab,
        y_scale,
        median_color,
        quartile_color,
        extreme_color,
        smooth_sigma,
        time_bins,
        width,
        height):
        """Private function generating density plots for all over_time functions"""
        self.logger.info ("\t\tComputing plot")

        lab1, dd1 = self.__over_time_data (df_level="all", field_name=field_name, smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab2, dd2 = self.__over_time_data (df_level="pass", field_name=field_name, smooth_sigma=smooth_sigma, time_bins=time_bins)

        # Plot initial data
        data= [
            go.Scatter(x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], mode="lines", line={"color":extreme_color}, connectgaps=True, legendgroup="Extreme"),
            go.Scatter(x=dd1["x"][1], y=dd1["y"][1], name=dd1["name"][1], mode="lines", fill="tonexty", line={"color":extreme_color}, connectgaps=True, legendgroup="Extreme"),
            go.Scatter(x=dd1["x"][2], y=dd1["y"][2], name=dd1["name"][2], mode="lines", line={"color":quartile_color}, connectgaps=True, legendgroup="Quartiles"),
            go.Scatter(x=dd1["x"][3], y=dd1["y"][3], name=dd1["name"][3], mode="lines", fill="tonexty", line={"color":quartile_color}, connectgaps=True, legendgroup="Quartiles"),
            go.Scatter(x=dd1["x"][4], y=dd1["y"][4], name=dd1["name"][4], mode="lines", line={"color":median_color}, connectgaps=True)]

        # Create update buttons
        updatemenus = [
            go.layout.Updatemenu (type="buttons", active=0, x=-0.07, y=0, xanchor='right', yanchor='bottom', buttons = [
                go.layout.updatemenu.Button (
                    label=lab1, method='restyle', args=[dd1]),
                go.layout.updatemenu.Button (
                    label=lab2, method='restyle', args=[dd2])])]

        # tweak plot layout
        layout = go.Layout (
            plot_bgcolor="whitesmoke",
            width = width,
            height = height,
            updatemenus = updatemenus,
            legend = {"x":-0.07, "y":1,"xanchor":'right',"yanchor":'top'},
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"},
            yaxis = {"title":y_lab, "zeroline":False, "type":y_scale, "showline":True, "rangemode":'nonnegative', "fixedrange":True},
            xaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":True, "rangemode":'nonnegative'})

        return go.Figure (data=data, layout=layout)

    def __over_time_data (self, df_level, field_name="read_len", smooth_sigma=1.5, time_bins=500):
        """Private function preparing data for qual_over_time"""
        self.logger.debug ("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        # get data
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        data = df[field_name].dropna().values

        # Bin data in categories
        t = (df["start_time"]/3600).values
        x = np.linspace (t.min(), t.max(), num=time_bins)
        t = np.digitize (t, bins=x, right=True)

        # List quality value per categories
        bin_dict = defaultdict (list)
        for bin_idx, val in zip (t, data) :
            bin = x[bin_idx]
            bin_dict[bin].append(val)

        # Aggregate values per category
        val_name = ["Min", "Max", "25%", "75%", "Median"]
        stat_dict = defaultdict(list)
        for bin in x:
            if bin in bin_dict:
                p = np.percentile (bin_dict[bin], [0, 100, 25, 75, 50])
            else:
                p = [np.nan,np.nan,np.nan,np.nan,np.nan]
            for val, stat in zip (val_name, p):
                stat_dict[val].append(stat)

        # Values smoothing
        if smooth_sigma:
            for val in val_name:
                stat_dict [val] = gaussian_filter1d (stat_dict [val], sigma=smooth_sigma)

        # make data dict
        data_dict = dict(
            x = [x,x,x,x,x],
            y = [stat_dict["Min"], stat_dict["Max"], stat_dict["25%"], stat_dict["75%"], stat_dict["Median"]],
            name = val_name)

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict)

    #~~~~~~~BARCODE_COUNT METHODS AND HELPER~~~~~~~#
    def barcode_counts (self,
        colors:list=["#f8bc9c", "#f6e9a1", "#f5f8f2", "#92d9f5", "#4f97ba"],
        width:int= None,
        height:int=500,
        plot_title:str="Percentage of reads per barcode"):
        """
        Plot a mean quality over time
        * colors
            List of colors (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError ("No barcode information available")
        self.logger.info ("\t\tComputing plot")

        # Prepare all data
        lab1, dd1 = self.__barcode_counts_data (df_level="all")
        lab2, dd2 = self.__barcode_counts_data (df_level="pass")

        # Plot initial data
        data= [go.Pie (labels=dd1["labels"][0] , values=dd1["values"][0] , sort=False, marker=dict(colors=colors))]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.2, y=0, xanchor='left', yanchor='bottom', buttons = [
                dict (label=lab1, method='restyle', args=[dd1]),
                dict (label=lab2, method='restyle', args=[dd2])])]

        # tweak plot layout
        layout = go.Layout (
            plot_bgcolor="whitesmoke",
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"})

        return go.Figure (data=data, layout=layout)

    def __barcode_counts_data (self, df_level):
        """Private function preparing data for barcode_counts"""
        self.logger.debug ("\t\tPreparing data for {} reads".format(df_level))

        # get data
        df = self.pass_df if df_level == "pass" else self.all_df
        counts = df["barcode"].value_counts()
        counts = counts.sort_index()

        # Extract label and values
        data_dict = dict (
            labels = [counts.index],
            values = [counts.values])

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict)

    #~~~~~~~BARCODE_COUNT METHODS AND HELPER~~~~~~~# ################################################################################ ADD TABLE AS IN ALIGNMENTS
    def channels_activity (self,
        colorscale:list = [[0.0,'rgba(255,255,255,0)'], [0.01,'rgb(255,255,200)'], [0.25,'rgb(255,200,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(0,0,0)']],
        smooth_sigma:float=1,
        time_bins:int=100,
        width:int=None,
        height:int=600,
        plot_title:str="Output per channel over experiment time"):
        """
        Plot a yield over time
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (y axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\t\tComputing plot")

        # Define maximal number of channels
        n_channels = 3000 if self.is_promethion else 512

        # Prepare all data
        lab1, dd1 = args=self.__channels_activity_data(df_level="all", count_level="reads", n_channels=n_channels, smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab2, dd2 = args=self.__channels_activity_data(df_level="pass", count_level="reads", n_channels=n_channels, smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab3, dd3 = args=self.__channels_activity_data(df_level="all", count_level="bases", n_channels=n_channels, smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab4, dd4 = args=self.__channels_activity_data(df_level="pass", count_level="bases", n_channels=n_channels, smooth_sigma=smooth_sigma, time_bins=time_bins)

        # Plot initial data
        data = [go.Heatmap(x=dd1["x"][0], y=dd1["y"][0], z=dd1["z"][0], xgap=0.5, colorscale=colorscale, hoverinfo="x+y+z")]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.06, y=0, xanchor='right', yanchor='bottom', buttons = [
                dict (label=lab1, method='restyle', args=[dd1]),
                dict (label=lab2, method='restyle', args=[dd2]),
                dict (label=lab3, method='restyle', args=[dd3]),
                dict (label=lab4, method='restyle', args=[dd4])])]

        # tweak plot layout
        layout = go.Layout (
            plot_bgcolor="whitesmoke",
            width = width,
            height = height,
            updatemenus = updatemenus,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"},
            xaxis = {"title":"Channel id", "zeroline":False, "showline":False, "nticks":20, "showgrid":False},
            yaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":False, "hoverformat":".2f", "fixedrange":True})

        return go.Figure (data=data, layout=layout)

    def __channels_activity_data (self, df_level, count_level="bases", n_channels=512, smooth_sigma=2, time_bins=150):
        """Private function preparing data for channels_activity"""
        self.logger.debug ("\t\tPreparing data for {} {}".format(df_level, count_level))

        # Get data and scaling factor
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        sf = self.pass_scaling_factor if df_level == "pass" else self.all_scaling_factor

        # Bin data in categories
        t = (df["start_time"]/3600).values
        bins = np.linspace (t.min(), t.max(), num=time_bins)
        t = np.digitize (t, bins=bins, right=True)

        # Count values per categories
        z = np.ones((len(bins), n_channels), dtype=np.int)
        if count_level == "bases":
            for t_idx, channel, n_bases in zip(t, df["channel"], df["read_len"]):
                z[t_idx][channel-1]+=n_bases
        elif count_level == "reads":
            for t_idx, channel in zip(t, df["channel"]):
                try:
                    z[t_idx][channel-1]+=1
                except IndexError:
                    print (t_idx, channel)
                    raise
        # Scale counts in case of downsampling
        z=z*sf

        # Time series smoothing
        if smooth_sigma:
            z = gaussian_filter1d (z.astype(np.float32), sigma=smooth_sigma, axis=0)

        # Define x and y axis
        x = ["c {}".format(i) for i in range(1, n_channels+1)]
        y = bins[1:]

        # Make data dict
        data_dict = dict (x=[x], y=[y], z=[z])

        label = "{} {}".format(df_level.capitalize(), count_level.capitalize())
        return (label, data_dict)

    #~~~~~~~ALIGNMENT_SUMMARY METHOD~~~~~~~#
    def alignment_summary (self,
        colors:list=["#f44f39","#fc8161","#fcaf94","#828282"],
        width:int= None,
        height:int=500,
        plot_title:str="Summary of reads alignment"):
        """
        Plot a basic alignment summary
        * colors
            List of colors (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_alignment:
            raise pycoQCError ("No Alignment information available")
        self.logger.info ("\t\tComputing plot")

        df = self.alignments_df
        # Create empty multiplot figure
        fig = make_subplots(rows=1, cols=2, column_widths=[0.4, 0.6], specs=[[{"type": "table"},{"type": "pie"}]])

        # plot Table
        data = go.Table(
            columnwidth = [3,2,2],
            header = {"values":list(df.columns), "align":"center", "fill_color":"grey", "font_size":14, "font_color":"white", "height":40},
            cells = {"values":df.values.T , "align":"center", "fill_color":"whitesmoke", "font_size":12, "height":30})
        fig.add_trace (data, row=1, col=1)

        # plot Pie plot
        data = go.Pie (
            labels=df["Alignments"],
            values=df["Counts"],
            sort=False,
            marker={"colors":colors},
            name="Pie plot",
            textinfo='label+percent')
        fig.add_trace (data, row=1, col=2)

        # Change the layout
        fig.update_layout(
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"})

        return fig

    #~~~~~~~ALIGNMENT RATE METHOD AND HELPER~~~~~~~#
    def alignment_rate (self,
            colors:list=["#fcaf94","#828282","#fc8161","#828282","#f44f39","#d52221","#828282","#828282","#828282","#828282"],
            width:int=None,
            height:int=600,
            plot_title:str="Bases alignment rate"):
        """
        Plot a basic alignment summary
        * colors
            List of colors (hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        # Verify that alignemnt information are available
        if not self.has_identity_freq:
            raise pycoQCError ("No identity frequency information available")
        self.logger.info ("\t\tComputing plot")

        # Extract Data
        bc_bases = self.all_df["read_len"].sum()
        s = self.all_df[[ "read_len", "align_len", "insertion", "deletion", "soft_clip", "mismatch"]].dropna().sum()
        total_error = s["insertion"]+s["deletion"]+s["mismatch"]
        matching = s["align_len"]-total_error
        unmapped = bc_bases-s["read_len"]
        ct = namedtuple("ct", ["Bases","Counts","Total_freq","Aligned_freq"])
        l = [
            ct("Basecalled", bc_bases, 1, 1),
            ct("Unmapped reads", unmapped, unmapped/bc_bases, 1),
            ct("Mapped reads", s["read_len"], s["read_len"]/bc_bases, 1),
            ct("Softclip", s["soft_clip"], s["soft_clip"]/bc_bases, 1),
            ct("Aligned", s["align_len"], s["align_len"]/bc_bases, 1),
            ct("Matching", matching, matching/bc_bases, matching/s["align_len"]),
            ct("Non-matching", total_error, total_error/bc_bases, total_error/s["align_len"]),
            ct("Insertions", s["insertion"], s["insertion"]/bc_bases, s["insertion"]/s["align_len"]),
            ct("Deletions", s["deletion"], s["deletion"]/bc_bases, s["deletion"]/s["align_len"]),
            ct("Mismatches", s["mismatch"], s["mismatch"]/bc_bases, s["mismatch"]/s["align_len"])]

        # Cast to df and compute percentage
        df = pd.DataFrame(l)

        # plot Table
        data1 = go.Table(
            columnwidth = [3,2,2,2],
            header = {"values":["Bases","Bases Count","% Total","% Aligned"], "align":"center", "fill_color":["grey"], "font_size":14, "font_color":"white", "height":40},
            cells = {"values":df.values.T , "format":["", ".3e", ".3p", ".3p"], "align":"center", "fill_color":"whitesmoke", "font_size":12, "height":30})

        data2 = go.Sankey(
            arrangement = "freeform",
            node = go.sankey.Node(
                pad = 20,
                thickness = 30,
                line = {"width":0},
                label = df["Bases"].values,
                x = [0,0.25,0.25,0.5,0.5,1,0.75],
                y = [0,1,0,0.9,0,0,0.8],
                color = colors),
            link = go.sankey.Link(
                source = [0, 0, 2, 2, 4, 4, 6, 6, 6],
                target = [1, 2, 3, 4, 5, 6, 7, 8, 9],
                value = df["Counts"].values[1:],
                hoverinfo = "none"))

        # Create multipanel figure
        fig = make_subplots(rows=1, cols=2, column_widths=[0.4, 0.6], specs=[[{"type": "table"},{"type": "sankey"}]])
        fig.add_trace (data1, row=1, col=1)
        fig.add_trace (data2, row=1, col=2)
        fig.update_layout(
            width = width,
            height = height,
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"})

        return fig

    #~~~~~~~ALIGNMENT COVERAGE METHOD AND HELPER~~~~~~~#
    def alignment_coverage (self,
        nbins:int=500,
        color:str='rgba(70,130,180,0.70)',
        smooth_sigma:int=1,
        width:int= None,
        height:int=500,
        plot_title:str="Coverage overview"):
        """
        Plot coverage over all the references
        * nbins
            Number of bins to divide the coverage into.
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignemnt information are available
        if not self.has_alignment:
            raise pycoQCError ("No Alignment information available")
        self.logger.info ("\t\tComputing plot")

        ref_offset_dict = self._ref_offset(self.ref_len_dict, "left", ret_type="dict")
        df = self.all_df[["ref_id", "ref_start", "ref_end", "align_len"]].dropna()
        steps = self.total_ref_len//nbins
        mean_cov = round(df["align_len"].sum()/self.total_ref_len, 2)

        # Compute coverage by interval
        l = []
        for line in df.itertuples():
            l.append(int(ref_offset_dict[line.ref_id]+line.ref_start))
        l = np.array(l)
        bins = np.arange(0, self.total_ref_len, steps)
        l = np.digitize(l,bins)
        y = np.bincount(l, weights=df["align_len"])/steps

        # Time series smoothing
        if smooth_sigma:
            y = gaussian_filter1d (y, sigma=smooth_sigma)

        # Plot coverage area
        data1 = go.Scatter (
            x=list(range(nbins+1)),
            y=y,
            name="Mean coverage",
            hoveron="points",
            hoverinfo="y",
            fill='tozeroy',
            fillcolor=color,
            mode='none',
            showlegend=True,
            connectgaps=True)

        # Plot mean coverage
        data2 = go.Scatter (
            x=[0,nbins],
            y=[mean_cov,mean_cov],
            name=f"Overall coverage<br>{mean_cov}X",
            mode="lines",
            hoverinfo="skip",
            line= {'color':'gray','width':2,'dash':'dot'})

        updatemenus = [
            dict (type="buttons", x=-0.2, y=0, xanchor='left', yanchor='bottom',  buttons = [
                dict (label="log", method='relayout', args=[{"yaxis":{"title":"Mean Coverage", "type":"log", "zeroline":False, "fixedrange":True}}]),
                dict (label="linear", method='relayout', args=[{"yaxis":{"title":"Mean Coverage", "type":"linear", "zeroline":False, "fixedrange":True}}])])]

        # Add chromosome shading and labels
        x_lab_coord = np.array(self._ref_offset(self.ref_len_dict, coordinates="middle", ret_type="list"))*nbins/self.total_ref_len
        x_lab = list(self.ref_len_dict.keys())
        shapes = []
        x_shape_coord = np.array(self._ref_offset(self.ref_len_dict, coordinates="left", ret_type="list")[1:])*nbins/self.total_ref_len
        for i in range(0, len(self.ref_len_dict)-2, 2):
            shapes.append(go.layout.Shape(type="rect",x0=x_shape_coord[i],x1=x_shape_coord[i+1], y0=0,y1=1, yref="paper", opacity=0.5, layer="below", fillcolor="lightgrey", line_width=0))

        # Tweak plot layout
        layout = go.Layout (
            width = width,
            height = height,
            plot_bgcolor="whitesmoke",
            updatemenus = updatemenus,
            shapes=shapes,
            hovermode = "closest",
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            xaxis = {"zeroline":False, "showline":True, "ticktext":x_lab, "tickvals":x_lab_coord, "tickangle":-45, "showgrid":False},
            yaxis = {"title":"Mean Coverage", "type":"log", "zeroline":False, "fixedrange":True},
            title = {"text":plot_title, "xref":"paper" ,"x":0.5, "xanchor":"center"})

        return go.Figure(data=[data1,data2], layout=layout)

    def _ref_offset (self, rlen, coordinates="left", ret_type="dict"):
        offset = [] if ret_type=="list" else OrderedDict()
        cumsum=0
        for ref, rlen in rlen.items ():
            # Define return val
            if coordinates =="left":
                v = cumsum
            elif coordinates =="middle":
                v = cumsum + rlen/2
            else:
                v = cumsum + rlen
            # Add to appropriate collection
            if ret_type =="list":
                offset.append(v)
            else:
                offset[ref]=v
            cumsum+=rlen
        return offset

    #~~~~~~~PRIVATE METHODS~~~~~~~#
    @staticmethod
    def _compute_percentiles (data):
        return list(np.quantile(data.dropna(), q=np.linspace(0,1,101)))

    @staticmethod
    def _compute_N50 (data):
        data = data.dropna().values
        data.sort()
        half_sum = data.sum()/2
        cum_sum = 0
        for v in data:
            cum_sum += v
            if cum_sum >= half_sum:
                return int(v)
