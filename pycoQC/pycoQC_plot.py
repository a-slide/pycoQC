# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import *
import warnings

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import plotly.graph_objs as go

# Local lib import
from pycoQC.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL SETTINGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set seed for deterministic random sampling
SEED = 42
np.random.RandomState(seed=SEED)

# Silence futurewarnings
warnings.filterwarnings("ignore", category=FutureWarning)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class pycoQC_plot ():

    def __init__ (self,
        df:pd.DataFrame,
        min_pass_qual:int=7,
        sample:int=100000,
        verbose:bool=False,
        quiet:bool=False):
        """
        * df
            pandas dataframe obtained with pycoQC_parse
        * min_pass_qual
            Minimum quality to consider a read as 'pass
        * sample
            If not None a n number of reads will be randomly selected instead of the entire dataset for ploting function (deterministic sampling)
        """

        # Set logging level
        self.logger = get_logger (name=__name__, verbose=verbose, quiet=quiet)
        self.logger.warning ("Loading plotting interface")

        # Save args to self values
        self.min_pass_qual = min_pass_qual
        self.sample = sample

        # Save df wiews and compute scaling factors
        self.all_df = df
        if sample and len(self.all_df)>sample:
            self.all_sample_df = self.all_df.sample(n=sample, random_state=SEED)
            self.all_scaling_factor = len(self.all_df)/sample
        else:
            self.all_sample_df = self.all_df
            self.all_scaling_factor = 1

        self.pass_df = df[df["mean_qscore"]>=min_pass_qual]
        if sample and len(self.pass_df)>sample:
            self.pass_sample_df = self.pass_df.sample(n=sample, random_state=SEED)
            self.pass_scaling_factor = len(self.pass_df)/sample
        else:
            self.pass_sample_df = self.pass_df
            self.pass_scaling_factor = 1

    def __str__(self):
        m = ""
        m+= "\tBarcode: {}\n".format(self.has_barcodes)
        m+= "\tAlignment: {}\n".format(self.has_alignment)
        m+= "\tPromethion: {}\n".format(self.is_promethion)
        m+= "\tAll reads: {:,}\n".format(len(self.all_df))
        m+= "\tAll bases: {:,}\n".format(int(self.all_df["num_bases"].sum()))
        m+= "\tAll median read length: {:,}\n".format(np.median(self.all_df["num_bases"]))
        m+= "\tPass reads: {:,}\n".format(len(self.pass_df))
        m+= "\tPass bases: {:,}\n".format(int(self.pass_df["num_bases"].sum()))
        m+= "\tPass median read length: {:,}\n".format(np.median(self.pass_df["num_bases"]))
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
    def is_promethion (self):
        return self.all_df["channel"].max() > 512

    #~~~~~~~SUMMARY_STATS_DICT METHOD AND HELPER~~~~~~~#

    def summary_stats_dict (self,
        barcode_split:bool=False,
        run_id_split:bool=False):
        """
        Return a dictionnary containing exhaustive information about the run.
        * barcode_split
            Add statistics split per barcode
        * run_id_split
            Add statistics split per run_id
        """
        self.logger.info ("\tCompute overall summary statistics")
        d = OrderedDict ()
        for df, lab in ((self.all_df, "All Reads"), (self.pass_df, "Pass Reads")):

            d[lab] = OrderedDict()
            d[lab]["General_stats"] = self._compute_stats(df)

            if run_id_split:
                d[lab]["run_id_stats"] = OrderedDict ()
                for id, sdf in df.groupby("run_id"):
                    d[lab]["run_id_stats"][id] = self._compute_stats(sdf)

            if self.has_barcodes and barcode_split:
                d[lab]["barcode_stats"] = OrderedDict ()
                for id, sdf in df.groupby("barcode"):
                    d[lab]["barcode_stats"][id] = self._compute_stats(sdf)
        return d

    def _compute_stats (self, df):
        d = OrderedDict ()
        d["Number of reads"] = len(df)
        d["Number of bases"] = int(df["num_bases"].sum())
        d["Quality score quantiles"] = self._compute_quantiles (df["mean_qscore"])
        d["Read length quantiles"] = self._compute_quantiles (df["num_bases"])
        d["N50"] = int(self._compute_N50(df["num_bases"]))
        d["Run Duration"] = float(np.ptp(df["start_time"])/3600)
        d["Active channels"] = int(df["channel"].nunique())
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
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        if groupby:
            self.logger.info ("\tPlotting reads summary by {}".format(groupby))
        else:
            self.logger.info ("\tPlotting overall reads summary")

        # Prepare all data
        lab1, dd1 = self.__summary_data (df_level="all", groupby=groupby)
        lab2, dd2 = self.__summary_data (df_level="pass", groupby=groupby)

        # Plot initial data
        data = [go.Table(header = dd1["header"][0], cells = dd1["cells"][0], columnwidth = [60, 20])]

        # Create update buttons
        updatemenus = [
            dict (type="buttons", active=0, x=-0.05, y=1, xanchor='right', yanchor='top', buttons = [
                dict (label=lab1, method='restyle', args=[dd1]),
                dict (label=lab2, method='restyle', args=[dd2])])]

        # Autodefine height depending on the numbers of run_ids
        if not height:
            height=300+(30*self.all_df[groupby].nunique()) if groupby else 300

        # tweak plot layout
        layout = go.Layout (updatemenus=updatemenus, width=width, height=height, title=plot_title)

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
            header = [{"values":header, "fill":{"color":"lightgrey"}, "align":"center", "font":{"color":'black', "size":12}, "height":40}],
            cells  = [{"values":cells, "fill":{"color":["white"]}, "align":"center", "font":{"color":'black', "size":12}, "height":30}])

        label = "{} Reads".format(df_level.capitalize())
        return (label, data_dict)

    def __df_to_cell (self, df, id=None):
        """Extract information from sub-dataframes and return a list of values"""
        l = []
        if id:
            l.append (id)
        l.append ("{:,}".format(len(df)))
        l.append ("{:,}".format(df["num_bases"].sum()))
        l.append ("{:,.2f}".format(df["num_bases"].median()))
        l.append ("{:,.2f}".format(self._compute_N50(df["num_bases"])))
        l.append ("{:,.2f}".format(df["mean_qscore"].median()))
        l.append ("{:,}".format(df["channel"].nunique()))
        l.append ("{:,.2f}".format(np.ptp(df["start_time"])/3600))
        return l

    #~~~~~~~1D DISTRIBUTION METHODS AND HELPER~~~~~~~#

    def reads_len_1D (self,
        color:str="lightsteelblue",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Distribution of read length"):
        """
        Plot a distribution of read length (log scale)
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting read length distribution")

        # Prepare all data
        lab1, dd1, ld1 = self.__reads_1D_data (df_level="all", field_name="num_bases", xscale="log", nbins=nbins, smooth_sigma=smooth_sigma)
        lab2, dd2, ld2 = self.__reads_1D_data (df_level="pass", field_name="num_bases", xscale="log", nbins=nbins, smooth_sigma=smooth_sigma)

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
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = plot_title,
            xaxis = {"title":"Read length (log scale)", "type":"log", "zeroline":False, "showline":True},
            yaxis = {"title":"Read density", "zeroline":False, "showline":True, "fixedrange":True, "range":ld1["yaxis.range"]})

        return go.Figure (data=data, layout=layout)

    def reads_qual_1D (self,
        color:str="salmon",
        nbins:int=200,
        smooth_sigma:float=2,
        width:int=None,
        height:int=500,
        plot_title:str="Distribution of read quality scores"):
        """
        Plot a distribution of quality scores
        * color
            Color of the area (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting read quality distribution")

        # Prepare all data
        lab1, dd1, ld1 = self.__reads_1D_data (df_level="all", field_name="mean_qscore", nbins=nbins, smooth_sigma=smooth_sigma)
        lab2, dd2, ld2 = self.__reads_1D_data (df_level="pass", field_name="mean_qscore", nbins=nbins, smooth_sigma=smooth_sigma)

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
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = plot_title,
            xaxis = {"title":"Read quality scores", "zeroline":False, "showline":True},
            yaxis = {"title":"Read density", "zeroline":False, "showline":True, "fixedrange":True, "range":ld1["yaxis.range"]})

        return go.Figure (data=data, layout=layout)

    def __reads_1D_data (self, df_level="all", field_name="num_bases", xscale="linear", nbins=200, smooth_sigma=2):
        """Private function preparing data for reads_len_1D and reads_qual_1D"""
        self.logger.debug ("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        # Get data
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        data = df[field_name].values

        # Count each categories in log or linear space
        min = np.nanmin(data)
        max = np.nanmax(data)
        if xscale == "log":
            count_y, bins = np.histogram (a=data, bins=np.logspace (np.log10(min), np.log10(max)+0.1, nbins))
        elif xscale == "linear":
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

    def reads_len_qual_2D (self,
        colorscale = [[0.0,'rgba(255,255,255,0)'], [0.1,'rgba(255,150,0,0)'], [0.25,'rgb(255,100,0)'], [0.5,'rgb(200,0,0)'], [0.75,'rgb(120,0,0)'], [1.0,'rgb(70,0,0)']],
        len_nbins:int=200,
        qual_nbins:int=75,
        smooth_sigma:float=2,
        width:int=None,
        height:int=600,
        plot_title:str="Mean read quality per sequence length"):
        """
        Plot a 2D distribution of quality scores vs length of the reads
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/ (Not recommanded to change)
        * len_nbins
            Number of bins to divide the read length values in (x axis)
        * qual_nbins
            Number of bins to divide the read quality values in (y axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting read length vs read quality 2D distribution")

        # Prepare all data
        lab1, dd1 = self.__reads_2D_data (df_level="all", len_nbins=len_nbins, qual_nbins=qual_nbins, smooth_sigma=smooth_sigma)
        lab2, dd2 = self.__reads_2D_data (df_level="pass", len_nbins=len_nbins, qual_nbins=qual_nbins, smooth_sigma=smooth_sigma)

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
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = plot_title,
            xaxis = {"title":"Estimated read length", "showgrid":True, "zeroline":False, "showline":True, "type":"log"},
            yaxis = {"title":"Read quality scores", "showgrid":True, "zeroline":False, "showline":True,})

        return go.Figure (data=data, layout=layout)

    def __reads_2D_data (self, df_level, len_nbins, qual_nbins, smooth_sigma=1.5):
        """ Private function preparing data for reads_len_qual_2D """
        self.logger.debug ("\t\tPreparing data for {} reads".format(df_level))

        # Extract data field from df
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df
        len_data = df["num_bases"]
        qual_data = df["mean_qscore"]

        len_min, len_med, len_max = np.percentile (len_data, (0, 50, 100))
        qual_min, qual_med, qual_max = np.percentile (qual_data, (0, 50, 100))

        len_bins = np.logspace (start=np.log10((len_min)), stop=np.log10(len_max)+0.1, num=len_nbins, base=10)
        qual_bins = np.linspace (start=qual_min, stop=qual_max, num=qual_nbins)
        z, y, x = np.histogram2d (x=qual_data, y=len_data, bins=[qual_bins, len_bins])

        if smooth_sigma:
            z = gaussian_filter(z, sigma=smooth_sigma)

        z_min, z_max = np.percentile (z, (0, 100))

        # Extract label and values
        data_dict = dict (
            x = [x, [len_med]],
            y = [y, [qual_med]],
            z = [z, None],
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
            Color of cumulative yield area (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * interval_color
            Color of interval yield line (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting sequencing output over experiment time")

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
            width = width,
            height = height,
            updatemenus = updatemenus,
            legend = {"x":-0.05, "y":1,"xanchor":'right',"yanchor":'top'},
            title = plot_title,
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
            y = np.bincount(t, weights=df["num_bases"].values)

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

    def len_over_time (self,
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
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting read length over experiment time")

        # Prepare all data
        lab1, dd1 = self.__over_time_data (df_level="all", field_name="num_bases", smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab2, dd2 = self.__over_time_data (df_level="pass", field_name="num_bases", smooth_sigma=smooth_sigma, time_bins=time_bins)

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
                go.layout.updatemenu.Button (label=lab1, method='restyle', args=[dd1]),
                go.layout.updatemenu.Button (label=lab2, method='restyle', args=[dd2])])]

        # tweak plot layout
        layout = go.Layout (
            width = width,
            height = height,
            updatemenus = updatemenus,
            legend = {"x":-0.07, "y":1,"xanchor":'right',"yanchor":'top'},
            title = plot_title,
            yaxis = {"title":"Read length (log scale)", "type":"log", "zeroline":False, "showline":True, "rangemode":'nonnegative', "fixedrange":True},
            xaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":True, "rangemode":'nonnegative'})

        return go.Figure (data=data, layout=layout)

    def qual_over_time (self,
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
            Color of median line color (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in (x axis)
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting read quality over experiment time")

        # Prepare all data
        lab1, dd1 = self.__over_time_data (df_level="all", field_name="mean_qscore", smooth_sigma=smooth_sigma, time_bins=time_bins)
        lab2, dd2 = self.__over_time_data (df_level="pass", field_name="mean_qscore", smooth_sigma=smooth_sigma, time_bins=time_bins)

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
            width = width,
            height = height,
            updatemenus = updatemenus,
            legend = {"x":-0.07, "y":1,"xanchor":'right',"yanchor":'top'},
            title = plot_title,
            yaxis = {"title":"Mean quality", "zeroline":False, "showline":True, "rangemode":'nonnegative', "fixedrange":True},
            xaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":True, "rangemode":'nonnegative'})

        return go.Figure (data=data, layout=layout)

    def __over_time_data (self, df_level, field_name="num_bases", smooth_sigma=1.5, time_bins=500, sample=100000):
        """Private function preparing data for qual_over_time"""
        self.logger.debug ("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        # get data
        df = self.pass_sample_df if df_level == "pass" else self.all_sample_df

        # Bin data in categories
        t = (df["start_time"]/3600).values
        x = np.linspace (t.min(), t.max(), num=time_bins)
        t = np.digitize (t, bins=x, right=True)

        # List quality value per categories
        bin_dict = defaultdict (list)
        for bin_idx, val in zip (t, df[field_name].values) :
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
            List of colors (hex, rgb, rgba, hsl, hsv or any CSV named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting barcode distribution")

        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError ("No barcode information available")

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
            legend = {"x":-0.2, "y":1,"xanchor":'left',"yanchor":'top'},
            updatemenus = updatemenus,
            width = width,
            height = height,
            title = plot_title)

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

    #~~~~~~~BARCODE_COUNT METHODS AND HELPER~~~~~~~#

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
            With of the ploting area in pixel
        * height
            height of the ploting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info ("\tPlotting channel activity")

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
            width = width,
            height = height,
            updatemenus = updatemenus,
            title = plot_title,
            xaxis = {"title":"Channel id", "zeroline":False, "showline":False, "nticks":20, "showgrid":False},
            yaxis = {"title":"Experiment time (h)", "zeroline":False, "showline":False, "hoverformat":".2f", "fixedrange":True})

        return go.Figure (data=data, layout=layout)

    def __channels_activity_data (self, df_level, count_level="bases", n_channels=512, smooth_sigma=2, time_bins=150, sample=100000):
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
            for t_idx, channel, n_bases in zip(t, df["channel"], df["num_bases"]):
                z[t_idx][channel-1]+=n_bases
        elif count_level == "reads":
            for t_idx, channel in zip(t, df["channel"]):
                z[t_idx][channel-1]+=1
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

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    @staticmethod
    def _compute_quantiles (data):
        d = OrderedDict ()
        quantil_lab = ("Min","C1","D1","Q1","Median","Q3","D9","C99","Max")
        quantile_val = np.quantile(data, q=[0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1])
        for lab, val in (zip(quantil_lab, quantile_val)):
            d[lab] = val
        return d

    @staticmethod
    def _compute_N50 (data):
        data = data.values.copy()
        data.sort()
        half_sum = data.sum()/2
        cum_sum = 0
        for v in data:
            cum_sum += v
            if cum_sum >= half_sum:
                return v
