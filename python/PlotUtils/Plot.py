"""
  Plot.py:
   Base classes for convenience plot wrappers.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    February 2012
"""

import array
import copy
import itertools
import operator
import os.path
import pprint
import sys
import uuid

import ROOT

import Histogram

# If all the plots to be made from Plot and OverlayPlot instances
# will be made using the same chain, you can set this module-level
# variable to that Chain for simplicity.  Then you don't need to
# specify the chain for each (Overlay)Plot instance.
# ( Set it from your top-level module after importing Plot:
#   Plot.CHAIN = FileManager.ChainManager( .... ) )
CHAIN=None

class Plot:
	""" Plot wrapper class.
	
	This class is designed to fulfill two purposes:
	- unify the interfaces for ROOT's basic plot classes (TH1s, TGraphs, etc.)
	- allow easy creation of plots from TChains

	The two sets of functionality can be used independently.
	To create a plot from a TChain, you'll typically call the constructor directly.
	For example, to plot a 2D histogram of 'var2' vs. 'var1' from a chain called
	'mychain', you might do:
	 plot = Plot(
	    plot_name = "myhisto",
	    plot_title = "My sample 2D histogram",
	    plot_nbins = [ 100, 25, ],
	    plot_ranges = [ (0,100), (0,25), ],
	    axis_labels = [ "var1", "var2", ],
	    draw_expr = "var2:var1",  # ROOT TTree->Draw() style here
	    chain = mychain,
	 )
	which you can then, for instance, Draw() or SaveToDisk().

	On the other hand, if you just want to wrap a ROOT TH1D that you already have,
	you might do
	 plot = Plot.BuildFromROOT(myhisto)
	which, again, you can Draw() and SaveToDisk() and such.  See the BuildFromROOT()
	documentation for more details.
	(The benefit of this latter functionality is much extended by the presence
	of the OverlayPlot, also within this module, whose documentation
	you may also wish to see.)
	"""
	         
	def __init__(
		self,
		plot_name,
		plot_nbins=None,
		plot_title="",
		plot_ranges=None,
		plot_bins=None,
		axis_labels=[],
		log_axes=[],
		plot_type=float,
		base_root_class=None,
		draw_options="",
		draw_overflow=False,
		draw_overflow_boundary=True,
		draw_expr=None,
		cuts=None,
		chain=None,
		show_stats=False,
		norm_variable=True,
		build_now=True
	):
		""" Plot constructor.
		
		Mandatory arguments:
		 - draw_expr -- the expression to plot from the TChain, in TTree->Draw() style
		 - plot_name -- name of the plot
		 - plot_nbins -- list of numbers of bins for each axis that is binned.  Note that this means
		                  that it should have one entry for a 1-D histogram, two for a 2-D histogram, etc.
		 - plot_ranges and plot_bins -- specifies the axis binning.  Either provide plot_ranges, which is
		                                  a list of tuples of (low_bin, high_bin) ranges, so as to provide
		                                  fixed bins; or provide plot_bins, which is a list of the low bin 
		                                  edges (and can have specify variable bin widths).

		Optional arguments:
		 - axis_labels -- list of axis labels you wish to use for your axes
		 - build_now -- create a TH?? object during initialization.  (You usually won't need this.)  Defaults to True.
		 - chain -- The TChain to use for plotting.  You may either provide this at initialization, or you
		            may set the module-level variable CHAIN instead (if all Plots made in this session will
		            be made from the same TChain).
		 - draw_options -- extra ROOT Draw() options you want passed when drawn (e.g., "p x0")
		 - draw_overflow -- draw the overflow (and/or underflow) as separate bin(s) (if content is nonzero).  False by default.
		 - log_axes -- a list of integers indicating which axes (counted from X-axis=1!) you want on a log scale
		 - plot_title -- plot title (printed across the top of the plot when displayed).  Can be left blank.
		 - plot_type -- what type of plot you want: float (which corresponds to TH?D) or int (TH?I).
		                 (Other ROOT histogram types aren't supported for now).  Defaults to float.
		 - base_root_class -- if you want to manually specify the ROOT class used internally
		                      (instead of allowing the PlotUtils machinery to try to guess intelligently)
		                      you can do it here.
		 - show_stats -- Do you want the ROOT stat box drawn?  Defaults to False.
		"""
		self._plot = None

		self.axis_labels = axis_labels
		self.cuts = cuts
		self.draw_options = draw_options
		self.draw_overflow = draw_overflow
		self.draw_expr = draw_expr
		self.log_axes = log_axes
		self.title = plot_title
		self.show_stats = show_stats
		self.draw_overflow_boundary = draw_overflow_boundary
		self.overflow_boundary = None
		self.norm_variable = norm_variable
		self.additional_graphics = []  # extra stuff to draw when drawing the plot
		
		self._plot_config = {
			"name": plot_name,
			"title": plot_title,
			"nbins": plot_nbins,
			"bins": plot_bins,
			"plot_type": plot_type,
			"ranges": plot_ranges,
			"base_root_class": base_root_class,
		}

		self.chain = chain or CHAIN
		
		self._built = False
		
		if build_now:
			self.Build()
			
	def __getattr__(self, attr):
		""" Check the underlying plot for properties that weren't found by normal attribute lookup.
		
		The Plot class is a wrapper around an internal '_plot' element, which contains the actual plot.
		For it to work transparently, method and attribute lookups that apply to the internal '_plot'
		need to be passed along to it.  This method does that.
		
		Note that __getattr__() is ONLY called when regular attribute lookup fails: so if this class
		happens to shadow a method in the '_plot' by giving it the same name, then this class's
		version is the only one that will be called.
		"""
		if hasattr(self._plot, attr):
			try:
				return getattr(self._plot, attr)
			# we want to return the attribute error for me,
			# the top-level object, not my Histogram object
			except AttributeError:
				pass
		
		# this is ok since __getattr__ is only called
		# when normal attribute lookup fails.
		raise AttributeError("'Plot' object has no attribute '%s'" % attr)

	def Build(self):
		""" Build the underlying plot. """ 
		if self._built:
			return
			
		if self._plot is None:
			config = copy.deepcopy(self._plot_config)
			if "plot_type" in config:
				config["histo_type"] = config["plot_type"]
				del config["plot_type"]
			self._plot = self._BuildPlot(config)
		
		self._built = True
	
	@staticmethod
	def BuildFromROOT(root_plot, **kwargs):
		""" Build a Plot from a ROOT object.
		
		Use this method if you wish to wrap a ROOT histogram/graph type with a Plot.
		All of the options one can supply to the constructor can be passed here,
		but generally you will only need a few besides the mandatory 'root_plot':
		 - plot_name -- if you wish to rename this plot (otherwise will use root_plot.GetName())
		 - plot_title -- if you want a different title than root_plot.GetTitle()
		 - axis_labels -- if you want different axis labels than those already on root_plot)
		 - axis_ranges -- if you want to resize the plot from root_plot
		 
		Please see the documentation for __init__ for more information on these options
		(including what format to specify them in).
		"""
		
		if isinstance(root_plot, ROOT.TH1):
			h = Histogram.BuildFromROOT(root_plot)
			if "plot_name" not in kwargs:
				kwargs["plot_name"] = h.name
			if "plot_title" not in kwargs:
				kwargs["plot_title"] = h.title
			if "axis_labels" not in kwargs:
				kwargs["axis_labels"] = []
				for axis_num in range(1, h.Dimension()+1):
					axis = h.GetAxis(axis_num)
					kwargs["axis_labels"].append(axis.GetTitle())
			kwargs["plot_ranges"] = h.ranges
			kwargs["plot_bins"] = h.bins
		else:
			h = root_plot

		kwargs["build_now"] = False

		plot = Plot(**kwargs)
		plot._plot = h
		plot._built = True
		return plot
	
	def _BuildPlot(self, plot_config):
		""" Internal method used to build the plot. """
		return Histogram.BuildNew(**plot_config)
		
	def Draw(self, canvas=None):
		""" Draw the plot, using the specified canvas if appropriate. """
		if isinstance(canvas, ROOT.TCanvas):
#			canvas.Clear()
			canvas.cd()
		else:
			canvas = ROOT.TCanvas(str(uuid.uuid4()))
			canvas.cd()

		# variable-width bins: normalize to smallest bin size.
		if self.norm_variable and (hasattr(self, "NormalizeVariableAxes") or hasattr(self._plot, "NormalizeVariableAxes")):
#			print "normalizing plot:", self.name
			self.NormalizeVariableAxes()
#		if self._plot_config["bins"] is not None and :
#			self._plot.NormalizeVariableAxes()

		if hasattr(self._plot, "SetStats"):
			self._plot.SetStats(self.show_stats)
		
		# only supported for 1D for now
		plot_to_draw = self._plot
		if self.draw_overflow:
			if not hasattr(self._plot, "IsBinOverflow"):
				raise Exception("Can only draw overflow bins for TH1 derivatives, not '%s'" % type(self._plot))
			if len(self._plot_config["nbins"]) != 1:
				raise Exception("Can only draw overflow bins for 1D plots")

			# need to retain this, or it disappears when the Python object goes out of scope
			self._alternate_plot = self.MakeOverflowBinPlot(self._plot)
			plot_to_draw = self._alternate_plot
			
			self.overflow_boundary = ROOT.TLine(
				plot_to_draw.GetBinLowEdge(plot_to_draw.GetNbinsX()),
				plot_to_draw.GetMinimum(),
				plot_to_draw.GetBinLowEdge(plot_to_draw.GetNbinsX()),
				plot_to_draw.GetMaximum(),
			)
			self.overflow_boundary.SetLineStyle(2) # dashed
			self.overflow_boundary.SetLineWidth(3)
		
		draw_options = self.draw_options
		if len(canvas.GetListOfPrimitives()) > 0 and "same" not in draw_options.lower():
			draw_options += " same"
		plot_to_draw.Draw(draw_options)
		
		# set any log axes requested
		axes = ("x", "y", "z")
		axis_names = dict( enumerate(axes) )
		for axis in self.log_axes:
			axis_name = axis if axis in axes else axis_names[axis]
			axis_name = axis_name.lower()
#			print "Setting logarithmic %s axis..." % axis_name
			getattr(canvas, "SetLog" + axis_name)(True)

		# force axis labels.  maybe not a great idea?
		self.SetLabels()
		
		canvas.Update()
		
		for obj in self.additional_graphics:
			if hasattr(obj, "Draw"):
				obj.Draw("same")
				
		if self.draw_overflow and self.draw_overflow_boundary:
			self.overflow_boundary.Draw("same")
				
		ROOT.gSystem.ProcessEvents()
		
		return canvas
	
	@staticmethod
	def MakeOverflowBinPlot(plot):
		if not hasattr(plot, "Clone"):
			raise Exception("Can't add an overflow bin to a plot that can't be cloned!")

		# no overflow if the plot hasn't been filled yet...
		if plot.GetXaxis().GetXbins().fN < 1:
			return plot
	
		# definitely don't want the replacement histogram
		# (with the same name!) added to the current directory
		adddir_status = ROOT.TH1.AddDirectoryStatus()
		ROOT.TH1.AddDirectory(False)
		newplot = plot.Clone()
		ROOT.TH1.AddDirectory(adddir_status)
		
		bins = newplot.GetXaxis().GetXbins()
		# can't slice because 'bins' is a PyROOT buffer...
		newbins = [bins[i] for i in range(newplot.GetNbinsX()+1)]  

		# 1 is the magic number here because
		# we always use the "width" option to TH1::Scale(),
		# which divides each bin content by its width
		# in the units of the axis.
		# thus, if the overflow bin is of unit width,
		# it will never be affected by this kind of renormalization.
		newbins.append( newbins[-1] + 1 )
		
		# this plunks the old overflow into the new "last" bin.
		newplot.SetBins(len(newbins)-1, array.array('d', newbins)) 
		
		return newplot
	
	def Plot(self, draw_expr="", cuts="", chain=None):
		""" Create the plot from the arguments.
		
		This method actually creates the Plot from a chain.
		(You shouldn't use it if your Plot was created from
		 an already-existing ROOT histogram.)
		 
		You can override the draw expression, cuts, and chain to use
		that were specified in the Plot constructor if you wish
		(see __init__()'s documentation for more details on those).
		"""
		
		chain = chain or self.chain
		if chain is None:
			raise ValueError("Must specify chain to use either in Plot constructor or as Plot() method parameter...")
	
		if draw_expr == "" and self.draw_expr is None:
			raise ValueError("Expression to draw must be specified either in plot constructor or as 'expr' argument to Plot() method...")
		elif draw_expr == "":
			draw_expr = self.draw_expr
			
		if cuts == "" and self.cuts is not None:
			cuts = self.cuts
			
#		self.Build()
		if hasattr(self._plot, "SetDirectory"):
			self._plot.SetDirectory(ROOT.gDirectory)

		adding_already = ROOT.TH1.AddDirectoryStatus()
		if not adding_already:
			ROOT.TH1.AddDirectory(True)
#		print "drawing expression '%s' with cuts '%s'" % (draw_expr, cuts)
		chain.Draw("%s>>%s" % (draw_expr, self._plot_config["name"]), cuts, "goff " + self.draw_options)
		
		# if the plot wasn't built to begin with,
		# then we can retrieve it from the ROOT namespace.
		if self._plot is None:
			self._plot = Histogram.BuildFromROOT(ROOT.gDirectory.Get(self._plot_config["name"]))

		if not adding_already:
			self._plot.SetDirectory(None)

		self.SetLabels()
		
	def SaveToDisk(self, filename, canvas=None):
		""" Save this plot to disk using the specified filename, and optionally, the given canvas. """
		
		if self._plot is None:
			raise Exception("Can't write a non-plot to disk!")
		
		canvas = self.Draw(canvas)
		canvas.Print(filename)
				
	def SetLabels(self):
		""" Force the axis labels to correspond to the configuration given to __init__(). """
		self._plot.SetTitle(self.title)
	
		for i in range(len(self.axis_labels)):
			axis = self._plot.GetAxis(i+1) if hasattr(self._plot, "GetAxis") else Histogram.GetAxis(self._plot, i+1)
#			print self.axis_labels[i]
			axis.SetTitle(self.axis_labels[i])
	

class OverlayPlot(Plot):
	""" Class for creating an overlay of several different plots in the same canvas.
	
	Because this class inherits from the Plot class of this module, its design
	is very similar: it can both create a new overlay from a TChain or wrap
	a collection of already-built ROOT plots.
	
	Please see the documentation for __init__ for more information on how to do this.
	    
	"""
	OVERLAY_DEFAULTS = {
		"colors": {},
		"fill_colors": {},
		"legend_map": {},
		"legend_style": {},
		"legend_position": [0.7, 0.7, 0.9, 0.9],
		"legend_header": None,
		"show_stats": False,
		"normalize_by": "self",
#		"norm_variable": True,
	}
	
	def __init__(
		self,
		keys,
		base_plot_name,
		event_list_expr=None,
		selection={},
		stack=False,
		**kwargs
	):
		""" Constructor.
		
		There are a few mandatory attributes that must be passed to the constructor
		if you want to plot from a TChain (the use case for calling this constructor directly):
		 - keys -- this is the list of names you plan to index your individual plots by.
		 - base_plot_name -- the base name for your plots ("_<key>" will be added to this).
		 - selection -- a dictionary containing the cut expressions that select the samples
			           for each key.  For instance: selection={"e": "pid==11", "mu": "pid==13"}
		Other options:
		 - draw_order -- a list specifying the keys in the order you want them drawn on a canvas.
		 - event_list_expr -- if you wish to use a ROOT TEventList, this is the expression
			                 used to generate it.
		 - stack -- set stack=True if you want your plots drawn stacked instead of overlaid.
		 - draw_options -- same idea as in the Plot() constructor.  Here, however, you can specify
		                 it in two ways:
		                  * as a regular string, in which case it will be applied to all the plots;
		                  * as a dictionary, in which case each item will be applied to the plot
		                    whose key matches the dictionary key.
		 - legend_map -- a dictionary of the legend text entries you want for each key
		 - legend_style -- a dictionary of the legend entry styles you want for each key
		 - normalize_by -- how to normalize the histograms relative to one another.  See Normalize()
		                   for more details (as well as how to specify this parameter).
		
		You can also pass most of the options that are valid for the Plot constructor as well;
		however, because of various bugs in the ROOT THStack, not all of them work correctly yet.
		"""

		self._keys = keys
		self.event_list_expr = event_list_expr
		self.event_list = None
		self.selection = selection or dict( zip( self._keys, ["",] * len(self._keys) ) )
		self._plots = {}
		self.stack = stack
		kwargs["plot_name"] = base_plot_name
		self.name = base_plot_name
		
		# save this.  we'll need it in a moment.
		build_now = True if "build_now" not in kwargs else kwargs["build_now"]
		
		self._autoranged = False
		self._normalized = False
		
		self._alternate_plots = dict.fromkeys(self._keys)

		for parameter in OverlayPlot.OVERLAY_DEFAULTS:
			if parameter in kwargs:
				# if the item is keyed, make sure all the keys exist!
				if hasattr(kwargs[parameter], "keys"):
					keepers = {}
					for k in kwargs[parameter]:
						if k in self._keys:
							keepers[k] = kwargs[parameter][k]
					kwargs[parameter] = keepers
				setattr(self, parameter, kwargs[parameter])
				del kwargs[parameter]
			else:
				setattr(self, parameter, OverlayPlot.OVERLAY_DEFAULTS[parameter])
		
		if "draw_options" in kwargs:
			if isinstance(kwargs["draw_options"], basestring):
				kwargs["draw_options"] = dict( [ (key, kwargs["draw_options"]) for key in self._keys ] )
#			else:
#				self.draw_options = kwargs["draw_options"]
#			del kwargs["draw_options"]
		else:
			kwargs["draw_options"] = {}
			for key in self._keys:
				kwargs["draw_options"][key] = "hist x"
				
		if "draw_order" in kwargs:
			self.draw_order = kwargs["draw_order"]

			# ensure there are only valid keys in the list			
			new_order = []
			for item in self.draw_order:
				if item in self._keys:
					new_order.append(item)
			self.draw_order = new_order
			# ensure all the keys are represented
			self.draw_order += list( set(self._keys) - set(self.draw_order) )
			del kwargs["draw_order"]
		else:
			self.draw_order = self._keys
		
		if self.legend_map is not None and len(self.legend_map) > 0:
			args = self.legend_position
			if self.legend_header is not None:
				args += [self.legend_header,]
			self.legend = ROOT.TLegend(*args)
		else:
			self.legend = None
			
		if "save_individually" in kwargs:
			self.save_individually = kwargs["save_individually"]
			del kwargs["save_individually"]
			
		# don't ever want the Plot superclass
		# to do the building for an OverlayPlot.
		# we'll do it below manually.
		kwargs["build_now"] = False
		Plot.__init__(self, **kwargs)
				
		del self._plot

		if not build_now:
			return

		self.Build()
		
	def __getattr__(self, attr):
		# want to override the Plot.__getattr__() method
		# because we have more than one plot here...
		if hasattr(self._plots, attr):
			return getattr(self._plots, attr)
		
		# this is sensible because __getattr__() is ONLY called
		# if the attribute can't be found via normal attribute lookup.
		raise AttributeError("'OverlayPlot' object has no attribute '%s'" % attr)
		
	# the next handful of methods are here to provide a dictionary-like interface.  for example:
	#  plot = OverlayPlot( ... )
	#  single_plot = plot["single_plot_name"]
	# etc.
	
	def __getitem__(self, key):
		""" Make this object function like a dictionary (get a plot by key from the OverlayPlot). """
		if key not in self._plots:
			raise KeyError("Plot not found in overlay: '%s'" % key)
		
		return self._plots[key]
		
	def __iter__(self):
		for key in self._plots:
			yield key
	
	def iteritems(self):
		for key, value in self._plots.iteritems():
			yield key, value
	
	iterkeys = __iter__
	
	def itervalues(self):
		for val in self._plots.itervalues():
			yield val
			
	def keys(self):
		return self._keys

	def ApplyROOTHistoMethod(self, method, *args, **kwargs):
		""" Apply a ROOT histogram method to all of the individual plots within this overlay.
		
		Useful for doing things like changing line thicknesses, for example:
		  overlay_plot.ApplyROOTHistoMethod("SetLineWidth", [2,])
		"""

		any_success = False				
		for h in self._plots.itervalues():
			if hasattr(h, method):
				getattr(h, method)(*args, **kwargs)
				any_success = True
				
		if not any_success:
			raise AttributeError("ROOT method '%s' doesn't exist in any contained plots" % method)
			
	def AutoRange(self, use_errors = True):
		""" Change the axes so that all the plots in an overlay are visible.
		
		Note that this method is only used if the plots are not compatible with
		a THStack (otherwise, the THStack does this for us).
		
		If use_errors is True (default), the autorange procedure will zoom
		the y-axis of a 1-D histogram such that the error bars are also visible.
		"""
		
		if self._autoranged:
			return
		
		if len(self._plots) <= 1:
			return

#		for h in self._plots.itervalues():
#			if not hasattr(h, "FindFirstBinAbove"):
##				print "don't have necessary method in", type(h)
#				return

		try:		
			plot_dimension = self.Dimension()
			if plot_dimension is None:
				return
			plot_dim_list = [0,] if plot_dimension == 1 else xrange(plot_dimension-1)
		except AttributeError:
			plot_dimension = None
			plot_dim_list = []
		
		#plot_dim_list = range(plot_dimension)

		if use_errors:
			plot_extremes = []
			for p in self._plots.itervalues():
				if hasattr(p, "GetBinError"):
					plot_extremes.append( (p.GetMinimum() - p.GetBinError(p.GetMinimumBin()), p.GetMaximum() - p.GetBinError(p.GetMaximumBin())) )
				elif hasattr(p, "GetY") and hasattr(p, "GetEYlow") and hasattr(p, "GetEYhigh"):
					N = int(p.GetN())
					y = p.GetY()
					yerr_low = p.GetEYlow()
					yerr_high = p.GetEYhigh()
					
					# there doesn't seem to be a way to check for a buffer
					# that's built from a NULL pointer...
					try:
						ylows = [y[i]-yerr_low[i] for i in range(N)] if yerr_low else [0.] * N
					except IndexError:
						ylows = [0.] * N
					try:
						yhighs = [y[i]+yerr_high[i] for i in range(N)] if yerr_high else [0.] * N
					except IndexError:
						yhighs = [0.] * N
				
					plot_extremes.append( (min(ylows), max(yhighs)) )
		else:
			plot_extremes = [ (self._plots[h].GetMinimum(), self._plots[h].GetMaximum()) for h in self._plots ]

		axis_extremes = [ [ [], [] ] for i in plot_dim_list ]
		for h in self._plots:
			if hasattr(h, "FindFirstBinAbove"):
				for axis in plot_dim_list:
					lowside = self._plots[h].FindFirstBinAbove(0, axis+1)
					if lowside >= 0:				
						axis_extremes[axis][0].append(lowside)
					highside = self._plots[h].FindLastBinAbove(0, axis+1)
					if highside >= 0:
						axis_extremes[axis][1].append(highside)
#		print plot_extremes
#		print axis_extremes
		plotmax = max(plot_extremes, key=operator.itemgetter(1))[1] * 1.1
		plotmin = min(plot_extremes, key=operator.itemgetter(0))[0] * 0.9
		
#		print plotmin, plotmax
		
		extreme_finder = (min, max)
		multiplier = (0.9, 1.1)
		
		axis_bin_limits = []
		for axis in axis_extremes:
			axis_bin_limits.append( [int(extreme_finder[side](axis[side]) * multiplier[side]) for side in range(2)] )
		
#		print axis_bin_limits

		for h in self._plots:
			for dim in plot_dim_list:
				print "setting axis %d to bin range: %s" % (dim+1, axis_bin_limits[dim])
				self._plots[h].GetAxis(dim+1).SetRange(*axis_bin_limits[dim])
				
#			print "setting axis %d to user range: %s" % (plot_dimension+1, [plotmin, plotmax])
#			self._plots[h].GetAxis(plot_dimension+1).SetRangeUser(0 if plotmin > 0 else plotmin, plotmax)
			self._plots[h].SetMinimum(0 if plotmin > 0 else plotmin)
			self._plots[h].SetMaximum(plotmax)
		
		self._autoranged = True
		
	def Build(self):
		""" Build the underlying plots. """ 
		if self._built:
			return
	
		title = "" if "title" not in self._plot_config else self._plot_config["title"]
		self._plots = Histogram.HistoStack(self.name, title)
		self._alternate_plots = Histogram.HistoStack(self.name, title)  # for drawing overflow bins

		for stack in (self._plots, self._alternate_plots):
			ROOT.SetOwnership(stack, False)  # don't let any Canvas assume ownership of the stack!
		
	
		for key in self.draw_order:
			if key not in self._plots:
				config = copy.deepcopy(self._plot_config)
				config["name"] = self._KeyedPlotName(config["name"], key)
#				print "making plot with name:", config["name"]
				if "plot_type" in config:
					config["histo_type"] = config["plot_type"]
					del config["plot_type"]
				
				opt = self.draw_options[key]
				if self.show_stats:
					opt += " sames"  # I hate you, ROOT
				self._plots.Add(self._BuildPlot(config), key, option=opt)
	
		self._built = True

	@staticmethod
	def BuildFromROOT(root_plots, **kwargs):
		""" Build an OverlayPlot from a collection of ROOT objects.
		
		Use this method if you wish to stuff a collection of ROOT histogram/graph types
		into an overlay.
		
		The mandatory argument 'root_plots' must be provided as a keyed dictionary:
		  root_plots = {
		  	"key1": plot1,
		  	"key2": plot2,
		  	...
		  }
		All of the options one can supply to the OverlayPlot constructor can be passed here,
		but generally you will only need a few besides 'root_plots':
		 - plot_name -- if you wish to specify the base name for this overlay (otherwise will use a UUID)
		 - axis_labels -- if you want different axis labels than those already on root_plots
		 - axis_ranges -- if you want to resize the plots from root_plots
		 
		Please see the documentation for __init__ for more information on these (and other) options.
		"""
		if "base_plot_name" not in kwargs:
			kwargs["base_plot_name"] = "overlay_" + str(uuid.uuid4())

		kwargs["keys"] = root_plots.keys()
		kwargs["build_now"] = False
		plot = OverlayPlot(**kwargs)

		# if these are all TH1s or Plot.Plots, we can use the wrapper.
		# otherwise, you're on your own...
		if all( [ isinstance(p, Plot) for p in root_plots.itervalues() ] ) and len(set([p.dimension for p in root_plots.itervalues()])) == 1:
			title = kwargs["plot_title"] if "plot_title" in kwargs else ""
			collection = Histogram.HistoStack(kwargs["base_plot_name"], title)
			ROOT.SetOwnership(collection, False)  # don't let any Canvas assume ownership of the stack!

			
			for key in plot.draw_order:
				opts = plot.draw_options[key]
				if plot.show_stats:
					opts = " ".join([opts, "sames"]).strip()

				collection.Add(root_plots[key]._plot, key, opts)
		elif all( [ isinstance(p, ROOT.TH1) for p in root_plots.itervalues() ] ) and len(set([isinstance(p, ROOT.TH3) for p in root_plots.itervalues()])) == 1 and len(set([isinstance(p, ROOT.TH2) for p in root_plots.itervalues()])) == 1:
			title = kwargs["plot_title"] if "plot_title" in kwargs else ""
			collection = Histogram.HistoStack(kwargs["base_plot_name"], title)
			ROOT.SetOwnership(collection, False)  # don't let any Canvas assume ownership of the stack!
			
			for key in plot.draw_order:
				if key in plot.draw_options:
					opts = plot.draw_options[key]
				else:
					opts = root_plots[key].GetOption()
				
				if plot.show_stats:
					opts = " ".join([opts, "sames"]).strip()
				
				collection.Add(Histogram.BuildFromROOT(root_plots[key]), key, opts)
		else:
#			print "using collection by itself"
			collection = root_plots
		
		plot._plots = collection
		plot._alternate_plots = type(collection)()
		plot._built = True
		
#		print plot._plots.values()
		
		return plot
		
	def Dimension(self):
		""" Determine the dimension of the plot. """
		if len(self._plots) == 0:
			return None
		
		try:
			dimensions = set( [h.Dimension() for h in self._plots.itervalues()] )
		except AttributeError:
			return None
			
		if len(dimensions) != 1:
			return None
		
		return list(dimensions)[0]
		
	def Divide(self, other_plot):
		""" Divide all the plots in this overlay by a single other_plot. """
		if isinstance(other_plot, OverlayPlot):
			for plot in other_plot._plots:
				if plot not in self._plots:
					print >> sys.stderr,  "warning: OverlayPlot contains no plot '%s' to divide into" % plot
					continue
				self._plots[plot].Divide(other_plot._plots[plot]._plot)
		else:
			for plot in self._plots.itervalues():
				plot.Divide(other_plot)
	
	def Draw(self, plot=None, canvas=None, clone_histos=False):
		""" Draw the overlay or one plot from it, using the specified canvas if appropriate.
		
		If 'plot' is left unspecified (or passed as None), this method will draw the whole
		overlay.  If, on the other hand, one of the keys from the set used to build the
		overlay is specified, only that plot will be drawn.
		
		If 'canvas' is given, that canvas will be used to draw.
		
		If 'clone_histos' is True, the ROOT histograms being drawn will first be
		cloned using the TObject::Clone method.  (This is a good idea if you're going
		to use the output of Draw() as the input to something else.)
		"""
		
		if isinstance(canvas, ROOT.TVirtualPad):
#			canvas.Clear()
			canvas.cd()
		else:
			canvas = ROOT.TCanvas(str(uuid.uuid4()))
			canvas.cd()
	
		if plot is not None:
			self.SetLabels(self._plots)
			if hasattr(self._plots[plot], "NormalizeVariableAxes") and self.norm_variable:
				self._plots[plot].NormalizeVariableAxes()
			if plot in self.draw_options and hasattr(self._plots[plot], "SetDrawOption"):
#				print "using draw option:", self.draw_options[plot]
				self._plots[plot].SetDrawOption(self.draw_options[plot])
			
			if clone_histos:
				plot_to_draw = self._plots[plot].Clone()
			else:
				plot_to_draw = self._plots[plot]
			if self.draw_overflow:
				self._alternate_plots[plot] = self.MakeOverflowBinPlot(self._plots[plot])
				plot_to_draw = self._alternate_plots[plot]
			plot_to_draw.Draw(self.draw_options[plot])
		else:
			if self.legend is not None and self.legend_map is not None:
				self.legend.Clear()
			have_one = False

#			print ("" if self.stack else "NOT"), "stacking plot:", self.GetName()
			legend_items = []
			for key in self.draw_order:
				self.Normalize(norm_axes=self.norm_variable)
				
				if key in self.colors:
					self._plots[key].SetLineColor(self.colors[key])
					self._plots[key].SetMarkerColor(self.colors[key])
				
				if key in self.fill_colors:
					self._plots[key].SetFillColor(self.fill_colors[key])

				if self.legend_map and key in self.legend_map:
					args = []
					if isinstance(self._plots[key], ROOT.TObject):
						args.append(self._plots[key])
						args.append(self.legend_map[key])
					elif hasattr(self._plots[key], "_plot") and isinstance(self._plots[key]._plot, ROOT.TObject):
						args.append(self._plots[key]._plot)
						args.append(self.legend_map[key])
					else:
						raise TypeError("I can't figure out what type your plot is...")
					
					if isinstance(self.legend_style, basestring):
						args.append(self.legend_style)
					elif key in self.legend_style:
						args.append(self.legend_style[key])
						
					legend_items.append(args)

#				print self.draw_options
				if key in self.draw_options and hasattr(self._plots[key], "SetDrawOption"):
#					print "setting option for plot '%s' to: '%s'" % (key, self.draw_options[key])
					self._plots[key].SetDrawOption(self.draw_options[key])
#					print "  ...option is: '%s'" % self._plots[key].GetDrawOption()

				# the remainder is only for when we're not using THStacks
				if isinstance(self._plots, ROOT.THStack):
					continue

				draw_opt = self.draw_options[key]
				if not have_one:
					if isinstance(self._plots[key], (ROOT.TEfficiency, ROOT.TGraph)):
						first_phrase = "a"
					else:
						first_phrase = ""
					draw_opt = " ".join( (draw_opt, first_phrase) )
					have_one = True
				else:
					if isinstance(self._plots[key], (ROOT.TEfficiency, ROOT.TGraph)):
						same_string = ""
					else:
						same_string = "same" if not self.show_stats else "sames"
					draw_opt = " ".join( (draw_opt, same_string ))

				draw_options = " ".join( (self.draw_options[key], draw_opt) ).strip()
				# set stats if asked
				if self.show_stats:
					if hasattr(self._plots[key], "SetStats"):
						self._plots[key].SetStats(True)
				else:
					if hasattr(self._plots[key], "SetStats"):
						self._plots[key].SetStats(False)
#				print "Drawing key '%s' with options: '%s'" % (key, draw_options)
#				print "key '%s' has stats box on:", not(self._plots[key].TestBit(self._plots[key].kNoStats))
				self._plots[key].Draw(draw_options)

			# legend items need to be added in reverse order when stacking
			if self.legend_map:
				for args in ( reversed(legend_items) if self.stack else legend_items ):
					self.legend.AddEntry(*args)					

			if self.draw_overflow:
				has_overflow = self.UpdateOverflowBinPlots()
			plots_to_draw = self._alternate_plots if self.draw_overflow else self._plots
				
			if isinstance(plots_to_draw, ROOT.THStack):
				# find the options which all the plots share...
				common_opts = set.intersection( *[set(opts.split(" ")) for opts in self.draw_options.itervalues()] )
				common_opts = " ".join(common_opts)
				options_used = (common_opts + (" nostack" if not self.stack else "")).strip()
#				print "drawing stack with options: '%s'" % options_used
				
				plots_to_draw.Draw( options_used )
				self.SetLabels(plots_to_draw)
				
				# the axis ranging must also be applied to the THStack.  sigh.
				if self._plot_config["ranges"] is not None:
					for axis, axis_range in enumerate(self._plot_config["ranges"]):
						if axis_range is not None:
							plots_to_draw.GetAxis(axis+1).SetRangeUser(*axis_range)

				if self.draw_overflow and has_overflow:
					template_plot = self._alternate_plots[self._alternate_plots.keys()[0]]
					self.overflow_boundary = ROOT.TLine(
						template_plot.GetBinLowEdge(template_plot.GetNbinsX()),
						self._alternate_plots.GetHistogram().GetMinimum(),
						template_plot.GetBinLowEdge(template_plot.GetNbinsX()),
						self._alternate_plots.GetHistogram().GetMaximum(),
					)
					self.overflow_boundary.SetLineStyle(2) # dashed
					self.overflow_boundary.SetLineWidth(3)

			else:
				if self.Dimension() is not None:	
					self.AutoRange()


			ROOT.gSystem.ProcessEvents()

			if self.show_stats:
				for i, key in enumerate(self.draw_order):	
#					print key
					statsbox = self._plots[key].FindObject("stats")
					if not statsbox:
						print "warning: couldn't find stats box for key: '%s'" % key
						continue
					if key in self.colors:
						statsbox.SetLineColor(self.colors[key])
						statsbox.SetTextColor(self.colors[key])
				
					height = min(0.2, 0.75/len(plots_to_draw))
					statsbox.SetY1NDC(0.85 - height * (i+1))
					statsbox.SetY2NDC(0.85 - height * i)

					statsbox.Draw()

			for obj in self.additional_graphics:
				if hasattr(obj, "Draw"):
					obj.Draw("same")

			if self.draw_overflow and self.draw_overflow_boundary and self.overflow_boundary:
				self.overflow_boundary.Draw("same")

			self.SetLabels(plots_to_draw, use_key=False)
			if self.legend is not None:
				self.legend.SetFillColor(ROOT.kWhite)
				self.legend.Draw()

		## end else (plot is not None)

		axes = ("x", "y", "z")
		axis_names = dict( enumerate(axes, start=1) )
		for axis in self.log_axes:
			if axis in axes:
				axis_name = axis
			elif axis in axis_names:
				axis_name = axis_names[axis]
			else:
				raise NameError("No axis named '%s'" % axis)
			axis_name = axis_name.lower()
#			print "Setting logarithmic %s axis..." % axis_name
			getattr(canvas, "SetLog" + axis_name)(True)
		canvas.Update()
		
		ROOT.gSystem.ProcessEvents()
		
		return canvas
		
	def GetHistogram(self):
		""" Override the ROOT THStack method here so as to get the right one
		    (depending on which was most recently drawn).
		"""
		for obj in (self._plots, self._alternate_plots):
			if hasattr(obj, "GetHistogram"):
				hist = obj.GetHistogram()
				if hist:
					return hist
		
		raise ReferenceError("GetHistogram() only works after a THStack has been drawn")
		
	def GetName(self):
		return self.name
	
	def _GetAxis(self, axis):
		""" Retrieve the OverlayPlot axis corresponding to the 'axis' identifier.
		
		You can specify either the axis name ('x', 'y', 'z'), or the axis number
		(1, 2, 3).
		"""
		if isinstance(self._plots, ROOT.THStack):
			return self._plots.GetAxis(axis)
		else:
			return getattr(self._plots[self.draw_order[0]], "Get%saxis" % axis)()

	
	def GetXaxis(self):
		""" Helper method to get the X axis directly. """
		return self._GetAxis("X")
	
	def GetYaxis(self):
		""" Helper method to get the Y axis directly. """
		return self._GetAxis("Y")
	
	def GetZaxis(self):
		""" Helper method to get the Z axis directly. """
		return self._GetAxis("Z")
	
	def _KeyedPlotName(self, name, key):
		return "%s_key=%s" % (name, key)

	def MapROOTHistoMethod(self, method, arg_map):
		""" Applies a ROOT histogram method to the histograms specified (as identified by their keys).
		    
		    This is similar to ApplyROOTHistoMethod(), except that this method allows you
		    to specify different arguments for the method for each key.
		    
		    Arguments to the method for each histogram should be passed as a dictionary:
		    arg_map = {
		      key1: {
		        "args": [arg1, arg2, ...],
		        "kwargs": { key1: val1, key2: val2, ... },
		      },
		      key2: {
		        "args": [arg1, ...]
		        "kwargs": { key1: val1, ... },
		      },
		      ...
		    } """

		if self.Dimension() is None:
			return
		
		if all( [hasattr(h, method) for h in self._plots.itervalues()] ):
			for key, arg_dict in arg_map.iteritems():
				if "args" not in arg_dict:
					arg_dict["args"] = []
				if "kwargs" not in arg_dict:
					arg_dict["kwargs"] = {}
				getattr(self._plots[key], method)(*arg_dict["args"], **arg_dict["kwargs"])
		else:
			raise AttributeError("ROOT histogram method '%s' doesn't exist" % method)
	
	def Normalize(self, norm_axes=True):
		""" Normalize the plots.
		
		There are two normalizations that happen here:
		 (1) Axes with variable bins are normalized to the smallest bin width.
		 (2) Histograms are normalized relative to each other.  How this happens is
		     governed by the value of the parameter 'normalize_by', specified at
		     initialization:
		      - if 'normalize_by' is None, no normalization is performed.
		      - if 'normalize_by' is the string "self", the histograms are all 
		        normalized to area 1.
		      - if 'normalize_by' is a keyed dictionary, the histograms are each
		        normalized by dividing by the number contained in normalize_by[key].
		"""
		if self._normalized:
			return

		self._normalized = True

		if norm_axes:
			for plot in self._plots.itervalues():
				if hasattr(plot, "NormalizeVariableAxes"):
					plot.NormalizeVariableAxes()

		if self.normalize_by is None:
			return
	
		for key in self._plots:
			# user can specify normalization by individual histogram...
			if hasattr(self._plots[key], "Normalize"):
				if hasattr(self.normalize_by, "keys") and key in self.normalize_by:
					self._plots[key].Normalize(self.normalize_by[key])
				else:
					self._plots[key].Normalize(self.normalize_by)
		

	def Plot(self, draw_expr="", selection=None, other_cuts="", chain=None):
		""" Create the plot from the arguments.
		
		See the documentation for Plot.Plot() for further details.
		"""
		chain = chain or self.chain
		if chain is None:
			raise ValueError("Must specify chain to use either in Plot constructor or as Plot() method parameter...")
	
		draw_expr = draw_expr or self.draw_expr
		if not draw_expr:
			raise ValueError("Expression to draw must be specified either in plot constructor or as 'draw_expr' argument to Plot() method...")
		
		other_cuts = other_cuts or self.cuts

		selection = selection or self.selection
		if isinstance(selection, basestring) and "%(key)s" in selection:
			selection = dict( [ (key, selection % {"key": key}) for key in self.keys() ] )

		if not all([key in selection for key in self.keys()]):
			raise Exception("selection must contain all the keys!")
			
		self.Build()
		
		adding_already = ROOT.TH1.AddDirectoryStatus()
		if not adding_already:
			ROOT.TH1.AddDirectory(True)
		
		if self.event_list_expr is not None and self.event_list is None:
			chain.Draw(">>%s_evtlist" % self.name, self.event_list_expr)
			self.event_list = ROOT.gDirectory.Get("%s_evtlist" % self.name)
			self.event_list.SetReapplyCut(True)
		
		if self.event_list is not None:
			chain.SetEventList(self.event_list)
			
		for key in self.keys():
			plot_name = self._KeyedPlotName(self._plot_config["name"], key)
			s = selection[key]
			if other_cuts:
				s = "(%s) && (%s)" % (s, other_cuts)
			
			opts = ("goff " + self.draw_options[key]).strip()
			
#			print "  Drawing plot with name: %s, '%s'" % (type(plot_name), plot_name)
#			print "    draw expression: %s, '%s'" % (type(draw_expr), draw_expr)
#			print "    draw selection criteria: %s, '%s'" % (type(s), s)
#			print "    draw options: %s, '%s'" % (type(opts), opts)
			
			# we need to associate the plot with the current TDirectory
			# so that the "Draw" expression used below writes to it,
			# instead of creating a new one with the same name.
			# we'll detach it from the directory afterwards.
			if hasattr(self._plots[key], "SetDirectory"):
#				print ROOT.gROOT.gDirectory
				self._plots[key].SetDirectory(ROOT.gROOT)
#			else:
#				print self._plots[key]
			
			N = chain.Draw(
				"%s>>%s" % (draw_expr, plot_name),
				s,
				opts
			)
#			print "  ... %d entries in this plot." % self._plots[key].GetEntries()
#			plot = chain.Get(plot_name)

			# detach the plot as explained above
			if not adding_already and hasattr(self._plots[key], "SetDirectory"):
				self._plots[key].SetDirectory(None)
			
		if not adding_already:
			ROOT.TH1.AddDirectory(False)
			

		# be sure to reset the chain's event list.
		# other plots might use this chain,
		# and they might not want to use the same event list.
		if self.event_list is not None:
			chain.SetEventList(None)
	
		self.SetLabels()
	
	def SaveToDisk(self, base_filename, canvas=None, save_individually=False):
		""" Save this overlay to disk using the specified filename, and optionally, the given canvas.
		
		This method currently saves each member of the overlay individually first, followed by
		the full overlay.
		"""
		if len(self._plots) == 0:
			raise Exception("No plots to write...")
		

		# save all the plots individually first (if requested)
		if save_individually or (hasattr(self, "save_individually") and self.save_individually):
			for key in self.keys():
				canvas = self.Draw(plot=key, canvas=canvas)
			
				d, f = os.path.split(base_filename)
				f, ext = os.path.splitext(f)
				f = "%s_key=%s%s" % (f, str(key), ext)
				canvas.Print(os.path.join(d, f))
		
		# now a composite:
		canvas = self.Draw(plot=None, canvas=canvas)  # 'plot=None' means 'draw them all!'
		canvas.Print(base_filename)

	def SetLabels(self, plots=None, use_key=True):
		""" Force the axis labels to correspond to the configuration given to __init__(). """
		plots = plots or self._plots
#		print pprint.pformat(plots._key_map)
		for key, plot in plots.iteritems():
			if use_key:
				selection = self.selection[key]
				title = "%s (%s)" % (self.title, selection)
			else:
				title = self.title
				
			if not plot:
				continue
				
			plot.SetTitle(title)

			for i in range(len(self.axis_labels)):
				axis = plot.GetAxis(i+1) if hasattr(plot, "GetAxis") else Histogram.GetAxis(plot, i+1)
				axis.SetTitle(self.axis_labels[i])
				axis.SetTitleColor(ROOT.kBlack)  # sometimes they get reset to transparent??
		
		# if we're working with a THStack, it's got its OWN axes too
		if isinstance(plots, ROOT.THStack):
			for i in range(len(self.axis_labels)):
				axis = plots.GetAxis(i+1)
				if axis:
#					print "setting title:", self.axis_labels[i]
					axis.SetTitle(self.axis_labels[i])
					axis.SetTitleColor(ROOT.kBlack)  # sometimes they get reset to transparent??

	def UpdateOverflowBinPlots(self):
		has_overflow = False
		for key in self.draw_order:
			if self[key].GetBinContent(self[key].GetNbinsX()) > 0:
				has_overflow = True
				break
			
		for key in self.draw_order:
			if has_overflow:
				self._alternate_plots[key] = self.MakeOverflowBinPlot(self._plots[key])
			else:
				self._alternate_plots[key] = self._plots[key].Clone()
		self._alternate_plots.Modified()
		
		return has_overflow
	
