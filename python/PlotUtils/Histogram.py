"""
  Histograms.py:
   Wrapper class for ROOT histograms that does some useful stuff.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    February 2012
"""

import array
import inspect
import itertools
import numbers
import operator
import pprint
import sys
import time
import traceback

import ROOT

# turn off the stupid ROOT "file directory" thing
# that always results in missing histograms when
# TFiles are closed (see the ROOT documentation
# for TH1, under the section "Creating Histograms"
# http://root.cern.ch/root/html/TH1.html)
ROOT.TH1.AddDirectory(False)

# also ALWAYS store the sum of the squares of the errors.
# this enables the histogram to maintain correct errors when scaled (etc.).
ROOT.TH1.SetDefaultSumw2(True)

# for now, we just use 'integer' and 'double' histograms
ROOT_HISTO_TYPES = {
	int: "I",
	float: "D",
}

AXIS_MAP = {
	1: "X",
	2: "Y",
	3: "Z",
}

def GetAxis(obj, identifier):
	""" Returns the axis specified by the identifier.
	
	Axes can be specified by letter ('x', 'y', 'z') or number (1, 2, 3).
	"""
	if identifier in AXIS_MAP:
		identifier = AXIS_MAP[identifier]

	if isinstance(identifier, basestring):
		try:
			return getattr(obj, "Get%saxis" % identifier.upper())()
		except AttributeError:
			raise TypeError("Can't get axis from %s (object: %s)" % (type(obj), str(obj)))
	else:
		raise ValueError("Don't understand axis identifier: '%s'" % str(identifier))

# use multiple inheritance
# to create a class that inherits
# from both _Histogram and the
# appropriate TH??.  return an instance
# of that class.
def BuildNew(name, title, histo_type, nbins, bins=None, ranges=None, base_root_class=None):
	""" Build a Histogram object that inherits from ROOT.TH1 and _Histogram.
	
	This method uses the magic of Python's dynamic type mechanism
	to create a custom Histogram class that inherits both from
	the correct ROOT.TH?? derivative and from the _Histogram
	class within this module.  By so doing, we need only write
	ONE wrapper for all of the various TH1, TH2, TH3 types that
	can imbue all of them with the extra convenience methods
	of _Histogram.
	"""
	
	if not hasattr(nbins, "__len__"):
		raise TypeError("nbins must be specified as a list.  For 1D plots, use the syntax [ nbins, ] (note the comma).")
	
	final_bins = [None,] * len(nbins)
	final_ranges = [None,] * len(nbins)

	if bins is None:
		final_ranges = ranges
	elif ranges is None:
		final_bins = bins
	else:
		for i in range(len(nbins)):
			if bins[i] is not None and ranges[i] is None:
				final_bins[i] = bins[i]
			elif bins[i] is None and ranges[i] is not None:
				final_ranges[i] = ranges[i]
			elif bins[i] is None and ranges[i] is None:
				raise ValueError("Histogram axis must specify either axis ranges (for fixed size bins) or bin lower limits (for variable bin sizes)...")
			else:
				print "bins:", bins
				print "ranges:", ranges
				raise ValueError("Histogram axis can't specify both axis ranges and bin lower limits...")

	if base_root_class is not None:
		root_class = base_root_class
	else:
		histo_obj = "TH"
		histo_obj += str(len(nbins))
		histo_obj += ROOT_HISTO_TYPES[histo_type]
	
		root_class = getattr(ROOT, histo_obj)
	
	class Histogram(root_class, _Histogram, _ROOTHistoMixin):
		def __init__(self):
			
			_Histogram.__init__(self, histo_type, nbins, final_ranges, final_bins)
			
			bin_parameters = []
			for i, n in enumerate(nbins):
				bin_parameters.append(n)
				if final_ranges[i] is None:
					try:
						bin_parameters.append( array.array('d', final_bins[i]) )
					except TypeError:
						print >> sys.stdout, "Couldn't understand your variable-bin specification for histogram '%s': '%s'.  Did you forget the brackets for a 1D histogram when specifying the bins?  (Syntax: [ [bin0, bin1, bin2, ...], ]  -- note the extra comma after the interior close bracket.)" % (name, final_bins[i])
						raise
				else:
					map(bin_parameters.append, final_ranges[i])
			root_class.__init__(self, name, title, *bin_parameters)

	return Histogram()

def BuildFromROOT(root_histogram):
	""" Wrap a ROOT histogram with the _Histogram extensions. """
	if not isinstance(root_histogram, ROOT.TH1):
		raise TypeError("Object of type %s is not a ROOT histogram!" % type(root_histogram))

	# if we've already wrapped this thing once, we don't need to do it again
	if isinstance(root_histogram, _Histogram):
		return root_histogram

	# determine the dimension.
	# (count backwards because ROOT's screwy inheritance model
	#  has TH2 and TH3 inheriting from TH1)
	for dim in range(3, 0, -1):
		if isinstance(root_histogram, getattr(ROOT, "TH%d" % dim)):
			dimension = dim
			break

	name = root_histogram.GetName()
	title = root_histogram.GetTitle()
	
	nbins = []
	bins = []
	ranges = []

	# determine the type.  it's a float unless it's an integer...
	histo_type = float
	integer_types = (ROOT.TH1I, ROOT.TH2I, ROOT.TH3I)
	if isinstance(root_histogram, integer_types):
		histo_type = int
	
	# discover the bin widths.
	for axis_num in range(1, dimension+1):
		axis = GetAxis(root_histogram, axis_num) #getattr(root_histogram, "Get%saxis" % AXIS_MAP[axis_num])()
		nbins.append(axis.GetNbins())
		if axis.IsVariableBinSize():
			ranges.append(None)
			bins.append(axis.GetXbins())
		else:
			bins.append(None)
			ranges.append( (axis.GetBinLowEdge(1), axis.GetBinUpEdge(axis.GetNbins())) )	# remember, the 0 bin is underflow and the N+1 bin is overflow...

	root_class = type(root_histogram)
	class Histogram(root_class, _Histogram, _ROOTHistoMixin):
		def __init__(self, base_histogram):
			# super() resolves left to right in the inheritance model.
			# for some reason that I haven't bothered to figure out,
			# 'root_class.__init__' doesn't work here.
			super(Histogram, self).__init__(base_histogram)	
			_Histogram.__init__(self, histo_type, nbins, ranges, bins)

##### the check below was crashing the RecoStabilityPlots...
	# we want to ensure that we're only wrapping a ROOT histogram.
	# check that the type reported is a ROOT type...
	# (for some reason you can't get TProfile in the list the normal way, though...)
#	if not (root_class in [c[1] for c in inspect.getmembers(ROOT, inspect.isclass)] + [ROOT.TProfile,]):
#		print inspect.getmembers(root_histogram)
#		raise Exception("Can't wrap object of type \"%s\" with PlotUtils.Histogram" % root_class)
			
	return Histogram(root_histogram)
			
class _ROOTHistoMixin(object):
	""" Wrapper properties around ROOT histogram methods.
	    No new functionality; just for convenience. """
	@property
	def name(self):
		return self.GetName()
	
	@property
	def title(self):
		return self.GetTitle()
	

class _Histogram(object):
	""" The customizations to be bolted on top of the ROOT histogram. """
	def __init__(self, type, nbins, ranges, bins):

		self.type = type
		self.nbins = nbins
		self.bins = bins
		self.ranges = ranges
		
		self._axes_normalized = False
		
	def __eq__(self, other):
		if isinstance(other, basestring):
			return self.name == other
		elif not hasattr(other, "name"):
			return False
			
		# this is risky...
		return self.name == other.name
		
	def __repr__(self):
		return "%d-D histogram, name '%s', title '%s'" % (self.Dimension(), self.name, self.title)
		
	def Dimension(self):
		return len(self.nbins)

	@property
	def dimension(self):
		return self.Dimension()

	# i.e., import the global GetAxis() function into this class's scope
	GetAxis = GetAxis
		
	def Normalize(self, normalize_by=None):
		""" Normalize the histogram.
		
		How this happens depends on 'normalize_by':
	      - if 'normalize_by' is None, no normalization is performed.
	      - if 'normalize_by' is the string "self", the histogram is 
	        normalized to area 1.
	      - if 'normalize_by' is a number, the histogram is normalized
	        by dividing by that number.
	      - furthermore, if the histogram is two-dimensional:
	        * if 'normalize_by' is the string "col",
	          each column will be normalized so it sums to 1
	        * if 'normalize_by' is the string "row",
	          each row will be normalized so that it sums to 1
		"""

		denominator = None
		sum_axis = None
		
		if isinstance(normalize_by, numbers.Number):
			denominator = normalize_by
		elif normalize_by is None:
			return
		elif normalize_by.lower() == "self":
			denominator = self.Integral()
		elif normalize_by.lower() == "row":
			if self.dimension != 2:
				raise TypeError("Can only row-normalize two-dimensional histograms.  This histogram is of dimension: %d" % self.dimension)
			sum_axis = "X"
		elif normalize_by.lower() == "col":
			if self.dimension != 2:
				raise TypeError("Can only column-normalize two-dimensional histograms.  This histogram is of dimension: %d" % self.dimension)
			sum_axis = "Y"
		else:
			print "I don't understand your normalization denominator: '%s'" % normalize_by
		
		if denominator is None and sum_axis is not None:
			other_axis = list(set(["X","Y"]) - set([sum_axis,]))[0]
#			print "sum_axis: %s; other_axis: %s" % (sum_axis, other_axis)
			for i in range(1, getattr(self, "Get%saxis"%other_axis)().GetNbins()+1):
				bin_contents = [self.GetBinContent( *([i,j] if sum_axis=="Y" else [j,i]) ) for j in range(1, getattr(self, "Get%saxis"%sum_axis)().GetNbins()+1)]
#				print normalize_by, i, " bin contents:", bin_contents 
				norm = sum(bin_contents)  # add together the contents of all the bins along the appropriate axis
				if norm == 0:
					continue
				for j in range(1, getattr(self, "Get%saxis"%sum_axis)().GetNbins()+1):
					coords = [i,j] if sum_axis=="Y" else [j,i]
					args = coords + [self.GetBinContent(*coords)/norm,]
					self.SetBinContent(*args)
		elif denominator > 0:
			self.Scale(1./denominator)
		
			
	def NormalizeVariableAxes(self):
		""" Scales bins in histograms with variable bins
		    so that they show the appropriate multiple of 
		    the minimum bin width times their entries.  """

		if self._axes_normalized:
			return
			
		# don't re-normalize non-variable bins...
		if not self.GetXaxis().IsVariableBinSize():
			return
		
		# TH1.Scale now handles this correctly?
		self.Scale(1.0, "width")
		self._axes_normalized = True

#		if self._axes_normalized:
#			return
#			
##		print "normalizing..."
#	
#		for axis_num, bins in enumerate(self.bins):
#			if bins is None:
#				continue
#			axis = self.GetAxis(axis_num+1)
#			min_bin_width = min( (axis.GetBinWidth(bin) for bin in range(axis.GetNbins())) )
#			
#			if min_bin_width <= 0:
#				raise ValueError("Histogram with nonsensical minimum bin width: %f" % min_bin_width)

#			axis_bins = [ self.GetAxis(i+1).GetNbins() for i in range(self.Dimension()) ]
##			print axis_bins
#			if len(axis_bins) > 1:
#				total_n_bins = 1
#				for n_bins in axis_bins:
#					total_n_bins *= n_bins
#			else:
#				total_n_bins = axis_bins[0]
#			for i in xrange(1,total_n_bins+1):
#				scale_factor = 1. / (self.GetBinWidth(i) / min_bin_width)
#				if scale_factor == 1:
#					continue
##					print scale_factor
#				self.SetBinContent(i, self.GetBinContent(i) * scale_factor)

#		self._axes_normalized = True

	SIXTY_EIGHT_PERCENT_COVERAGE = [0.159, 0.841]
	def ProfileQuantile(self, axis, graph_name=None, quantile_divisions=SIXTY_EIGHT_PERCENT_COVERAGE):
		""" Since a TProfile is a TH1 under the hood, it can't handle asymmetric errors.
		    This is a self-rolled equivalent which does.
		    
		    The default quantile divisions get you the 68% confidence interval around the mean.
		    
		    If you don't specify the output graph_name,
		    the default is to append "_quantile_errors" to this histogram's name. """
		
		if self.dimension != 2:
			raise TypeError("Can't profile a histogram with other than 2 dimensions.  Yours is a %d-D histogram..." % self.dimension)
			
		if len(quantile_divisions) != 2:
			raise ValueError( "quantile_divisions must be an iterable of length 2.  You gave me a %s of length %d..." % (type(quantile_divisions), len(quantile_divisions)) )

		# note: we project onto the *other* axis so we can get the statistical properties of its distribution and use them for mean and error bars
		other_axis = list( set(["x","y"]) - set(axis.lower()) )[0]
		projection_method = "Projection%s" % other_axis.capitalize()  
		if not hasattr(self, projection_method):
			raise TypeError("This histogram doesn't have a '%s' method, and I need it to profile the %s axis." % (projection_method, axis))

		axis = self.GetAxis(axis)
		other_axis = self.GetAxis(other_axis)
		
		points = []
		errors_x = []
		errors_y = []
		quantiles = array.array('d', [0,]*len(quantile_divisions))
		quantile_divisions = array.array('d', quantile_divisions)
		for bin_num in range(1, axis.GetNbins()+1):  # remember, 0 is underflow and Nbins+1 is overflow...
			proj = getattr(self, projection_method)("%s_tmp_%d" % (self.name, time.time()), bin_num, bin_num)  # random name, just to be safe
			
			x = axis.GetBinCenter(bin_num)
			y = proj.GetMean()
			points.append( (x,y) )
		
			# for now, make the errors for an empty projection 0
			if proj.Integral() == 0:
				n_qs = 2
				quantiles[0] = quantiles[1] = 0
			else:
				n_qs = proj.GetQuantiles(len(quantiles), quantiles, quantile_divisions)
			if n_qs != 2:
				raise RuntimeError("GetQuantiles() didn't compute all the quantiles in PlotUtils.Histogram._Histogram.ProfileQuantile().  I requested 2 quantiles, but got %d of them..." % n_qs)
			
			bin_width = axis.GetBinWidth(bin_num)
			errors_x.append( (bin_width/2., bin_width/2.) )
			errors_y.append( (y-quantiles[0], quantiles[1]-y) )
			
#			print points
#			print errors_x
#			print errors_y

		out_graph = ROOT.TGraphAsymmErrors(
			len(points),
			array.array('d', [operator.getitem(i, 0) for i in points]), # x-values
			array.array('d', [operator.getitem(i, 1) for i in points]), # y-values
			array.array('d', [operator.getitem(i, 0) for i in errors_x]), # low-sides of x errors
			array.array('d', [operator.getitem(i, 1) for i in errors_x]), # high-sides of x errors
			array.array('d', [operator.getitem(i, 0) for i in errors_y]), # low-sides of y errors
			array.array('d', [operator.getitem(i, 1) for i in errors_y]), # high-sides of y errors
		)
		name = graph_name or "%s_quantile_errors" % self.name
		out_graph.SetName(name)
			
		return out_graph			

	def ProfileQuantileX(self, graph_name=None, quantile_divisions=SIXTY_EIGHT_PERCENT_COVERAGE):
		return self.ProfileQuantile("x", graph_name, quantile_divisions)

	def ProfileQuantileY(self, graph_name=None, quantile_divisions=SIXTY_EIGHT_PERCENT_COVERAGE):
		return self.ProfileQuantile("y", graph_name, quantile_divisions)
		
class HistoStack(ROOT.THStack):
	""" A wrapper around ROOT.THStack which adds some dictionary-like functionality.
	
	You can use this class to grab a histogram from a THStack as if
	it were a dictionary:
	  stack["histo"]
	will return the histogram with name=="stack" from the stack.
	"""
	def __init__(self, name="", title=""):
		super(HistoStack, self).__init__(name, title)
		
		self._key_map = {}
		
	def __contains__(self, key):
		if len(self._key_map) > 0:
			return key in self._key_map
		
		hists = self.GetHists()
		if hists:
			return bool(hists.FindObject(key))
		
		return False
		
#	def __del__(self):
#		if self:
#			idx = ROOT.gROOT.GetListOfCleanups().IndexOf(self)
#			if not idx < 0:
#				ROOT.gROOT.GetListOfCleanups().RemoveAt(idx)
#			while len(self.GetHists()) > 0:
#				self.GetHists().RemoveLast()
#			print self.GetHists().GetSize()
#			print [i for i in self.GetHists()]
	
	def __delitem__(self, key):
		if key not in self:
			raise KeyError(key)
		
		obj = self._key_map[key].pop()
		self.GetHists().Remove(obj)
	
	def __iter__(self):
		if len(self._key_map) > 0:
			for h in self._key_map:
				yield h
		else:
			hists = self.GetHists()
			if not hists:
				hists = {}
			for h in hists:
				yield h.GetName()
	
	def iteritems(self):
		# if the user has added items by key,
		# use those.
		if len(self._key_map) > 0:
			for k, h in self._key_map.iteritems():
				yield k, h
		else:
			hists = self.GetHists()
			if not hists:
				hists = {}

			# otherwise, use the histogram names.
			for h in hists:
				yield h.GetName(), h
	
	iterkeys = __iter__
	
	def itervalues(self):
		if len(self._key_map) > 0:
			for h in self._key_map.itervalues():
				yield h 
		else:
			hists = self.GetHists()
			if not hists:
				hists = {}

			for h in hists:
				yield h
	
	def __getitem__(self, key):
		if key in self._key_map:
			return self._key_map[key]
			
		hists = self.GetHists()
		if hists:
			h = self.GetHists().FindObject(key)
			if h:
				return h
		
		raise KeyError("HistoStack has no histogram '%s'" % key)
		
	def keys(self):
		return [k for k in self.iterkeys()]

	def __len__(self):
		if len(self._key_map) > 0:
			return len(self._key_map)
		
		return len(self.GetHists())
		
	def __setitem__(self, key, h):
#		print "Want to set key '%s' with obj:" % key, h
#		print " Current key list:"
#		print pprint.pformat(self._key_map)
#		traceback.print_stack()
		
		# need the inverse key map so that we can rebuild later
		inv_key_map = dict((v,k) for k, v in self._key_map.iteritems())
		
		objs_to_replace = set()
		if key in self._key_map:
			objs_to_replace.add(self._key_map[key])
		if self.GetHists() and h in list(self.GetHists()):
			objs_to_replace.add(h)
		
#		print objs_to_replace
		
		added = False
		new_objs = []
		if len(objs_to_replace) > 0:
			# this complicated replacement procedure is necessary
			# because it is the only way to insert a new histogram
			# into the middle of the TList while retaining the option 
			# associated with the old one. boo.
			current = ROOT.TIter(self.GetHists())
			obj = current()
			while obj:
#				print "considering obj:", obj

				if key in self._key_map:
					del self._key_map[key]
			
				if obj in objs_to_replace:
#					print "replacing obj:", obj
					# replace the old one of this key...
					if obj in inv_key_map:
						del inv_key_map[obj]
					#self.Add(h, key=key, option=current.GetOption())
					new_objs.append( {"obj": h, "key": key, "opt": current.GetOption()} )
					added = True
				else:
					# ... but put the other objects back in as they were
					k = inv_key_map.get(obj, obj.GetName() if hasattr(obj, "GetName") else None)
					#self.Add(obj, key=key, option=current.GetOption())
					new_objs.append( {"obj": obj, "key": k, "opt": current.GetOption()} )

				obj = current()
		else:
			self.Add(h, key=key, option=h.GetOption() if hasattr(h, "GetOption") else "")
			return

		if self.GetHists():
			self.GetHists().Clear()
		self._key_map.clear()
		for h_info in new_objs:
			self.Add(h_info["obj"], key=h_info["key"], option=h_info["opt"])
	
	def values(self):
		return [v for v in self.itervalues()]
	
	def Add(self, histogram, key=None, option=""):
#		print "Adding histogram with key '%s':" % key, histogram
		if len(self._key_map) > 0 and key is None:
			raise ValueError("You have already specified a key for a histogram.  Therefore, you must specify keys for ALL histograms in this stack.")
		
		# don't let a histogram be added twice.  this will foul up __setitem__.
		if histogram in self._key_map.values() or histogram in (list(self.GetHists()) if self.GetHists() else []):
#			print self._key_map.values()
			raise ValueError("Histogram '%s' is already in this stack.  Can't add again!" % histogram)
		
		if len(option) == 0 and hasattr(histogram, "GetOption"):
			option = histogram.GetOption()
		
		if key is not None:
			self._key_map[key] = histogram
		
#		print "adding histogram with option: '%s'" % option
		super(HistoStack, self).Add(histogram, option)

	def GetAxis(self, identifier):
		if identifier in AXIS_MAP:
			identifier = AXIS_MAP[identifier]

		if isinstance(identifier, basestring):
			try:
#				print "Getting %s axis" % identifier
				return getattr(self, "Get%saxis" % identifier.upper())()
			except AttributeError:
				raise TypeError("Can't get axis from %s (requested axis: %s)" % (type(self), str(identifier)))
		else:
			raise ValueError("Don't understand identifier: '%s'" % str(identifier))
