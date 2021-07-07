"""
  Plotter.py:
   Plot manager class.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    June 2012
"""

import re
import uuid
import pprint
import os.path

import ROOT

import PlotUtils.LoadPlotUtilsLib

from PlotUtils.FileList import FileList
from PlotUtils import FileManager
from PlotUtils import Plot


class Plotter(object):
	""" Class to manage plotting.
	
	    This class is designed to be subclassed.
	    Subclasses should implement a method called 'Draw<name>Plots'
	    (accepting at least one parameter besides 'self': the plot collection)
	    for each collection of plots they are drawing; then, to
	    draw them, instantiate an instance of the subclass,
	    and call 'MakePlots("name")'.  'Draw<name>Plots' should save
	    any plots to the collection it was passed:
	      collection['plot_name'] = plot
	     or
	      collection.Add(plot, name='plot_name') (where the 'name' argument is optional).
	    This mechanism allows for easy modularization of the plots to
	    be drawn (so that they can be drawn selectively with minimal effort).
	    
	    You can save individual plots to disk by calling SavePlotToDisk('method/name'),
	    where 'method' is the '<name>' from above, and 'name' is 'plot_name'
	    from above.  You can also save an entire collection by simply
	    using SavePlotCollectionToDisk('collection').  The type of file(s)
	    to be written out is specified by setting the extensions in
	    plotter.plot_types (the default is '.eps' files only).
	    
	    Plotter objects can also be used in conjunction with PlotUtils.FileManagers
	    to quickly make plots from TChains (hence most of the parameters
	    to the constructor).  Have a look at what FileManager can do to
	    get an idea.
	
	"""
	def __init__(self, tree_names=None, file_list=None, output_location=None, name=None, friend_trees=True, testing=False):
		""" Constructor.
		
		All of the arguments are optional.  They perform as follows:
		  - tree_names -- A list of the names of the ROOT TChains/TTrees you want the child FileManager to manage.
		  - file_list -- The name of a text file containing a list of input files you want
		                 the child FileManager to load.
		  - output_location -- The directory you want the child FileManager to write output files
		                       (including plots) to.
		  - name -- this Plotter's name.  A UUID will be used if none is supplied
		  - friend_trees -- whether to attach extra trees to the first one as 'friends'.
		  - testing -- Testing mode (True/False).  Passed along to the FileManager
		               (see the PlotUtils.FileManager documentation).
		"""
		self.testing = testing
		
		if file_list is not None:
			self.SetInputFileList(file_list, tree_names, output_location, friend_trees)
			self._file_mgr.output_location = output_location
		else:
			self._file_mgr = FileManager.FileManager(output_location=output_location)

		self.name = name or uuid.uuid4()
		
		self.plot_collections = {}
		
		self.plot_types = [ "eps", ]
		
	def __getitem__(self, key):
		""" Simulate dictionary lookup.
		
		Plots can be retrieved using the syntax 'method/plotname',
		where 'method' corresponds to the <name> in a 'Draw<name>Plots' method.
		"""
		try:
			method, plotname = key.split("/", 1)
		except ValueError:
			raise Exception("Plot must be specified as 'method/plotname'")

		return self.plot_collections[method][plotname]
	
	def AddPlot(self, collection_name, plot_name, plot):
		""" Directly adds a plot to the collection 'collection_name'
		    with name 'plot_name'. """
		     
		if collection_name not in self.plot_collections:
			self.plot_collections[collection_name] = {}
		
		self.plot_collections[collection_name][plot_name] = plot
	
	def CreatePlot(self, *args, **kwargs):
		if "chain" not in "kwargs" and self._file_mgr.GetChain() is not None:
			kwargs["chain"] = self._file_mgr.GetChain()
		
		return Plot.Plot(*args, **kwargs)
	
	def CreateOverlayPlot(self, *args, **kwargs):
		if "chain" not in "kwargs" and self._file_mgr is not None:
			kwargs["chain"] = self._file_mgr.GetChain()
		
		return Plot.OverlayPlot(*args, **kwargs)
	
	@classmethod
	def GetDrawMethods(cls):
		""" Return the list of 'draw' methods
		    ( which can be called using 'MakePlots("<name>") )'
		     that are defined for this class. """
		
		matched_methods = map(lambda x: re.match(r"Draw(\w+)Plots", x), dir(cls))
		return filter(None, (m.group(1) if m else m for m in matched_methods))
	
#	draw_methods = property(draw_methods)	

	def GetChain(self):
		""" Return the 'master' chain used by default for drawing plots. """
		return self._file_mgr.GetChain()
		
	def MakePlots(self, method, *args, **kwargs):
		""" Make (i.e., draw) the plots in a plot collection.
		
		    The mandatory 'method' argument must correspond
		    to the <name> in a method Draw<name>Plots,
		    defined when subclassing this class.
		    See the class documentation for more.
		    
		    You can pass extra arguments to the method by
		    simply including them in the call to MakePlots.
		"""
		
		if method not in self.plot_collections:
			self.plot_collections[method] = PlotCollection()
		
#		try:
		getattr(self, "Draw%sPlots" % method)(self.plot_collections[method], *args, **kwargs)
#		except AttributeError:
#			raise Exception("You need to define plotting method 'Draw%sPlots' before calling MakePlots with '%s' as an argument." % (method, method))
#		except TypeError:
#			raise Exception("Method 'Draw%sPlots' must be defined as accepting two arguments: 'self' and the plot collection it will store its plots in." % method)

	def SetInputFileList(self, file_list, tree_names=None, output_location=None, friend_trees=True):
		""" Set the list of files that should be used for input. 
		
		The list can be specified in 3 ways:
		 (1) A FileList object (see the documentation in FileList.py)
		 (2) A path to a file containing a list of filenames
		     that are in the format to be accepted by FileList
		 (3) A Python list of file names.
		"""

		files = []
		if isinstance(file_list, FileList):
			files = file_list.GetFiles(testing=self.testing)
		if isinstance(file_list, basestring):
			files = FileList(file_list).GetFiles(testing=self.testing)
		elif hasattr(file_list, "__iter__"):
			files = file_list
		else:
			raise TypeError("file_list must be a FileList, a file path (string) or a list of file names...")
		
		self._file_mgr = FileManager.FileManager(tree_names, files, output_location, friend_trees)

	def SaveCanvasToDisk(self, canvas, name=None, plot_types=None, file_mgr=None):
		""" Save the contents of a ROOT canvas to disk.
		
		    If the 'name' argument is omitted, canvas.GetName() will be used.
		    If 'plot_types' is omitted, the Plotter's defaults will be used.
		"""
		
		canvas.Update()
		ROOT.gSystem.ProcessEvents()
		
		file_mgr = file_mgr or self._file_mgr
		plot_types = plot_types or self.plot_types
		name = name or canvas.GetName()
		for plot_type in plot_types:
			filename = self._file_mgr.GetOutputPlotName(name, plot_type)
			d, f = os.path.split(filename)
			f, ext = os.path.splitext(f)
			f = "%s%s" % (f, ext)
			canvas.Print(os.path.join(d, f))

	
	def SavePlotToDisk(self, plot, name=None, plot_types=None, file_mgr=None):
		""" Save an image of a plot (either PlotUtils.Plot or a TH1 derivative) to disk.
		
		    If the 'name' argument is omitted, plot.GetName() will be used.
		    If 'plot_types' is omitted, the Plotter's defaults will be used.
		    If 'file_mgr' (a PlotUtils.FileManager instance determines
		      where to write the output) is omitted, the Plotter's default will be used.
		"""
	
		file_mgr = file_mgr or self._file_mgr
		plot_types = plot_types or self.plot_types

		if isinstance(plot, basestring):
			try:
				method, plotname = plot.split("/", 1)
				plot = self.plot_collections[method][plotname]
			except ValueError:
				raise Exception("Plot must be specified as 'method/plotname', not as '%s'" % plot)
			except KeyError:
				raise Exception("Can't find plot: '%s'" % plot)

		name = name or plot.GetName()

		if hasattr(plot, "SaveToDisk"):
			for plot_type in plot_types:
				filename = file_mgr.GetOutputPlotName(name, plot_type)
				plot.SaveToDisk(filename)
		elif hasattr(plot, "Draw"):
			canvas = ROOT.TCanvas()
			canvas.cd()
			plot.Draw()
			self.SaveCanvasToDisk(canvas, name, plot_types)
			
		else:
			print "Not sure what to do with object '%s' of %s: skipping." % (name, type(plot))
				
			
	def SavePlotCollectionToDisk(self, collection, collection_name=None, file_mgr=None, subdir=False, *extra_args, **extra_kwargs):
		""" Save a collection of plots (corresponding to a 'Draw<name>Plots' method) to disk. 
		
		   If the parameter 'file_mgr' is omitted, the Plotter's default will be used.
		   The parameter 'subdir' can be used to specify whether the filename
		     should contain a subdirectory corresponding to the plot collection's name.
		     
		   Any other parameters passed will be handed on to MakePlots, if it is called.
		"""

#		print pprint.pformat(self.plot_collections)
#		print pprint.pformat(collection)
		
		if isinstance(collection, basestring):
			collection_name = collection
			if collection not in self.plot_collections:
				if collection in self.GetDrawMethods():
					self.MakePlots(collection, *extra_args, **extra_kwargs)
				else:
					raise ValueError("Unknown plot collection: '%s'" % collection)
			collection = self.plot_collections[collection]
				
		for plotname, plot in collection.iteritems():
#			print plotname, plot, plot._histos.values()

			# if this collection is composed of further collections,
			# recurse on down...
			if isinstance(plot, dict):
				self.SavePlotCollectionToDisk(plot, collection_name=collection_name, file_mgr=file_mgr, subdir=subdir, *extra_args, **extra_kwargs)
			else:
				if subdir and collection_name is not None:
					plotname = "%s/%s" % (collection_name.lower(), plotname)

#				print "Saving plot '%s' (%s)" % (plotname, plot)
#				print "  collection_name:", collection_name, "; subdir:", subdir

				# we allow canvases in here too
				# in case there's a complicated construction
				if isinstance(plot, ROOT.TCanvas):
					self.SaveCanvasToDisk(plot, name=plotname, file_mgr=file_mgr)
				else:
					self.SavePlotToDisk(plot, name=plotname, file_mgr=file_mgr)
				
	def SavePlotCollectionsToDisk(self):
		""" Save all known collections of plots to disk.  (See SavePlotCollectionToDisk().)"""

		for collection in self.plot_collections:
			self.SavePlotCollectionToDisk(collection)
	

class PlotCollection(dict):
	""" A helper class that, in addition to normal dictionary function,
	    allows objects to be added without a key (which is then inferred
	    from the object's GetName() method). 
	    
	    This object also can keep a list of objects to be retained when
	    given an object that has a "GetListOfPrimitives" (like ROOT.TCanvas).
	    This helps to prevent ROOT and the Python interpreter trying to
	    delete the same object twice when it goes out of scope in Python... """
	
	def __init__(self):
		dict.__init__(self)
		
		self._retentions = {}
		
	def __delitem__(self, key):
		if key not in self:
			raise KeyError("Cannot find key '%s'" % key)
	
		if key in self._retentions:
			del self._retentions[key]
		
		super(PlotCollection, self).__delitem__(self, name)
		
	def __setitem__(self, name, plot):
		super(PlotCollection, self).__setitem__(name, plot)	
		
		self._RetainRelatedObjects(plot, name)
		
	def update(self, other_dict):
		""" Though it pains me, I need to override the dict update()
		    method since it bypasses __setitem__(). """
		
		if not hasattr(other_dict, "iteritems"):
			raise AttributeError("'other_dict' object passed to update() has no method 'iteritems'")
		
		for k, v in other_dict:
			self.__setitem__(k, v)
	    
	def Add(self, plot, name=None):
		if not name and not hasattr(plot, "GetName"):
			raise NameError("Object '%s' has no GetName() method.  You must pass 'name' to Add() as well.")
		name = name or plot.GetName()
		self.__setitem__(name, plot)	
		
	def Retain(self, obj):
		""" Add an object to the retention list directly.
		
		    The 'top name' in this case is None. """
		
		if None not in self._retentions:
			self._retentions[None] = []
		self._retentions[None].append(obj)
		
	def _RetainRelatedObjects(self, obj, topname):
		if hasattr(obj, "GetListOfPrimitives"):
			if topname not in self._retentions:
				self._retentions[topname] = []
			for relobj in obj.GetListOfPrimitives():
				if hasattr(relobj, "GetListOfPrimitives"):
					self._RetainRelatedObjects(relobj, topname)
				
				if relobj not in self._retentions[topname]:
#					print "Retaining object: '%s' under top name '%s'" % (relobj, topname)
					self._retentions[topname].append(relobj)	
