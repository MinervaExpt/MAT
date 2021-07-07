"""
  Plotter.py:
   Plot manager class.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    June 2013
"""

import array
import datetime
import inspect
import os.path
import pprint
import re
import sys
import time
import traceback
import uuid

import ROOT

import PlotUtils.LoadPlotUtilsLib

from PlotUtils.FileList import FileList
from PlotUtils import FileManager
from PlotUtils import Plot


class Plotter(ROOT.TPySelector):
	def __init__(self, tree_names, file_list, output_location=None, name=None, testing=False, max_testing_files=25):
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

		self.name = uuid.uuid4() if name is None else name
		self.mangling_prefix = "%s_" % self.name  # see _DeManglePlotname() and _ManglePlotname().
		
		# --- methods used during the processing of the TTree/TChain.
		self._initialize_methods = []
		self._process_methods = []
		self._finalize_methods = []

		# --- which branches contain the variables you want to use?
		self._branches_to_enable = []
		
		self.testing = testing
		
		self.plot_collections = {}
		
		self.plot_types = [ "eps", ]
		
		self.SetInputFileList(file_list, tree_names, output_location, max_testing_files=max_testing_files)
		
		super(Plotter, self).__init__(self._file_mgr.GetChain(), self)
		self._file_mgr.output_location = output_location

		# the "option" needs to be set with the name of the Python module
		# containing the class derived from TPySelector for TPySelector
		# to be happy.  this is completely bizarre, but I can't find
		# any other way around it.
		self.SetOption(self.__module__)
		
		self._total_entries = None
		
	def __contains__(self, key):
		""" Simulate the 'in' operator.  See __getitem__() for syntax. """
		
		try:
			return bool(self._ParseKey(key))
		except KeyError:
			return False
		
	def __getitem__(self, key):
		""" Simulate dictionary lookup.
		
		Plots can be retrieved using the syntax 'collection/plotname'.
		
		Entire collections can be retrieved with just 'collection'.
		"""
		
		return self._ParseKey(key)
	
	def _DeManglePlotname(self, name):
		# de-mangle...
		return name.replace(self.mangling_prefix, "")
		
	def _ManglePlotname(self, name):
		if name.find(self.mangling_prefix) < 0:
			return self.mangling_prefix + name
		
	
	def _ParseKey(self, key):
		""" Parse and return a plot or plot collection
		    using a key of the form 'collection/plotname' or 'collection'. """
		
		split_key = key.split("/", 1)
		if len(split_key) == 1:
			if key in self.plot_collections:
				return self.plot_collections[key]
				
			raise KeyError("Plot must be specified as 'method/plotname'")
		
		split_key[1] = self._DeManglePlotname(split_key[1])
		return self.plot_collections[split_key[0]][split_key[1]]
		
	@property
	def output_location(self):
		return self._file_mgr.output_location

	def CreatePlot(self, *args, **kwargs):
		if "chain" not in "kwargs" and self._file_mgr.GetChain() is not None:
			kwargs["chain"] = self._file_mgr.GetChain()
			
		# a bit of name-mangling since every name is global in ROOT...
		if "plot_name" in kwargs and str(self.name) not in kwargs["plot_name"]:
			kwargs["plot_name"] = self._ManglePlotname(kwargs["plot_name"])
		
		return Plot.Plot(*args, **kwargs)
	
	def CreateOverlayPlot(self, *args, **kwargs):
		if "chain" not in "kwargs" and self._file_mgr is not None:
			kwargs["chain"] = self._file_mgr.GetChain()
		
		return Plot.OverlayPlot(*args, **kwargs)
		
	def LoopTree(self, tree=None):
		tree = tree or self._file_mgr.GetChain()
		
		# count the total number of entries NOW
		# so that we don't have to do it repeatedly later
		self._total_entries = int(tree.GetEntries())
		self._total_entries_div10 = int(self._total_entries/10.)
		
		# store start time to use for duration estimates
		self._start_time = time.time()
		self._last_update = 0

		# use TTree.Process with this object as the TSelector
		tree.Process(self)
		
	def _NotifyException(self, e):
		""" TPySelector's C++ backend doesn't print out the
		    entire stack trace when an exception happens.
		    This method does.
		    
		    Don't call this method unless an exception is
		    currently being handled, or it will crash the
		    program with a different exception (because
		    sys.exc_info() will return (None, None, None)). """
		    
		print "Traceback (most recent call last):"
		traceback.print_tb(sys.exc_info()[2])
		print "%s: %s" % (e.__class__.__name__, e)
	
	def RegisterMethods(self):
		""" Derived classes should override this method.
		    Add all the methods you want called to
		    self._initialize_methods, self._process_methods,
		    self._finalize_methods as appropriate. """
		pass
	
	def RegisterBranches(self):
		""" Derived classes should override this method.
		    Add all the branches in the ntuple that contain
		    the variables that you want to use. """
		pass
	
	#---------------
	# Methods overridden from TSelector
	#---------------
	def SlaveBegin(self, tree):
		# we want to ONLY use the relevant branches.
		# otherwise, ntuple loops are SLOOOOOOOW...
		self.RegisterBranches()
		tree.SetBranchStatus("*", 0)
		for branch in self._branches_to_enable:
#			print "enabling branch:", branch
			if "*" in branch or hasattr(tree, branch):
				tree.SetBranchStatus(branch, 1)
			else:
				print "Warning: selection '%s' matched no branches..."  % branch
#		map(tree.SetBranchStatus, self._branches_to_enable, len(self._branches_to_enable) * [1,])

		self.RegisterMethods()
		
		try:
			for method in self._initialize_methods:
				method(tree)
		except Exception as e:
			self._NotifyException(e)
			sys.exit(1)
			
		print "Looping over %d events:" % self._total_entries
		print " 0%",
		sys.stdout.flush()
	
	def Process(self, event):
		if self.GetEntry(event, 0) <= 0:
			return False

		N = self._total_entries
		n = self.fChain.GetChainEntryNumber(event)
		print '\r',
		
		now = time.time()
		if N < 10 or (n+1) % self._total_entries_div10 == 0 or now - self._last_update > 0.99:
			secs_elapsed = now - self._start_time
			frac_done = float(n) / N
			if N < 10:
				for i in range(n+2):
					print " %d/%d" % (n+1, N),
				sys.stdout.flush()
			else:
				curr_perc = int(float(n+1)/N * 100 + 0.5)  # the '+0.5' so as to get rounding to nearest int
				for i in range(0,curr_perc+10,10):
					print " %d%%" % i,
				sys.stdout.flush()
		
			if frac_done > 0:
				remain = datetime.timedelta(seconds = int(secs_elapsed / frac_done - secs_elapsed))
				print "   (Estimated completion in: %s)" % remain,
				sys.stdout.flush()
				
			self._last_update = now
		
		try:
			for method in self._process_methods:
#				print "executing method:", method
				method()
		
			return True
		except Exception as e:
			self._NotifyException(e)
			sys.exit(1)  # re-raising doesn't stop immediately
		
	def SlaveTerminate(self):
		print
		print "Finalizing:"
		
		try:
			for method in self._finalize_methods:
				method()
		except Exception as e:
			self._NotifyException(e)
			sys.exit(1)
		
		print "Finished."
	
	#---------------
	# This object's interface extensions
	#---------------
	
	def AddCollection(self, collection_name):
		""" Directly create a new plot collection with name 'collection_name'. """

		if collection_name not in self.plot_collections:
			self.plot_collections[collection_name] = PlotCollection()
		
		return self.plot_collections[collection_name]
		
	
	def AddPlot(self, collection_name, plot, plot_name=None):
		""" Directly adds a plot to the collection 'collection_name'
		    with name 'plot_name'. """
		     
		if collection_name not in self.plot_collections:
			self.plot_collections[collection_name] = PlotCollection()
		
		# de-mangle for easier access
		if hasattr(plot, "GetName") and plot_name is None:
			plot_name = self._DeManglePlotname(plot.GetName())
			
		self.plot_collections[collection_name].Add(plot, plot_name)
	
	def SetInputFileList(self, file_list, tree_names=None, output_location=None, friend_trees=True, max_testing_files=None):
		""" Set the list of files that should be used for input. 
		
		The list can be specified in 3 ways:
		 (1) A FileList object (see the documentation in FileList.py)
		 (2) A path to a file containing a list of filenames
		     that are in the format to be accepted by FileList
		 (3) A Python list of file names.
		"""

		files = []
		other_args = {"max_testing_files": max_testing_files} if max_testing_files else {}
		if isinstance(file_list, FileList):
			files = file_list.GetFiles(testing=self.testing, **other_args)
		if isinstance(file_list, basestring):
			files = FileList(file_list, ignore_missing=True).GetFiles(testing=self.testing, **other_args)
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
#			print "      (%s)" % f
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
#				print "      (%s)" % filename
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

		print "Saving %d plots in collection: '%s'" % (len(collection), collection_name)
				
		for plotname, plot in collection.iteritems():
#			print plotname, plot, plot._histos.values()
			# if this collection is composed of further collections,
			# recurse on down...
			if isinstance(plot, dict):
				self.SavePlotCollectionToDisk(plot, collection_name=collection_name, file_mgr=file_mgr, subdir=subdir, *extra_args, **extra_kwargs)
			else:
				print "  Saving plot '%s'..." % plotname
				
				if subdir and collection_name is not None:
					dirname = collection_name.lower()
					dirpath = os.path.join( (file_mgr or self._file_mgr).output_location, dirname)
					if not os.path.isdir(dirpath):
						print "Warning: output directory '%s' does not exist.  Attempting to create it." % dirpath
						os.mkdir(dirpath, 0755 )
						
					plotname = "%s/%s" % (dirname, plotname)

#				print "Saving plot '%s' (%s)" % (plotname, plot)
#				print "  collection_name:", collection_name, "; subdir:", subdir

				# we allow canvases in here too
				# in case there's a complicated construction
				if isinstance(plot, ROOT.TCanvas):
					self.SaveCanvasToDisk(plot, name=plotname, file_mgr=file_mgr)
				else:
					self.SavePlotToDisk(plot, name=plotname, file_mgr=file_mgr)
				
	def SavePlotCollectionsToDisk(self, subdir=False):
		""" Save all known collections of plots to disk.  (See SavePlotCollectionToDisk().)"""

		for collection in self.plot_collections:
			self.SavePlotCollectionToDisk(collection, subdir=subdir)
	

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
		
		super(PlotCollection, self).__delitem__(key)
		
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
