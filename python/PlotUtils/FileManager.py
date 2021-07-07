"""
  FileManager.py:
   Contains classes for managing input and output files,
   with preference towards ROOT files.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    February 2012
"""

import os.path

import ROOT

class FileManager(object):
	"""  Tool for managing input files and output file locations.
	
	The FileManager can be used for a few different purposes:
	 (1) It can take an input file list of ROOT ntuples and bundle
	     them together into a ROOT TChain.  (A convient tool to use
	     for making this list is PlotUtils.FileList, on which see
	     its own documentation.)
	 (2) It can be configured with an output file location and used
	     to generate file paths for files to be written to that location.
	These mechanisms can be used both separately and in conjunction with each other.
	
	Note: this design is unfortunate and needs to be revisited.
	"""
	def __init__(self, tree_names=None, ntuples=None, output_location=None, friend_trees = True, master_tree = None):
		""" Constructor.
		
		Options:
		 - tree_names -- List of the names of the TTrees within the files you're loading that
		                 you want access to.
		 - ntuples -- List of the paths to the ntuples you want to load into a TChain.
		 - output_location -- The directory where you would like any output files written.
		 - friend_trees -- Whether to attach extra trees to the master tree (see next option)
		                   as friend trees (see the ROOT TTree documentation).
		 - master_tree -- When multiple tree_names are given, one of them is loaded first,
		                  and the others are added as "friend" chains.  (See the TChain
		                  documentation on the ROOT web site.)  Which one should be the
		                  first can be specified here.  If not specified, the first item
		                  in 'tree_names' is used.
		
		If you wish to use the ROOT ntuple list -> TChain functionality, you need to
		specify at least 'tree_names' (a list) and 'ntuples'.
		
		If you only wish to use output filename construction, you need only specify
		'output_location'.
		"""
		self.chains = {}
		self._chain_loaded = False
		self._master_chain = None
		self._use_friends = friend_trees
		
		self.tree_names = tree_names
		self.ntuples = ntuples
		self.output_location = output_location
		
		if output_location is not None and not os.path.isdir(output_location):
			raise OSError("Specified output path '%s' is not accessible" % output_location)
		
		if master_tree is not None:
			self.master_tree = master_tree
		elif self.tree_names is not None:
			self.master_tree = self.tree_names[0]
		else:
			self.master_tree = None

	def GetChain(self, tree_names=None, overwrite=False, force_load=False):
		""" Retrieve the master chain. Load the ntuples if necessary.
		
		You can override the tree_names specified in the constructor if you wish.
		If the parameter 'overwrite' is False (the default), duplicate chain names
		will result in an exception being raised.
		"""
		
		if self._chain_loaded and not force_load:
			return self._master_chain
		
		tree_names = tree_names or self.tree_names
		if tree_names is None:
			return None
			
		for chain in tree_names:
			if (chain not in self.chains) or overwrite:
				self.chains[chain] = ROOT.TChain(chain)
				
				for filename in self.ntuples:
#					print "adding file to chain:", filename
					self.chains[chain].AddFile(filename)
			else:
				raise Exception("Not overwriting chain '%s'" % chain)
		
		# clear the list of friends so we can start fresh	
		self._master_chain = self.chains[self.master_tree]
		friends = self._master_chain.GetListOfFriends()
		if friends:
			for f in friends:
				self._master_chain.RemoveFriend(f)
		
		# now add all the other chains as friends if requested
		if self._use_friends:
			for chain in self.chains:
				if self.chains[chain].GetName() == self._master_chain.GetName():
					continue
			
				self._master_chain.AddFriend(chain)
		
		self._chain_loaded = True
		
		return self._master_chain
	
	def GetOutputPlotName(self, name, ext="eps"):
		""" Create a file path for a plot with given extension.
		
		Uses the output_location defined in the constructor (on which
		see documentation) to generate a fully-qualified path to a
		plot file, using the given extension.
		"""
		
		return self.GetFileName("%s.%s" % (name, ext))
	
	def GetFileName(self, name):
		""" Create a file path for a file with given name.
		
		Uses the output_location defined in the constructor (on which
		see documentation) to generate a fully-qualified path to a
		file with the specified name.
		"""
	
		return os.path.join(self.output_location, name)
