"""
  FileList.py:
   Class to load and store a list of files from a source text file.
   
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2012
"""

import os.path

class FileList:
	""" Class for loading a list of files from a text file index. 
	
	The text file must contain one file path per line.
	Leading and trailing spaces, as well as blank lines
	and lines beginning with '#' are ignored."""
	def __init__(self, source_file, load_now=True, ignore_missing=False):
		self._source_fname = source_file
		self._ignore_missing = ignore_missing
		
		self._filenames = None
		self._loaded = False
		
		self.max_testing_files = 25
		
		if load_now:
			self.Load()
	
	def GetFiles(self, testing=False, max_testing_files=None):
		""" Build the file list.
		
		If 'testing' is True, then no more than 'max_testing_files'
		will be loaded, regardless of how many are listed in the source file.
		"""
		self.Load()
		
		max_testing_files = max_testing_files or self.max_testing_files
		
		# make an ordered list so as to be determinate.
		l = sorted(self._filenames)
		
		if testing and max_testing_files is not None and len(l) > max_testing_files:
			return l[:max_testing_files]
		else:
			return l
	
	def Load(self):
		""" Load and parse the file index. """
		if self._loaded:
			return
			
		self._filenames = set()

#		print "opening from file:", self._source_fname
		with open(self._source_fname, "r") as source_file:
			for line in source_file:
				# comments, blank lines
				line = line.strip()
				if len(line) == 0 or line[0] in ('#', '\n'):
					continue

				line = os.path.expandvars(line)
				
				if not os.path.exists(line):
					if not self._ignore_missing:
						raise IOError("Cannot find file: %s" % line)
					else:
						print "Warning: could not find file:", line
						continue
				
#				print "  adding file to index:", line
				self._filenames.add(line)
		
		self._loaded = True
