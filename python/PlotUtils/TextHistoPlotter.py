"""
 module: TextHistoPlotter
 author: J. Wolcott
 
  Plot the values stored in a CSV-style text file
  in a (shudder) ROOT histogram (TH1F).
"""

import sys
import os.path
import optparse

import ROOT

COLUMN = 2
OUT_FILE = "text_histos.root"

class TextHistoPlotter:
	""" Class for plotting the contents of a text .csv file in a ROOT histogram. """
	def __init__(self, col=COLUMN, sep=None, overwrite=True):
		self.col = col
		self.sep = sep
		
		self.overwrite = overwrite
		
	def _ParseLine(self, line, line_num, col):
		""" Internal method for parsing a line of a .csv file. """
		line = line.strip()
		if len(line) == 0:
			return None

		# comments
		if line[0] in ('#', '"'):
			return None
		
		columns = line.split(self.sep)
		
		if len(columns) < col:
			raise Exception("Not enough columns in line %d!" % line_num)
		
		return columns
		
	def PlotFile(self, in_file, col=None, out_file=None, histo_title=None, axis_labels=None, histogram=None):
		""" Plot the file.
		
		Arguments:
		 - in_file (mandatory) -- Which file to plot.
		 - col -- Override the default 'col' parameter (specified in the constructor,
		on which see documentation).
		 - histo_title -- Title for the ROOT histogram.
		 - axis_labels -- List of titles for the axis labels of the ROOT histogram.
		 - histogram -- If you pass in a ROOT histogram via this parameter,
		                it will be Fill()ed with the contents of the .csv file.
		                Otherwise, a new histogram will be created.
		"""
		title = histo_title if histo_title is not None else in_file
		col = col or self.col
		
		if histogram is None:
			histo_name = os.path.basename(in_file)
			histogram = ROOT.TH1F(histo_name, title, 100, 0, 0)   # request automatic binning

			if axis_labels is not None:
				axis_getters = {
					0: histogram.GetXaxis,
					1: histogram.GetYaxis,
				}
				for i, label in enumerate(axis_labels.split(";")):
					axis_getters[i]().SetTitle(label)
		
		with open(in_file, "r") as text_file:
			for line_num, line in enumerate(text_file):
				columns = self._ParseLine(line, line_num, col)
				if columns is None:
					continue
					
				histogram.Fill(float(columns[col-1]))
		
		if out_file is None and __name__ == "__main__":
			out_file = OUT_FILE
		
		if out_file is not None:
			out_file = ROOT.TFile(out_file, "RECREATE" if self.overwrite else "UPDATE")
			histogram.Write(histogram.GetName(), ROOT.TObject.kOverwrite)
			out_file.Close()
		else:
			return histogram
			
	def GetVals(self, in_file, col=None):
		""" Gets the values in the selected column as a GENERATOR. """
		
		col = col or self.col
		
		with open(in_file, "r") as text_file:
			for line_num, line in enumerate(text_file):
				columns = self._ParseLine(line, line_num, col)
				if columns is None:
					continue
					
				yield float(columns[col-1])

###############################
# handle command-line arguments
if __name__ == "__main__":
	ROOT.gROOT.SetBatch(True)
	
	help_preamble = """
%prog: Plot values from CSV-style text file(s) in ROOT histogram.

Usage: %prog <options> FILE1 [FILE2 ...] <options>
"""

	parser = optparse.OptionParser(usage=help_preamble)
	parser.add_option(
		"-o", "--out-file",
		dest="out_file",
		action="append",
		help="Save histogram to file 'FILE'.  Default: save all histograms to '%s' in current directory." % OUT_FILE,
		metavar="FILE",
	)
	parser.add_option(
		"", "--col",
		dest="col",
		action="store",
		help="Use column number COL (indexed from 1) as the column containing the data.  Default: %default.  (See also '--sep'.)",
		metavar="COL",
		default=COLUMN
	)
	parser.add_option(
		"", "--sep",
		dest="sep",
		action="store",
		help="Use 'SEP' as the column separator character(s).  Default: any whitespace character.  (See also '--col'.)",
		metavar="SEP",
	)
	parser.add_option(
		"", "--histo_title",
		dest="histo_title",
		help="Histogram title.  Default: use input filename.",
		action="store"
	)
	parser.add_option(
		"", "--axis_labels",
		dest="axis_labels",
		help="Axis label(s).  Specify in format: 'axis1;axis2'.",
		action="store"
	)
	 
	options, args = parser.parse_args()
	
	if len(args) == 0:
		parser.print_help()
		print "Must specify at least one file to read in."
		sys.exit(1)
	
	histo_opt_names = ["out_file", "histo_title", "axis_labels"]
	
	for attr in histo_opt_names:
		if hasattr(options, attr) and len(getattr(options, attr)) != len(args):
			parser.print_help()
			print "If option '--%s' is specified, you must specify one instance for each input file." % attr
			sys.exit(1)
	

	# go!
	plotter_args = {}
	for arg in ("sep", "col"):
		if hasattr(options, arg):
			plotter_args[arg] = options[arg]
	plotter = TextHistoPlotter(**plotter_args)

	for i, filename in enumerate(args):
		plotter_args = {"in_file": filename}
		for opt in histo_opt_names:
			if len(getattr(options, opt)) > 0:
				plotter_args[opt] = getattr(options, opt)[i]
		
		plotter.PlotFile(**plotter_args)
