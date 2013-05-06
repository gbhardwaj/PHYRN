#!/usr/bin/env python
# encoding: utf-8
desc = """Created by Gaurav Bhardwaj on 02/06/2011.
Copyright (c) 2012 Department of Biochemistry and Molecular Medicine, University of California, Davis. All rights reserved.
"""


import sys
import Bio
import datetime
from optparse import OptionParser
from Bio import SeqIO
from datetime import datetime

parser = OptionParser(description = "Chops fasta sequences at specified positions with optional additional buffer.", usage = "Usage: python %prog -i full_length_file -b boundary text file -o output file", epilog = desc)

parser.add_option("-i", action = "store", type = "string", dest = "Input_File", help = "Full-Length FASTA Sequence File")
parser.add_option("-b", action = "store", type = "string", dest = "Boundary_File", help = "List of Boundaries")
parser.add_option("-o", action = "store", type = "string", dest = "Output_file", help = "Chopped_file (Default: chopped_output.fa)", default = "chopped_output.fa")
parser.add_option("-x", action = "store", type = "int", dest = "Buffer", help = "Buffer Region (Default: 0)", default = 0)
(options, args) = parser.parse_args()
if len(sys.argv) < 2 :
	parser.print_help()
	sys.exit(0)

prog_start = datetime.now()
buffer = options.Buffer			#buffer size
seq_start = {}
seq_stop = {}

fout = open(options.Output_file,"w")	
fob = open(options.Boundary_File,"rU")
for line in fob:
	line = line.replace("\n","")
	items = line.split("\t")
	seq_start[items[0]] = int(items[1])
	seq_stop[items[0]] = int(items[2])
	
fob.close()
handle = open(options.Input_File,"rU")
for rdrp in Bio.SeqIO.parse(handle, "fasta"):
	#seq_dir[rdrp.id] = rdrp.seq
	seq_len = len(rdrp.seq)
	check = seq_start.get(rdrp.id)
	if check == None:
		print rdrp.id + " does not have boundaries defined"
	
	else:
		if seq_start[rdrp.id]-buffer >= 0:
			x = seq_start[rdrp.id]-buffer
		else:
			x = 0
		if seq_stop[rdrp.id]+buffer > seq_len:
			y = seq_len
		else:
			y = seq_stop[rdrp.id]+buffer
		print >> fout,">"+rdrp.id+"_chopped"
		print >> fout,rdrp.seq[x:y]

fout.close()
handle.close()
prog_stop = datetime.now()
print "Successfully Done! Program finished with duration of ", prog_stop - prog_start