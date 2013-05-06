#!/usr/bin/env python
# encoding: utf-8

desc = """
Created by Gaurav Bhardwaj on 02/06/2012.
Copyright (c) 2012 Department of Biochemistry and Molecular medicine, University of California, Davis. All rights reserved.

By the way, Universe is nothing but a giant pineapple!!! 

"""

import sys
import Bio
import os
from optparse import OptionParser

from Bio import SeqIO
from os import mkdir

parser = OptionParser(description = "Submits job array to make PSSMs", usage = "USAGE: python %prog -i input_fasta_file -d full_path_nr -n new_db_name -e e_value", epilog = desc)

parser.add_option("-i", action = "store", type = "string", dest = "Input_File", help = "Input FASTA Sequence File")
parser.add_option("-d", action = "store", type = "string", dest = "database", help = "Database path to search (Default: nr database)", default = "/user/local/bin/db/nr/nr")
parser.add_option("-n", action = "store", type = "string", dest = "newdb", help = "Name of New database")
parser.add_option("-e", action = "store", type = "string", dest = "evalue", help = "E-value (Default: 1e-6)", default = 1e-6)
parser.add_option("-r", action = "store", type = "int", dest = "iter", help = "Number of iterations (default: 6)", default = 6)
(options, args) = parser.parse_args()

if len(sys.argv) < 2 :
	parser.print_help()
	sys.exit(0)

db_name = options.newdb
e_value = options.evalue
db_desc = str(options.iter)+" iterations at "+str(options.evalue) #change here
dirname = db_name
os.mkdir("./"+ dirname+"/")
os.mkdir("./"+ dirname+"/"+ 'pssm')
os.mkdir("./"+ dirname+"/"+ 'cdd')
num = 00001
handle = open(options.Input_File, "rU") #change here
keytx = open("./"+ dirname+ "/" + "pssm_source.txt","w")
masters_name = db_name+"_masters.fa"
masters = open("./"+ dirname+ "/" + masters_name,"w")

for sequence in Bio.SeqIO.parse(handle,"fasta"):
	fname = db_name+str(num)+".fa"
	fob = open("./"+ dirname+"/"+ 'pssm'+ "/" + fname,"w")
	print >> fob,">"+db_name+str(num)+" "+sequence.id+" "+db_desc
	print >> fob, sequence.seq
	print >> keytx, db_name+str(num)+"\t"+ sequence.id
	print >> masters,">"+db_name+str(num)+" "+sequence.id+" "+db_desc
	print >> masters, sequence.seq
	num = num+1
	fob.close()

masters.close()
keytx.close()		
handle.close()
os.chdir("./"+ dirname+"/"+ 'pssm'+'/')

#Generate job file#

job_file = open("job_file.sh","w")
blastpgp_path = '/user/local/bin/blast-2.2.25/bin/blastpgp'
nr_path = options.database

print >> job_file,"#!/bin/bash"
print >> job_file,"#$"+" "+"-N"+" "+db_name+"_pssm_genration"
print >> job_file,"#$ -o psiblast.log" 
print >> job_file,"#$ -j y"
print >> job_file,"#$ -cwd"
print >> job_file,blastpgp_path+" "+"-d"+" "+nr_path+" "+"-i"+" "+db_name+"${SGE_TASK_ID}.fa"+" "+"-o"+" "+db_name+"${SGE_TASK_ID}.out"+" "+"-C"+" "+db_name+"${SGE_TASK_ID}.asnt"+" "+"-e"+" "+"1"+" "+"-m"+" "+"0"+" "+"-s"+" "+"T"+" "+"-a"+" "+"1"+" "+"-h"+" "+str(e_value)+" "+"-j"+" "+str(options.iter)+" "+"-M"+" "+"BLOSUM62"+" "+"-F"+" "+"F"+" "+"-u"+" "+"1"+" "+"-J"+" "+"T"

job_file.close()

##Submit job File##
lim = num-1
job_run_cmd = "qsub -p -5 -t 1:"+str(lim)+" job_file.sh"
os.system(job_run_cmd)
