#!/usr/bin/env python
# encoding: utf-8
"""
Created by Zhenhai Zhang on 2009-07-09.
Copyright (c) 2009 Bioinformatics and Genomics@PSU. All rights reserved.
Usage: cscore_2_matrix.py -p pssm_file -q query_file -suffix suffix -f folder
"""
#import modules: 
import sys, csv, os
from Bio import SeqIO
from datetime import datetime
# variable definition
sep = "\t"

#
# argument process functions from here
#
def processParas(para_list, **keywords):
	# process the parameter information, all parameter start with "-paraname"
	# return a dictionary (paraName, paraValue)
	# remove the first parameter which is the program name
	para_list = para_list[1 :]
	kwgs, values = para_list[ :: 2], para_list[1 :: 2]
	
	if len(kwgs) != len(values):
		print "number of keywords and values does not equal"
		sys.exit(0)
	
	kwgs = map(lambda x : keywords[x[1 :]], kwgs)
	values = map(evalValues, values)
	
	return dict(zip(kwgs,values))
	
def evalValues(v):
	# Evaluate strings and return a value corresponding to its real type (int, float, list, tuple)
	try:	return eval(v)
	except:	return v

def getParas(my_dict, *args):
		if len(args) == 1:	return my_dict[args[0]]
		else:	return (my_dict[arg] for arg in args)
		

#
# files access functions start here
#

def getFiles(align_folder, suffix):
	all_files = os.listdir(align_folder)
	return [x for x in all_files if x.endswith(suffix)]


def get_all_ids(f):
	
	# process fasta format file and return all ids as a lit
	print "getting ids from file:", f
	
	total, result = 0, []
	for entry in SeqIO.parse(open(f, "rU"), "fasta"):
		result.append(entry.id)
		total += 1
	
	print "finished...%d in total..." %total
	return result

def init_score_dict():
	
	global pssm_list
	result = dict()
	
	for pssm in pssm_list:
		result[pssm] = 0.0
		
	return result

def swapext(ext):
	
	if ext == ".txt":
		return ".tab"
	elif ext == ".tab":
		return ".txt"
	
def get_aline_from_dict(score_dict, pssm_list):
	return [score_dict[x] for x in pssm_list]
	
def get_prompt_no(total):
	return int(round(float(total) / 10, 0))		# prompt every 

def check_unique(l, desc):
	print "checking uniqueness of %s" %desc
	for x, y in zip(l, l[ 1 : ]):
		if x == y:
			print "ERROR! non-unique %s found %s" %(desc, x)
			sys.exit(0)
def get_file_name_only(s):
	if s.find("/") > -1:
		s = s[s.rindex("/") + 1 :]
	return s
		
if __name__ == '__main__':
	start = datetime.now()
	print "program started at:", start
	#get arguments
	if len(sys.argv) < 9:
		print __doc__
		sys.exit(0)
	
	dict_args = processParas(sys.argv, p="pssm_file", q="query_file", suffix="suffix", f="folder")
	pssm_file, query_file, suffix, align_folder = getParas(dict_args, "pssm_file", "query_file", "suffix", "folder")
	pssm_list, query_list, align_files, alert = get_all_ids(pssm_file), get_all_ids(query_file), getFiles(align_folder, suffix), False
	
	# sort pssms and queries coz closer names probably are closer to each other in evolution perspective
	pssm_list.sort()
	query_list.sort()
	
	check_unique(pssm_list, "PSSM")
	check_unique(query_list, "query")

	
	pssm_no, query_no, align_no = len(pssm_list), len(query_list), len(align_files)	

	prompt_no = get_prompt_no(align_no)
	
	
	# original script in following line:
	#outfile, total_file = pssm_file + "_" + query_file + "_score" + swapext(suffix), 0
	
	# zhenhai modified to incorporate folder names together with file names
	outfile, total_file = get_file_name_only(pssm_file) + "_" + get_file_name_only(query_file) + "_score" + swapext(suffix), 0
	print outfile
	
	writer = csv.writer(open(outfile, "w"), delimiter = sep)

	writer.writerow(["query"] + pssm_list)
	
	print "Total %d PSSMs, %d query sequences, %d alignment files" %(len(pssm_list), len(query_list), len(align_files))
	for align_file in align_files:
		
		print "calculating score in file:", align_file
		score_dict = init_score_dict()
		
		if not align_folder.endswith("/"):
			align_file = align_folder + "/" + align_file
		else:
			align_file = align_folder + align_file
		
		reader, total = csv.reader(open(align_file, "rU"), delimiter = sep), 0
		reader.next()		# skip title line
		for row in reader:
			query, pssm, identity, query_coverage = row[0], row[2], float(row[9]), float(row[12])
			score = identity * query_coverage
			if score > score_dict[pssm]:
				score_dict[pssm] = score
			
			total += 1
		
		aline = get_aline_from_dict(score_dict, pssm_list)

		writer.writerow([query] + aline)
		
		try:
			query_list.remove(query)
		except:
			print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"
			print "============ ATTENTION           |           ATTENTION ============"
			print "============ ATTENTION          ---          ATTENTION ============"
			
			print "REPEAT ALIGNMENT of query %s in file %s" %(query, align_file)

			print "============ ATTENTION          ---          ATTENTION ============"
			print "============ ATTENTION           |           ATTENTION ============"
			print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"

			alert = True
			
		if total > pssm_no:
			alert = True
			print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"
			print "============ ATTENTION           |           ATTENTION ============"
			print "============ ATTENTION          ---          ATTENTION ============"
			
			print "Alignment number: %d EXCEEDED total PSSM number: %d in file %s" %(total, pssm_no, align_file)
			
			print "============ ATTENTION          ---          ATTENTION ============"
			print "============ ATTENTION           |           ATTENTION ============"
			print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"
			
			
			
		print "finished processing file %s; %d hits found from total %d pssms" %(align_file, total, pssm_no)
		
		total_file += 1
		if total_file % prompt_no == 0:
			print "%d%s of total task finished... " %(total_file / prompt_no * 10, "%")
	
	print "%d file finished... %d queries have no hits with any pssms" %(total_file, len(query_list))
	
	for query in query_list:
		aline = [query] + [0.0] * pssm_no
		writer.writerow(aline)
	
	if not alert:
		print "Congrats! All done....... result saved in ", outfile
	else:
		print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"
		print "============ ATTENTION           |           ATTENTION ============"
		print "============ ATTENTION          ---          ATTENTION ============"
		
		print "WARNING: Please check SCREEN OUTPUT or LOG FILE, possibly NON-UNIQUE QUERY or MULTIPLE HITS FOR SINGLE QUERY found!!!"
		
		print "============ ATTENTION          ---          ATTENTION ============"
		print "============ ATTENTION           |           ATTENTION ============"
		print "============ ATTENTION      =8-(o_o)-8=      ATTENTION ============"
		
	
	end = datetime.now()
	print "program ended at:", end
	print "program duration:", end - start	
		
	
	
