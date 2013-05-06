#!/usr/bin/env python
# encoding: utf-8
"""
Created by Zhenhai Zhang on 2009-07-10.
Copyright (c) 2009 Bioinformatics and Genomics@PSU. All rights reserved.
Usage: score_matrix_2_mega_eu.py matrix_file
"""
#import modules: do not swap the sequential import; sep in os.path would overide the one in mytools
import sys, csv, os
from datetime import datetime
from numpy import array, sqrt

sep = "\t"
linesep = os.linesep
def EuclideanDistance(p, q):
	""" compute euclidean distance of two numpy array """
	return sqrt(sum((p - q) * (p - q)))

def retrieve_query_score_list(f):
	
	print "retrieving scores from file:", f
	reader, q_list, q_score, total, pssm_no = csv.reader(open(f, "rU"), delimiter = sep), [], dict(), 0, 0
	reader.next()	# skip title line
	
	for row in reader:
		query, score_list = row[0], array(map(float, row[1 :]))
		q_list.append(query)
		q_score[query] = score_list
		
		total += 1
	pssm_no = len(score_list)
	print " %d records found..." %(total)
	
	return q_list, q_score, pssm_no
	
if __name__ == '__main__':
	
	#get arguments
	if len(sys.argv) < 2:
		print __doc__
		sys.exit(0)
	start = datetime.now()
	print "program started at ", start
	matrix_file = sys.argv[1].strip()
	
	query_list, query_score_dict, pssm_no = retrieve_query_score_list(matrix_file)
	query_no = len(query_list)
	
	
	head, ext = os.path.splitext(matrix_file)
	outfile = head + "_eu.meg"
	
	handle = open(outfile, "w")
	
	# writing headers
	handle.write("#mega" + linesep)
	handle.write("!Title " + outfile + ";" + linesep)
	handle.write("!Format DataType=Distance DataFormat=LowerLeft NTaxa=" + str(query_no) + ";" + linesep)
	handle.write("!Description" + linesep)
	handle.write("  No. of Taxa : " + str(query_no) + linesep)
	handle.write("  Created on " + datetime.now().strftime("%A %d %B %Y %I:%M%p") + linesep)
	handle.write("  " + str(query_no) + " x " + str(query_no) + " pairwise distance matrix (from" + str(query_no) + " x " + str(pssm_no) + " data matrix)" + linesep)
	handle.write("  Score: cpro (Product score assigned to a profile by multiplying %identity and %coverage of the alignment produced by that profile.)" + linesep)
	handle.write(";" + linesep + linesep)
	
	# writing query names
	for ind, query in enumerate(query_list):
		aline = "[" + str(ind + 1) + "] #" + query + linesep
		handle.write(aline)
	handle.write(linesep)		# blank line
	
	handle.write("[" + " ".join(map(str, range(1, query_no + 1))) + "]" + linesep)
	for i in range(query_no):
		aline = ["[" + str(i + 1) + "]"]
		
		for j in range(0, i):
			aline.append(str(EuclideanDistance(query_score_dict[query_list[i]], query_score_dict[query_list[j]])))
		
		handle.write(" ".join(aline) + linesep)
		
		if i % 100 == 0:
			print "%d record finished..." %i
		
	handle.close()
	end = datetime.now()
	print "program ended at ", start
	print "program finished with duration of ", end - start	
	
	
	
	
		
	


