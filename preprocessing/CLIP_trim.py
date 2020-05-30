#!/usr/bin/env python

# Caleb Lareau, Broad Institute

##### IMPORT MODULES #####
import os
import re
import regex
import sys
import gzip

from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"

opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing")

opts.add_option("-n", "--nreads", default = 500000000, help="Number of reads in each split output file")
opts.add_option("-o", "--output", help="Output sample convention")

options, arguments = opts.parse_args()

print(options)

def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


def chunk_writer_gzip(filename, what):
	'''
	Basic function to write a chunk of a fastq file
	to a gzipped file
	'''
	with gzip.open(filename, 'wt') as out_write:
				out_write.writelines(what)
	return(filename)		

def formatRead(title, sequence, quality):
	"""
	Takes three components of fastq file and stiches them together in a string
	"""
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))


# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### INPUTS #####
a = options.fastq1
b = options.fastq2
outname = options.output
o = options.output

cpu = int(options.ncores)
n = int(options.nreads)

# Parse input files
extension = a.split('.')[-1]
if extension == "fastq" or extension == "fq":
	sys.exist("Quitting... GZIP your .fastq files!")
elif extension == "gz":
	print("Found supplied .fastq.gz files")
else:
	sys.exit("ERROR! The input files (-a , -b) a *.fastq.gz")


#------------------------------
	
def debarcode_v1(duo):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = duo[0]; listRead2 = duo[1]
	
	# parameters to return
	fq1 = ""
	fq2 = ""
	
	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	
	UMI = sequence2[0:10]
	sequence2 = sequence2[10:]
	quality2 = quality2[10:]
	
	# Return the barcode with underscores + the biological sequence learned	
	fq1 = formatRead(UMI + "_" + title1, sequence1, quality1)
	fq2 = formatRead(UMI + "_" + title2, sequence2, quality2)
	return(fq1, fq2)



with gzip.open(a, "rt") as f1:
	with gzip.open(b, "rt") as f2:
		
		# Establish iterators
		it1 = batch_iterator(FastqGeneralIterator(f1), n)
		it2 = batch_iterator(FastqGeneralIterator(f2), n)
		
		# iterate over batches of length n
		for i, batch1 in enumerate(it1):
			batch2 = it2.__next__()
			output = o 
			
			# parallel process the barcode processing and accounting of failures.
			pool = Pool(processes=cpu)
			pm = pool.map(debarcode_v1, zip(batch1, batch2))
			pool.close()
			
			# Aggregate output
			fastq1 = [item[0] for item in pm]
			fastq2 = [item[1] for item in pm]

			# Export one chunk in parallel
			filename1 = output +'_1.fastq.gz'
			filename2 = output +'_2.fastq.gz'
			
			pool = Pool(processes=2)
			toke = pool.starmap(chunk_writer_gzip, [(filename1, fastq1), (filename2, fastq2)])
			pool.close()
