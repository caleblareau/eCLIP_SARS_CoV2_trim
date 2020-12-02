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
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes by parsing out the UMI and trimming"

opts.add_option("-a", "--fastq1", help="<Read1> as fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> as fastq.gz")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing")

opts.add_option("-n", "--nreads", default = 10000000, help="Number of reads in each split output file")
opts.add_option("-o", "--output", default = "clip_trim", help="Output sample convention")


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




if __name__ == "__main__":	

	options, arguments = opts.parse_args()

	print(options)

	# return usage information if no argvs given
	if len(sys.argv)==1:
		os.system(sys.argv[0]+" --help")
		sys.exit()

	##### INPUTS #####
	a = options.fastq1
	b = options.fastq2
	outname = options.output

	cpu = int(options.ncores)
	n = int(options.nreads)
	print(outname)
	# File handling hell
	with gzip.open(a, "rt") as f1:
		with gzip.open(b, "rt") as f2:
			with gzip.open(outname + '_1.fastq.gz', "wt") as out_f1:
				with gzip.open(outname + '_2.fastq.gz', "wt") as out_f2:
				
					# Establish iterators
					it1 = batch_iterator(FastqGeneralIterator(f1), n)
					it2 = batch_iterator(FastqGeneralIterator(f2), n)
		
					# iterate over batches of length n
					for i, batch1 in enumerate(it1):
						batch2 = it2.__next__()			
						# parallel process the barcode processing and accounting of failures.
						pool = Pool(processes=cpu)
						pm = pool.map(debarcode_v1, zip(batch1, batch2))
						pool.close()
			
						# Aggregate output
						fq_data = list(map(''.join, zip(*[item for item in pm])))
						out_f1.writelines(fq_data[0])
						out_f2.writelines(fq_data[1])
