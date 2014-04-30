#!/usr/bin/python

######################################
# E. Roberson                        #
# Created 2014-01-29                 #
# Deindexes single-end illumina data #
######################################

###########
# Imports #
###########
import argparse
from collections import namedtuple
from datetime import datetime
import os
import re
import sys

####################
# Version and name #
####################
SCRIPT_PATH=sys.argv[0]
SCRIPT_NAME=SCRIPT_PATH.split('/')[-1].split('\\')[-1]
VERSION = '1.0.0'

#######################
# Constants & globals #
#######################
COMP_BASE = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
MUT_ARRAY = {'A':('T','C','G','N'), 'T':('A','C','G','N'), 'C':('A','G','T','N'), 'G':('A','C','T','N')}
OUTPUTEXTENSION = '.fqs'
JNKEXT = '.fqs'

############
# Counters #
############
readCount = 0
ambigCount = 0
matchCount = 0
junkReadCount = 0

########
# Fxns #
########
def nonZero( count ):
	return count if count > 0 else 1

def errMsg( s ):
	print s
	sys.exit(1)

def setQualScale( filename ):
	try:
		INFH = open(filename, 'r')
	except:
		errMsg( "Error when trying to open file [%s] to determine quality scale" % (filename) )
	
	# phred		 0 to 93, ASCII 33 to 126; offset 33; phred scale
	# sol1.0	-5 to 62, ASCII 59 to 126; offset 64; log odds
	# sol1.3+	 0 to 62, ASCII 64 to 126; offset 64; phred scale
	# Qsanger = -10 * log10(p)
	# Qsol1.0 = -10 * log10( p / 1-p )
	# 33 - 58 == phred
	# 59 - 63 == sol1.0
	# sol1.3 if no evidence for either
	
	maxReads = 10000
	readCount = 0
	minQual = 256
	localIlmn1p0Min = 59
	localIlmn1p3Min = 64
	localPhredMin = 33
	
	buffer = ''
	
	while readCount < maxReads and INFH:
		buffer = INFH.readline() # @ identifier
		buffer = INFH.readline() # seq
		buffer = INFH.readline() # + identifier
		buffer = INFH.readline().rstrip("\r\n") # quality
		
		readCount += 1
		
		for qScore in buffer:
			if ord(qScore) < minQual:
				minQual = ord(qScore)
				
			if minQual < localIlmn1p0Min:
				readCounts = maxReads
				break
	
	INFH.close()
	
	if minQual < localIlmn1p0Min:
		return "phred"
	elif minQual < localIlmn1p3Min:
		return "illumina1.0"
	else:
		return "illumina1.3"
		
def seqList( seq, mutationDictionary, mismatch=True ):
	#########################
	# Clean up the sequence #
	#########################
	seq = seq.upper().rstrip("\r\n")
	
	outputList = [seq]
	
	if mismatch == True:
		for base in mutationDictionary[seq[-1]]:
			outputList.append( seq[:-1] + base )
			
	return outputList

def boolString( boolVal ):
	if boolVal:
		return "True"
	return "False"

def getSequenceAmount( totalSeqBp ):
	if totalSeqBp > 1E15:
		ext = "Pbp"
		adjBp = totalSeqBp / 1E15
	elif totalSeqBp > 1E12:
		ext = "Tbp"
		adjBp = totalSeqBp / 1E12
	elif totalSeqBp > 1E9:
		ext= "Gbp"
		adjBp = totalSeqBp / 1E9
	elif totalSeqBp > 1E6:
		ext = "Mbp"
		adjBp = totalSeqBp / 1E6
	elif totalSeqBp > 1E3:
		ext = "Kbp"
		adjBp = totalSeqBp / 1E3
	else:
		ext = "bp"
		adjBp = totalSeqBp
	return SequenceInterpretation( totalSeqBp, adjBp, ext )
	
#####################
# class Definitions #
#####################
class Index:	
	def __init__(self, index, outBase, outExt, mutDict):
		self.index = index.upper().rstrip("\r\n")
		self.matches = 0
		self.filename = outBase + "_" + self.index + outExt
		
		try:
			self.filehandle = open( self.filename, 'w' )
		except:
			errMsg( "Trouble opening filehandle [%s] in Index object creation" % (self.filename) )
		
		self.regEx = re.compile( '^' + '|^'.join( seqList( self.index, mutDict ) ) )
	
	def output(self, identifier, sequence, quality, index):
		self.matches += 1
		self.filehandle.write( "%s%s\n%s\n+\n%s\n" % (identifier, index, sequence, quality) )
			
	def __del__(self):
		self.filehandle.close()
		print( "%s %s" % (self.index, self.matches) )

############
# argparse #
############
parser = argparse.ArgumentParser(prog=SCRIPT_NAME, epilog="%s v%s" % (SCRIPT_NAME, VERSION))

parser.add_argument('seqFastq', help="Single-end reads from amplicon library (FASTQ)")
parser.add_argument('indexFastq', help="Index reads (FASTQ)")
parser.add_argument('outputBase', help="Basename for output files")
parser.add_argument('indexString', help="Comma-separated list of indexes")
parser.add_argument('--noMismatch', help="Disallows single terminal mismatch", action='store_false', required=False, default=True, dest='mismatch')
parser.add_argument('--Ncutoff', help="Percent N base cutoff before sequence junked", required=False, default=0.85, dest='n_perc_cutoff', type=float)

args = parser.parse_args()

INDEXES = [val for val in args.indexString.upper().rstrip("\r\n").split(',')]

################
# Input scale? #
################
qScale = setQualScale( args.seqFastq )

if qScale != 'phred':
	errMsg( "WHOA. %s quality scale detected. Convert to phred first. Bailing." % (qScale) )

#########################
# print settings to log #
#########################
print "%s v%s" % (SCRIPT_PATH, VERSION)
print "Options set"
print "==========="
print "Input-scale: %s" % (qScale)
print "Indexes: %s" % (', '.join(INDEXES))
print "Mismatches allowed: %s" % (boolString(args.mismatch))
print "Sequence file: %s" % (args.seqFastq)
print "Index file: %s" % (args.indexFastq)
print "Output base: %s" % (args.outputBase)
print "Output extension: %s" % (OUTPUTEXTENSION)
sys.stdout.flush() # forces writing immediately. otherwise sits in buffer until the program finishes

##########################
# Build index dictionary #
##########################
indexDict = {}
for indexSeq in INDEXES:
	#############################
	# Die if index listed twice #
	#############################
	if indexSeq in indexDict:
		errMsg( "%s listed twice - aborting" % (indexSeq) )
	
	######################
	# Populate the index #
	######################
	indexDict[indexSeq] = Index(indexSeq, args.outputBase, OUTPUTEXTENSION, MUT_ARRAY)

######################
# Set up ambig files #
######################
try:
	ambigSeqFH = open( args.outputBase + "_seq_" + "noindex" + JNKEXT, 'w' )
except:
	errMsg( "Trouble opening file for ambiguous sequences" )
	
try:
	ambigNdxFH = open( args.outputBase + "_ndx_" + "noindex" + JNKEXT, 'w' )
except:
	errMsg( "Trouble opening file for ambiguous indexes" )
	
######################
# Set up junk  files #
######################
try:
	junkSeqFH = open( args.outputBase + "_seq_" + "junk" + JNKEXT, 'w' )
except:
	errMsg( "Trouble opening file for junk sequences" )
	
try:
	junkNdxFH = open( args.outputBase + "_ndx_" + "junk" + JNKEXT, 'w' )
except:
	errMsg( "Trouble opening file for junk indexes" )

###################
# start the clock #
###################
analysisStartTime = datetime.now()

#####################
# process the files #
#####################
try:
	seqFH = open( args.seqFastq, 'r' )
except:
	errMsg( "Couldn't open sequence input [%s]" % ( args.seqFastq ) )

try:
	ndxFH = open( args.indexFastq, 'r' )
except:
	errMsg( "Couldn't open index input [%s]" % ( args.indexFastq ) )

while seqFH and ndxFH:
	seqId = seqFH.readline().rstrip('\r\n')
	seqBases = seqFH.readline().rstrip('\r\n')
	seqRepId = seqFH.readline().rstrip('\r\n')
	seqQual = seqFH.readline().rstrip('\r\n')
	
	ndxId = ndxFH.readline().rstrip('\r\n')
	ndxBases = ndxFH.readline().rstrip('\r\n')
	ndxRepId = ndxFH.readline().rstrip('\r\n')
	ndxQual = ndxFH.readline().rstrip('\r\n')
	
	if len(seqId) == 0 and len(ndxId) == 0:
		break
	
	readCount += 1
	
	if seqId.split(' ')[0] != ndxId.split(' ')[0]:
		errMsg( "Reads out of sync!!!\nR1: %s\nR2: %s\n" % (seqId, ndxId) )
	
	# check for high N seq, write to junk
	nperc = float( seqBases.count('N') ) / float( len(seqBases) )
	
	if nperc >= args.n_perc_cutoff:
		junkReadCount += 1
		junkSeqFH.write( "%s\n%s\n+\n%s\n" % (seqId, seqBases, seqQual) )
		junkNdxFH.write( "%s\n%s\n+\n%s\n" % (ndxId, ndxBases, ndxQual) )
	else:
		matched = False
		for indexSeq in indexDict:
			if matched:
				break
				
			if indexDict[indexSeq].regEx.search( ndxBases ):
				matched = True
				matchCount += 1
				indexDict[indexSeq].output( seqId, seqBases, seqQual, indexSeq )
			
		if not matched:
			ambigCount += 1
			ambigSeqFH.write( "%s\n%s\n+\n%s\n" % (seqId, seqBases, seqQual) )
			ambigNdxFH.write( "%s\n%s\n+\n%s\n" % (ndxId, ndxBases, ndxQual) )
seqFH.close()
ndxFH.close()

#############################
# Print summary information #
#############################
print """
Deindexing completed in %s seconds

###########
# Summary #
###########
""" % (datetime.now() - analysisStartTime)

print "%s total reads processed" % (readCount)
print "%s reads (%.1f%%) ignored due to high 'N' content" % (junkReadCount, float(junkReadCount)/float( nonZero(readCount) )*100.0)
print "%s/%s (%.1f%%) reads matched specified indexes" % (matchCount, readCount, float(matchCount)/float( nonZero(readCount) )*100.0)

# indexes and counts are printed on exit as part of __del__

ambigSeqFH.close()
ambigNdxFH.close()
junkSeqFH.close()
junkNdxFH.close()
