#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Name: plasmodiumAfrica
# Description: Reads a GenBank files with multiple entries corresponding to 'circumsporozoite protein plasmodium falciparum'. Tries to retrieve geographical information from as many entries as possible. Returns some stats.
# Author: Denis Jacob Machado (dmachado[at]uncc[dot]edu)
# Date: April 18, 2022

##
# LIBRARIES
##

import argparse, re, sys
from Bio import SeqIO

##
# FUNCTIONS
##

def READ_GEOGRAPHY(filePath):
	selected_geography = []
	with open(filePath, "r") as handle: # This is an example while loop to read a FASTA file using BioPython
		for line in handle.readlines():
			line = line.strip()
			if line:
				line = re.sub(r"[\s\n\r\t ]+", " ", line.upper())
				if not line in selected_geography:
					selected_geography.append(line)
	return selected_geography

def READ_GENBANK(filePath):
	entries_geography = {}
	with open(filePath, "r") as handle: # This is an example while loop to read a FASTA file using BioPython
		for record in SeqIO.parse(handle, "genbank"):
			for feature in record.features:
				try:
					country = "/".join(feature.qualifiers["country"])
					country = re.sub(r"[\s\n\r\t ]+", " ", country.upper())
					country = re.sub(",", "/", country)
				except:
					pass
				else:
					entries_geography[record.id] = country
					break
			if not record.id in entries_geography:
				entries_geography[record.id] = "?"
	return entries_geography

def FIND_MATCHES_IN_SELECTION(selected_geography, entries_geography):
	sys.stdout.write("Record ID, Country, Match\n")
	missing = 0
	total = len(entries_geography)
	matches = 0
	for recid in entries_geography:
		geo = entries_geography[recid]
		if geo == "?":
			missing += 1
			sys.stdout.write("{}, {}, {}\n".format(recid, geo, "?"))
		else:
			match = False
			for query in selected_geography:
				if query in geo:
					match = True
					matches += 1
					sys.stdout.write("{}, {}, {}\n".format(recid, geo, 1))
					break
			if match == False:
				sys.stdout.write("{}, {}, {}\n".format(recid, geo, 0))
	stats = """COUNTRY SEARCH:
- Entries = {}
- Missing localities = {} ({}%)
- Entries in selected localities = {} ({}%)
- Entries in other localities = {} ({}%)
---
""".format(total, missing, (missing/total)*100.00, matches, (matches/total)*100, total - (missing + matches), ((total - (missing + matches))/total)*100.00)
	sys.stderr.write(stats)
	return

def main(args): # This is the beginning of the definition of the main function
	selected_geography = READ_GEOGRAPHY(args.geography)
	entries_geography = READ_GENBANK(args.sequences)
	FIND_MATCHES_IN_SELECTION(selected_geography, entries_geography)
	return

##
# INITIALIZATION
##

if __name__ == '__main__':
	# The arguments are defined in the next few lines
	parser=argparse.ArgumentParser()
	parser.add_argument("-s", "--sequences", help = "Sequence data in GenBank format", type = str, default = True)
	parser.add_argument("-g", "--geography", help = "Text file with geographical search terms. Case insensitive. One line per entry.", type = str, default = True)
	args = parser.parse_args()
	main(args) # This will call the main fuction

exit() # Quit this script
