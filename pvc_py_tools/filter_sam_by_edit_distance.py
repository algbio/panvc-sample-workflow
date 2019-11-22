#!/usr/bin/python3

# Copyright (c) Tuukka Norri 2019
# Licenced under the MIT licence.

import argparse
import re
import sys


def find_optional_start(line):
	"""Find the starting position of the optional fields on the given line."""
	pos = 0
	for i in range(0, 11): 
		pos = line.find("\t", pos) 
		if -1 == pos:
			return -1
		pos += 1
	return pos

def check_edit_distance(fields, max_ed):
	"""Check the optional SAM fields for edit distance."""
	for field in fields:
		match = ed_re.match(field)
		if match:
			ed = int(match.group(1))
			return (ed <= max_ed)
	return False


parser = argparse.ArgumentParser(description = 'Filter SAM records by edit distance.')
parser.add_argument('--max-ed', metavar = 'N', type = int, required = True, help = "Maximum edit distance")

ed_re = re.compile(r"^NM:i:(\d+)$")

args = parser.parse_args()
for line in sys.stdin:
	if line.startswith("@"):
		# Output header lines.
		sys.stdout.write(line)
		continue
	
	# Find the end of the mandatory fields.
	pos = find_optional_start(line)
	if -1 == pos:
		continue
	
	# Check for edit distance and output.
	optional_fields_part = line[pos:-1]
	optional_fields = optional_fields_part.split("\t")
	if check_edit_distance(optional_fields, args.max_ed):
		sys.stdout.write(line)
