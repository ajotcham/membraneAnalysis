# assume a log file only has a single column format here across different runs here 

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="write lammps log output to a csv/data file")
parser.add_argument('-i', action="store", dest="input_file")
parser.add_argument('-o', action="store", dest="output_format")
args = parser.parse_args()

def purge(tokens): return [t for t in tokens if len(t) >= 1] # purge empty strings from token lists

with open(args.input_file,'r') as f: lines = f.readlines()
lines = [l.strip('\n') for l in lines]

mpi_line_idxs = []
looptime_line_idxs = []
error_line_idxs = []
linecounter = 0
for line in lines:
	if line.startswith("Per MPI rank memory allocation"): mpi_line_idxs.append(linecounter)
	if line.startswith("Loop time"): looptime_line_idxs.append(linecounter)
	if line.startswith("ERROR"): error_line_idxs.append(linecounter)
	linecounter += 1

assert len(error_line_idxs) <= 1
assert len(mpi_line_idxs) - len(looptime_line_idxs) <= 1

header_idxs = []
for idx in mpi_line_idxs:
	header_idxs.append(idx+1)
headers = []
for idx in header_idxs:
	headers.append(purge(lines[idx].split(' ')))
# make sure headers for all runs have the same number of columns
headlengths = [len(header) for header in headers]
headlengths = list(set(headlengths))
assert len(headlengths) == 1
header = headers[0] # assume same thermo_style across runs

thermo_lines = []
if len(mpi_line_idxs) == len(looptime_line_idxs): # all runs ran fine
	for i in range(len(mpi_line_idxs)):
		thermo_lines.extend(lines[mpi_line_idxs[i]+2:looptime_line_idxs[i]])
elif len(mpi_line_idxs) - len(looptime_line_idxs) == 1: # errored out or didn't finish
	for i in range(len(mpi_line_idxs)-1):
		thermo_lines.extend(lines[mpi_line_idxs[i]+2:looptime_line_idxs[i]])
	thermo_lines.extend(lines[mpi_line_idxs[-1]+2:]) # last run that errored/timed out
thermo_lines = [purge(line.split(' ')) for line in thermo_lines]

# don't include last line if sim timed out
if (len(mpi_line_idxs) - len(looptime_line_idxs)) and (len(error_line_idxs) == 0):
	thermo_lines = thermo_lines[:-1]

# make sure data is numeric
#thermo_data = []
#for line in thermo_lines:
#	thermo_data.append([float(t) for t in line])
thermo_data = thermo_lines
#print(thermo_data)

out_file = args.input_file[:-3] + args.output_format
df = pd.DataFrame(thermo_data, columns=header)
if args.output_format == 'csv': df.to_csv(out_file)
elif args.output_format == 'dat': df.to_csv(out_file, sep=" ")

