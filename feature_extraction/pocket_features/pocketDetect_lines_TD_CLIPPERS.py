#!/software/miniconda3/envs/pyrosetta2/bin/python
import sys

# Script based on tm3 file from the protein-pocket-detecting program CLIPPERS ("Protein pockets: inventory, shape and comparison" - Coleman & Sharp 2009)
# Detects main pocket by finding the deepest one and getting the biggest possible volume without extending to much to the surface.
# By Benjamin Basanta - Sept. 2015

# Give tm3 file as 1st argument
tm3_file = open(sys.argv[1],'r')

def GetParent(line):
    cols = line.split()
    parent = cols[-2]
    return parent

def GetPocketStats(line):
    cols = line.split()
    stats= cols[0]+' '+cols[1]+' '+cols[2]+' '+cols[6]+' '+cols[8]+' '+cols[9]+' '+cols[4]+' '+cols[5]
    return stats

lines = tm3_file.readlines()
pocket = '1'
#print 'Name Surface_area Vol height mouths area_of_biggest_mouth mean_TD max_TD'
lowest_mean_TD = 0
pocket_tree_lines_raw = []
pocket_tree_lines = []
# Take the first 50 pockets in the tree starting for the deepest (group 1) and get the lowest mean_TD.
# Pockets further than 50 usually extend too far out in the surface, at least for the NTF2 cases.
for n in range(0,50):
    for line in lines:
        cols = line.split()
        if cols[0] == pocket:
            pocket = GetParent(line)
            pocket_tree_lines_raw.append(line)
            pocket_tree_lines.append(GetPocketStats(line))
            lowest_mean_TD = float(GetPocketStats(line).split()[6])

# Using the lowest mean_TD, calculate the mean_TD cutoff to call that a pocket: cutoff = max_TD - (max_TD - lowest_mean_TD)/2

mean_TD_adjust = 0.75 # Changed default to 0.75 beacuse that gave better results: pocket detection without native information overlapped better with pocket detection with that information.
if len(sys.argv) == 3: # Optional, change the cutoff in commandline. This is only useful if you want to "calibrate" this script to detect similar positions as something else.
    mean_TD_adjust = float(sys.argv[2])

max_TD = float(pocket_tree_lines[0].split()[-1])
cutoff = max_TD - (max_TD - lowest_mean_TD) * mean_TD_adjust
key_pocket = ''
key_pocket_raw = ''
filename = sys.argv[1][:-19]
min_mean_TD_deltas = {}

# Get the main pocket using the cutoff:
# Walk through the pocket "tree" and get the mean_TD, store the difference between
# mean_TD and the cutoff, so you can retreave the pocket closest to the cutoff
for i,pock in enumerate(pocket_tree_lines):
    cols = pock.split()
    if int(cols[0])<120:
        mean_TD = float(cols[-2])
        min_mean_TD_deltas[mean_TD - cutoff] = i

# Sorting to get the pocket closest to the cutoff that is smaller than cutoff:
delta_subset = sorted(i for i in min_mean_TD_deltas.keys() if i > 0)

# Print the pocket with the smaller difference to the intended cutoff:
key_pocket_raw = filename + " " + pocket_tree_lines_raw[min_mean_TD_deltas[delta_subset[0]]][:-2]

print key_pocket_raw
