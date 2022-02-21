from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import array, dot, set_printoptions
import argparse

parser = argparse.ArgumentParser(description='Compare protein structures through RMSD (root mean squared deviation) using quaternion')
parser.add_argument('--tar', type=str,
	metavar='<pdb>', help='path to target (standard) protein structure')
parser.add_argument('--mod', type=str,
	metavar='<pdb>', help='path to model (predicted) protein structure')
arg = parser.parse_args()

# Matrices initialization
tar_vecs = []
with open(arg.tar, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if ' CA ' in line: # process out the non-CAs
			x = float(line.split()[6]) # acquire x position
			y = float(line.split()[7]) # acquire y position
			z = float(line.split()[8]) # acquire z position
			tar_vecs.append([x,y,z])

mod_vecs = []
with open(arg.mod, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if ' CA ' in line:
			x = float(line.split()[6])
			y = float(line.split()[7])
			z = float(line.split()[8])
			mod_vecs.append([x,y,z])

tar_vecs = array(tar_vecs)
mod_vecs = array(mod_vecs)

sup = SVDSuperimposer()
sup.set(tar_vecs, mod_vecs)
sup.run()
rms = sup.get_rms()
rot, tran = sup.get_rotran()
y_on_x1 = dot(mod_vecs, rot) + tran
y_on_x2 = sup.get_transformed

print(rms)
