import argparse
import math

parser = argparse.ArgumentParser(description='Compare protein structures')
parser.add_argument('--tar', type=str,
	metavar='<pdb>', help='path to target protein')
parser.add_argument('--mod', type=str,
	metavar='<pdb>', help='path to model protein')
arg = parser.parse_args()

def rmsd(v1, v2):
	sum_edistance_sqd = 0
	for i in range(len(v1)):
		edistance = math.dist(v1[i] , v2[i]) # obtain eucladian distance between the repective CAs
		sum_edistance_sqd += math.pow(edistance,2) # obtain the sum of educladian distance squared
	
	return math.sqrt(sum_edistance_sqd/len(v1))	# the RMSD

# Vectors initialization
tar_vec = []
with open(arg.tar, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if 'CA' in line: # process out the non-CAs
			x = float(line.split()[6]) # acquire x position
			y = float(line.split()[7]) # acquire y position
			z = float(line.split()[8]) # acquire z position
			tar_vec.append([x,y,z])

mod_vec = []
with open(arg.mod, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if 'CA' in line:
			x = float(line.split()[6])
			y = float(line.split()[7])
			z = float(line.split()[8])
			mod_vec.append([x,y,z])

if len(tar_vec) != len(mod_vec):
	raise IndexError('Proteins must be the same size!')
else:
	prolen = len(tar_vec)
	
preopt_rmsd = rmsd(tar_vec, mod_vec)
print(f'pre-optimization score: {preopt_rmsd:.3f}')

# Get barycenters of target and model vectors.... Can put into a function
tar_bc = [0, 0, 0]
for position in tar_vec:
	tar_bc[0] += position[0] # sum up the weight of CAs in x direction
	tar_bc[1] += position[1] # sum up the weight of CAs in y direction
	tar_bc[2] += position[2] # sum up the weight of CAs in z direction
	
mod_bc = [0, 0, 0]
for position in mod_vec:
	mod_bc[0] += position[0]
	mod_bc[1] += position[1]
	mod_bc[2] += position[2]
	
for i in range(3):
	tar_bc[i] = tar_bc[i]/prolen # Obtain the barycenter, or center of mass, of the vector
	mod_bc[i] = mod_bc[i]/prolen

# Linearly translate the vectors
for i in range(prolen):
	for j in range(3):
		tar_vec[i][j] = tar_vec[i][j] - tar_bc[j]
		mod_vec[i][j] = mod_vec[i][j] - mod_bc[j]

post_translation_rmsd = rmsd(tar_vec, mod_vec)
print(f'post-translation score: {post_translation_rmsd:.3f}')
