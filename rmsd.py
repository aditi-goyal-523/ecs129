import argparse
import math
import numpy as np

parser = argparse.ArgumentParser(description='Compare protein structures through RMSD (root mean squared deviation) using quaternion')
parser.add_argument('--tar', type=str,
	metavar='<pdb>', help='path to target (standard) protein structure')
parser.add_argument('--mod', type=str,
	metavar='<pdb>', help='path to model (predicted) protein structure')
arg = parser.parse_args()

def rmsd(v1, v2):
	sum_edistance_sqd = 0
	for i in range(len(v1)):
		edistance = math.dist(v1[i] , v2[i]) # obtain eucladian distance between the repective CAs
		sum_edistance_sqd += math.pow(edistance,2) # obtain the sum of educladian distance squared
	
	return math.sqrt(sum_edistance_sqd/len(v1))	# the RMSD

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

if len(tar_vecs) != len(mod_vecs):
	raise IndexError('Proteins must be the same size!')
else:
	prolen = len(tar_vecs)

tar_vecs = np.array(tar_vecs)
mod_vecs = np.array(mod_vecs)


preopt_rmsd = rmsd(tar_vecs, mod_vecs)
print(f'pre-optimization RMSD: {preopt_rmsd:.3f}')

# Get barycenters of target and model vectors.... Can put into a function
tar_bc = np.mean(tar_vecs, axis=0)
mod_bc = np.mean(mod_vecs, axis=0)

# Linearly translate the vectors
tar_vecs = tar_vecs - tar_bc
mod_vecs = mod_vecs - mod_bc

post_translation_rmsd = rmsd(tar_vecs, mod_vecs)
print(f'post-translation RMSD: {post_translation_rmsd:.3f}')

# get R from SVD
product = np.transpose(mod_vecs) @ tar_vecs
u, s, vh = np.linalg.svd(product)
R = vh.T * u.T

# get fancy F in terms of the matrix elements of fancy R for quaternion algorithm
F = [[0 for i in range(4)] for j in range(4)]
F[0][0] = R[0][0] + R[1][1] + R[2][2]
F[0][1] = R[1][2] - R[2][1]
F[0][2] = R[2][0] - R[0][2]
F[0][3] = R[0][1] - R[1][0]
F[1][0] = R[1][2] - R[2][1]
F[1][1] = R[0][0] - R[1][1] - R[2][2]
F[1][2] = R[0][1] + R[1][0]
F[1][3] = R[0][2] + R[2][0]
F[2][0] = R[2][0] - R[0][2]
F[2][1] = R[0][1] + R[1][0]
F[2][2] =-R[0][0] + R[1][1] - R[2][2]
F[2][3] = R[1][2] + R[2][1]
F[3][0] = R[0][1] - R[1][0]
F[3][1] = R[0][2] + R[2][0]
F[3][2] = R[1][2] + R[2][1]
F[3][3] =-R[0][0] - R[1][1] + R[2][2]

F = np.array(F)
eigen_val, eigen_vec = np.linalg.eig(F)
eigen_max = max(eigen_val)
	
e_sqrd = 0
for i in range(prolen):
	e_sqrd = np.dot(tar_vecs[i],tar_vecs[i]) + np.dot(mod_vecs[i],mod_vecs[i]) - (2 * eigen_max)

e_sqrd = e_sqrd/prolen

e = math.sqrt(e_sqrd) 

print(f'post-transformation RMSD: {e:.3f}') 


