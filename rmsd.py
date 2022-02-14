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

def AL(vec):
	vec = [0, vec[0], vec[1], vec[2]]
	al = [[0 for i in range(4)] for j in range(4)]
	
	al[0][0] = vec[0]; al[0][1] = -vec[1]; al[0][2] = -vec[2]; al[0][3] = -vec[3]
	al[1][0] = vec[1]; al[1][1] =  vec[0]; al[1][2] = -vec[3]; al[1][3] =  vec[2]
	al[2][0] = vec[2]; al[2][1] =  vec[3]; al[2][2] =  vec[0]; al[2][3] = -vec[1]
	al[3][0] = vec[3]; al[3][1] = -vec[2]; al[3][2] =  vec[1]; al[3][3] =  vec[0]
	
	return np.array(al)
	
def AR(vec):
	vec = [0, vec[0], vec[1], vec[2]]
	ar = [[0 for i in range(4)] for j in range(4)]
	
	ar[0][0] = vec[0]; ar[0][1] = -vec[1]; ar[0][2] = -vec[2]; ar[0][3] = -vec[3]
	ar[1][0] = vec[1]; ar[1][1] =  vec[0]; ar[1][2] =  vec[3]; ar[1][3] = -vec[2]
	ar[2][0] = vec[2]; ar[2][1] = -vec[3]; ar[2][2] =  vec[0]; ar[2][3] =  vec[1]
	ar[3][0] = vec[3]; ar[3][1] =  vec[2]; ar[3][2] = -vec[1]; ar[3][3] =  vec[0]
	
	return np.array(ar)

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

# get F
F = np.array([[0.0 for i in range(4)] for j in range(4)])

for i in range(prolen):
	F -= np.matmul(AL(tar_vecs[i]),AR(mod_vecs[i])) / prolen

# get RMSD
eigen_val, eigen_vec = np.linalg.eig(F) # Get eigenvalues and eigenvectors of F
eigen_max = max(eigen_val) # Get maximum eigenvalue of F
	
e_sqrd = 0
for i in range(prolen):
	e_sqrd = np.dot(tar_vecs[i],tar_vecs[i]) + np.dot(mod_vecs[i],mod_vecs[i]) - (2 * eigen_max)

e_sqrd = e_sqrd/prolen
#print(e_sqrd)

e = math.sqrt(e_sqrd) 

print(f'post-transformation RMSD: {e:.3f}')

