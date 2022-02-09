import argparse

parser = argparse.ArgumentParser(description='Compare protein structures')
parser.add_argument('--target', type=str,
	metavar='<pdb>', help='path to target protein')
parser.add_argument('--model', type=str,
	metavar='<pdb>', help='path to model protein')
arg = parser.parse_args()

# Vector initialization

target_vector = []
with open(arg.target, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if 'CA' in line: # processing out the non-CA atoms
			x = line.split()[6] # acquiring x position
			y = line.split()[7] # acquiring y position
			z = line.split()[8] # acquiring z position
			target_vector.append([x,y,z])

model_vector = []
with open(arg.model, 'r') as fh:
	lines = fh.readlines()
	for line in lines:
		if 'CA' in line:
			x = line.split()[6]
			y = line.split()[7]
			z = line.split()[8]
			model_vector.append([x,y,z])

