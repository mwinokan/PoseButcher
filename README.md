# PoseButcher

Pose butcher segments a ligand into categories:

- GOOD:

	* fragment space: in catalytic fragment space
	* pocket X: in desirable pocket X

- BAD:
	
	* solvent space: Heading out of the protein/crystal
	* protein clash: Clashing with the protein

## Usage

1. Create the butcher:

`butcher = Butcher(protein, hits, pockets)`

2. Chop up a ligand pose:

`result = butcher.chop(pose)`
