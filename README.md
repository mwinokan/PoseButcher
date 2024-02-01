
# PoseButcher

> "A good butcher always trims the fat"

PoseButcher is a tool for categorising and segmenting virtual hits with reference to experimental protein structures and (fragment) hits.

Ligand atoms are tagged with categories:

	- GOOD:

		* fragment space: within the fragment bolus
		* pocket X: in a specified catalytic/allosteric pocket X

	- BAD:
		
		* protein clash: Clashing with the protein
		* solvent space: Heading out of the protein/crystal

## Usage at a glance

	1. Create the butcher (see PoseButcher.__init__):

		from posebutcher import PoseButcher
		butcher = PoseButcher(protein, hits, pockets)

	2. Chop up a posed virtual hit (rdkit.ROMol with a conformer):

		result = butcher.chop(mol)

	3. Tag a compound based on its pocket occupancy and clashes:

		tags = butcher.tag(mol)

	4. Explore the expansion opportunities from a given atom in a virtual hit

		result = butcher.explore(mol)

	5. (Coming soon) Trim a parts of a compound that clash with a protein or leave the crystal

		mol = butcher.trim(mol)

	6. (Coming soon) Score how well a virtual hit recapitulates shape and colour of the fragment bolus

		score: float = butcher.score(mol)

## Sample output

### Protein, ligand, and catalytic site pockets with Open3D:

![butcher_3d](https://github.com/mwinokan/PoseButcher/assets/36866506/19c7895d-328b-4bdb-aca1-539f59fbea28)

### A chopped up ligand:

![butcher_2d](https://github.com/mwinokan/PoseButcher/assets/36866506/f4947703-83f9-45f5-8e70-108918582938)

### Exploring expansion vector into a desirable pocket

![butcher_vector](https://github.com/mwinokan/PoseButcher/assets/36866506/af34e286-3dfe-4f5a-871e-ff8d28292a3d)

## Installation

1. Install PoseButcher:

`pip install --upgrade posebutcher`

*N.B. open3d version 0.18 or greater is required*

## Examples

PoseButcher ships with some open access test data from the XChem group at Diamond Light Source, funded by the ASAP consortium.

To run the "example.ipynb":

```
git clone git@github.com:mwinokan/PoseButcher
cd PoseButcher
jupyter lab
```
