# PoseButcher

Pose butcher segments a ligand into categories:

- GOOD:

	* fragment space: in catalytic fragment space
	* pocket X: in desirable pocket X

- BAD:
	
	* solvent space: Heading out of the protein/crystal
	* protein clash: Clashing with the protein

## Usage

1. Install:

`$ pip install git+https://github.com/mwinokan/PoseButcher.git#egg=posebutcher`

2. Import:

`from posebutcher import Butcher`

2. Create the butcher:

`butcher = Butcher(protein, hits, pockets)`

3. Chop up a ligand pose:

`result = butcher.chop(pose)`

## Examples

```
from posebutcher import Butcher
protein = Path / mp.System
pose = rdkit.Mol / mp.AtomGroup / ase.Atoms
pockets = {'P1':'r101,r122,...'} / or mp.Primitives?
butcher = Butcher(protein, pockets)
result : dict = butcher(pose, draw = ['2d,3d,debug'])
```
