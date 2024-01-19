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

protein = 'test_data/template.pdb'
hits = 'test_data/filtered.sdf'
pockets = {}

butcher = Butcher(protein, hits, pockets)

mol_df = PandasTools.LoadSDF('test_data/BBS_AMU_products_fragalysis.sdf')

pose = mol_df['ROMol'].values[2]

result = butcher.chop(pose)
```
