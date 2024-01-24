# PoseButcher

*A good butcher always trims the fat*

Pose butcher segments a ligand into categories:

- GOOD:

	* fragment space: in catalytic fragment space
	* pocket X: in desirable pocket X

- BAD:
	
	* solvent space: Heading out of the protein/crystal
	* protein clash: Clashing with the protein

## Usage

0. Dependencies:

	* MolParse `pip install molparse`
	* Open3d `pip install open3d`

1. Install PoseButcher:

`pip install git+https://github.com/mwinokan/PoseButcher.git#egg=posebutcher`

2. Import:

`from posebutcher import Butcher`

2. Create the butcher:

`butcher = Butcher(protein, hits, pockets)`

3. Chop up a ligand pose:

`result = butcher.chop(pose)`

## Example

Run this from the repository root directory, use a jupyter notebook 

```
from posebutcher import Butcher

protein = 'test_data/hits/A71EV2A-x0310_0A_bound.pdb'
hits = 'test_data/filtered.sdf'

pockets = {
    "P1":  dict(type='sphere', atoms=['GLY 127 O', 'PRO 107 CG', 'VAL 124 CG1'], radius='mean'),
    "P1'": dict(type='sphere', atoms=['VAL 84 CG1', 'TYR 90 CD2', 'SER 87 CB'], radius='mean'),
}

butcher = Butcher(protein, hits, pockets)

df = PandasTools.LoadSDF('test_data/BBS_AMU_products_fragalysis.sdf')

mol = df['ROMol'].values[45]

result = butcher.chop(mol)

```
