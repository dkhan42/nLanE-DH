# nLanE-DH
**PLEASE USE PYSCF 2.8.0 OR HIGHER**

Non-linear and non-empirical double hybrid (nLanE-DH) density functional derived through the adiabatic connection.

Usage example He2+ molecule:
```
from pyscf import gto
from nLanE import return_energy

geom = "He 0 0 0; He 0 0 1.2" 
mol = gto.M(atom=geom, basis='cc-pvqz', unit="Angstrom", spin=1, charge=1, verbose=0)
energy = return_energy(mol)
```

Or supply xyz file:
```
from pyscf import gto
from nLanE import return_energy

geom = "path_to_xyz.xyz"
mol = gto.M(atom=geom, basis='cc-pvqz', unit="Angstrom")
energy = return_energy(mol)
```

To Do : Density fitting, RI-MP2
