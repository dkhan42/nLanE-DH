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


For fast calculation use density-fitted SCF and MP2 calculations along with a simpler functional for the SCF:
```
from pyscf import gto
from nLanE_DF import return_energy #This script only does density-fitted calculations

geom = "path_to_xyz.xyz"
mol = gto.M(atom=geom, basis='cc-pvqz', unit="Angstrom")
energy = return_energy(mol, xcscf = 'R2SCAN') #Can use PBE for greater efficiency
```



To re-produce results from the paper use these settings with PySCF 2.8.0:
```
from pyscf import gto
from nLanE import return_energy

geom = "path_to_xyz.xyz"
mol = gto.M(atom=geom, basis='def2-qzvpp', unit="Angstrom")
energy = return_energy(mol, xcscf = 'SCAN', W1x = 'SCAN', W1c = 'SCAN') #Everything with SCAN and no density-fitting
```
