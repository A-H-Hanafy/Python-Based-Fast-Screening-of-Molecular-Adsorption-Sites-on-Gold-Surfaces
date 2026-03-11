PbO Adsorption on Au(111) Surface
=================================

Description
-----------

This code calculates the optimized lattice constant (``a``) for bulk Gold (Au) and constructs
a four-layer (111) slab with a ``3 × 3`` supercell in order to determine the surface energy.

The script also calculates the energy of an isolated horizontal PbO molecule in vacuum
as a reference state. The PbO molecule is then placed on several high-symmetry adsorption
sites on the Au(111) surface in order to simulate the adsorption process.

Finally, geometry optimization is performed for each configuration and the adsorption
energy (:math:`E_{ads}`) is calculated to determine the most stable binding site.


Adsorption Results
------------------

CHGNet v0.3.0 initialized with 412,525 parameters
CHGNet will run on cpu
CHGNet will run on cpu
The lattice constant a =  4.170016067634946
The Energy of Gold Surface Au(111) =  -113.38893
The Energy of PbO (Isolated) =  -8.688824

+-----------+-------------+--------------+-------------+
| Site      | Orientation | E_total (eV) | E_ads (eV)  |
+===========+=============+==============+=============+
| ontop     | Flat        |  -123.2050   |  -1.1272    |
+-----------+-------------+--------------+-------------+
| bridge    | Flat        |  -123.2222   |  -1.1444    |
+-----------+-------------+--------------+-------------+
| fcc       | Flat        |  -123.3423   |  -1.2645    |
+-----------+-------------+--------------+-------------+
| hcp       | Flat        |  -123.3230   |  -1.2453    |
+-----------+-------------+--------------+-------------+

Most Stable: fcc (Flat) at -1.2645 eV
