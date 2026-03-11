The code calculates the optimized lattice constant ($a$) for bulk Gold (Au) and constructs a four-layer (111) slab with a $3 \times 3$ supercell to determine its surface energy. 
It also calculates the energy of an isolated, horizontal PbO molecule in a vacuum as a reference. 
The PbO molecule is then placed at various high-symmetry sites on the Gold surface to simulate the adsorption process. 
Finally, the script performs a geometry optimization for each configuration, calculates the total system energy, and determines the adsorption energy ($E_{ads}$) to identify the most stable 
binding site.

CHGNet v0.3.0 initialized with 412,525 parameters
CHGNet will run on cpu
CHGNet will run on cpu
The lattice constant a =  4.170016067634946
The Energy of Gold Surface Au(111) =  -113.38893
The Energy of PbO (Isolated) =  -8.688824

Site       | Orientation  | E_total (eV) | E_ads (eV)
----------------------------------------
ontop      | Flat         |  -123.2050 |    -1.1272
bridge     | Flat         |  -123.2222 |    -1.1444
fcc        | Flat         |  -123.3423 |    -1.2645
hcp        | Flat         |  -123.3230 |    -1.2453
----------------------------------------
Most Stable: fcc (Flat) at -1.2645 eV
