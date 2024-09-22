# OPTMAN optical model code
This repository contains the optical model code OPTMAN. The code performs coupled channels optical 
model calculations for incident nucleons using rigid rotor, spherical soft rotor and deformed soft 
rotor formalisms and allows fit of the potential. Soft rotor model is implemented within the [SHEMMAN]
code.


[SHEMMAN]: https://github.com/IAEA-NDS/SHEMMAN

The code implements Lane-consistent dispersive optical potential with relativistic corrections and 
OpenMP parallelization by incident energies.
The obtained coupled-channels optical model potentials allow
for the description of nucleon induced reactions up to 200 MeV, including (p,n) reactions with excitation 
of isobaric analog states.

## Manuals

Programs OPTMAN and SHEMMAN version 8 (2004) [link1]
Supplement to OPTMAN code, manual version 10 (2008) [link2]
Program OPTMAN Version 14 (2013), User's Guide. [link3]

[link1]: https://inis.iaea.org/collection/NCLCollectionStore/_Public/36/116/36116793.pdf?r=1
[link2]: https://jopss.jaea.go.jp/pdfdata/JAEA-Data-Code-2008-025.pdf
[link3]: https://inis.iaea.org/collection/NCLCollectionStore/_Public/44/117/44117922.pdf?r=1

There are no updates for the manual after 2014, however main model developments can be found in [paper1], 
[paper2] and [paper3]

[paper1]: http://link.aps.org/doi/10.1103/PhysRevC.94.064605
[paper2]: https://link.aps.org/doi/10.1103/PhysRevC.102.059901
[paper3]: https://www.epj-conferences.org/10.1051/epjconf/202023903003

Updated version of the input file format is in `docs` folder.

## Known issues

The code may produce incorrect values if compiled without OpenMP and several incident energies are provided 
in an input file.

## Feedback

This code is currently maintained by Dmitry Martyanov and Roberto Capote 

