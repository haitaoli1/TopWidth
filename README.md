# TopWidth
Mathematica Package to calculate the top decay width with NNLO corrections  in QCD and NLO corrections in EW. 

# Requirement 
The HPL package is required to generate the numerics of the harmonic polylogarithm, which can be download from https://krone.physik.uzh.ch/data/HPL/ .  

HPL is supposed to be initialized through "\<\<HPL`". If not plase  set the path "$HPLPath="the:\path\of\the\installation".

# Download
Download the package through 

git clone https://github.com/haitaoli1/TopWidth.git

go to the directory "TopWidth", run the example notebook "example.nb". 

# Functions 
TopWidth[order, mbCorr, WwidthCorr, EWcorr, mu] returns the top quark width in the Standard Model;  
order = 0, 1, or 2; The QCD corrections are calculated accordingly;  
mbCorr=0 or 1; if mbCorr=1, the b mass effects included at LO in QCD for order = 0, NLO in QCD for order=1;  mb effects in NNLO QCD is not available;  
WwidthCorr=0 or 1; When WwidthCorr=0  the on-shell W is produced in top quark decay; WwidthCorr=1, the off-shell effects included up to NNLO in QCD according to order;  
EWcorr=0 or 1; For EWcorr=1 and order>=1, the NLO EWcorrections are included.;  

SetParameters:: SetParameters[mtv,mbv,mwv,wwidthv,mz, GFv]  set the parameters for mt, mw, mb, wwidth, mz, and GF.

# License
TopWidth is covered by the GNU General Public License 3.

TopWidth is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# Contacts
If you have any problem, contact Haitao Li by haitao.li@sdu.edu.cn. 
