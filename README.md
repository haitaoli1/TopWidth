# TopWidth
Mathematica Package to calculate the top decay width. The QCD corrections are included up to N3LO leading color contribution. The NLO EW corrections are implemented. The b mass effects at NLO QCD corrections are considered. The finite W width effects are included for all the QCD corrections.   

# Requirement 
The HPL package is required to generate the numerics of the harmonic polylogarithm, which can be downloaded from https://krone.physik.uzh.ch/data/HPL/ .  

HPL is supposed to be initialized through "\<\<HPL`". If not please  set the path "$HPLPath="the:\path\of\the\installation".

# Download
Download the package through 

git clone https://github.com/haitaoli1/TopWidth.git

go to the directory "TopWidth", run the example notebook "example.nb". 

# Functions 
TopWidth[order, mbCorr, WwidthCorr, EWcorr, mu] returns the top quark width in the Standard Model;  
order = 0, 1, 2, or 3; The QCD corrections are calculated accordingly;  
mbCorr=0 or 1; if mbCorr=1, the b mass effects included at LO in QCD for order = 0, NLO in QCD for order=1;  mb effects in NNLO QCD is not available;  
WwidthCorr=0 or 1; When WwidthCorr=0  the on-shell W is produced in top quark decay; WwidthCorr=1, the off-shell effects included up to NNLO in QCD according to order;  
EWcorr=0 or 1; For EWcorr=1 and order>=1, the NLO EWcorrections are included.;  

SetParameters:: SetParameters[mtv,mbv,mwv,wwidthv,mz, GFv]  set the parameters for mt, mw, mb, wwidth, mz, and GF.

# Reference 

If TopWidth is used in your research please cite our paper   
[1] "Analytic result for the top-quark width at next-to-next-to-leading order in QCD", by Long-Bin Chen, Hai Tao Li, Jian Wang, Yefan Wang, arXiv:2212.06341.
[2] "Analytic three-loop QCD corrections to top-quark and semileptonic b->u decays", by Long-Bin Chen, Hai Tao Li, Zhao Li, Jian Wang, Yefan Wang, Quan-Feng Wu, arXiv:2309.00762.

# License
TopWidth is covered by the GNU General Public License 3.

TopWidth is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 or any later version of the License.

Copyright (C) 2023-2024 Long-Bin Chen, Hai Tao Li, Zhao Li, Jian Wang, Yefan Wang, Quan-Feng Wu

# Contacts
If you have any question about this program, contact Haitao Li by haitao.li@sdu.edu.cn. 
