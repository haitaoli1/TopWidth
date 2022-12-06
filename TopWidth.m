(* ::Package:: *)

(* ::Title:: *)
(*TopWidth*)


Print["         (****** TopWidth-0.0 ******)  "];
Print["    Authors: Long-Bin Chen, Hai Tao Li, Jian Wang, YeFan Wang"];
Print["    TopWidth[QCDorder, mbCorr, WwidthCorr, EWcorr, mu] is provided for top width calculations "];
Print["    Please cite the paper for reference: arXiv:2211:XXXXX"];
Print[""];


BeginPackage["TopWidth`",{"HPL`"}]

(*topDec::usage = "topDec.m is a package that calcualtes top quark width."*)



TopWidth::usage="TopWidth Function:: TopWidth[order, mbCorr, WwidthCorr, EWcorr, mu] returns the top quark width in the Standard Model;
order = 0, 1, or 2; The QCD corrections are calculated accordingly;
mbCorr=0 or 1; if mbCorr=1, the b mass effects included at LO in QCD for order = 0, NLO in QCD for order=1;  mb effects in NNLO QCD is not available;
WwidthCorr=0 or 1; When WwidthCorr=0  the on-shell W is produced in top quark decay; WwidthCorr=1, the off-shell effects included up to NNLO in QCD according to order;
EWcorr=0 or 1; For EWcorr=1 and order>=1, the NLO EWcorrections are included.";
SetParameters::usage="TopWidth Function::SetParameters[mtv,mbv,mwv,wwidthv,GFv];
set the parameters for mt, mw, mb, wwidth and GF.
";



Begin["`Private`"]
(*Symbol["HPL`HPL"];
HPL=Symbol["HPL`HPL"];*)
(*mt=Symbol["Global`mt"];
mw=Symbol["Global`mw"];
mb=Symbol["Global`mb"];
\[CapitalGamma]w=Symbol["Global`\[CapitalGamma]w"];
GF=Symbol["Global`GF"];
rep2Num=Symbol["Global`rep2Num"];*)
(*Protected[repbeta, alsGeneralb2,alphasb2, X0, X1, Xl, Xh, XF, XA, X2, \[CapitalGamma]0, Li2, gammat,mbFun0,mbFun1,mbgammat ];*)
(* default parameter settings *)
Vtb=1;
SetParameters[mtv_,mbv_,mwv_,wwidthv_,GFv_]:=Module[{},
rep2Num={mt->mtv,mb->mbv,mw->mwv, \[CapitalGamma]w-> wwidthv,GF-> GFv};
]
rep2Num ={mt-> 17269/100,GF-> 11663788 10^-12, mb-> 478/100, mw->80377/1000 , \[CapitalGamma]w-> 2085/1000};
TopWidth[QCDorder_, mbCorr_, WwidthCorr_, EWcorr_, mu_]:=Module[{gt},
gt=0;
If[QCDorder<0 || QCDorder>2, Print["Please check the input of QCDorder: QCDorder should be 0,1,2;" Return[0]]];
If[mbCorr!=0 && mbCorr!=1, Print["Please check the input of mbCorr: mbCorr should be 0 or 1;" Return[0]]];
If[WwidthCorr!=0 && WwidthCorr!=1, Print["Please check the input of WwidthCorr: WwidthCorr should be 0 or 1;" Return[0]]];
If[EWcorr!=0 && EWcorr!=1, Print["Please check the input of EWcorr: EWcorr should be 0 or 1;" Return[0]]];
If[WwidthCorr==0, gt=tWidth0[QCDorder,mu];];
If[mbCorr==1, gt=gt+tWdithmbcorr[QCDorder,mu]];
If[EWcorr==0,gt=gt+tWidthEW[mu]];
If[WwidthCorr==1,gt=gt+tWidthGammaW[QCDorder,mu]];
gt//N
];

tWidth0[QCDorder_,mu_]:=Module[{als, x0,x1,x2,threshold},
threshold=0;
If[mu>= (mt/.rep2Num), threshold=1, threshold=0];
x0=X0/.w-> (mw/mt)^2/.rep2Num;
x1=0;
x2=0;
If[QCDorder>=1, x1=X1/.w-> (mw/mt)^2/.rep2Num];
If[QCDorder==2, x2=X2/.w-> (mw/mt)^2/.rep2Num/.nh-> threshold ];
als = alphasb2[mu^2];
Return[((GF*mt^3Vtb^2)/(8Sqrt[2]\[Pi])/.rep2Num)(x0 + als/\[Pi] x1 + If[x2!= 0, (als/\[Pi])^2 (x2 + 1/4 x1 (11-(2 (5+threshold))/3)Log[mu^2/(mt^2/.rep2Num)] ),0] )//N];
]

tWdithmbcorr[QCDorder_,mu_]:=Module[{bcorr,als},
bcorr=0;
als = alphasb2[mu^2];
If[QCDorder>= 0,  bcorr = mbFun0/2-X0/.w-> (mw/mt)^2/.rep2Num;];
If[QCDorder>= 1,  bcorr = bcorr+ als/\[Pi] (-1/3 mbFun1-X1)/.w-> (mw/mt)^2/.rep2Num;];
bcorr((GF*mt^3Vtb^2)/(8Sqrt[2]\[Pi])/.rep2Num)
]

tWidthGammaW[QCDorder_,mu_]:=Module[{als,x0, x1, x2,threshold},
threshold=0;
If[mu>= (mt/.rep2Num), threshold=1, threshold=0];
x1=0;
x2=0;
x0=nIntegrate[X0,w];
If[QCDorder>=1,x1=nIntegrate[X1,w];];
If[QCDorder==2,x2=nIntegrate[X2/.nh-> threshold,w];];
als = alphasb2[mu^2];
Return[((GF*mt^3Vtb^2)/(8Sqrt[2]\[Pi])/.rep2Num)(x0 + als/\[Pi] x1 + If[x2!= 0, (als/\[Pi])^2 (x2 + 1/4 x1 (11-(2 (5+threshold))/3)Log[mu^2/(mt^2/.rep2Num)] ),0] )//N];
];
nIntegrate[f_,w_]:=Module[{x},
1/\[Pi] NIntegrate[(f/.w-> x) ((mw/mt)^2 (\[CapitalGamma]w/mw))/((x-(mw/mt)^2)^2+(mw/mt)^4 (\[CapitalGamma]w/mw)^2)/.rep2Num,{x,0,1}]
];


tWidthEW[mu_]:=0;

(* running alphas-- approxiedmated 3loop results *)
repbeta={\[Beta]0-> 11-(2 nf)/3, \[Beta]1-> 102-(38 nf)/3,\[Beta]2-> 2857/2-(5033 nf)/18+(325 nf^2)/54};
alsGeneralb2[Q2_]:=  -((4 \[Pi] \[Beta]0)/\[Beta]1) 1/(1-(\[Beta]2 \[Beta]0)/\[Beta]1^2+ProductLog[-1,-(\[Beta]0/\[Beta]1 ) (lam2General/Q2)^(\[Beta]0^2/\[Beta]1) Exp[(\[Beta]2 \[Beta]0)/\[Beta]1^2-1]])/.repbeta//Simplify;

(* \[CapitalLambda]^2 from alphas[(*mz*) 91.1876] == 0.118  *)
lamb2Nf5={lam2General->0.16668689669525621`};
(* when mu>mt, nf=5, new \[CapitalLambda]^2 is provided *)
lamb2Nf6={lam2General->0.021763628679063205`};
alphasb2[mu2_]:=If[ mu2>= (mt^2/.rep2Num), alsGeneralb2[mu2]/.lamb2Nf6/.nf-> 6, 
alsGeneralb2[mu2]/.lamb2Nf5/.nf-> 5
]
Li2[x_]:=PolyLog[2,x];
(* analytical expressions up to NNLO with mb=0 and W width \[CapitalGamma]w= 0 *)
X0=2w^3-3w^2+1;
Xl=-((-1+w)*(-22-120*w+99*w^2))/36+(Pi^2*(23-12*w-111*w^2+124*w^3))/108-(w*(-86-25*w+106*w^2)*HPL[{0},w])/36+(Pi^2*(-1+w)^2*(1+2*w)*HPL[{1},w])/3-((-1+w)*(6-143*w-35*w^2+124*w^3)*HPL[{1},w])/(36*w)+((5+15*w-39*w^2+7*w^3)*HPL[{2},w])/9+((-1+w)^2*(1+2*w)*HPL[{3},w])/3-((-1+w)*(-37-55*w+38*w^2)*HPL[{1,0},w])/18-((-1+w)^2*(5+4*w)*HPL[{1,1},w])/3-((-1+w)^2*(1+2*w)*HPL[{2,0},w])/3+(2*(-1+w)^2*(1+2*w)*HPL[{2,1},w])/3-(2*(-1+w)^2*(1+2*w)*HPL[{1,1,0},w])/3+(-1+w)^2*(1+2*w)*Zeta[3];
Xh=(12775-12528*w-9237*w^2+15902*w^3)/1296-(Pi^2*(-53+29*w+99*w^2-137*w^3+14*w^4))/(54*(-1+w))+((9+344*w-498*w^2+168*w^3+265*w^4)*HPL[{1},w])/(54*w)-((23-8*w-18*w^2+32*w^3+19*w^4)*HPL[{2},w])/(9*(-1+w))+((1-12*w-3*w^2+2*w^3)*HPL[{3},w])/3-((1-12*w-3*w^2+2*w^3)*Zeta[3])/3;
XF=((-1+w)*(-86-106*w+177*w^2))/16+(Pi^4*(-11-328*w-191*w^2+42*w^3))/720-(Pi^2*(119-120*w+177*w^2+120*w^3))/48-(Pi^2*(1+4*w-3*w^2+2*w^3)*HPL[{-2},w])/6+(Pi^2*(-1+w)*(3+8*w+5*w^2)*HPL[{-1},w])/6+(Pi^2*w*(4-13*w+16*w^2)*HPL[{0},w])/48+((-4-196*w-385*w^2+8*w^3)*HPL[{0},w])/16-(Pi^2*(-1+w)*(37+66*w+15*w^2)*HPL[{1},w])/12+((-1+w)*(-2-175*w-449*w^2+34*w^3)*HPL[{1},w])/(16*w)+(Pi^2*(15+76*w-3*w^2+18*w^3)*HPL[{2},w])/12-((2+w-220*w^2-159*w^3+80*w^4)*HPL[{2},w])/(8*w)-((3-4*w-9*w^2+18*w^3)*HPL[{3},w])/2-((3+4*w-2*w^2+4*w^3)*HPL[{4},w])/2-2*(1+4*w-3*w^2+2*w^3)*HPL[{-2,2},w]+2*(-1+w)*(3+8*w+5*w^2)*HPL[{-1,2},w]-(Pi^2*(-1+w)^2*(1+2*w)*HPL[{1,0},w])/12+((-1+w)*(-2-63*w-99*w^2+22*w^3)*HPL[{1,0},w])/(8*w)-((-1+w)^2*(-1+24*w+29*w^2)*HPL[{1,1},w])/(4*w)+(3*(-1+w)^2*HPL[{1,2},w])/2-3*(-1+w)^2*(1+2*w)*HPL[{1,3},w]+((12+10*w-15*w^2+2*w^3)*HPL[{2,0},w])/4+((-54+4*w+57*w^2+4*w^3)*HPL[{2,1},w])/4-(-1+w)^2*(1+2*w)*HPL[{2,2},w]+((3+16*w-2*w^2+4*w^3)*HPL[{3,0},w])/2+w*(-16-7*w+2*w^2)*HPL[{3,1},w]+((-1+w)*(-26-25*w+w^2)*HPL[{1,1,0},w])/2+(3*(-1+w)^2*(1+2*w)*HPL[{1,2,0},w])/2-((-1-28*w-11*w^2+2*w^3)*HPL[{2,1,0},w])/2-(3*Pi^2*(-12+w^2)*(HPL[{-1},1-w]-Log[2]))/16-(Pi^2*(-76+27*w^2+16*w^3)*Log[2])/16-((-12+w^2)*(HPL[{-1,0,0},1-w]-(3*Zeta[3])/4))/4+((-12+w^2)*((Pi^2*HPL[{-1},1-w])/6-HPL[{-1,2},1-w]-(5*Zeta[3])/8))/4+((-212-192*w+199*w^2+400*w^3)*Zeta[3])/32+6*(-1+w)^2*(1+2*w)*HPL[{1},w]*Zeta[3];
XA=((-1+w)*(-785-3381*w+1188*w^2))/576+(Pi^4*(11-312*w-385*w^2+86*w^3))/1440-(Pi^2*(-505-1272*w+1005*w^2+2614*w^3))/864+(Pi^2*(1+4*w-3*w^2+2*w^3)*HPL[{-2},w])/12-(Pi^2*(-1+w)*(3+8*w+5*w^2)*HPL[{-1},w])/12-(Pi^2*w*(-4+w+16*w^2)*HPL[{0},w])/32+(w*(-2420-2317*w+2542*w^2)*HPL[{0},w])/288-(Pi^2*(-1+w)*(-16+45*w+57*w^2)*HPL[{1},w])/24+((-1+w)*(66-1873*w-1387*w^2+1352*w^3)*HPL[{1},w])/(144*w)+(Pi^2*(11+44*w+33*w^2+10*w^3)*HPL[{2},w])/24+((-224+354*w+759*w^2+32*w^3)*HPL[{2},w])/72+((-13-24*w+27*w^2+16*w^3)*HPL[{3},w])/6-((1+16*w+8*w^2)*HPL[{4},w])/4+(1+4*w-3*w^2+2*w^3)*HPL[{-2,2},w]-(-1+w)*(3+8*w+5*w^2)*HPL[{-1,2},w]-(Pi^2*(-1+w)^2*(1+2*w)*HPL[{1,0},w])/8+((-1+w)*(-485-827*w+466*w^2)*HPL[{1,0},w])/72+((-1+w)^2*(65+166*w)*HPL[{1,1},w])/24+((-1+w)^2*(1+2*w)*HPL[{1,2},w])/2+(-1+w)^2*(1+2*w)*HPL[{1,3},w]+((34+42*w-117*w^2+38*w^3)*HPL[{2,0},w])/24-((68+72*w-231*w^2+124*w^3)*HPL[{2,1},w])/24+((1+4*w+8*w^2)*HPL[{3,0},w])/4-((1+4*w+8*w^2)*HPL[{3,1},w])/2+((-1+w)*(-22-73*w+41*w^2)*HPL[{1,1,0},w])/12-((-1+w)^2*(1+2*w)*HPL[{1,2,0},w])/4+((3+12*w+13*w^2+2*w^3)*HPL[{2,1,0},w])/4+(3*Pi^2*(-12+w^2)*(HPL[{-1},1-w]-Log[2]))/32+(Pi^2*(-76+27*w^2+16*w^3)*Log[2])/32+((-12+w^2)*(HPL[{-1,0,0},1-w]-(3*Zeta[3])/4))/8-((-12+w^2)*((Pi^2*HPL[{-1},1-w])/6-HPL[{-1,2},1-w]-(5*Zeta[3])/8))/8-((-36-96*w-425*w^2+560*w^3)*Zeta[3])/64-(3*(-1+w)^2*(1+2*w)*HPL[{1},w]*Zeta[3])/2;
X2=CF(TR*nl*Xl+TR*nh*Xh+CF*XF+CA*XA)/.{CF-> 4/3,CA-> 3, nl-> 5,TR-> 1/2};
X1=CF w (-((\[Pi]^2 (-1+w)^2 (1+2 w))/(3 w))+1/4 (4+5/w-15 w+6 w^2)-((-1+w)^2 (5+4 w) Log[1-w])/(2 w)+(-1+w+2 w^2) Log[w]-((-1+w)^2 (1+2 w) Log[1-w] Log[w])/w-(2 (-1+w)^2 (1+2 w) PolyLog[2,w])/w)/.{CF-> 4/3,CA-> 3, nl-> 5,TR-> 1/2};
\[CapitalGamma]0=(GF*mt^3Vtb^2)/(8Sqrt[2]\[Pi]);
gammat=\[CapitalGamma]0(X0+\[Alpha]s/\[Pi] X1  +(\[Alpha]s/\[Pi])^2*X2  );

(* analytical expressions for finite mb *)
mbFun0=2 ((1-mb^2/mt^2)^2+((1+mb^2/mt^2) mw^2)/mt^2-(2 mw^4)/mt^4) Sqrt[1+mb^4/mt^4+mw^4/mt^4-2 (mb^2/mt^2+(mb^2 mw^2)/mt^4+mw^2/mt^2)];
mbFun1=-(1/mt^6)(mt^2 (5 mb^4-22 mb^2 mt^2+5 mt^4+9 (mb^2+mt^2) mw^2-6 mw^4) Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]-6 mt^2 (-mb^2+mt^2) (mb^2+mt^2-mw^2) Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4] Log[mb/mt]+4 mt^2 ((mb^2-mt^2)^2+(mb^2+mt^2) mw^2-2 mw^4) Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4] (3 Log[mb/mt]+2 Log[mw/mt]-2 Log[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])+(mb^6-mb^4 (11 mt^2+2 mw^2)-3 (mt^2-mw^2)^2 (mt^2+4 mw^2)+mb^2 (mt^4+12 mt^2 mw^2+5 mw^4)) Log[(mb^2+mt^2-mw^2-mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])]+4 (mb^2-mt^2) ((mb^2-mt^2)^2+(mb^2+mt^2) mw^2-4 mw^4) Log[-((mb^2-mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]-1/2 (mb^2+mt^2-mw^2) (mb^4+mt^4+mt^2 mw^2-2 mw^4+mb^2 (-2 mt^2+mw^2)) (4 \[Pi]^2-16 Li2[(mb^2+mt^2-mw^2-mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])]-16 Li2[(mb^2-mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2-mt^2+mw^2-mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])]-8 Li2[(2 mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])]+8 Li2[-((mb^2-mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]-4 Log[-((2 mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2-mt^2+mw^2-mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]^2+4 Log[(2 mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])] Log[(2 mt^4 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mw^2 (mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]-4 Log[(8 mt^6 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/((mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]) (-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])^2)] Log[-((mb^2-mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]-8 Log[(mb^2+mt^2-mw^2-mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])/(mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])] Log[(8 mt^2 (mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2)))/((mb^2+mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])^2 (-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))]+Log[-((mw^2 (-mb^2+mt^2+mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4]))/(mt^2 (mb^2-mt^2-mw^2+mt^2 Sqrt[(mb^4+(mt^2-mw^2)^2-2 mb^2 (mt^2+mw^2))/mt^4])))]^2));
mbgammat=(GF mt^3 Vtb^2)/(16 Sqrt[2] \[Pi]) ( matF0 -2 /3 \[Alpha]s/\[Pi] matF1  );

End[]

EndPackage[]



