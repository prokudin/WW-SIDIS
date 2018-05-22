(* ::Package:: *)

Print["Package WW-SIDIS contains the set of TMDs calculated with WW approximation and SIDIS structure functions"];
Print["Copyright: Alexei Prokudin (PSU Berks), Kemal Tezgin (UConn), Version 1 (05/21/2018)"];
Print["e-mail: prokudin@jlab.org"];
Print["https://github.com/prokudin/WW-SIDIS"];
Print["___________________________________________________________________________"];
Print["Contains the following functions: "];
BeginPackage["wwsidis`"];
f1u::usage="f1u[x_,Q2_ ] is the unpolarised collinear PDF for u quark";
f1d::usage="f1d[x_,Q2_ ] is the unpolarised collinear PDF for d quark";
?f1u;
?f1d;
h1Lu::usage="h1Lu[x_,Q2_ ] is the unpolarised collinear PDF for u quark";
h1Ld::usage="f1d[x_,Q2_ ] is the unpolarised collinear PDF for d quark";
g2::usage="f1u[x_,Q2_ ] is the unpolarised collinear PDF for u quark";
g2n::usage="f1d[x_,Q2_ ] is the unpolarised collinear PDF for d quark";
g1Tperpu::usage="f1u[x_,Q2_ ] is the unpolarised collinear PDF for u quark";
g1Tperpd::usage="f1d[x_,Q2_ ] is the unpolarised collinear PDF for d quark";
ALL::usage="ALL[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AUTSivers::usage="AUTSivers[pion_, x_, z_, Q2_, PT_] is the Sivers asymmetry";
AUTCollins::usage="AUTCollins[pion_, x_, z_, Q2_, PT_] is the Collins asymmetry";
Begin["`Private`"];
DSShplus= ReadList["./Grids/fragmentationpiplus.dat",Real,RecordLists-> True];
DSShminus= ReadList["./Grids/fragmentationpiminus.dat",Real,RecordLists-> True];
uhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,3]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}]
dhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,4]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
shplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,5]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
ubhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,6]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
dbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,7]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
sbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,8]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
uhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,3]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}]
dhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,4]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
shminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,5]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
ubhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,6]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
dbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,7]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
sbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,8]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
<<mstwpdf.m;
prefix="Grids/mstw2008lo";
Timing[ReadPDFGrid[prefix,0]];
upv[x_,Q2_]:=xf[0,x,Sqrt[Q2],8]/x;
dnv[x_,Q2_]:=xf[0,x,Sqrt[Q2],7]/x;
usea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-2]/x;
dsea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-1]/x;
str[x_,Q2_]:=xf[0,x,Sqrt[Q2],3]/x;
sbar[x_,Q2_]:=xf[0,x,Sqrt[Q2],-3]/x;
up[x_,Q2_]:=upv[x,Sqrt[Q2]]+usea[x,Sqrt[Q2]];
dn[x_,Q2_]:=dnv[x,Sqrt[Q2]]+dsea[x,Sqrt[Q2]];
upbar[x_,Q2_]:=usea[x,Sqrt[Q2]];
dnbar[x_,Q2_]:=dsea[x,Sqrt[Q2]];
g1param= ReadList["./Grids/g1.dat",Real,RecordLists-> True];
g1u=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,3]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1d=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,4]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1s=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,5]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1ubar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,6]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1dbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,7]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1sbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,8]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
sb= ReadList["./Grids/SofferBound.dat",Real,RecordLists-> True];
sb1u=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,3]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1d=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,4]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1s=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,5]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1ubar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,6]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1dbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,7]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1sbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,8]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
(*2005 fit*)
Clear[avk];
avk=0.25;
(* distribution*)
f1u[x_,Q2_ ]:= up[x,Q2];
f1d[x_,Q2_ ]:= dn[x,Q2] ;
f1ubar[x_,Q2_ ]:= upbar[x,Q2];  
f1dbar[x_,Q2_ ]:= dnbar[x,Q2] ; 
f1s[x_,Q2_ ]:= str[x,Q2];
f1sdbar[x_,Q2_ ]:= sbar[x,Q2];


f1uTMD[x_,Q2_ ,kt_]:= up[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];
f1dTMD[x_,Q2_,kt_ ]:= dn[x,Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];
f1ubarTMD[x_,Q2_,kt_ ]:= upbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];
f1dbarTMD[x_,Q2_ ,kt_]:= dnbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];
f1sTMD[x_,Q2_,kt_ ]:= str[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];
f1sdbarTMD[x_,Q2_,kt_ ]:= sbar[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];

(*2005 fit*)
Clear[avp];
avp = 0.2;
(* fragmentation*)
 
D1uTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", uhplus[z, Q2],If[pion == "pi-", uhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1dTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", dhplus[z, Q2],If[pion == "pi-", dhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1ubarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", ubhplus[z, Q2],If[pion == "pi-", ubhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1dbarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", dbhplus[z, Q2],If[pion == "pi-", dbhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1sTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", shplus[z, Q2],If[pion == "pi-", shminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1sbarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", sbhplus[z, Q2],If[pion == "pi-", sbhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
 
 
D1u[pion_, z_, Q2_ ] := If[pion == "pi+", uhplus[z, Q2], If[pion == "pi-", uhminus[z, Q2]]]
D1d[pion_, z_, Q2_ ] := If[pion == "pi+", dhplus[z, Q2], If[pion == "pi-", dhminus[z, Q2]]]
D1ubar[pion_, z_, Q2_ ] := If[pion == "pi+", ubhplus[z, Q2],If[pion == "pi-", ubhminus[z, Q2]]]
D1dbar[pion_, z_, Q2_ ] := If[pion == "pi+", dbhplus[z, Q2],If[pion == "pi-", dbhminus[z, Q2]]]
D1s[pion_, z_, Q2_ ] := If[pion == "pi+", shplus[z, Q2],If[pion == "pi-", shminus[z, Q2]]]
D1sbar[pion_, z_, Q2_ ] := If[pion == "pi+", sbhplus[z, Q2], If[pion == "pi-", sbhminus[z, Q2]]]
 
(*Lattice fit*)
Clear[avkg];
avkg = avk 0.76 ;


(* Helicity functions with DGLAP Evolution*)


g1uTMD[x_, Q2_, kt_] := g1u[x, Q2] 1/(\[Pi] avkg) Exp[-kt^2/avkg];
g1dTMD[x_, Q2_, kt_] := g1d[x, Q2] 1/(\[Pi] avkg) Exp[-kt^2/avkg];
 
(*2013 fit*)
Clear[NuT, NdT, alphaT, betaT];
NuT = 0.46;
NdT = -1.000;
alphaT = 1.11;
betaT = 3.64;

(* Transversity function with DGLAP Evolution*)

h1u[x_, Q2_] := (NuT x^alphaT (1 - x)^betaT (alphaT + betaT)^(alphaT + betaT)/((alphaT^alphaT) (betaT^betaT))) sb1u[x, Q2];
h1d[x_, Q2_] := (NdT x^alphaT (1 - x)^betaT (alphaT + betaT)^(alphaT + betaT)/((alphaT^alphaT) (betaT^betaT))) sb1d[x, Q2];

h1uTMD[x_, Q2_, kt_] := h1u[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];
h1dTMD[x_, Q2_, kt_] := h1d[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];

(*2013 fit*)
Clear[MC, NCfav, NCunf, alphaC, betaC, Mh];
MC = Sqrt[1.5];
NCfav = 0.49;
NCunf = -1.000;
alphaC = 1.06;
betaC = 0.07;
Mh = 0.135;

(*Collins function with DGLAP Evolution*)

H1perpFavHalfMoment[x_, Q2_] := Sqrt[2 E] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2]/4 Sqrt[(\[Pi] avp)/(MC^2 + avp)^3];
H1perpUnfHalfMoment[x_, Q2_] := Sqrt[2 E] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2]/4 Sqrt[(\[Pi] avp)/(MC^2 + avp)^3];

H1perpFavFirstMoment[x_, Q2_] := Sqrt[E/2] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2] (MC^3 avp)/(Mh x (MC^2 + avp)^2);
H1perpUnfFirstMoment[x_, Q2_] := Sqrt[E/2] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2] (MC^3 avp)/(Mh x (MC^2 + avp)^2);

H1perpFavTMD[x_, Q2_, kt_] := (x Mh)/MC Exp[-kt^2/MC^2] Sqrt[2 E] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2] 1/(\[Pi] avp) Exp[-kt^2/avp];
H1perpUnfTMD[x_, Q2_, kt_] := (x Mh)/MC Exp[-kt^2/MC^2] Sqrt[2 E] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2] 1/(\[Pi] avp) Exp[-kt^2/avp];

H1perpFirstMoment[quark_, pion_, x_, Q2_] := If[(quark == "u" && pion == "pi+") || (quark == "bard" && pion == "pi+") || (quark == "d" && pion == "pi-") || (quark == "baru" && pion == "pi-"),H1perpFavFirstMoment[x, Q2], H1perpUnfFirstMoment[x, Q2]]
 
(*2015 fit*)
Clear[MC5, NCfav5, NCunf5, alphaC5,  Mh5, avp5];
MC5 = Sqrt[0.28];
NCfav5 = 0.9;
NCunf5 = -0.37;
alphaC5 = 2.02;

Mh5 = 0.135;
avp5 = 0.12;

MC5*MC5*avp5 /(avp5 + MC5*MC5)

(*Collins function with DGLAP Evolution*)

H1perpFavFirstMoment5[x_, Q2_] := Sqrt[E/2] (NCfav5 x^alphaC5) D1u["pi+", x, Q2] (MC5^3 avp5)/(Mh5 x (MC5^2 + avp5)^2);
H1perpUnfFirstMoment5[x_, Q2_] := Sqrt[ E/2] (NCunf5) D1d["pi+", x, Q2] (MC5^3 avp5)/(Mh5 x (MC5^2 + avp5)^2);
H1perpFirstMoment5[quark_, pion_, x_, Q2_] := If[(quark == "u" && pion == "pi+") || (quark == "bard" && pion == "pi+") || (quark == "d" && pion == "pi-") || (quark == "baru" && pion == "pi-"),H1perpFavFirstMoment5[x, Q2], H1perpUnfFirstMoment5[x, Q2]]
 
(*2011 fit*)

Clear[Ms, avks, Nu, Nusea, Nd, Ndsea, Nst, Nstbar, alphauv, alphadv, asea, beta, Mp];
Ms = Sqrt[0.19];
avks = avk Ms^2/(avk + Ms^2);
Nu = 0.4;
Nusea = 0.;
Nd = -0.97;
Ndsea = 0.;
Nst = 0.;
Nstbar = 0.;
alphauv = 0.35;
alphadv = 0.44;
asea = 1.;
beta = 3.46;
betauv = 2.6;
betadv = 0.9;
avPTs[z_] := Sqrt[avp^2 + avks^2 z^2]
Mp = 0.938;

(* Sivers function with DGLAP Evolution*)

usiv[x_, Q2_] := (Nu x^alphauv (1 - x)^betauv (alphauv + betauv)^(alphauv + betauv)/((alphauv^alphauv) ( betauv^betauv))) up[x, Q2]
dsiv[x_, Q2_] := (Nd x^alphadv (1 - x)^betadv (alphadv + betadv)^(alphadv + betadv)/((alphadv^alphadv) ( betadv^betadv))) dn[x, Q2]
ubarsiv[x_, Q2_] := (Nusea x^asea (1 - x)^beta (asea + beta)^(asea + beta)/((asea^asea) (beta^beta))) upbar[x, Q2]
dbarsiv[x_, Q2_] := (Ndsea x^asea (1 - x)^beta (asea + beta)^(asea + beta)/((asea^asea) (beta^beta))) dnbar[x, Q2]
ssiv[x_, Q2_] := (Nst x^asea (1 - x)^beta (asea + beta)^(asea + beta)/((asea^asea) ( beta^beta))) str[x, Q2]
sbarsiv[x_, Q2_] := (Nstbar x^asea (1 - x)^beta (asea + beta)^(asea + beta)/((asea^asea) (beta^beta))) sbar[x, Q2]
 
(*Sivers function with DGLAP Evolution*)

f1TperpuFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) usiv[x, Q2] avks^2/avk ;
f1TperpdFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) dsiv[x, Q2] avks^2/avk;
f1TperpubarFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) ubarsiv[x, Q2] avks^2/avk ;
f1TperpdbarFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) dbarsiv[x, Q2] avks^2/avk;
f1TperpsFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) ssiv[x, Q2] avks^2/avk ;
f1TperpsbFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp Ms) sbarsiv[x, Q2] avks^2/avk;


f1TperpuTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] usiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpdTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] dsiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpubarTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] ubarsiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpdbarTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] dbarsiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpsTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] ssiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpsbTMD[x_, Q2_, kt_] := -(Mp/Ms) Sqrt[2 E] sbarsiv[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avks];
 
 
(*2015 fit*)
Clear[MTT, avkTT, NuTT, NdTT, alphaTT, betaTT];
MTT = Sqrt[0.18];
avkTT = avk MTT^2/(avk + MTT^2);
NuTT = 1;
NdTT = -1;
alphaTT = 2.5;
betaTT = 2.;
(* Pretzelosity function with DGLAP Evolution*)

huTT[x_, Q2_] := E (NuTT x^alphaTT (1 - x)^betaTT (alphaTT + betaTT)^(alphaTT + betaTT)/((alphaTT^alphaTT) (betaTT^betaTT))) (up[x, Q2] - g1u[x, Q2])
hdTT[x_, Q2_] := E (NdTT x^alphaTT (1 - x)^betaTT (alphaTT + betaTT)^(alphaTT + betaTT)/((alphaTT^alphaTT) (betaTT^betaTT))) (dn[x, Q2] - g1d[x, Q2])


(*Pretzelosity function with DGLAP Evolution*)

h1TperpuFirstMoment[x_, Q2_] := 1/(2 MTT^2) huTT[x, Q2] avkTT^2/avk ;
h1TperpdFirstMoment[x_, Q2_] := 1/(2 MTT^2) hdTT[x, Q2] avkTT^2/avk;



h1TperpuTMD[x_, Q2_, kt_] := Mp^2/MTT^2  huTT[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avkTT];
h1TperpdTMD[x_, Q2_, kt_] := Mp^2/MTT^2  hdTT[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avkTT];


h1TperpuSecondMoment[x_, Q2_] := 1/(2 Mp^2 MTT^2) huTT[x, Q2] avkTT^3/avk ;
h1TperpdSecondMoment[x_, Q2_] := 1/(2 Mp^2 MTT^2) hdTT[x, Q2] avkTT^3/avk;
 
(*2015 fit*)
Clear[MBM, avkBM, NuBM, NdBM];
MBM = Sqrt[0.1];
avkBM = avk MBM^2/(avk + MBM^2);
NuBM = -0.49;
NdBM = -1;

(* Boer-Mulders function with DGLAP Evolution*)

h1perpu[x_, Q2_] := NuBM up[x, Q2]
h1perpd[x_, Q2_] := NdBM dn[x, Q2]


(*Boer-Mulders function with DGLAP Evolution*)

h1perpuFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp MBM) h1perpu[x, Q2] avkBM^2/avk;
h1perpdFirstMoment[x_, Q2_] := -Sqrt[ E/2] 1/(Mp MBM) h1perpd[x, Q2] avkBM^2/avk;


h1perpuTMD[x_, Q2_, kt_] := -Sqrt[ 2 E] Mp/MBM Exp[-kt^2/avkBM] h1perpu[x, Q2]/(\[Pi] avk);
h1perpdTMD[x_, Q2_, kt_] := -Sqrt[ 2 E] Mp/MBM Exp[-kt^2/avkBM] h1perpd[x, Q2]/(\[Pi] avk);

(*WW-type relations*)
 
h1Lu[x_, Q_] := -x^2 NIntegrate[h1u[y, Q]/y^2, {y, x, 1.}];
h1Ld[x_, Q_] := -x^2  NIntegrate[h1d[y, Q]/y^2, {y, x, 1.}];

gTu[x_, Q2_] := NIntegrate[g1u[y, Q2]/y, {y, x, 1.}];
gTd[x_, Q2_] := NIntegrate[g1d[y, Q2]/y, {y, x, 1.}];
gTubar[x_, Q2_] := NIntegrate[g1ubar[y, Q2]/y, {y, x, 1.}];
gTdbar[x_, Q2_] := NIntegrate[g1dbar[y, Q2]/y, {y, x, 1.}];
gTs[x_, Q2_] := NIntegrate[g1s[y, Q2]/y, {y, x, 1.}];
gTsbar[x_, Q2_] := NIntegrate[g1sbar[y, Q2]/y, {y, x, 1.}];

g2[x_, Q2_] := 1/2 * (4./9.  gTu[x, Q2] + 1./9. gTd[x, Q2] + 4./9. gTubar[x, Q2] + 1./9. gTdbar[x, Q2] + 1./9. gTs[x, Q2] + 1./9. gTsbar[x, Q2]) - 1/2 * (4./9. g1u[x, Q2] +  1./9. g1d[x, Q2] + 4./9. g1ubar[x, Q2] + 1./9. g1dbar[x, Q2] + 1./9. g1s[x, Q2] + 1./9. g1sbar[x, Q2]);

g2n[x_, Q2_] := 1/2 * (4./9.  gTd[x, Q2] + 1./9. gTu[x, Q2] + 4./9. gTdbar[x, Q2] + 1./9. gTubar[x, Q2] + 1./9. gTs[x, Q2] + 1./9. gTsbar[x, Q2]) - 1/2 * (4./9. g1d[x, Q2] +  1./9. g1u[x, Q2] + 4./9. g1dbar[x, Q2] + 1./9. g1ubar[x, Q2] + 1./9. g1s[x, Q2] + 1./9. g1sbar[x, Q2]);

g1Tperpu[x_, Q2_] := x NIntegrate[g1u[y, Q2]/y, {y, x, 1.}];
g1Tperpd[x_, Q2_] := x NIntegrate[g1d[y, Q2]/y, {y, x, 1.}];
g1Tperpubar[x_, Q2_] := x NIntegrate[g1ubar[y, Q2]/y, {y, x, 1.}];
g1Tperpdbar[x_, Q2_] := x NIntegrate[g1dbar[y, Q2]/y, {y, x, 1.}];
g1Tperps[x_, Q2_] := x NIntegrate[g1s[y, Q2]/y, {y, x, 1.}];
g1Tperpsbar[x_, Q2_] := x NIntegrate[g1sbar[y, Q2]/y, {y, x, 1.}];
 
(*F_{UU}*)
avPT[z_] := avp + avk z^2;


FUU[pion_, x_, z_, Q2_, PT_] := ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sdbar[x, Q2] D1sbar[pion, z, Q2]) Exp[(-PT^2/avPT[z])]/(\[Pi] avPT[z])

FUUIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sdbar[x, Q2] D1sbar[pion, z, Q2])

F2[x_, Q2_] := (4.0/9.0) f1u[x, Q2]  + (1.0/9.0) f1d[x, Q2] + (1.0/9.0) f1s[x, Q2]   + (4.0/9.0) f1ubar[x, Q2] + (1.0/9.0) f1dbar[x, Q2]  + (1.0/9.0) f1sdbar[x, Q2]

MULT[pion_, x_, z_, Q2_, PT_] := 2 \[Pi] PT FUU[pion, x, z, Q2, PT]/F2[x, Q2] ;
 
(*F_{LL}*)
avLLPT[z_] := avp + avkg z^2;


FLL[pion_, x_, z_, Q2_, PT_] := ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) Exp[(-PT^2/avLLPT[z])]/(\[Pi] avLLPT[z])


FLLIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2])

sjlab = Mp^2 + 2 5.7 Mp;
y[Q2_, x_] :=  Q2/(sjlab x);

ALL[pion_, x_, z_, Q2_, PT_] := (y[Q2, x] (2 - y[Q2, x] ) FLL[pion, x, z, Q2, PT])/((1 + (1 - y[Q2, x])^2) FUU[pion, x, z, Q2, PT]);

ALLIntegrated[pion_, x_, z_, Q2_] := (y[Q2, x] (2 - y[Q2, x] ) FLLIntegrated[pion, x, z, Q2])/((1 + (1 - y[Q2, x])^2) FUUIntegrated[pion, x, z, Q2]);

(*F_{UT}^{Sin\phi_h-Sin\phi_S}*)
avUTPTf1tperp[z_] := avp + avks z^2;

GUTf1tperp[z_, PT_] := -(Mp^2/(\[Pi] avUTPTf1tperp[z]^2)) Exp[(-PT^2/avUTPTf1tperp[z])]
AUTf1tperp[z_] := -(( Mp z Sqrt[\[Pi]])/ Sqrt[avUTPTf1tperp[z]])

FUTf1tperp[pion_, x_, z_, Q2_, PT_] := (2 z PT)/Mp ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) GUTf1tperp[z, PT]


FUTf1tperpIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) AUTf1tperp[z]

AUTSivers[pion_, x_, z_, Q2_, PT_] := FUTf1tperp[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTSiversIntegrated[pion_, x_, z_, Q2_] := FUTf1tperpIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];

(*F_{UT}^{Sin\phi_h+Sin\phi_S}*)
avUTPTh1[z_] := avp MC^2/(avp + MC^2) + avk z^2;

GUTh1[z_, PT_] := Mp^2/(\[Pi] avUTPTh1[z]^2) Exp[(-PT^2/avUTPTh1[z])]
AUTh1[z_] := ( Mh z Sqrt[\[Pi]])/ Sqrt[avUTPTh1[z]]

FUTh1[pion_, x_, z_, Q2_, PT_] := (2 z PT Mh)/Mp^2 ((4.0/9.0) h1u[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1d[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GUTh1[z, PT]


FUTh1Integrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1u[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1d[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AUTh1[z]

AUTCollins[pion_, x_, z_, Q2_, PT_] := FUTh1[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTCollinsIntegrated[pion_, x_, z_, Q2_] := FUTh1Integrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)
 (*F_{UU}*)

 
 
End[];
EndPackage[];
