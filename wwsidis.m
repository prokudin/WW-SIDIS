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
AUTh1tp::usage="AUTh1tp[pion_, x_, z_, Q2_, PT_] is the Pretzelosity asymmetry";
AUUcos2phi::usage="AUUcos2phi[pion_, x_, z_, Q2_, PT_] is the Boer-Mulders asymmetry";
ALT::usage="ALT[pion_, x_, z_, Q2_, PT_] is the Double Spin asymmetry";
AULsin2phi::usage="AULsin2phi[pion_, x_, z_, Q2_, PT_] is the Kotzinian-Mulders asymmetry";
ALTcosPhi::usage="ALTcosPhi[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AULsinphi::usage="AULsinphi[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
ALLcosphi::usage="ALLcosphi[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
ALTcos2PhiminusPhiS::usage="ALTcos2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AUUcosphi::usage="AUUcosphi[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AUTsinphiS::usage="AUTsinphiS[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AUTsin2PhiminusPhiS::usage="AUTsin2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] is the ... asymmetry";
AULsin2phiIntegrated::usage="AULsin2phiIntegrated[pion_, x_, z_, Q2_] is the P_T integrated Kotzinian-Mulders asymmetry";
AUTsin2PhiminusPhiSIntegrated::usage="AUTsin2PhiminusPhiSIntegrated[pion_, x_, z_, Q2_] is the P_T integrated ... asymmetry";
CrossSectionUU::usage="CrossSectionUU[pion_,x_,z_,Q2_,PT_,energy_,phi_] calculates the differential cross section of unpolarized beam and target";

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
 
(*Leading Structure Functions and Asymmetries*)
 
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
 
 
(*F_{UT}^{Sin[3\phi_h-\phi_S]}*)
avkP = avk MTT^2/(avk + MTT^2);
avUTPTh1tp[z_] := avp MC^2/(avp + MC^2) + avkP z^2;

GUTh1tp[z_, PT_] := (Mp^4 Mh)/(\[Pi] avUTPTh1tp[z]^4) Exp[(-PT^2/avUTPTh1tp[z])]
AUTh1tp[z_] := ( 3 Mh Mp^2 z^3 Sqrt[\[Pi]])/( 2 Sqrt[avUTPTh1tp[z]^3])


FUTh1tp[pion_, x_, z_, Q2_, PT_] := (2 (z^3) (PT^3) )/Mp^2 ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GUTh1tp[z, PT]


FUTh1tpIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AUTh1tp[z]

AUTh1tp[pion_, x_, z_, Q2_, PT_] := FUTh1tp[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTh1tpIntegrated[pion_, x_, z_, Q2_] := FUTh1tpIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UU}^{Cos[2\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avUUcos2phiPT[z_] := avc + avkBM z^2;

GfactorUUcos2phi[z_, PT_] := (4 Mp^3 Mh)/(\[Pi] avUUcos2phiPT[z]^3) Exp[(-PT^2/avUUcos2phiPT[z])]
AfactorUUcos2phi[z_] := ( 4 Mh Mp (z^2)  )/ avUUcos2phiPT[z]

FUUcos2phi[pion_, x_, z_, Q2_, PT_] := (z^2 PT^2)/Mp^2 ((4.0/9.0) h1perpuFirstMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1perpdFirstMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorUUcos2phi[z, PT]


FUUcos2phiIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1perpuFirstMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1perpdFirstMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorUUcos2phi[z]

AUUcos2phi[pion_, x_, z_, Q2_, PT_] := FUUcos2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUUcos2phiIntegrated[pion_, x_, z_, Q2_] := FUUcos2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
(*F_{LT}^{Cos[\phi_h-\phi_S]}*)
avkg1T = avkg; (*assumption*)
avLLPT[z_] := avp + avkg1T z^2;

GLT[z_, PT_] := (2 Mp^2)/(\[Pi] avLLPT[z]^2) Exp[(-PT^2/avLLPT[z])]
ALT[z_] := ( Mp z Sqrt[\[Pi]])/ Sqrt[avLLPT[z]]

FLT[pion_, x_, z_, Q2_, PT_] := (z PT)/Mp ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) GLT[z, PT]


FLTIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) ALT[z]

ALT[pion_, x_, z_, Q2_, PT_] := FLT[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTIntegrated[pion_, x_, z_, Q2_] := FLTIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
(*F_{UL}^{Sin[2\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avkh1l = avk; (* assumption *)
avULsin2phiPT[z_] := avc + avkh1l z^2;

GfactorULsin2phi[z_, PT_] := (4 Mp^3 Mh)/(\[Pi] avULsin2phiPT[z]^3) Exp[(-PT^2/avULsin2phiPT[z])]
AfactorULsin2phi[z_] := ( 4 Mh Mp (z^2)  )/ avULsin2phiPT[z]

FULsin2phi[pion_, x_, z_, Q2_, PT_] := (z^2 PT^2)/Mp^2 ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorULsin2phi[z, PT]


FULsin2phiIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorULsin2phi[z]

AULsin2phi[pion_, x_, z_, Q2_, PT_] := FULsin2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AULsin2phiIntegrated[pion_, x_, z_, Q2_] := FULsin2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*Sub-leading Structure Functions and Asymmetries*)

(*F_{LT}^{Cos[\phi_S]}*)
avLLPT[z_] := avp + avkg z^2;

GLTCOSPhi[z_, PT_] := 1 /(\[Pi] avLLPT[z]) Exp[(-PT^2/avLLPT[z])]
ALTCOSPhi[z_] := 1.

FLTcosPhi[pion_, x_, z_, Q2_, PT_] := (-2 Mp x)/Sqrt[Q2] ((4.0/9.0) gTu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) gTd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) gTs[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) gTubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) gTdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) gTsbar[x, Q2] D1sbar[pion, z, Q2] ) GLTCOSPhi[z, PT];


FLTcosPhiIntegrated[pion_, x_, z_, Q2_] := (-2 Mp x)/Sqrt[Q2] ((4.0/9.0) gTu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) gTd[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) gTs[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) gTubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) gTdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) gTsbar[x, Q2] D1sbar[pion, z, Q2] ) ALTCOSPhi[z];

ALTcosPhi[pion_, x_, z_, Q2_, PT_] := FLTcosPhi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTcosPhiIntegrated[pion_, x_, z_, Q2_] := FLTcosPhiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UL}^{Sin[\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avkh1l = avk; (* assumption *)
avULsinphiPT[z_] := avc + avkh1l z^2;
(* x hL = h_1L^perp(1) *)

GfactorULsinphi[z_, PT_] := (2 (Mp^2) )/(\[Pi] avULsinphiPT[z]^2) Exp[(-PT^2/avULsinphiPT[z])]
AfactorULsinphi[z_] := ( Mh z Sqrt[Pi]  )/ (avULsinphiPT[z]^(1/2))

FULsinphi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (z PT Mh)/Mp^2 (-2) ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorULsinphi[z, PT]


FULsinphiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) (-2) ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorULsinphi[z]

AULsinphi[pion_, x_, z_, Q2_, PT_] := FULsinphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AULsinphiIntegrated[pion_, x_, z_, Q2_] := FULsinphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{LL}^{Cos[\phi_h]}*)
avLLcosphiPT[z_] := avp + avkg z^2;

GfactorLLcosphi[z_, PT_] := avkg/(\[Pi] avLLcosphiPT[z]^2) Exp[(-PT^2/avLLcosphiPT[z])]
AfactorLLcosphi[z_] := ( z Sqrt[Pi] avkg)/(2 Mp Sqrt[avLLcosphiPT[z]]) ;

FLLcosphi[pion_, x_, z_, Q2_, PT_] := (-((2 Mp)/Sqrt[Q2])) (z PT)/Mp ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) GfactorLLcosphi[z, PT]


FLLcosphiIntegrated[pion_, x_, z_, Q2_] := (-((2 Mp)/Sqrt[Q2])) ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) AfactorLLcosphi[z]

ALLcosphi[pion_, x_, z_, Q2_, PT_] := FLLcosphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALLcosphiIntegrated[pion_, x_, z_, Q2_] := FLLcosphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{LT}^{Cos[2\phi_h-\phi_S]}*)
avkg1T = avkg; (*assumption*)
avLTcos2phiPT[z_] := avp + avkg1T z^2;

GLTcos2phi[z_, PT_] := (Mp^2 avkg1T)/(\[Pi] avLTcos2phiPT[z]^3) Exp[(-PT^2/avLTcos2phiPT[z])]
AfactorLTcos2phi[z_] := -(( z^2 avkg1T)/ avLTcos2phiPT[z])

FLTcos2phi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (-((z^2 PT^2)/Mp^2)) ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) GLTcos2phi[z, PT]


FLTcos2phiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) AfactorLTcos2phi[z]

ALTcos2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] := FLTcos2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTcos2PhiminusPhiSIntegrated[pion_, x_, z_, Q2_] := FLTcos2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UU}^{cos\phi_h}*)
avUUcosphi2PT[z_] := avp + avk z^2;
GfactorUUcosphi2[z_, PT_] := 1/(\[Pi] avUUcosphi2PT[z]) Exp[(-PT^2/avUUcosphi2PT[z])]
AfactorUUcosphi2[z_] := (  z Sqrt[Pi] avk)/(2 Mp  Sqrt[avUUcosphi2PT[z]]) ;

FUUcosphi[pion_, x_, z_, Q2_, PT_] := (-((2 Mp)/Sqrt[Q2])) (z PT 2 Mp)/avUUcosphi2PT[z] avk/(2 Mp^2) ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sdbar[x, Q2] D1sbar[pion, z, Q2]) GfactorUUcosphi2[z, PT]


FUUcosphiIntegrated[pion_, x_, z_, Q2_] := (-((2 Mp)/Sqrt[Q2])) ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sdbar[x, Q2] D1sbar[pion, z, Q2]) AfactorUUcosphi2[z]

AUUcosphi[pion_, x_, z_, Q2_, PT_] := FUUcosphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUUcosphiIntegrated[pion_, x_, z_, Q2_] := FUUcosphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UT}^{Sin[\phi_S]}*)
avkh1 = avk; (*asumption*)
avc = MC^2 avp/(MC^2 + avp);
avUTsinphi2PT[z_] := avc + avkh1 z^2;
GfactorUTsinphi2[z_, PT_] :=  (1 - PT^2/avUTsinphi2PT[z])/(\[Pi] avUTsinphi2PT[z]) Exp[(-PT^2/avUTsinphi2PT[z])]
AfactorUTsinphi2[z_] := 0. ;

FUTsinphiS[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (4 z^2 Mh Mp)/avUTsinphi2PT[z] avkh1/(2 Mp^2) ((4.0/9.0) h1u[x, Q2] If[pion == "pi+", H1perpFavFirstMoment[z, Q2], H1perpUnfFirstMoment[z, Q2]] + (1.0/9.0) h1d[x, Q2] If[pion == "pi+", H1perpUnfFirstMoment[z, Q2], H1perpFavFirstMoment[z, Q2]]) GfactorUTsinphi2[z, PT]


FUTsinphiSIntegrated[pion_, x_, z_, Q2_] := 0

AUTsinphiS[pion_, x_, z_, Q2_, PT_] := FUTsinphiS[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTsinphiSIntegrated[pion_, x_, z_, Q2_] := FUTsinphiSIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UT}^{Sin[2\phi_h-\phi_S]}*)
avUTsin2phi1PT[z_] := avp + avks z^2;

GfactorUTsin2phi1[z_, PT_] := (Mp^2)  /(\[Pi] avUTsin2phi1PT[z]) Exp[(-PT^2/avUTsin2phi1PT[z])]
AfactorUTsin2phi1[z_] := (Mp^2 z^2)/ avUTsin2phi1PT[z] ;

avkP = avk MTT^2/(avk + MTT^2);
avc = MC^2 avp/(MC^2 + avp);
avUTsin2phi2PT[z_] := avc + avkP z^2;
GfactorUTsin2phi2[z_, PT_] := 1/(\[Pi] avUTsin2phi2PT[z]) Exp[(-PT^2/avUTsin2phi2PT[z])]
AfactorUTsin2phi2[z_] := (4 Mh Mp z^2)/ avUTsin2phi1PT[z] ;

FUTsin2phi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (z^2 PT^2)/avUTsin2phi1PT[z]^2 avks/Mp^2 ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) GfactorUTsin2phi1[z, PT] + ((2 Mp)/Sqrt[Q2]) (z^2 PT^2 4 Mp Mh)/avUTsin2phi2PT[z]^2 (-1) ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorUTsin2phi2[z, PT]


FUTsin2phiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) avks/Mp^2 ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) AfactorUTsin2phi1[z] + ((2 Mp)/Sqrt[Q2]) (-1) ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorUTsin2phi2[z]

AUTsin2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] := FUTsin2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTsin2PhiminusPhiSIntegrated[pion_, x_, z_, Q2_] := FUTsin2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(* Calculation of Unpolarized Cross Section *)
p1[y_]:=(1-y)/(1-y+y^2/2);
p2[y_]:=y (1-y/2)/(1-y+y^2/2);
p3[y_]:=(2-y) Sqrt[1-y]/(1-y+y^2/2);
p4[y_]:=y Sqrt[1-y]/(1-y+y^2/2);
coupling=1/137;
y[x_,Q2_,energy_]:=Q2/(2*0.938*energy*x);
CrossSectionUU[pion_,x_,z_,Q2_,PT_,energy_,phi_]:=coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2)*(FUU[pion,x,z,Q2,PT]+Cos[phi] p3[y[x,Q2,energy]] FUUcosphi[pion,x,z,Q2,PT]+Cos[2 phi] p1[y[x,Q2,energy]] FUUcos2phi[pion,x,z,Q2,PT])

End[];
EndPackage[];
