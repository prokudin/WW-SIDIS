(* ::Package:: *)

Print["Package WW-SIDIS contains the set of TMDs calculated with WW approximation and SIDIS structure functions"];
Print["Contains the following functions: "];
Print["Uses Mathematica package for MSTW PDFs by Graeme Watt <Graeme.Watt(at)cern.ch>"];
Print["Copyright: Alexei Prokudin (PSU Berks), Kemal Tezgin (UConn), Version 1 (05/21/2018)"];
Print["e-mail: prokudin@jlab.org"];
Print["https://github.com/prokudin/WW-SIDIS"];
Print["___________________________________________________________________________"];


(*BeginPackage["wwsidis`",{"mstwpdf`"}];*)
BeginPackage["wwsidis`"];
f1u::usage="f1u[x_,Q2_ ] is the unpolarised collinear PDF for u quark";


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Read and interpolate DSS collinear fragmentation functions.*)


(* ::Input:: *)
(*DSShplus= ReadList["./Grids/fragmentationpiplus.dat",Real,RecordLists-> True];*)
(*DSShminus= ReadList["./Grids/fragmentationpiminus.dat",Real,RecordLists-> True];*)


(* ::Input:: *)
(*uhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,3]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}]*)
(*dhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,4]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];*)
(*shplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,5]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];*)
(*ubhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,6]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];*)
(*dbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,7]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];*)
(*sbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,8]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];*)


(* ::Input:: *)
(*uhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,3]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}]*)
(*dhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,4]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];*)
(*shminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,5]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];*)
(*ubhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,6]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];*)
(*dbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,7]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];*)
(*sbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,8]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];*)


(* ::Input:: *)
(*Plot[{z uhplus[z,3.],z dhplus[z,3.]},{z,0.01,1}]*)


(* ::Section::Closed:: *)
(*Read and interpolate MSTW collinear distribution functions. *)


(* ::Input:: *)
(*<<mstwpdf.m*)


(* ::Input:: *)
(*?ReadPDFGrid*)


(* ::Input:: *)
(*?xf*)


(* ::Input:: *)
(*Use of the central PDF set*)


(* ::Input:: *)
(*prefix="Grids/mstw2008lo";*)
(*Timing[ReadPDFGrid[prefix,0]]*)


(* ::Input:: *)
(*Clear[upv,dnv,usea,dsea,str,sbar];*)


(* ::Input:: *)
(*upv[x_,Q2_]:=xf[0,x,Sqrt[Q2],8]/x*)
(*dnv[x_,Q2_]:=xf[0,x,Sqrt[Q2],7]/x*)
(*usea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-2]/x*)
(*dsea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-1]/x*)
(*str[x_,Q2_]:=xf[0,x,Sqrt[Q2],3]/x*)
(*sbar[x_,Q2_]:=xf[0,x,Sqrt[Q2],-3]/x*)
(*up[x_,Q2_]:=upv[x,Sqrt[Q2]]+usea[x,Sqrt[Q2]]*)
(*dn[x_,Q2_]:=dnv[x,Sqrt[Q2]]+dsea[x,Sqrt[Q2]]*)
(*upbar[x_,Q2_]:=usea[x,Sqrt[Q2]]*)
(*dnbar[x_,Q2_]:=dsea[x,Sqrt[Q2]]*)
(**)


(* ::Input:: *)
(*Plot[{x up[x,1.],x dn[x,1.]},{x,0.01,0.99}]*)


(* ::Section::Closed:: *)
(*Read and interpolate helicity collinear distribution functions. Gluck:1998xa*)


(* ::Input:: *)
(*g1param= ReadList["./Grids/g1.dat",Real,RecordLists-> True];*)
(* *)


(* ::Input:: *)
(*g1u=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,3]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)
(*g1d=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,4]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)
(*g1s=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,5]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)
(*g1ubar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,6]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)
(*g1dbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,7]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)
(*g1sbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,8]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];*)


(* ::Input:: *)
(*Plot[{x*g1u[x,1.],x*g1d[x,1.]},{x,0.001,1}]*)


(* ::Section::Closed:: *)
(*Read and interpolate Soffer Bound collinear distribution functions.*)
(**)


(* ::Input:: *)
(*sb= ReadList["./Grids/SofferBound.dat",Real,RecordLists-> True];*)
(* *)


(* ::Input:: *)
(*sb1u=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,3]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)
(*sb1d=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,4]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)
(*sb1s=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,5]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)
(*sb1ubar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,6]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)
(*sb1dbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,7]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)
(*sb1sbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,8]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];*)


(* ::Input:: *)
(*Plot[{x*sb1u[x,1.],x*sb1d[x,1.]},{x,0.01,1}]*)


(* ::Chapter:: *)
(*8 "basis" functions Subscript[f, 1], Subscript[D, 1], Subscript[h, 1], Subscript[g, 1], *)
(*\!\(\*SubsuperscriptBox[\(H\), \(1\), \(\[UpTee]\)]\), \!\( *)
(*\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)], *)
(*\*SubsuperscriptBox[\(h\), \(1  T\), \(\[UpTee]\)]\), *)
(*\!\(\*SubsuperscriptBox[\(h\), \(1\), \(\[UpTee]\)]\)*)


(* ::Section::Closed:: *)
(*Subscript[f, 1]*)


(* ::Input:: *)
(*(*2005 fit*)*)
(*Clear[avk];*)
(*avk=0.25;*)
(*(* distribution*)*)
(**)
(**)
(*f1u[x_,Q2_ ]:= up[x,Q2];*)
(*f1d[x_,Q2_ ]:= dn[x,Q2] ;*)
(*f1ubar[x_,Q2_ ]:= upbar[x,Q2];  *)
(*f1dbar[x_,Q2_ ]:= dnbar[x,Q2] ; *)
(*f1s[x_,Q2_ ]:= str[x,Q2];*)
(*f1sdbar[x_,Q2_ ]:= sbar[x,Q2];*)
(**)
(**)
(*f1uTMD[x_,Q2_ ,kt_]:= up[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];*)
(*f1dTMD[x_,Q2_,kt_ ]:= dn[x,Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];*)
(*f1ubarTMD[x_,Q2_,kt_ ]:= upbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];*)
(*f1dbarTMD[x_,Q2_ ,kt_]:= dnbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];*)
(*f1sTMD[x_,Q2_,kt_ ]:= str[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];*)
(*f1sdbarTMD[x_,Q2_,kt_ ]:= sbar[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];*)
(**)
(**)


(* ::Input:: *)
(*Plot[{x f1u[x,1.7], x f1d[x,1.7]},{x,0.001,0.99}]*)


(* ::Input:: *)
(*Plot[{f1uTMD[0.1,1.,kt],f1dTMD[0.1,1.,kt]},{kt,0,1}]*)


(* ::Input:: *)
(**)
(*f1plot =Plot[{x f1u[x,2.4],x f1d[x,2.4]},{x,0.01,1},Frame-> True,FrameTicks-> {{Automatic,None},{Automatic,None}},PlotRange->{{0,1},{0,1}},PlotStyle->{{Thick,Red},{Thick,Blue},{Thick,Red,Dashed},{Thick,Blue,Dashed}}, FrameLabel->{"x","\!\(\*SubsuperscriptBox[\(xf\), \(1\), \(q\)]\)(x)"},BaseStyle->{FontSize->12,FontFamily->"Arial"},PlotLegends->Placed[{"u","d","\!\(\*OverscriptBox[\(u\), \(_\)]\)","\!\(\*OverscriptBox[\(d\), \(_\)]\)"},{0.7,0.7}]]*)


(* ::Section::Closed:: *)
(* End *)
(**)


End[];
EndPackage[];
