#!/usr/bin/env wolframscript

nGenerations = 6;

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 3, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]]; r = ToExpression[argv[[3]]],
	p = 5; q = 3; r = 4;
];

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm];
dataFolder = "data";
epsilon = 0.00000001;
Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

phi = (1 + Sqrt[5])/2;
HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], Sin[t]},
  {0, -Sin[t], Cos[t]}
};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);


dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

piOverHpq = ArcCos[Sqrt[Cos[Pi/p]^2 + Cos[Pi/q]^2]]//Simplify;
piOverHqr = ArcCos[Sqrt[Cos[Pi/q]^2 + Cos[Pi/r]^2]]//Simplify;

(* chi: p0 and p2 *)
coshchi = Cos[Pi/p] Cos[Pi/q] Cos[Pi/r]/Sin[piOverHpq]/Sin[piOverHqr]//Simplify;
inner02 = coshchi;

(* psi: p1 and p2 *)
coshpsi = Sin[Pi/p] Cos[Pi/r]/Sin[piOverHpq]//Simplify;
inner12 = coshpsi;

(* phi: p0 and p1, 2phi: from p0 to p1 to another p0. edge length *)
coshphi = Cos[Pi/p] Sin[Pi/r]/Sin[piOverHqr]//Simplify;
inner01 = coshphi;
cosh2phi = 2 coshphi^2 - 1//Simplify;
sinh2phi = Sqrt[cosh2phi^2 - 1]//Simplify;

Print[{piOverHpq, piOverHqr}];
Print[{coshphi, coshchi, coshpsi}];

vertexFigureCoordinates = {};

If[q == 3 && r == 4,
	vertexFigureCoordinates = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, -1}
	};
];

If[q == 3 && r == 5,
	vertexFigureCoordinates = {
		{phi, 1, 0},
		{phi, -1, 0},
		{1, 0, phi},
		{1, 0, -phi}
	};
];

If[q == 5 && r == 3,
	vertexFigureCoordinates = {
		{phi, 1/phi, 0},
		{phi, -1/phi, 0},
		{1, 1, 1},
		{1, 1, -1}
	};
];

normalizedVertexFigureCoordinates = Map[#/Norm[#]&, vertexFigureCoordinates];
Print[normalizedVertexFigureCoordinates//N];

(* vertex first *)
vertex1 = {1, 0, 0, 0};

seedVertices = Map[{cosh2phi, sinh2phi #[[1]], sinh2phi #[[2]], sinh2phi #[[3]]}&, normalizedVertexFigureCoordinates];
(* Print["seedVertices"]; *)
Print[seedVertices//N];
mirrors = {
	seedVertices[[1]] - vertex1,
	seedVertices[[1]] - seedVertices[[2]],
	seedVertices[[2]] - seedVertices[[3]],
	seedVertices[[3]] - seedVertices[[4]]
};

Print[mirrors//N];

face1 = {vertex1}//N;
finished = False;
For[i = 0, i < 20 && finished == False, i++,
	newVertex = FullSimplify[HDoubleReflect[face1[[-1]], mirrors[[1]], mirrors[[2]]]];
	If[!(ApproxSamePoint[newVertex, face1[[1]]]),
		AppendTo[face1, newVertex],
		finished = True;
	]
];

Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

Print[face1];

Exit[];



(* chi: p0 and p2 *)
inner02 = Cot[Pi/p]*Cot[Pi/q];
(* psi: p1 and p2 *)
inner12 = Cos[Pi/q]/Sin[Pi/p];
(* phi: p0 and p1 *)
inner01 = Cos[Pi/p]/Sin[Pi/q];

ch2psi = 2 * inner12 ^ 2 - 1;

p2 = {1, 0, 0};
p0 = {inner02, Sqrt[inner02^2-1], 0}//Simplify;
p1 = {inner12, Sqrt[inner12^2-1] Cos[Pi/p], Sqrt[inner12^2-1] Sin[Pi/p]}//Simplify;
(* other p2 such that the midpoint of p2 and otherp2 is p1 *)
otherp2 = {ch2psi, Sqrt[ch2psi^2-1] Cos[Pi/p], Sqrt[ch2psi^2-1] Sin[Pi/p]};

vertex1 = p0;
mirror1 = {0, 0, 1};
mirror2 = {0, -Sin[Pi/p], Cos[Pi/p]};
mirror3 = otherp2 - p2//Simplify;

mirrors = {mirror1, mirror2, mirror3}//N;

face1 = {vertex1}//N;
finished = False;
For[i = 0, i < 20 && finished == False, i++,
	newVertex = FullSimplify[HDoubleReflect[face1[[-1]], mirrors[[1]], mirrors[[2]]]];
	If[!(ApproxSamePoint[newVertex, face1[[1]]]),
		AppendTo[face1, newVertex],
		finished = True;
	]
];

Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

faces = {face1};
countStart = -1;
For[i = 0, i <= nGenerations && countStart != Length[faces], i++,
	Print["i = " <> IntegerString[i]];
	countStart = Length[faces];
	For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
		newFaces = Map[HReflect[#, mirrors[[mirrorIndex]] ]&, faces, {2}];
		faces = Union[faces, newFaces, SameTest -> sameCenter];
	];
];

Print["Making it more symmetric"];

(* Making it symmetric *)
centeredMirrors = Select[mirrors, #[[1]] == 0 &];

countStart = -1;
For[i = 0, i < 20 && countStart != Length[faces], i++,
	Print["i = " <> IntegerString[i]];
	countStart = Length[faces];
	For[mirrorIndex = Length[centeredMirrors], mirrorIndex >= 1, mirrorIndex--,
		newFaces = Map[HReflect[#, centeredMirrors[[mirrorIndex]] ]&, faces, {2}];
		faces = Union[faces, newFaces, SameTest -> sameCenter];
	];
];
faces
