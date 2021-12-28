#!/usr/bin/env wolframscript

exportToPov = True;
highlightCenter = False;
topView = False;
stellateFaces = True;

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 7; q = 3; (* default *)
];

stellationLevel = If[Length[argv] >= 3, ToExpression[argv[[3]]], 1];

nGenerations = If[Length[argv] >= 4, ToExpression[argv[[4]]], 3];


ratio = Tan[stellationLevel Pi/p]/Tan[Pi/p];
Print["Stellation ratio: " <> ToString[N[ratio]]];

shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm];
colors = <|
	3 -> Cyan,
	5 -> Blue,
	7 -> Red
|>;

color = If[KeyExistsQ[colors, p], colors[p], Red];

splitEdgeParts = 8;

dataFolder = "data";

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], Sin[t]},
  {0, -Sin[t], Cos[t]}
};
epsilon = 0.00000001;v
HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
getHyperboloid[v_] := {-v[[1]], v[[2]], v[[3]]};

scale[edge_, ratio_] := Module[
	{center},
	center = Mean[edge];
	Map[(# - center) * ratio + center&, edge]
];

scaleFace[face_, ratio_, stellationLevel_] := Module[
	{edges, edgeOrders},
	edges = Map[scale[#, ratio]&, getEdgesFromFace[face]];
	edgeOrders = Table[Mod[stellationLevel * k, p]+1, {k, 0, p-1}];
	(* edges[[edgeOrders]] *)
	Map[#[[1]]&, edges[[edgeOrders]]]
];


generateFaces[p_, q_] := Module[
	{},

	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

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
]

plotShape3D[faces_] := Module[
	{allEdges, edges, allVertices, vertices},

	(* selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &]; *)
	selectedFaces = faces;
	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];
	If[stellateFaces, 
		displayFaces = Map[scaleFace[#, ratio, stellationLevel]&, displayFaces]
	];

	Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces."];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges , SameTest -> sameCenter];
	Print["There are " <> IntegerString[Length[edges]] <> " edges in total."];
	(* norm = HNorm[edges[[1]][[1]]]; *)

	(* splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1]; *)
	(* splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."]; *)
	displayEdges = Map[getHyperboloid, edges, {2}];
	displayEdges = Map[scale[#, ratio] &, displayEdges];

	allVertices = Flatten[displayEdges, 1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];


	centeredFaces = Select[faces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	displayCenteredFaces = Map[getHyperboloid, centeredFaces, {2}];
	displayCenteredFaces = Map[(# + {0.001, 0, 0})&, displayCenteredFaces, {2}];

	centeredEdges = Flatten[Map[getEdgesFromFace, centeredFaces], 1];
	

	(* splitCenteredEdges = Map[HNormalize[#, norm]&, splitCenteredEdges, {2}]; *)
	Print["There are " <> IntegerString[Length[centeredEdges]] <> " centered edges after split."];
	displayCenteredEdges = Map[getHyperboloid, centeredEdges, {2}];

	(* centeredVertices = displayCenteredFaces[[1]]; *)

	rangeZoom = Sqrt[ratio];
	lighting = {{"Point", White, {10, -10, 10}}};
	range = If[topView, 1.6, 3] * rangeZoom;
	thickness = range / 80.0 / rangeZoom;
	viewCenter = {-5, 0, 0};

	elements = {
			{
				Opacity[1], FaceForm[color], Map[Polygon, displayFaces]
			},
			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[Purple], Map[Sphere[#, thickness * rangeZoom]&, vertices]
			}
	};

	If[highlightCenter,
		highlightColor = White;
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[highlightColor], Map[Polygon, displayCenteredFaces]
				},
				{
					Opacity[1], FaceForm[highlightColor], Map[Cylinder[#, thickness * 2.2]&, displayCenteredEdges]
				}
			}
		];
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> If[topView, {10, 0, 0}, {3, -3, 3}],
		ViewVertical -> If[topView, {0, 1, 0}, {1, 0, 0}],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];
	graph
];


Print[shape];

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
Print["Reading data for " <> shape <> " from " <> dataFileName];

(* faces = Get[dataFileName]; *)
faces = generateFaces[p, q];

graph = plotShape3D[faces];

outputFileName = If[highlightCenter, "highlightCenter_", ""] <> If[topView, "topView_", ""] <> If[stellateFaces, "stellateFaces_", ""] <> shape <> If[stellationLevel > 1, "_" <> IntegerString[stellationLevel], ""] <> "_" <> IntegerString[nGenerations] <> ".png";

imageSize = {600, 450};

imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, shape}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
outputFileName = FileNameJoin[{outputFolder, outputFileName}];

cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

If[exportToPov,
  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
  Run[cropCommand];
  ,
  Export[outputFileName, graph];
];
Print["Saved an image to " <> outputFileName];
