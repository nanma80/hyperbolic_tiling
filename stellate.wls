#!/usr/bin/env wolframscript

exportToPov = False;
highlightCenter = False;
topView = False;

exportToPov = True;
(* highlightCenter = True; *)
(* topView = True; *)

shapes = {
	"h_7_3"
	(* "h_3_7", *)
	(* "h_72_7", *)
	(* "h_7_72" *)
};

stellationLevel = 3;

ratio = Tan[stellationLevel Pi/7]/Tan[Pi/7];

splitEdgeParts = 8;
colors = {Red, Blue, Purple};
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

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001;

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
	edgeOrders = Table[Mod[stellationLevel * k, 7]+1, {k, 0, 6}];
	(* edges[[edgeOrders]] *)
	Map[#[[1]]&, edges[[edgeOrders]]]
];

plotShape3D[faces_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];
	color = colors[[If[highlightCenter, 1, stellationLevel]]];

	(* selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &]; *)
	selectedFaces = faces;
	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	stellatedFaces = Map[scaleFace[#, ratio, stellationLevel]&, selectedFaces];
	displayStellatedFaces = Map[getHyperboloid, stellatedFaces, {2}];
	If[!highlightCenter, displayFaces = displayStellatedFaces];

	Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces."];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges, SameTest -> sameCenter];
	If[!highlightCenter, edges = Map[scale[#, ratio]&, edges]];
	
	(* normSquares = Map[HNormSquare, Flatten[edges, 1]] is always constant, no matter what ratio is *)

	Print["There are " <> IntegerString[Length[edges]] <> " edges in total."];
	(* norm = HNorm[edges[[1]][[1]]]; *)

	(* splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1]; *)
	(* splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."]; *)
	displayEdges = Map[getHyperboloid, edges, {2}];

	(* allVertices = Flatten[displayFaces,1]; *)
	allVertices = Flatten[displayEdges, 1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	facesPerVertex = Table[
		Length[Select[displayStellatedFaces, (Length[Intersection[#, {vertices[[k]]}, SameTest -> ApproxSamePoint]] > 0)&]], {k, 1, Length[vertices]}
	];
	Print["Each vertex belongs to up to n faces: n = " <> IntegerString[Max[facesPerVertex]]];
	(* 
		stellation level = 1: each vertex belongs to up to 3 faces
		stellation level = 2: each vertex belongs to up to 7 faces
		stellation level = 3: each vertex belongs to up to 8, 9, 10, (infinite?) faces
	 *)

	centeredFaces = Select[stellatedFaces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	displayCenteredFaces = Map[getHyperboloid, centeredFaces, {2}];
	displayCenteredFaces = Map[(# + {0.001, 0, 0})&, displayCenteredFaces, {2}];

	centeredEdges = Flatten[Map[getEdgesFromFace, centeredFaces], 1];

	(* splitCenteredEdges = Map[HNormalize[#, norm]&, splitCenteredEdges, {2}]; *)
	Print["There are " <> IntegerString[Length[centeredEdges]] <> " centered edges after split."];
	displayCenteredEdges = Map[getHyperboloid, centeredEdges, {2}];

	(* centeredVertices = displayCenteredFaces[[1]]; *)

	rangeZoom = If[highlightCenter, Sqrt[stellationLevel], stellationLevel];
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
				Opacity[1], FaceForm[Blue], Map[Sphere[#, thickness * 2]&, vertices]
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


shapeIndex = 1;

(* For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++, *)
shape = shapes[[shapeIndex]];

Print[shape];

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
Print["Reading data for " <> shape <> " from " <> dataFileName];

faces = Get[dataFileName];

graph = plotShape3D[faces];

outputFileName = If[highlightCenter, "highlightCenter_", ""] <> If[topView, "topView_", ""] <> IntegerString[stellationLevel] <> "_hyperboloid.png";
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

(* ] *)
