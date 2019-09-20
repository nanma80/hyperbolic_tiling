#!/usr/bin/env wolframscript

exportToPov = True;
highlightCenter = True;
topView = False;
level = 2;

shapes = {
	(* "h_7_3", *)
	"h_3_7"
	(* "h_72_7", *)
	(* "h_7_72" *)
};

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

plotShape3D[faces_, highlightedFace1_, highlightedFace2_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];
	color = configs[shape][[3]];

	selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &];
	(* selectedFaces = faces; *)
	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces."];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges , SameTest -> sameCenter];
	Print["There are " <> IntegerString[Length[edges]] <> " edges in total."];
	(* norm = HNorm[edges[[1]][[1]]]; *)

	(* splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1]; *)
	(* splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."]; *)
	displayEdges = Map[getHyperboloid, edges, {2}];

	allVertices = Flatten[displayFaces,1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];


	centeredFaces1 = {highlightedFace1};
	centeredEdges1 = Flatten[Map[getEdgesFromFace, centeredFaces1], 1];
	displayCenteredFaces1 = Map[getHyperboloid, centeredFaces1, {2}];
	displayCenteredEdges1 = Map[getHyperboloid, centeredEdges1, {2}];
	Print["There are " <> IntegerString[Length[displayCenteredEdges1]] <> " displayCenteredEdges1 in total."];

	centeredFaces2 = {highlightedFace2};
	centeredEdges2 = Flatten[Map[getEdgesFromFace, centeredFaces2], 1];
	displayCenteredFaces2 = Map[getHyperboloid, centeredFaces2, {2}];
	displayCenteredEdges2 = Map[getHyperboloid, centeredEdges2, {2}];
	Print["There are " <> IntegerString[Length[displayCenteredEdges2]] <> " displayCenteredEdges2 in total."];

	(* centeredVertices = displayCenteredFaces[[1]]; *)

	lighting = {{"Point", White, {10, -10, 10}}};
	range = If[topView, 1.6, 2];
	thickness = range / 80.0;
	viewCenter = {-4, 0, 0};

	elements = {
			{
				Opacity[0.5], FaceForm[color], Map[Polygon, displayFaces]
			},
			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 2]&, vertices]
			}
	};



	If[level==1,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[Yellow], Map[Polygon, displayCenteredFaces1]
				},
				{
					Opacity[1], FaceForm[Yellow], Map[Cylinder[#, thickness * 2]&, displayCenteredEdges1]
				}
			}
		];
	];

	If[level==2,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[Yellow], Map[Cylinder[#, thickness * 2]&, displayCenteredEdges1]
				},
				{
					Opacity[1], FaceForm[Darker[Blue]], Map[Cylinder[#, thickness * 2]&, displayCenteredEdges2]
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




For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
	shape = shapes[[shapeIndex]];

	Print[shape];

	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
	Print["Reading data for " <> shape <> " from " <> dataFileName];

	faces = Get[dataFileName];

	w0 = Root[-27 - 756*#1^2 - 2688*#1^4 + 448*#1^6 & , 2, 0];
	r0 = Sqrt[3/(2 - 2*Sin[(3*Pi)/14])];
	
	seedNewFace1 = N[Table[{w0, r0 Cos[k 2 Pi / 7], r0 Sin[k 2 Pi / 7]}, {k, 0, 6}]];
	seedNewFace2 = N[Table[{w0, r0 Cos[k 4 Pi / 7], r0 Sin[k 4 Pi / 7]}, {k, 0, 6}]];

	graph = plotShape3D[faces, seedNewFace1, seedNewFace2];

	outputFileName = "facet_" <> IntegerString[level] <> "_hyperboloid.png";

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
]
