#!/usr/bin/env wolframscript

exportToPov = True;

shapes = {
	(* "h_7_3" *)
	"h_3_7"
	(* "h_72_7", *)
	(* "h_7_72" *)
};

splitEdgeParts = 8;
level = 2; (* 1 or 2 *)

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

plotShape3D[faces_, highlightedFace1_, highlightedFace2_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];
	color = configs[shape][[3]];

	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getPoincare, faces, {2}];

	Print["There are " <> IntegerString[Length[faces]] <> " faces."];

	allEdges = Flatten[Map[getEdgesFromFace, faces], 1];
	edges = Union[allEdges , SameTest -> sameCenter];
	Print["There are " <> IntegerString[Length[edges]] <> " edges in total."];
	norm = HNorm[edges[[1]][[1]]];

	splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1];
	splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
	Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."];
	displayEdges = Map[getPoincare, splitEdges, {2}];

	allVertices = Flatten[displayFaces,1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];


	centeredFaces1 = {highlightedFace1};
	displayCenteredFaces1 = Map[getPoincare, centeredFaces1, {2}];

	centeredEdges1 = Flatten[Map[getEdgesFromFace, centeredFaces1], 1];
	splitCenteredEdges1 = Flatten[ Map[splitEdge[#, splitEdgeParts]&, centeredEdges1], 1];
	splitCenteredEdges1 = Map[HNormalize[#, norm]&, splitCenteredEdges1, {2}];
	Print["There are " <> IntegerString[Length[splitCenteredEdges1]] <> " centered edges after split."];
	displayCenteredEdges1 = Map[getPoincare, splitCenteredEdges1, {2}];

	centeredVertices1 = displayCenteredFaces1[[1]];


	centeredFaces2 = {highlightedFace2};
	displayCenteredFaces2 = Map[getPoincare, centeredFaces2, {2}];

	centeredEdges2 = Flatten[Map[getEdgesFromFace, centeredFaces2], 1];
	splitCenteredEdges2 = Flatten[ Map[splitEdge[#, splitEdgeParts]&, centeredEdges2], 1];
	splitCenteredEdges2 = Map[HNormalize[#, norm]&, splitCenteredEdges2, {2}];
	Print["There are " <> IntegerString[Length[splitCenteredEdges2]] <> " centered edges after split."];
	displayCenteredEdges2 = Map[getPoincare, splitCenteredEdges2, {2}];

	centeredVertices2 = displayCenteredFaces2[[1]];


	thickness = 0.01;
	lighting = {{"Point", White, {10, -10, 10}}};
	range = 0.8;

	elements = {
			{
				Opacity[1], FaceForm[color], Map[Polygon, displayFaces]
			},
			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 2]&, vertices]
			},

			{
				Opacity[1], FaceForm[Yellow], Map[Cylinder[#, thickness * 1.2]&, displayCenteredEdges2]
			},
			{
				Opacity[1], FaceForm[Yellow], Map[Sphere[#, thickness * 2.3]&, centeredVertices2]
			}
		};

	If[level==2,
		elements = Join[elements,

		{

			{
				Opacity[1], FaceForm[Blue], Map[Cylinder[#, thickness * 1.1]&, displayCenteredEdges1]
			},
			{
				Opacity[1], FaceForm[Blue], Map[Sphere[#, thickness * 2.2]&, centeredVertices1]
			}

		}
		]
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> {100, 0, 0},
		ViewVertical -> {0, 0, 1},
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
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
	seedNewFace1 = N[Table[{w0, r0 Cos[k 4 Pi / 7], r0 Sin[k 4 Pi / 7]}, {k, 0, 6}]];
	seedNewFace2 = N[Table[{w0, r0 Cos[k 2 Pi / 7], r0 Sin[k 2 Pi / 7]}, {k, 0, 6}]];

	faces = Join[faces, { seedNewFace1, seedNewFace2 }];

	graph = plotShape3D[faces, seedNewFace1, seedNewFace2];

	outputFileName = "facet_" <> IntegerString[level] <> ".png";
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
