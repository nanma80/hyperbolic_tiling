#!/usr/bin/env wolframscript

exportToPov = True;

shapes = {
	"h_7_3",
	"h_3_7",
	"h_72_7",
	"h_7_72"
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

plotShape3D[faces_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];
	color = configs[shape][[3]];

	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getKlein, faces, {2}];

	Print["There are " <> IntegerString[Length[faces]] <> " faces."];

	allEdges = Flatten[Map[getEdgesFromFace, faces], 1];
	edges = Union[allEdges , SameTest -> sameCenter];
	Print["There are " <> IntegerString[Length[edges]] <> " edges in total."];
	norm = HNorm[edges[[1]][[1]]];

	splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1];
	splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
	Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."];
	displayEdges = Map[getKlein, splitEdges, {2}];

	allVertices = Flatten[displayFaces,1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];


	centeredFaces = Select[faces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	displayCenteredFaces = Map[getKlein, centeredFaces, {2}];

	centeredEdges = Flatten[Map[getEdgesFromFace, centeredFaces], 1];
	splitCenteredEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, centeredEdges], 1];
	splitCenteredEdges = Map[HNormalize[#, norm]&, splitCenteredEdges, {2}];
	Print["There are " <> IntegerString[Length[splitCenteredEdges]] <> " centered edges after split."];
	displayCenteredEdges = Map[getKlein, splitCenteredEdges, {2}];

	centeredVertices = displayCenteredFaces[[1]];

	thickness = 0.01;
	lighting = {{"Point", White, {10, -10, 10}}};
	range = 0.8;

	graph = Graphics3D[
		{
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
				Opacity[1], FaceForm[Magenta], Map[Cylinder[#, thickness * 1.1]&, displayCenteredEdges]
			},
			{
				Opacity[1], FaceForm[Magenta], Map[Sphere[#, thickness * 2.2]&, centeredVertices]
			}
		},
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

	graph = plotShape3D[faces];

	outputFileName = "klein.png";
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
