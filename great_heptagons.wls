#!/usr/bin/env wolframscript

exportToPov = False;
topView = False;
plotFaces = False;
plotEdges = False;
plotVertices = False;
highlightCommonFaces = False;
vertexFigure = False;
vertexFigureTopology = False;


exportToPov = True;
(* topView = True; *)
plotVertices = True;
(* plotFaces = True; *)
plotEdges = True;
highlightCommonFaces = True;
(* vertexFigure = True; *)
(* vertexFigureTopology = True; *)

(* 
	8: vertex figure has 9 lines. 
	58: 17 lines. All data from the precomputed data file 
*)
truncateThreshold = 8;
(* truncateThreshold = 58; *)


(* vfRatio = 1/2: connecting midpoints. = 1: connecting neighboring points of topVertex *)
vfRatio = 1/2;
vfRatio = 1;

(* always use convexhull *)
convexHull = True;

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
colorsConvexHull = {Red, Green, Orange};
dataFolder = "data";
epsilon = 0.00000001;
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

Boost[t_] :=
{
  {Cosh[t], Sinh[t], 0},
  {Sinh[t], Cosh[t], 0},
  {0,       0,       1}
};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001;
sameEdge[e1_, e2_] := sameCenter[e1, e2] && ((ApproxSamePoint[e1[[1]], e2[[1]]] && ApproxSamePoint[e1[[2]], e2[[2]]] ) || (ApproxSamePoint[e1[[1]], e2[[2]]] && ApproxSamePoint[e1[[2]], e2[[1]]] ));
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
	color = 
		If[convexHull,
			colorsConvexHull[[stellationLevel]],
			colors[[stellationLevel]]
		];

	selectedFaces = faces;

	selectedFaces = Select[selectedFaces, Abs[Mean[#][[1]]] < truncateThreshold &];
	selectedFaces = Union[selectedFaces, -selectedFaces];
	(* panning = 0; *)
	(* selectedFaces = Map[Boost[panning].# &, selectedFaces, {2}]; *)


	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	stellatedFaces = Map[scaleFace[#, ratio, stellationLevel]&, selectedFaces];

	convexHullOrder = Table[Mod[- 6 / stellationLevel * k, 7]+1, {k, 0, 6}];
	If[convexHull, stellatedFaces = Map[#[[convexHullOrder]]&, stellatedFaces]];

	displayStellatedFaces = Map[getHyperboloid, stellatedFaces, {2}];
	displayFaces = displayStellatedFaces;
	selectedFaces = stellatedFaces;

	Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces. Contributing " <> IntegerString[Length[selectedFaces] * 7 ] <> " edges before dedup."];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges, SameTest -> sameEdge];

	(* normSquares = Map[HNormSquare, Flatten[edges, 1]] is always constant, no matter what ratio is *)

	Print["There are " <> IntegerString[Length[edges]] <> " distinct edges in total."];
	(* norm = HNorm[edges[[1]][[1]]]; *)

	(* splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1]; *)
	(* splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."]; *)
	displayEdges = Map[getHyperboloid, edges, {2}];

	(* allVertices = Flatten[displayFaces,1]; *)
	allVertices = Flatten[displayEdges, 1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
(* 	Print["Tallying the z value of vertices"];
	Print[Sort[Tally[Map[#[[1]]&, vertices], Abs[#1 - #2]<0.00001&], (#1[[1]] > #2[[1]])& ] ]; *)

	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	facesPerVertex = Table[
		Length[Select[displayStellatedFaces, (Length[Intersection[#, {vertices[[k]]}, SameTest -> ApproxSamePoint]] > 0)&]], {k, 1, Length[vertices]}
	];
	Print["Each vertex belongs to up to " <> IntegerString[Max[facesPerVertex]] <> " faces."];
	topVertexIndices = Select[Range[Length[facesPerVertex]], facesPerVertex[[#]]==Max[facesPerVertex]&];

	topVertices = vertices[[topVertexIndices]];
	topVertex = SortBy[topVertices, Norm[# - {-3, -3, 3}]&][[1]];
	Print["topVertex"];
	Print[topVertex];

(* 	facesPerEdge = Table[
		Length[Select[displayStellatedFaces, (Length[Intersection[getEdgesFromFace[#], {displayEdges[[k]]}, SameTest -> sameEdge]] > 0)&]], {k, 1, Length[displayEdges]}
	];

	Print["Each edge belongs to up to " <> IntegerString[Max[facesPerEdge]] <> " faces."];
	Print[Tally[facesPerEdge]];
 *)	
	(* 
	Each edge belongs to up to 2 faces. tally: {{1, 910}, {2, 826}}
	*)

	centeredFaces = Select[stellatedFaces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	
	displayCenteredFaces = Select[displayStellatedFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] > 0)&];

	displayCenteredEdges = Flatten[Map[getEdgesFromFace, displayCenteredFaces], 1];

	truncatedEdges = {};

	For[fid = 1, fid <= Length[displayCenteredFaces], fid++,
		faceToTruncate = displayCenteredFaces[[fid]];
		vid = Select[Range[Length[faceToTruncate]], ApproxSamePoint[faceToTruncate[[#]], topVertex]&][[1]];
		faceToTruncate = RotateLeft[faceToTruncate, vid - 1];
		AppendTo[truncatedEdges, 
			{
				faceToTruncate[[1]] (1 - vfRatio) + faceToTruncate[[ Length[faceToTruncate] ]] vfRatio,
				faceToTruncate[[1]] (1 - vfRatio) + faceToTruncate[[2]] vfRatio
			}
		];
	];

	Print["Vertex figure has " <> IntegerString[Length[truncatedEdges]] <> " line segments"];
	verticesInVf = Flatten[truncatedEdges, 1];
	tallyVertices = Tally[verticesInVf, ApproxSamePoint[#1, #2]&];
	Print["tally of tallyVertices"];
	Print[Tally[Map[#[[2]]&, tallyVertices]]];

	roundedVertexFigure = Round[truncatedEdges, epsilon];
	vfTopologyGraph = GraphPlot3D[
		Map[(#[[1]] -> #[[2]])&, roundedVertexFigure],
		EdgeShapeFunction -> ({FaceForm[Red], Cylinder[#1, 0.03]} &), 
		VertexShapeFunction -> ({FaceForm[Blue], Sphere[#, 0.06]} &),
		(* GraphLayout -> "SpringElectricalEmbedding", *)
		GraphLayout -> "CircularEmbedding",
		ViewPoint -> {0, 0, 10}
	];

	rangeZoom = stellationLevel * 1.5;
	lighting = {{"Point", White, {50, -50, 50}}};
	range = 3 * rangeZoom;
	thickness = range / 80.0 / rangeZoom;
	viewCenter = {-2, 0, 0};

	elements = {};
	If[plotVertices,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[Blue], Map[Sphere[#, thickness * 4]&, vertices]
				}
			}
		];	
	];

	If[plotEdges,
		elements = Union[
			elements,
			{
				{
					Opacity[If[highlightCommonFaces, 1, 1]], FaceForm[If[plotFaces, White, color]], Map[Cylinder[#, thickness * (1 + 0.00001/(1+ Abs[#[[1]][[1]] - #[[2]][[1]]] )) ]&, displayEdges]
				}
			}
		];
	];

	If[highlightCommonFaces,
		highlightColor = White;
		elements = Union[
			elements,
			{
				{
					Opacity[If[plotFaces, 1, 0]], FaceForm[highlightColor], Map[Polygon, displayCenteredFaces]
				},
				{
					Opacity[1], FaceForm[highlightColor], Map[Cylinder[#, thickness * 2.2]&, displayCenteredEdges]
				}
			}
		];
	];

	If[vertexFigure,
		vfColor = Red;
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[vfColor], Map[Cylinder[#, thickness * 2.3]&, truncatedEdges]
				}
			}
		];
	];

	If[plotFaces,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[color], Map[Polygon, displayFaces]
				}
			}
		];
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> If[topView, {1000, 0, 0}, {3, -3, 3}],
		ViewVertical -> If[topView, {0, 1, 0}, {1, 0, 0}],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];

	If[vertexFigureTopology, vfTopologyGraph, graph]
];

shape = "h_7_3";

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
Print["Reading data for " <> shape <> " from " <> dataFileName];

faces = Get[dataFileName];

graph = plotShape3D[faces];

outputFileName = If[vertexFigure, "vertexFigure_", ""] <> If[vertexFigureTopology, "topology_", ""]<> If[highlightCommonFaces, "highlightCommonFaces_", ""] <> If[topView, "topView_", ""] <> If[plotFaces, "f", ""] <> If[plotEdges, "e", ""] <> If[plotVertices, "v", ""] <> If[plotFaces || plotEdges || plotVertices, "_", ""] <> "great_heptagons_" <> IntegerString[truncateThreshold] <> ".png";

(* imageSize = {600, 450}; *)
imageSize = {4, 3}/3 * 720;

imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, "great_heptagons"}];
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

Print["Saved image to " <> outputFileName];

