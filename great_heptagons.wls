#!/usr/bin/env wolframscript

exportToPov = False;
plotFaces = False;
plotEdges = False;
plotVertices = False;
highlightCommonFaces = False;
vertexFigure = False;
vertexFigureTopology = False;
edgeTopology = False;

exportToPov = True;
(* plotVertices = True; *)
plotFaces = True;
(* plotEdges = True; *)
(* highlightCommonFaces = True; *)
(* vertexFigure = True; *)
(* vertexFigureTopology = True; *)
(* edgeTopology = True; *)


(* viewAngle = "top"; *)
(* viewAngle = "side"; *)
viewAngle = "";

(* 
	4: 5 lines
	8: vertex figure has 9 lines. 
	28: 13 lines
	58: 17 lines. All data from the precomputed data file 
*)
(* truncateThreshold = 2; *)
(* truncateThreshold = 3; *)
(* truncateThreshold = 4; *)
(* truncateThreshold = 6; *)
truncateThreshold = 8;
(* truncateThreshold = 18; *)
(* truncateThreshold = 28; *)
(* truncateThreshold = 38; *)
(* truncateThreshold = 58; *)


(* vfRatio = 1/2: connecting midpoints. = 1: connecting neighboring points of topVertex *)
(* vfRatio = 1/4; *)
vfRatio = 1/2;
(* vfRatio = 1; *)

(* always use convexhull *)
convexHull = True;

If[edgeTopology, vertexFigureTopology = False];

stellationLevel = 3;
ratio = Tan[stellationLevel Pi/7]/Tan[Pi/7];

splitEdgeParts = 8;
dataFolder = "data";
epsilon = 0.00000001;
Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

RoundEpsilon[v_] := Round[v, epsilon];
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
getColor[center_] := Hue[ Mod[Arg[ center[[2]] + I center[[3]] ], 2 Pi]/2/Pi, (Tanh[center[[1]]/10 + 0.5 ] + 1)/2 , 1];


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

	selectedFaces = faces;

	selectedFaces = Select[selectedFaces, Abs[Mean[#][[1]]] < truncateThreshold &];
	selectedFaces = Union[selectedFaces, -selectedFaces];

(* 	heights = Map[Abs[Mean[#][[1]]]&, selectedFaces];
	heightsTally = Sort[Tally[heights, ApproxSamePoint[#1, #2]&]];
	Print[heightsTally];
 *)
(* 
for truncateThreshold = 28, we have 170 faces, in 8 levels

{{1.812922503760635, 2}, {3.0021462851076754, 14}, {5.674307866050231, 14}, {7.817215073975985, 28}, {13.821507644191334, 28}, {16.493669225133896, 28}, {21.308738014002206, 28}, {22.49796179534925, 28}}

 *)
	(* panning = 0; *)
	(* selectedFaces = Map[Boost[panning].# &, selectedFaces, {2}]; *)

	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	stellatedFaces = Map[scaleFace[#, ratio, stellationLevel]&, selectedFaces];

	convexHullOrder = Table[Mod[- 6 / stellationLevel * k, 7]+1, {k, 0, 6}];
	stellatedFaces = Map[#[[convexHullOrder]]&, stellatedFaces];

	displayStellatedFaces = Map[getHyperboloid, stellatedFaces, {2}];
	displayFaces = displayStellatedFaces;
	selectedFaces = stellatedFaces;

	Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces. Contributing " <> IntegerString[Length[selectedFaces] * 7 ] <> " edges before dedup."];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges, SameTest -> sameEdge];

	Print["There are " <> IntegerString[Length[edges]] <> " distinct edges in total."];
	displayEdges = Map[getHyperboloid, edges, {2}];
	allVertices = Flatten[displayEdges, 1];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];

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
(* 	singleEdges = {};
	For[eid = 1, eid <= Length[displayEdges], eid++,
		If[facesPerEdge[[eid]] == 1,
			AppendTo[singleEdges, displayEdges[[eid]]];
		];
	];

	singleEdgeLengths = Map[Norm[#[[1]] - #[[2]]]&, singleEdges];
	singleEdgeCenters = Map[Mean, singleEdges];
	singleEdgeEnds = Flatten[singleEdges, 1];

	Print["min length, min center norm, min center z, min ends norm, min ends z"];
	Print[Min[singleEdgeLengths]];
	Print[Min[Map[Norm, singleEdgeCenters]]];
	Print[Min[Map[Abs[#[[1]]]&, singleEdgeCenters]]];
	Print[Min[Map[Norm, singleEdgeEnds]]];
	Print[Min[Map[Abs[#[[1]]]&, singleEdgeEnds]]];

	facesWithSingleEdges = {};
	For[fid = 1, fid<=Length[displayFaces], fid++,
		edgesInFace = getEdgesFromFace[displayFaces[[fid]]];
		If[Length[Intersection[edgesInFace, singleEdges, SameTest -> sameEdge]] > 0,
			AppendTo[facesWithSingleEdges, displayFaces[[fid]]];
		];
	];

	Print["faces with single edges"];
	Print[Length[facesWithSingleEdges]];
	heights = Map[Abs[Mean[#][[1]]]&, facesWithSingleEdges];
	heightsTally = Sort[Tally[heights, ApproxSamePoint[#1, #2]&]];
	Print[heightsTally];

	Exit[]; *)

	displayCenteredFaces = Select[displayStellatedFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] > 0)&];
	otherFaces = Select[displayStellatedFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] == 0)&];
	allFacesAfterTruncation = otherFaces;

	displayCenteredEdges = Flatten[Map[getEdgesFromFace, displayCenteredFaces], 1];

	truncatedEdges = {};

	For[fid = 1, fid <= Length[displayCenteredFaces], fid++,
		faceToTruncate = displayCenteredFaces[[fid]];
		vid = Select[Range[Length[faceToTruncate]], ApproxSamePoint[faceToTruncate[[#]], topVertex]&][[1]];
		faceToTruncate = RotateLeft[faceToTruncate, vid - 1];
		truncatedCenteredFace = Join[
			{faceToTruncate[[1]] (1 - vfRatio) + vfRatio faceToTruncate[[ Length[faceToTruncate] ]]},
			{faceToTruncate[[1]] (1 - vfRatio) + vfRatio faceToTruncate[[2]]},
			faceToTruncate[[Range[2, Length[faceToTruncate] ]]]
		];
		AppendTo[allFacesAfterTruncation, truncatedCenteredFace];
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

(* 	For[k=1, k<=Length[truncatedEdges], k++,
		te = truncatedEdges[[k]];
		(* Print[{k, HInner[te[[1]], te[[2]]], HNorm[te[[1]] - te[[2]]]}]; *)
		For[k2=k+1, k2<=Length[truncatedEdges], k2++,
			te2 = truncatedEdges[[k2]];
			int = Intersection[ te, te2, SameTest -> ApproxSamePoint ];
			If[Length[int] == 1,
				comp1 = Complement[te, int, SameTest -> ApproxSamePoint][[1]];
				comp2 = Complement[te2, int, SameTest -> ApproxSamePoint][[1]];
				d12 = comp1 - comp2;
				d11 = comp1 - int[[1]];
				d22 = comp2 - int[[1]];
				normalizedInner1 = HInner[d12, d11]/HNorm[d12]/HNorm[d11];
				Print[{k, 
					k2, 
					HInner[comp1, comp2]//RoundEpsilon, 
					HNorm[comp1 - comp2]//RoundEpsilon,
					HInner[comp1 - int[[1]], comp2 - int[[1]]]//RoundEpsilon,
					normalizedInner1,
					ArcCos[normalizedInner1]//RoundEpsilon,
					(Pi/ArcCos[normalizedInner1])//RoundEpsilon
				}];
			];
		];
	];

	Exit[]; *)

	If[vertexFigure,
		displayFaces = allFacesAfterTruncation;
		displayEdges = Flatten[Map[getEdgesFromFace, displayFaces], 1];
	];

	roundedVertexFigure = RoundEpsilon[truncatedEdges];
	vfTopologyGraph = GraphPlot3D[
		Map[(#[[1]] -> #[[2]])&, roundedVertexFigure],
		EdgeShapeFunction -> ({FaceForm[Red], Cylinder[#1, 0.03]} &), 
		VertexShapeFunction -> ({FaceForm[Blue], Sphere[#, 0.06]} &),
		GraphLayout -> "CircularEmbedding",
		ViewPoint -> {0, 0, 10}
	];

	roundedEdges = RoundEpsilon[displayEdges];
	edgeTopologyGraph = GraphPlot3D[
		Map[(#[[1]] -> #[[2]])&, roundedEdges],		
		EdgeShapeFunction -> ({FaceForm[Red], Cylinder[#1, 0.03]} &), 
		VertexShapeFunction -> ({FaceForm[Blue], Sphere[#, 0.06]} &),
		GraphLayout -> "SpringElectricalEmbedding"
		(* GraphLayout -> "CircularEmbedding", *)
		(* ViewPoint -> {0, 0, 10} *)
	];

	rangeZoom = stellationLevel * 1.5;
	lighting = {{"Point", White, {50, -50, 50}}};
	range = If[viewAngle == "side", 0.8, 3] * rangeZoom;
	thickness = range / 80.0 / rangeZoom;
	viewCenter = {-2, 0, 0};

	elements = {};
	If[plotVertices,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 4]&, vertices]
				}
			}
		];	
	];

	If[plotEdges,
		elements = Union[
			elements,
			Table[
			{
					Opacity[If[highlightCommonFaces, 1, 1]], 
					FaceForm[If[plotFaces, White, getColor[Mean[displayEdges[[k]]]]]], 
					Cylinder[
						displayEdges[[k]], 
						thickness * (1 + 0.00001/(1+ Abs[displayEdges[[k]][[1]][[1]] - displayEdges[[k]][[2]][[1]]] )) 
					]
			},
			{k, 1, Length[displayEdges]}]
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
			Table[{Opacity[1], FaceForm[getColor[Mean[displayFaces[[k]]]]], Polygon[displayFaces[[k]]]}, {k, 1, Length[displayFaces]}]
		];
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> If[viewAngle == "top", {1000, 0, 0}, If[viewAngle == "side", {0, -2.52, 3.165}*100, {50, -50, 50}]],
		ViewVertical -> If[viewAngle == "top", {0, 1, 0}, {1, 0, 0}],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];

	If[vertexFigureTopology, vfTopologyGraph, 
		If[edgeTopology, edgeTopologyGraph, graph]
	]
];

shape = "h_7_3";

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
Print["Reading data for " <> shape <> " from " <> dataFileName];

faces = Get[dataFileName];

graph = plotShape3D[faces];

outputFileName = If[vertexFigure, "vertexFigure_", ""] <> If[vertexFigureTopology, "topology_", ""] <> If[edgeTopology, "edge_topology_", ""] <> If[highlightCommonFaces, "highlightCommonFaces_", ""] <> If[viewAngle == "", "", viewAngle <> "_"] <> If[plotFaces, "f", ""] <> If[plotEdges, "e", ""] <> If[plotVertices, "v", ""] <> If[plotFaces || plotEdges || plotVertices, "_", ""] <> "great_heptagons_" <> IntegerString[truncateThreshold] <> ".png";

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

