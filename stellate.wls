#!/usr/bin/env wolframscript

exportToPov = False;
highlightCenter = False;
topView = False;
convexHull = False;
plotFaces = False;
highlightCommonFaces = False;
antipodal = False;

exportToPov = True;
highlightCenter = True;
(* topView = True; *)
convexHull = True;
(* plotFaces = True; *)
highlightCommonFaces = True;
antipodal = True;

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

	(* selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &]; *)
	selectedFaces = faces;
	selectedFaces = Select[selectedFaces, Norm[Mean[#][[{2, 3}]]] < 8 &];
	(* norms of 2,3: 0, 2.239, 5.37, 7.6, 13.7... *)
	If[antipodal,
		selectedFaces = Union[selectedFaces, -selectedFaces];
	];
	(* panning = 0; *)
	(* selectedFaces = Map[Boost[panning].# &, selectedFaces, {2}]; *)


	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	stellatedFaces = Map[scaleFace[#, ratio, stellationLevel]&, selectedFaces];

	convexHullOrder = Table[Mod[- 6 / stellationLevel * k, 7]+1, {k, 0, 6}];
	If[convexHull, stellatedFaces = Map[#[[convexHullOrder]]&, stellatedFaces]];

	displayStellatedFaces = Map[getHyperboloid, stellatedFaces, {2}];
(* 	If[!highlightCenter, 
		displayFaces = displayStellatedFaces;
		selectedFaces = stellatedFaces];
 *)
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
	topVertex = SortBy[topVertices, Norm[# - {3, -3, 3}]&][[1]];

	(* 
		stellation level = 1: each vertex belongs to up to 3 faces
		stellation level = 2: each vertex belongs to up to 7 faces
		stellation level = 3: each vertex belongs to up to 8, 9, 10, (infinite?) faces
	 *)

(* 	facesPerEdge = Table[
		Length[Select[displayStellatedFaces, (Length[Intersection[getEdgesFromFace[#], {displayEdges[[k]]}, SameTest -> sameEdge]] > 0)&]], {k, 1, Length[displayEdges]}
	];
	Print["Each edge belongs to up to " <> IntegerString[Max[facesPerEdge]] <> " faces."];
 *)
	(* 
	facesPerEdge can be 2 for the {7/3} polyhedron. Appears to be always 1 for its convex hull. If we include -faces, {7/3} still has 2 faces per edge. The convex hull also has 2 faces per edge
	*)

	centeredFaces = Select[stellatedFaces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	
	displayCenteredFaces = If[highlightCommonFaces,
		Select[displayStellatedFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] > 0)&],
		Map[getHyperboloid, centeredFaces, {2}]
	];

	displayCenteredFaces = Map[(# + {0.001, 0, 0})&, displayCenteredFaces, {2}];

	displayCenteredEdges = Flatten[Map[getEdgesFromFace, displayCenteredFaces], 1];

	(* splitCenteredEdges = Map[HNormalize[#, norm]&, splitCenteredEdges, {2}]; *)
	Print["There are " <> IntegerString[Length[displayCenteredEdges]] <> " centered edges after split."];
	(* displayCenteredEdges = Map[getHyperboloid, centeredEdges, {2}]; *)

	(* centeredVertices = displayCenteredFaces[[1]]; *)

	rangeZoom = stellationLevel * If[highlightCommonFaces, 1.5, 1];
	lighting = {{"Point", White, {10, -10, 10}}};
	range = If[topView, 3, 3] * rangeZoom;
	thickness = range / 80.0 / rangeZoom;
	viewCenter = {-5, 0, 0};

	elements = 
		{
			{
				Opacity[If[highlightCommonFaces, 0.25, 1]], FaceForm[If[plotFaces, White, color]], Map[Cylinder[#, thickness * (1 + 0.00001/(1+ Abs[#[[1]][[1]] - #[[2]][[1]]] )) ]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[Blue], Map[Sphere[#, thickness * 4]&, vertices]
			}
		};

	If[highlightCenter || highlightCommonFaces,
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
	graph
];


shapeIndex = 1;

(* For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++, *)
shape = shapes[[shapeIndex]];

Print[shape];
Print[stellationLevel];

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
Print["Reading data for " <> shape <> " from " <> dataFileName];

faces = Get[dataFileName];

graph = plotShape3D[faces];

outputFileName = If[antipodal, "antipodal_", ""] <> If[highlightCenter, "highlightCenter_", ""] <> If[highlightCommonFaces, "highlightCommonFaces_", ""] <> If[topView, "topView_", ""] <> If[convexHull, "convexHull_", ""] <> If[!plotFaces, "skeleton_", ""] <> IntegerString[stellationLevel] <> "_hyperboloid.png";

(* imageSize = {600, 450}; *)
imageSize = {800, 600};

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

(* 
stellationLevel == 3: 455 vertices. 

Select[vertices, #[[1]]<0&]//Length
322

Select[vertices, #[[1]]>0&]//Length
133

Union[vertices, -vertices, SameTest -> ApproxSamePoint]
644

Union[Select[vertices, #[[1]]<0&], -Select[vertices, #[[1]]>0&], SameTest -> ApproxSamePoint]//Length

322

133 vertices with #[[1]] > 0 are antipodes of a subset of vertices with #[[1]] < 0
This is true even if we boost the shape. Boost keeps the antipodal relation.

stellationLevel == 3: 1281 edges.

epsilon = 0.000001;
Select[displayEdges, Abs[Mean[#][[1]]] < epsilon & ]//Length
42

Select[displayEdges, Mean[#][[1]] < -epsilon & ]//Length
847

Select[displayEdges, Mean[#][[1]] > epsilon & ]//Length
392

pnunion = Union[displayEdges, -displayEdges, SameTest -> sameEdge];
pnunion//Length
1736

Select[pnunion, Abs[Mean[#][[1]]] < epsilon & ]//Length
42

Select[pnunion, Mean[#][[1]] < -epsilon & ]//Length
847

Select[pnunion, Mean[#][[1]] > epsilon & ]//Length
847

Need to include edges around zero from both side, to be complete
negaUnion = Union[Select[displayEdges, Mean[#][[1]] < epsilon & ], -Select[displayEdges, Mean[#][[1]] > -epsilon & ], SameTest -> sameEdge];
negaUnion//Length
889

Select[negaUnion, Abs[Mean[#][[1]]] < epsilon & ]//Length
42

Select[negaUnion, Mean[#][[1]] < -epsilon & ]//Length
847

Select[negaUnion, Mean[#][[1]] > epsilon & ]//Length
0

So, the edges strictly positive is antipodal to some edges strictly negative. Some edges around zero.
*)