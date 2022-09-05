#!/usr/bin/env wolframscript

exportToPov = False;
plotFaces = False;
plotEdges = False;
plotVertices = False;
highlightCommonFaces = False;
vertexFigure = False;
vertexFigureTopology = False;
edgeTopology = False;
showHyperboloid = False;

exportToPov = True;
plotVertices = True;
plotFaces = True;
plotEdges = True;
(* highlightCommonFaces = True; *)
(* vertexFigure = True; *)
(* vertexFigureTopology = True; *)
(* edgeTopology = True; *)
showHyperboloid = True;

(* viewAngle = "top"; *)
(* viewAngle = "side"; *)
viewAngle = "";

(* 
3: 14 faces
4: 28 faces
6: 42 faces
18: 112 faces
 *)
truncateThreshold = 3;
(* truncateThreshold = 4; *)
(* truncateThreshold = 6; *)
(* truncateThreshold = 8; *)
(* truncateThreshold = 18; *)
(* truncateThreshold = 28; *)
(* truncateThreshold = 38; *)
(* truncateThreshold = 58; *)


(* vfRatio = 1/2: connecting midpoints. = 1: connecting neighboring points of topVertex *)
(* vfRatio = 1/4; *)
(* vfRatio = 1/2; *)
vfRatio = 1;

(* for reconstruction of faces *)
nGenerations = 5;

(* always use convexhull *)
convexHull = True;

If[edgeTopology, vertexFigureTopology = False];

stellationLevel = 3;
ratio = Tan[stellationLevel Pi/7]/Tan[Pi/7];
rangeZoom = stellationLevel * 1.5;
lighting = {{"Point", White, {50, -50, 50}}};
range = If[viewAngle == "side", 0.8, 4.5] * rangeZoom;

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

(* order vertices in a face in a loop *)
orderFace[f_] := Module[
	{vfCenter, verticesInVf, proj, p2d, args},
	vfCenter = f[[1]][[1]];
	verticesInVf = f[[2]];
	proj = Orthogonalize[{{1, 2, 3}, {1.23, 4.34, 0}}];
	p2d = Map[(# - vfCenter) . Transpose[proj] &, verticesInVf];
	args = Map[Arg[#[[1]] + #[[2]] I] &, p2d];
	{ {vfCenter}, verticesInVf[[Ordering[args]]]}
];

faceOverlap[f1_,f2_] := ApproxSamePoint[f1[[1]][[1]], f2[[1]][[1]]];
mergeFace[f1_,f2_] := {f1[[1]], Union[f1[[2]], f2[[2]], SameTest -> ApproxSamePoint]};

mergeToFaces[f1sOriginal_, f2s_] := Module[
	{fid, f2},
	f1s = f1sOriginal;
	For[f2id=1,f2id<=Length[f2s],f2id++,
		f2 = f2s[[f2id]];
		overlappingWithExisting = False;
		For[f1id=1,f1id<=Length[f1s],f1id++,
			f1 = f1s[[f1id]];
			If[faceOverlap[f1, f2],
				f1s[[f1id]] = mergeFace[f1, f2];
				overlappingWithExisting = True;
			];
		];
		If[overlappingWithExisting == False,
			AppendTo[f1s, f2];
		];
	];
	f1s
];

limitFace[f_, zRange_] := Module[
	{orderedFace, vfCenter, vertices},
	orderedFace = orderFace[f];
	vfCenter = orderedFace[[1]];
	vertices = orderedFace[[2]];

	inRangeVertices = Select[vertices, Abs[#[[1]]] < zRange& ];
	edges = getEdgesFromFace[vertices];
	For[eid=1,eid<=Length[edges],eid++,
		edge = edges[[eid]];
		If[ (Abs[edge[[1]][[1]]] - zRange) * (Abs[edge[[2]][[1]]] - zRange) < 0,
			solutions = Solve[Abs[lambda edge[[1]][[1]] + (1 - lambda) edge[[2]][[1]]] == zRange && lambda > 0 && lambda < 1, {lambda}];
			(* Print[Length[solutions]]; *)
			newVertices = (lambda edge[[1]] + (1 - lambda) edge[[2]])/.solutions;
			inRangeVertices = Join[inRangeVertices, newVertices];
		];
	];

	{vfCenter, inRangeVertices }
];

plotShape3D[faces_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];

	faceCenters = Map[#[[1]][[1]]&, faces];
	displayFaces = Map[#[[2]]&, faces];
	Print["There are " <> IntegerString[Length[displayFaces]] <> " faces"];
	displayEdges = Flatten[Map[getEdgesFromFace, displayFaces], 1];
	Print["There are " <> IntegerString[Length[displayEdges]] <> " edges in total."];
	vertices = Flatten[displayEdges, 1];
	
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	facesPerVertex = Table[
		Length[Select[displayFaces, (Length[Intersection[#, {vertices[[k]]}, SameTest -> ApproxSamePoint]] > 0)&]], {k, 1, Length[vertices]}
	];
	Print["Each vertex belongs to up to " <> IntegerString[Max[facesPerVertex]] <> " faces."];
	topVertexIndices = Select[Range[Length[facesPerVertex]], facesPerVertex[[#]]==Max[facesPerVertex]&];

	topVertices = vertices[[topVertexIndices]];
	topVertex = SortBy[topVertices, Norm[# - {-3, -3, 3}]&][[1]];
	Print["topVertex"];
	Print[topVertex];
	topVertexNorm = RoundEpsilon[HNorm[topVertex]];


	(* faces per edge analysis. Very slow  *)
	facesPerEdge = Table[
		Length[Select[displayFaces, (Length[Intersection[getEdgesFromFace[#], {displayEdges[[k]]}, SameTest -> sameEdge]] > 0)&]], {k, 1, Length[displayEdges]}
	];

	Print["Each edge belongs to up to " <> IntegerString[Max[facesPerEdge]] <> " faces."];
	Print[Tally[facesPerEdge]];
	
	(* 
	Each edge belongs to up to 2 faces.
	tally by threshold 
	4: {{2, 140}, {1, 224}}
	6: {{2, 252}, {1, 280}}
	18: {{2, 728}, {1, 392}}
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

	displayCenteredFaces = Select[displayFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] > 0)&];
	otherFaces = Select[displayFaces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] == 0)&];

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

(* issue: centers and these faces don't match *)
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

	thickness = range / 80.0 / rangeZoom;
	viewCenter = {0, 0, 0};

	elements = {};
	If[plotVertices,
		displayVertices = Select[vertices, RoundEpsilon[HNorm[#]] == topVertexNorm& ];
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 4]&, displayVertices]
				}
			}
		];
	];

	If[plotEdges,
		displayEdges = Select[displayEdges, RoundEpsilon[HNorm[#[[1]]]] == topVertexNorm || RoundEpsilon[HNorm[#[[2]]]] == topVertexNorm& ];
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
		splitFaces = {};
		For[fid=1,fid<=Length[displayFaces],fid++,
			splitEdges = getEdgesFromFace[displayFaces[[fid]]];
			AppendTo[splitFaces, Map[{#[[1]], #[[2]], faceCenters[[fid]] }&, splitEdges]];
		];

		opacity = 1;
		If[showHyperboloid, opacity = 0.1];

		elements = Union[
			elements,
			Table[
				{
					Opacity[opacity], 
					FaceForm[getColor[ faceCenters[[k]] ]], 
					Polygon[ splitFaces[[k]] ]
				}, 
			{k, 1, Length[displayFaces]}]
		];
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> If[viewAngle == "top", {1000, 0, 0}, If[viewAngle == "side", {0, -2.52, 3.165}*100, {30, -50, 50}]],
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
faces = Get[dataFileName];
selectedFaces = Select[faces, Abs[Mean[#][[1]]] < truncateThreshold + 2 &];
selectedFaces = Union[selectedFaces, -selectedFaces];
displayFaces = Map[getHyperboloid, selectedFaces, {2}];
stellatedFaces = Map[scaleFace[#, ratio, stellationLevel]&, selectedFaces];
convexHullOrder = Table[Mod[- 6 / stellationLevel * k, 7]+1, {k, 0, 6}];
stellatedFaces = Map[#[[convexHullOrder]]&, stellatedFaces];
faces = Map[getHyperboloid, stellatedFaces, {2}];

heights = Map[RoundEpsilon[Mean[#][[1]]]&, faces];
heightsTally = Tally[heights];
uniqueHeights = Map[#[[1]]&, heightsTally];
sortedHeights = Sort[uniqueHeights, (Abs[#1] < Abs[#2]) || (Abs[#1] == Abs[#2] && #1 < #2) &];
getLevelByHeight = <||>;
getFaceCountByHeight = <||>;
For[hid = 1, hid <= Length[sortedHeights], hid++,
	getLevelByHeight[sortedHeights[[hid]]] = hid;
	getFaceCountByHeight[heightsTally[[hid]][[1]] ] = heightsTally[[hid]][[2]];
];

Print[faces//Length];
Print[sortedHeights];

p = 7; 
q = 3;
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

topVertexFaceId = Select[Range[Length[faces]], Mean[faces[[#]]][[1]] < -1 && Mean[ faces[[#]] ][[1]] > -2 &][[1]];

Print[topVertexFaceId];

edgesByFace = Map[getEdgesFromFace, faces];

faceOverlapOnEdge[id1_, id2_] := Length[Intersection[edgesByFace[[id1]], edgesByFace[[id2]], SameTest -> sameEdge]] == 1;
neighborFaceIds = Select[Range[Length[faces]], faceOverlapOnEdge[#, topVertexFaceId]& ];
Print[neighborFaceIds];

neighborCenter = Mean[faces[[neighborFaceIds[[1]]  ]]];
topFaceCenter = Mean[faces[[topVertexFaceId]]];
newMirror3 = neighborCenter - topFaceCenter;

(* replace old mirror3 by newMirror3 *)
mirrors = {mirror1, mirror2, newMirror3};
Print["mirrors"];
Print[mirrors];

allVertices = Flatten[faces, 1];
vertices = Union[allVertices , SameTest -> ApproxSamePoint];
Print[{Length[allVertices], Length[vertices]}];

facesPerVertex = Table[
	Length[Select[faces, (Length[Intersection[#, {vertices[[k]]}, SameTest -> ApproxSamePoint]] > 0)&]], {k, 1, Length[vertices]}
];
Print["Each vertex belongs to up to " <> IntegerString[Max[facesPerVertex]] <> " faces."];
topVertexIndices = Select[Range[Length[facesPerVertex]], facesPerVertex[[#]]==Max[facesPerVertex]&];

topVertices = vertices[[topVertexIndices]];
topVertex = SortBy[topVertices, Norm[# - {-3, -3, 3}]&][[1]];
Print["topVertex"];
Print[topVertex];

topFaces = Select[faces, (Length[Intersection[#, {topVertex}, SameTest -> ApproxSamePoint]] > 0)&];
truncatedEdges = {};

For[fid = 1, fid <= Length[topFaces], fid++,
	faceToTruncate = topFaces[[fid]];
	vid = Select[Range[Length[faceToTruncate]], ApproxSamePoint[faceToTruncate[[#]], topVertex]&][[1]];
	faceToTruncate = RotateLeft[faceToTruncate, vid - 1];
	truncatedCenteredFace = Join[
		{faceToTruncate[[1]] (1 - vfRatio) + vfRatio faceToTruncate[[ Length[faceToTruncate] ]]},
		{faceToTruncate[[1]] (1 - vfRatio) + vfRatio faceToTruncate[[2]]},
		faceToTruncate[[Range[2, Length[faceToTruncate] ]]]
	];
	AppendTo[truncatedEdges, 
		{
			faceToTruncate[[1]] (1 - vfRatio) + faceToTruncate[[ Length[faceToTruncate] ]] vfRatio,
			faceToTruncate[[1]] (1 - vfRatio) + faceToTruncate[[2]] vfRatio
		}
	];
];

Print["Vertex figure has " <> IntegerString[Length[truncatedEdges]] <> " line segments"];
verticesInVf = Union[Flatten[truncatedEdges, 1], SameTest -> ApproxSamePoint];
Print[Length[verticesInVf]];

faceCenter = topVertex/HInner[topVertex, topVertex]*HInner[topVertex, verticesInVf[[1]]];
reconstructFaces = {{{faceCenter}, verticesInVf}};

countStart = -1;

For[i = 0, i <= nGenerations && countStart != Length[reconstructFaces], i++,
	Print["i = " <> IntegerString[i]];
	countStart = Length[reconstructFaces];
	For[mirrorIndex = Length[mirrors], mirrorIndex >= 1, mirrorIndex--,
		newFaces = Map[HReflect[#, mirrors[[mirrorIndex]] ]&, reconstructFaces, {3}];
		(* reconstructFaces = Union[reconstructFaces, newFaces, SameTest -> faceOverlap]; *)
		reconstructFaces = mergeToFaces[reconstructFaces, newFaces];
	];
];
Print[Length[reconstructFaces]];

invertedFaces = Map[(-#)&, reconstructFaces, {2}];
int = Intersection[reconstructFaces, invertedFaces, SameTest -> faceOverlap];
Print["Length[int]"];
Print[Length[int]];
reconstructFaces = mergeToFaces[reconstructFaces, invertedFaces];

Print["before invert"];
Print[Length[reconstructFaces]];

Print["Making it more symmetric"];

(* Making it symmetric *)
centeredMirrors = Select[mirrors, #[[1]] == 0 &];

countStart = -1;
For[i = 0, i < 20 && countStart != Length[reconstructFaces], i++,
	Print["i = " <> IntegerString[i]];
	countStart = Length[reconstructFaces];
	For[mirrorIndex = Length[centeredMirrors], mirrorIndex >= 1, mirrorIndex--,
		newFaces = Map[HReflect[#, centeredMirrors[[mirrorIndex]] ]&, reconstructFaces, {3}];
		(* reconstructFaces = Union[reconstructFaces, newFaces, SameTest -> faceOverlap]; *)
		reconstructFaces = mergeToFaces[reconstructFaces, newFaces];
	];
];


Print["after merging with invert"];
Print[Length[reconstructFaces]];
faceLengths = Map[Length[#[[2]]]&, reconstructFaces];
Print[Tally[faceLengths]];

newCenterTs = Map[Abs[#[[1]][[1]][[1]] ]&, reconstructFaces];
newCenterTsTally = Tally[newCenterTs, ApproxSamePoint[#1, #2]& ];
Print[newCenterTsTally];

reconstructFaces = Select[reconstructFaces, Abs[#[[1]][[1]][[1]] ] < truncateThreshold &];
Print["after filtering by threshold"];
Print[Length[reconstructFaces]];
faceLengths = Map[Length[#[[2]]]&, reconstructFaces];
Print[Tally[faceLengths]];


Print["limit face"];
reconstructFaces = Map[limitFace[#, range]&, reconstructFaces];

Print["ordering faces"];
reconstructFaces = Map[orderFace, reconstructFaces];
Print[Length[reconstructFaces]];
(* Print[reconstructFaces]; *)

If[showHyperboloid, 
	reconstructFaces = reconstructFaces[[{5}]];
	faceCenter1 = reconstructFaces[[1]][[1]][[1]];
];

graph = plotShape3D[reconstructFaces];

If[showHyperboloid, 
	plotPoints = 20;
	maxRecursion = 2;

	hyperboloidGraphSide1 = ContourPlot3D[
		x^2 - y^2 - z^2 == topVertexNorm^2,
		{x, -range, range}, 
		{y, -1.1 range, 1.1 range}, 
		{z, -1.1 range, 1.1 range},
		Axes -> None,
		BoundaryStyle -> None,
		MaxRecursion -> maxRecursion,
		PlotPoints -> plotPoints,
		MeshFunctions -> {0 &}, 
		Boxed -> False,
		MeshShading -> { { Lighter[Yellow] , Opacity[0.7] } },
		RegionFunction -> Function[{x, y, z}, HInner[faceCenter1, {x,y,z}] < 1.1 HInner[faceCenter1, faceCenter1]],
		RegionBoundaryStyle -> None
	];

	hyperboloidGraphSide2 = ContourPlot3D[
		x^2 - y^2 - z^2 == topVertexNorm^2,
		{x, -range, range}, 
		{y, -1.1 range, 1.1 range}, 
		{z, -1.1 range, 1.1 range},
		Axes -> None,
		BoundaryStyle -> None,
		MaxRecursion -> maxRecursion,
		PlotPoints -> plotPoints,
		MeshFunctions -> {0 &}, 
		Boxed -> False,
		MeshShading -> { { Lighter[Yellow] , Opacity[0.2] } },
		RegionFunction -> Function[{x, y, z}, HInner[faceCenter1, {x,y,z}] > 0.9 HInner[faceCenter1, faceCenter1]],
		RegionBoundaryStyle -> None
	];

	graph = Show[graph, hyperboloidGraphSide1, hyperboloidGraphSide2]
	(* graph = hyperboloidGraph *)
];


outputFileName = If[vertexFigure, "vertexFigure_", ""] <> If[vertexFigureTopology, "topology_", ""] <> If[edgeTopology, "edge_topology_", ""] <> If[highlightCommonFaces, "highlightCommonFaces_", ""] <> If[viewAngle == "", "", viewAngle <> "_"] <> If[plotFaces, "f", ""] <> If[plotEdges, "e", ""] <> If[plotVertices, "v", ""] <> If[plotFaces || plotEdges || plotVertices, "_", ""] <> If[showHyperboloid, "hyperboloid_", ""] <> "great_heptagons_vfch_" <> IntegerString[truncateThreshold] <> ".png";

(* imageSize = {4, 3}/3 * 720; *)
imageSize = {4, 3}/3 * 1080;

imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, "great_heptagons_vfch"}];
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

