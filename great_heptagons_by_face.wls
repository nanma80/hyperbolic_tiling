#!/usr/bin/env wolframscript

exportToPov = False;
plotFaces = False;
plotEdges = False;
plotVertices = False;
vertexFigure = False;
vertexFigureTopology = False;

frameCountPerLevel = 100;

exportToPov = True;
(* plotVertices = True; *)
plotFaces = True;
(* plotEdges = True; *)
(* vertexFigure = True; *)
(* vertexFigureTopology = True; *)

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
(* truncateThreshold = 8; *)
truncateThreshold = 8;
(* truncateThreshold = 28; *)
(* truncateThreshold = 58; *)


(* vfRatio = 1/2: connecting midpoints. = 1: connecting neighboring points of topVertex *)
(* vfRatio = 1/4; *)
vfRatio = 1/2;
(* vfRatio = 1; *)

(* always use convexhull *)
convexHull = True;

stellationLevel = 3;
ratio = Tan[stellationLevel Pi/7]/Tan[Pi/7];

splitEdgeParts = 8;
dataFolder = "data";
epsilon = 0.00000001;
Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

imageSize = {4, 3}/3 * 720;
imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, "great_heptagons"}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

shapeWithParams = If[vertexFigure, "vertexFigure_", ""] <> If[vertexFigureTopology, "topology_", ""] <> If[viewAngle == "", "", viewAngle <> "_"] <> If[plotFaces, "f", ""] <> If[plotEdges, "e", ""] <> If[plotVertices, "v", ""] <> If[plotFaces || plotEdges || plotVertices, "_", ""] <> "great_heptagons_" <> IntegerString[truncateThreshold];


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
getColor[center_] := Hue[ Mod[Arg[ center[[2]] + I center[[3]] ], 2 Pi]/2/Pi, (Tanh[center[[1]]/10 + 0.5 ] + 1)/2, 1];
(* getColor[center_] := Red; *)

scale[edge_, ratio_] := Module[
	{center},
	center = Mean[edge];
	Map[(# - center) * ratio + center&, edge]
];

scaleFace[face_, ratio_, stellationLevel_] := Module[
	{edges, edgeOrders},
	edges = Map[scale[#, ratio]&, getEdgesFromFace[face]];
	edgeOrders = Table[Mod[stellationLevel * k, 7]+1, {k, 0, 6}];
	Map[#[[1]]&, edges[[edgeOrders]]]
];

plotShape3D[faces_, heightLevel_, progressInLevel_, frameIndex_] := Module[
	{displayFaces, allEdges, edges, allVertices, vertices},

	viewPointRotate = RotationMatrix[2 Pi/frameCount * frameIndex, {1, 0, 0}];
	displayFaces = faces;
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

	rangeZoom = stellationLevel * 1.5;
	lighting = {{"Point", White, {50, -50, 50}}};
	(* range = If[viewAngle == "side", 0.8, 3] * rangeZoom; *)
	range = If[viewAngle == "side", 0.8, 3] * rangeZoom/2;
	thickness = range / 80.0 / rangeZoom;
	viewCenter = {-2, 0, 0};

	elements = {};
	If[plotVertices,
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 4]&, vertices.viewPointRotate]
				}
			}
		];	
	];

	If[plotEdges,
		elements = Union[
			elements,
			Table[
			{
					Opacity[1],
					FaceForm[If[plotFaces, White, getColor[Mean[displayEdges[[k]]]]]], 
					Cylinder[
						displayEdges[[k]].viewPointRotate, 
						thickness * (1 + 0.00001/(1+ Abs[displayEdges[[k]][[1]][[1]] - displayEdges[[k]][[2]][[1]]] )) 
					]
			},
			{k, 1, Length[displayEdges]}]
		];
	];


	If[vertexFigure,
		vfColor = Red;
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[vfColor], Map[Cylinder[#, thickness * 2.3]&, truncatedEdges.viewPointRotate]
				}
			}
		];
	];

	If[plotFaces,
		faceElements = {};
		For[fid = 1, fid <= Length[displayFaces], fid++,
			displayFace = displayFaces[[fid]];
			faceCenter = Mean[displayFace];
			height = RoundEpsilon[faceCenter[[1]]];
			faceLevel = getLevelByHeight[height];
			faceCountSameHeight = getFaceCountByHeight[height];
			normalizedArg = Mod[Arg[ faceCenter[[2]] - I faceCenter[[3]] ] + 2 Pi/7, 2 Pi]/2/Pi * 0.5;

			opacity = If[
				faceLevel < heightLevel,
				1,
				If[
					faceLevel > heightLevel,
					0,
					If[progressInLevel < normalizedArg, 
						0,
						If[progressInLevel < normalizedArg + 0.3, 
							(progressInLevel - normalizedArg)/0.3,
							1
						]
					]
				]
			];

			If[opacity > 0,
				color = getColor[faceCenter];
				AppendTo[faceElements, 
					{
						Opacity[opacity], 
						FaceForm[color],
						Polygon[displayFace.viewPointRotate]
					}
				];
			];
		];

		elements = Union[
			elements,
			faceElements
		];
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> If[viewAngle == "top", {1000, 0, 0}, If[viewAngle == "side", {0, -2.52, 3.165}*1000, {30, -50, 50}*80]],
		ViewVertical -> If[viewAngle == "top", {0, 1, 0}, {1, 0, 0}],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];

	If[vertexFigureTopology, vfTopologyGraph, graph]
];

exportFrame[faces_, heightLevel_, progressInLevel_, frameIndex_, outputFolderName_] := Module[
	{outputFrameFolder},
	Print[DateString["ISODateTime"] <> " Generating frame for:"];
	Print[{frameIndex, heightLevel, progressInLevel}];
	(* Return[]; *)
	(* If[frameIndex <= 175, Return[]]; *)
	(* If[Mod[frameIndex, 2] != frameIndexMod, Return[]]; *)
	outputFrameFolder = FileNameJoin[{outputFolder, outputFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	outputFileName = "frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFrameFolder, outputFileName}];

	graph = plotShape3D[faces, heightLevel, progressInLevel, frameIndex];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];

	Print["Saved image to " <> outputFileName];
];

exportAnimation[] := Module[
	{
		gifFileName,
		gifFolder,
		frameIndex,
		explosionFactor,
		gifCommand
	},

	outputFrameFolderName = shapeWithParams <> "_frames";

	outputFrameFolder = FileNameJoin[{outputFolder, outputFrameFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	gifFileName = shapeWithParams;
	gifFolder = outputFolder;
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	shape = "h_7_3";
	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
	faces = Get[dataFileName];
	selectedFaces = Select[faces, Abs[Mean[#][[1]]] < truncateThreshold &];
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

	frameCount = (Length[sortedHeights] - 2) * frameCountPerLevel;

	For[levelId = 3, levelId <= Length[sortedHeights], levelId++,
		For[progressId = 1, progressId <= frameCountPerLevel, progressId++,
			frameIndex = (levelId - 3) * frameCountPerLevel + progressId;
			Print["exportFrame"];
			Print[{levelId, progressId, N[progressId/frameCountPerLevel], frameIndex, frameCount}];

			exportFrame[faces, levelId, progressId/frameCountPerLevel, frameIndex, outputFrameFolderName];
		];
	];
	
	delay = If[frameCount >= 50, 3, Floor[500/frameCount]];

	(* gifCommand = "magick convert -duplicate 1,-2-1 -quiet +repage -alpha remove -loop 0 -delay " *)
	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " 
		<> IntegerString[delay] <> " -dispose 2 " <> outputFrameFolder <> "\\frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

	Print[gifCommand];
	Run[gifCommand];
	Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];
	
  mp4Command = "ffmpeg -i " <> gifFolder <> "\\" <> gifFileName <> ".gif" 
      <> " -movflags faststart -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" "
      <> gifFolder <> "\\" <> gifFileName <> ".mp4";

  mp4FileName = gifFolder <> "\\" <> gifFileName <> ".mp4";
  If[FileExistsQ[mp4FileName], DeleteFile[mp4FileName]];
  Print[mp4Command];
  Run[mp4Command];
  Print["Saved the video to " <> gifFolder <> "\\" <> gifFileName <> ".mp4"];
];


exportAnimation[];

(* exportFrame[faces, 1, 0, 1, "test"];
exportFrame[faces, 2, 0, 2, "test"];
exportFrame[faces, 3, 0, 3, "test"];
exportFrame[faces, 4, 0, 4, "test"];
exportFrame[faces, 5, 0, 5, "test"];
exportFrame[faces, 6, 0, 6, "test"];
exportFrame[faces, 7, 0, 7, "test"];
exportFrame[faces, 8, 0, 8, "test"];
 *)

(* shape = "h_7_3";
dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
faces = Get[dataFileName];
selectedFaces = Select[faces, Abs[Mean[#][[1]]] < truncateThreshold &];
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
Table[exportFrame[faces, 5, k, 200 + 100 k, "test"], {k, 0/100, 5/100, 1/100}];
 *)
