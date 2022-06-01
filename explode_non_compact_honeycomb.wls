#!/usr/bin/env wolframscript

maxExplosionFactor = 10.0;
frameCount = 150;

(* maxExplosionFactor = 8.0;
frameCount = 6;
 *)
cellThreshold = 10;
(* cellThreshold = 40; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] < 3,
	Print["Usage: wolframscript.exe explode_honeycomb.wls <p> <q> <r> <cellThreshold>"];
];

If[Length[argv] >= 3, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]]; r = ToExpression[argv[[3]]],
	p = 6; q = 3; r = 3;
];

If[Length[argv] >= 4, 
	cellThreshold = ToExpression[argv[[4]]]
];

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "explode_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm];

dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 1.25;
colors = Join[{Red, Blue, Green, Yellow, Magenta, LightCyan, Brown, Orange, Pink, Darker[Gray, 0.2], Cyan}, RandomColor[80]];

outputFolder = FileNameJoin[{imageFolder, "explode_honeycombs"}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
exportToPov = True;
splitEdgeParts = 15;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

phi = (1 + Sqrt[5])/2;
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
splitFaceToTriangles[face_] := Table[{face[[1]], face[[k-1]], face[[k]]}, {k, 3, Length[face]}];
HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], Sin[t]},
  {0, -Sin[t], Cos[t]}
};

cellCenterToOrigin = ArrayFlatten[{{1, 0}, {0, RotationMatrix[{-{1, 1, 1}, {1, 0, 0}}]}}];
Boost[boost_] := {{Cosh[boost], Sinh[boost], 0, 0}, {Sinh[boost], Cosh[boost], 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
centralReflect[point_, mirror_] := - HReflect[point, mirror];
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];
ApproxSamePointLoose[point1_, point2_] := Round[point1, 10 epsilon] == Round[point2, 10 epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Mean[N[set1]], Mean[N[set2]]];
sameCellCenter[cell1_, cell2_] := sameCenter[Flatten[cell1, 1], Flatten[cell2, 1]];
getCellCenter[cell_] := Total[Flatten[cell, 1]]/Length[Flatten[cell, 1]];
getCellHeight[cell_] := getCellCenter[cell][[1]];
opaqueFactor = 2/3;
getOpacityByExplosion[factor_] := 
  If[factor < maxExplosionFactor*opaqueFactor, 1, 
   1/(1 - opaqueFactor) - 1/(1 - opaqueFactor)/maxExplosionFactor*factor];

getNonCompactCellCenter[cell_] := Module[
	{samplePoints, face1Points, face2id, centerSolution},
	
	face1Points = cell[[1]][[Range[3]]];
	For[face2id = 1, face2id <= Length[cell[[2]]], face2id++,
		samplePoints = Join[face1Points, {cell[[2]][[face2id]]}];
		If[Abs[Det[samplePoints]] > epsilon, 
			Break[]
		];
	];

	If[Abs[Det[samplePoints]] <= epsilon,
		Print["Cannot construct samplePoints"];
		Return[]
	];

	centerSolution = 
		Solve[HInner[{w, x, y, z}, samplePoints[[1]]] == 
		   HInner[{w, x, y, z}, samplePoints[[2]]] && 
		  HInner[{w, x, y, z}, samplePoints[[1]]] == 
		   HInner[{w, x, y, z}, samplePoints[[3]]] && 
		  HInner[{w, x, y, z}, samplePoints[[1]]] == 
		   HInner[{w, x, y, z}, samplePoints[[4]]] && 
		   w == 1, 
		   {w, x, y, z}
	   ];
	{w, x, y, z} /. centerSolution[[1]]
];

getSplitEdges[edges_] := Module[
	{splitEdges, eid, edge, norm, splitParts},
	splitEdges = {};
	For[eid = 1, eid <= Length[edges], eid++,
		edge = edges[[eid]];
		norm = HNorm[edge[[1]]];
		splitParts = splitEdge[edge, splitEdgeParts];
		splitParts = Map[HNormalize[#, norm]&, splitParts, {2}];
		splitEdges = Join[splitEdges, splitParts];
	];
	splitEdges
];

getHyperboloid[v_] := {v[[2]], v[[3]], v[[4]]};
getKlein[v_] := {v[[2]], v[[3]], v[[4]]}/v[[1]];
getPoincare[v_] := {v[[2]], v[[3]], v[[4]]}/(1+v[[1]]);
getPoincareExploded[v_, norm_] := norm {v[[2]], v[[3]], v[[4]]}/(norm+v[[1]]);
explodedCell[cell_, explosionFactor_] := Module[
	{cellCenter},
	cellCenter = getNonCompactCellCenter[cell];
	Map[
			(# + cellCenter * explosionFactor)&, 
			cell, {2}
	]
];

explodeCellsByStep[cells_, refCellCenters_, explosionFactor_, explosionStep_, heightsMap_] := Module[
	{inactiveCells, activeCells},

	If[explosionStep == 0,
			inactiveCells = {};
			activeCells = cells,
		If[explosionStep == 2,
			inactiveCells = Select[cells, Length[Intersection[{getNonCompactCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &];
			activeCells = Select[cells, Length[Intersection[{getNonCompactCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]==0 &],
			(* explosionStep == 1 *)
			inactiveCells = Select[cells, Norm[getNonCompactCellCenter[#][[Range[2,4]]]] < epsilon &];
			activeCells = Select[cells, Norm[getNonCompactCellCenter[#][[Range[2,4]]]] > epsilon && Length[Intersection[{getNonCompactCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &]
			(* inactiveCells = {}; *)
			(* activeCells = Select[cells, Length[Intersection[{getNonCompactCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &]; *)
		];
	];

	For[acId = 1, acId <= Length[activeCells], acId++,
		cell = activeCells[[acId]];
		height = getCellHeight[cell];
		layer = heightsMap[Round[height, epsilon]];
		nonlinearFactor = getNonlinearFactor[explosionFactor, layer];
		Print["{acId, nonlinearFactor, layer}"];
		Print[{acId, nonlinearFactor, layer}];
		activeCells[[acId]] = explodedCell[cell, nonlinearFactor];
	];
	(* activeCells = Map[explodedCell[#, explosionFactor]&, activeCells]; *)
	{activeCells, inactiveCells}
];


filterCells[cells_, cellThreshold_] := Module[
	{cellCenterHeights, heightTally, tallyCounts, includeCount, heightThreshold, levelIndex},

	cellCenterHeights = Map[getCellHeight, cells];
	heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
	tallyCounts = Map[#[[2]]&, heightTally];
	Print[tallyCounts];

	includeCount = 0;
	heightThreshold = -1;
	For[levelIndex = 1, levelIndex <= Length[heightTally], levelIndex++,
		If[includeCount < cellThreshold,
			includeCount += heightTally[[levelIndex]][[2]];
			If[includeCount >= cellThreshold,
				heightThreshold = heightTally[[levelIndex]][[1]] + epsilon
			];
		];
	];

	If[heightThreshold == -1,
		heightThreshold = heightTally[[Length[heightTally]]][[1]] + epsilon
	];

	filteredCells = Select[cells, getCellCenter[#][[1]] <= heightThreshold&];
	Print["Filter from " <> IntegerString[Length[cells]] <> " cells down to " <> IntegerString[Length[filteredCells]] <> " cells."];

	filteredCells
];

getNonlinearFactor[linear_, layer_] := Module[
	{layerOffset, linearOffset, returnValue},

	minFactor = 0.004;
	layerOffset = If[layer > 1, layer - 1, 0];
	linearOffset = linear + layerOffset*delayFactor;
	returnValue = linearOffset^2/maxExplosionFactor;
	(* If[returnValue > maxExplosionFactor, returnValue = maxExplosionFactor]; *)
	If[linearOffset < 0 || returnValue < minFactor, returnValue = minFactor];
	returnValue
];

dataFileName = FileNameJoin[{dataFolder, "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm] <> "_" <> IntegerString[cellThreshold] <> ".wl"}];

If[FileExistsQ[dataFileName], 
    Print["Reading from from " <> dataFileName];
    data = Get[dataFileName],
    (* file do not exist *)
    Print["data file "<> dataFileName <> " does not exist"];
    Print["cell threshold can only be 100, 300, or 1000"];
    Exit[];
];

cells = data["cells"];
Print[cells//Length];

(* 
353:
{1, 20, 60, 12, 120}
1, 21, 81, 93, 213

535:
{1, 12, 60, 120}
1, 13, 73, 193
 *)

newCellThreshold = cellThreshold;
layers = 1;

If[{p,q,r} == {6,3,3} && cellThreshold == 10 && False,
	(* 1 + 20 + 60 + 12 *)
	newCellThreshold = 4;
	layers = 1;
];

(* for testing, 1 + 20 *)
(* newCellThreshold = 21; *)

cellCenterHeights = Map[getCellCenter[#][[1]]&, cells];
heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
tallyCounts = Map[#[[2]]&, heightTally];
Print["heights"];
Print[heightTally];
Print[tallyCounts];
layers = Length[tallyCounts];
Print["Setting layers to: " <> IntegerString[layers]];

delayFactor = If[layers >= 4, 2, 3];

cellCenters = Map[getNonCompactCellCenter, cells];
cellCenters = Sort[cellCenters];
Map[Print, cellCenters];
Map[Print[HNorm[#]]&, cellCenters];
(* Exit[]; *)

exportFrame[explosionFactor_, frameIndex_, outputFolderName_, copyLast_] := Module[
	{},
	Print[DateString["ISODateTime"] <> " Generating frame for:"];
	Print[{explosionFactor, frameIndex}];
	If[frameIndex < 145, Return[]];
	outputFrameFolder = FileNameJoin[{outputFolder, outputFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	outputFileName = "frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFrameFolder, outputFileName}];

	If[copyLast,
		lastOutputFileName = "frame_" <> IntegerString[frameIndex-1, 10, 4] <> ".png";
		lastOutputFileName = FileNameJoin[{outputFrameFolder, lastOutputFileName}];
		Print["Copying last file from " <> lastOutputFileName <> " to " <> outputFileName];
		If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
		CopyFile[lastOutputFileName, outputFileName];
		Return[];
	];

	cells = data["cells"];
	cells = filterCells[cells, newCellThreshold];

	cellCenterHeights = Map[getCellCenter[#][[1]]&, cells];
	heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
	tallyCounts = Map[#[[2]]&, heightTally];
	Print["heights"];
	Print[heightTally];
	Print[tallyCounts];
	heightsMap = <||>;
	For[htIndex = 1,htIndex<=Length[heightTally],htIndex++,
		heightsMap[Round[heightTally[[htIndex]][[1]], epsilon]] = htIndex;
	];
	Print[heightsMap];

	colorMap = <|
		1 -> Green,
		2 -> Red,
		3 -> Blue
	|>; 
	range = 1.0;
	cellsCenterHeights = Map[Round[getCellHeight[#], epsilon]&, cells];
	colorsByCells = Map[colorMap[heightsMap[#]]&, cellsCenterHeights];
	Print[colorsByCells];

	cellsActiveAndInactive = explodeCellsByStep[cells, {}, explosionFactor, 0, heightsMap];
	cellsActive = cellsActiveAndInactive[[1]];
	cellsInactive = cellsActiveAndInactive[[2]];
	Print[Map[Length, {cellsActive, cellsInactive}]];

	facesActive = Join[Flatten[cellsActive, 1]];	
	edgesActive = Join[ Flatten[Map[getEdgesFromFace, facesActive], 1]];
	facesInactive = Join[Flatten[cellsInactive, 1]];	
	edgesInactive = Join[ Flatten[Map[getEdgesFromFace, facesInactive], 1]];

	Print["The counts of active and inactive cells, faces, edges"];
	Print[Map[Length, {cellsActive, cellsInactive}]];
	Print[Map[Length, {facesActive, facesInactive}]];
	Print[Map[Length, {edgesActive, edgesInactive}]];

	projectTo3d[v_] := getPoincareExploded[v, HNorm[v]];

	edgesActive = getSplitEdges[edgesActive];
	edgesInactive = getSplitEdges[edgesInactive];

	surfaces = {};
	cellsCenters = Map[getCellCenter, Join[cellsActive, cellsInactive]];

	For[activeFlag = 1, activeFlag <= 2, activeFlag ++,
		cells = cellsActiveAndInactive[[activeFlag]];
		opacity = If[activeFlag == 1, 
			getOpacityByExplosion[explosionFactor],
			1
		];
		Print["{opacity, Length[cells]}"];
		Print[{opacity, Length[cells]}];

		For[cid = 1, cid <= Length[cells], cid++,
			(* If[Mod[cid, 3] == 1, *)
			Print[DateString["ISODateTime"] <> " Generating surface for cell id " <> IntegerString[cid]];
			(* ]; *)
			sampleCell = cells[[cid]];
			isCenterCell = Norm[getCellCenter[sampleCell][[Range[2,4]]] ] < epsilon;
			cellOpacity = If[isCenterCell, 1, opacity];
			cellHeight = getCellHeight[sampleCell];
			triangles = Flatten[Map[splitFaceToTriangles, sampleCell], 1];
			color = colorsByCells[[cid]];

			For[triangleIndex = 1, triangleIndex <= Length[triangles], triangleIndex++,
				triangle = triangles[[triangleIndex]];
				norm = HNorm[triangle[[1]]];
				(* problem: cell center is not real center. different faces have different norms *)
				(* Print[norm]; *)

				planarPoints[l1_, l2_] := 
					l2 (l1 triangle[[1]] + (1 - l1) triangle[[2]]) +
					(1 - l2) (l1 triangle[[1]] + (1 - l1) triangle[[3]]);
				normalizedPoints[l1_, l2_] := HNormalize[planarPoints[l1, l2], norm];
				poincarePoints[l1_, l2_] := getPoincareExploded[normalizedPoints[l1, l2], norm];
				surface = ParametricPlot3D[
					poincarePoints[l1, l2], {l1, 0, 1}, {l2, 0, 1},
					PlotStyle -> {color, Opacity[cellOpacity]}, Mesh -> None
				];
				AppendTo[surfaces, surface];
			];

			If[isCenterCell, 
				range = Norm[poincarePoints[0, 0]] * rangeFactor;
			];
		];
	];

	Print[DateString["ISODateTime"] <> " Populated surfaces. Count: " <> IntegerString[Length[surfaces]]];

	Print["Range of plot: " <> ToString[range]];
	thickness = range * 0.008;
	thicknessDecayExponent = 1.0;

	graphElements = {};

	Print[DateString["ISODateTime"] <> " Generated graph elements"];


	graph = Graphics3D[
		graphElements,
		ViewPoint -> {2, -4, 1}*30,
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
		Lighting -> lighting
	];

	If[Length[surfaces] > 1,
		graph = Apply[Show, Join[{graph}, surfaces]];
	];

	Print[DateString["ISODateTime"] <> " Defined graph"];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  Print[DateString["ISODateTime"] <> " Exporting to POV ray"];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];
	Print[DateString["ISODateTime"] <> " Saved an image to " <> outputFileName];
];

exportAnimation[] := Module[
	{
		gifFileName,
		gifFolder,
		frameIndex,
		explosionFactor,
		gifCommand
	},

	shapeWithParams = shape <> "_" <> IntegerString[newCellThreshold];
	outputFrameFolderName = shapeWithParams <> "_frames";

	outputFrameFolder = FileNameJoin[{outputFolder, outputFrameFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	gifFileName = shapeWithParams;
	gifFolder = outputFolder;
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	frameIndex = 1;
	explosionFactors = Table[k * 1.0, 
		{k, maxExplosionFactor, (-layers + 1)*delayFactor, 
			-(maxExplosionFactor + (layers - 1)*delayFactor)/(frameCount - 1)}];

	Print["explosionFactors"];
	Print[explosionFactors];

	expId = 1;
	explosionFactor = explosionFactors[[expId]];
	(* frameIndex++; *)
	exportFrame[explosionFactor, frameIndex++, outputFrameFolderName, False];

	For[fid = 1, fid <= Ceiling[frameCount / 6], fid++,
		exportFrame[0, frameIndex++, outputFrameFolderName, True]
	];

	For[expId = 2, expId <= Length[explosionFactors], expId++,
		explosionFactor = explosionFactors[[expId]];
		exportFrame[explosionFactor, frameIndex++, outputFrameFolderName, False];
	];

	For[fid = 1, fid <= Ceiling[frameCount / 6], fid++,
		exportFrame[0, frameIndex++, outputFrameFolderName, True]
	];

	delay = If[frameCount >= 50, 3, Floor[500/frameCount]];

	gifCommand = "magick convert -duplicate 1,-2-1 -quiet +repage -alpha remove -loop 0 -delay "
	(* gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay "  *)
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

(* exportFrame[0.0, 11, "test", False]; *)
(* exportFrame[-2.0, 12, "test", False]; *)

(* exportFrame[4.0, 10, "test", False]; *)
(* exportFrame[1.0, 20, "test", False]; *)
(* exportFrame[0.3, 21, "test", False]; *)
(* exportFrame[0.0, 22, "test", False]; *)
(* exportFrame[-2.0, 23, "test", False]; *)
(* exportFrame[-3.0, 24, "test", False]; *)

(* exportFrame[0.5, 2, "test", False]; *)
(* exportFrame[1.0, 3, "test", False]; *)
(* exportFrame[4.0, 4, "test", False]; *)

(* 
for ($num = 118 ; $num -le 134 ; $num++) { cp frame_0117.png "frame_0$num.png" }
 *)