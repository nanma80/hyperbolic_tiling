#!/usr/bin/env wolframscript

(* Objective: visualize 2{5, 3, 4}[5{4, 3, 5}] *)

maxExplosionFactor = 10.0;
frameCount = 80;

(* maxExplosionFactor = 8.0; *)
(* frameCount = 6; *)

(* models = {"hyperboloid", "klein", "poincare"}; *)
(* models = {"klein", "poincare"}; *)
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
(* models = {"poincare"}; *)

(* 
explosion steps: 
0: explode all cells
1: hide the moving cells in step 2, and explode the other cells in 534 and 435
2: move away the cells of 435 that don't share centers with any 534 cells
 *)

(* 2: 1 cell. 4: 7 cells per component. 5: 31 cells per component *)
(* heightThreshold = 2; *)
(* heightThreshold = 4; *)
(* heightThreshold = 5; *)
(* heightThreshold = 7; *)

showVertices = False;

argv = Rest @ $ScriptCommandLine;
If[Length[argv] < 1,
	Print["Usage: wolframscript.exe transition_compound_534.wls <displayCompIndex>"];
];

displayCompIndex = If[Length[argv] >= 1, 
	ToExpression[argv[[1]]],
	0
];

explosionStep = 1;
show534 = True;

(* explosionFactor = explosionPercent / 100.0; *)

(* Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"]; *)
(* shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm]; *)
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 1.25;
colors = {Red, Green, Blue, Purple, Yellow};
color534 = Gray;
colorVertices = Cyan;

outputFolder = FileNameJoin[{imageFolder, "explode_534"}];
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
getOpacityByExplosion[factor_] := If[factor < maxExplosionFactor/2, 1, 2 - 2/maxExplosionFactor * factor];
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
explodedCell[cell_, explosionFactor_] := Map[
	(# + Mean[Map[Mean,cell]] * (HNorm[First[First[cell]]]/HNorm[Mean[Map[Mean, cell]]])^1.5 * explosionFactor)&, 
	cell, {2}];

explodeCellsByStep[cells_, refCellCenters_, explosionFactor_, explosionStep_] := Module[
	{inactiveCells, activeCells},

	If[explosionStep == 0,
			inactiveCells = {};
			activeCells = cells,
		If[explosionStep == 2,
			inactiveCells = Select[cells, Length[Intersection[{getCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &];
			activeCells = Select[cells, Length[Intersection[{getCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]==0 &],
			(* explosionStep == 1 *)
			inactiveCells = Select[cells, Norm[getCellCenter[#][[Range[2,4]]]] < epsilon &];
			activeCells = Select[cells, Norm[getCellCenter[#][[Range[2,4]]]] > epsilon && Length[Intersection[{getCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &]
			(* inactiveCells = {}; *)
			(* activeCells = Select[cells, Length[Intersection[{getCellCenter[#]}, refCellCenters, SameTest -> ApproxSamePoint ]]>0 &]; *)
		];
	];
	activeCells = Map[explodedCell[#, explosionFactor]&, activeCells];
	{activeCells, inactiveCells}
];

getAllComponentsThroughCache[heightThreshold_] := Module[
	{data, subsetDataFileName},
	data = <||>;
	subsetDataFileName = FileNameJoin[{dataFolder, "compound_534_435_cell_first_" <> IntegerString[heightThreshold] <> ".wl"}];
	If[FileExistsQ[subsetDataFileName], 
      Print["Reading from from " <> subsetDataFileName];
      data = Get[subsetDataFileName],
      (* file do not exist *)
      Print["Run compound_534_cell_first.wls first with the same heightThreshold to create cache data"];
      Exit[];
  ];
  data
];

heightThreshold = If[show534, 7, 5];
(* heightThreshold = 5; *)

data = getAllComponentsThroughCache[heightThreshold];

exportFrame[faceOpacity_, frameIndex_, outputFolderName_, copyLast_] := Module[
	{},
	(* Return[]; *)
	explosionFactor = 0;
	componentCells = data["componentCells"];
	cells534 = data["convexHullCells"];
	(* vertices534 = data["convexHullVertices"]; *)
	range = 0.6;
	componentCellCenterHeights = Apply[Join, Map[getCellHeight, componentCells, {2}]];
	cells534CenterHeights = Map[getCellHeight, cells534];
	componentHeightsTally = Sort[Tally[componentCellCenterHeights, ApproxSamePoint]];
	cells534HeightsTally = Sort[Tally[cells534CenterHeights, ApproxSamePoint]];
	Print["Heights of components and 534 before explosion"];
	Print[componentHeightsTally];
	Print[cells534HeightsTally];

	(* 
	explode the cells: 
	when height threshold = 5:
	{{1.851229586821858, 5}, {2.995352392457251, 30}, {4.84658197927913, 120}}
	{{1.8512295868219169, 1}, {4.846581979279199, 12}}

	explosion steps: 
	0: explode all cells
	1: hide the moving cells in step 2, and explode the other cells in 534 and 435
	2: move away the cells of 435 that don't share centers with any 534 cells
	*)

	componentCellCenters = Map[getCellCenter, componentCells, {2}];
	cells534Centers = Map[getCellCenter, cells534];

	cells534ActiveAndInactive = explodeCellsByStep[cells534, cells534Centers, explosionFactor, explosionStep];
	cells534Active = cells534ActiveAndInactive[[1]];
	cells534Inactive = cells534ActiveAndInactive[[2]];
	componentCellsActiveAndInactive = Map[explodeCellsByStep[#, cells534Centers, explosionFactor, explosionStep]&, componentCells];
	componentCellsActive = Map[#[[1]]&, componentCellsActiveAndInactive];
	componentCellsInactive = Map[#[[2]]&, componentCellsActiveAndInactive];

	(* vertexNorm = HNorm[cells534[[1]][[1]][[1]]]; *)
	(* Print["vertexNorm"]; *)
	(* Print[vertexNorm]; *)

	(* componentFaces = Map[Join[Flatten[#, 1]] &, componentCells]; *)
	(* componentEdges = Map[Join[Flatten[Map[getEdgesFromFace, #], 1]]&, componentFaces]; *)
	(* componentVertices = Map[Join[Flatten[#, 1]]&, componentFaces]; *)

	faces534Active = Join[Flatten[cells534Active, 1]];	
	edges534Active = Join[ Flatten[Map[getEdgesFromFace, faces534Active], 1]];
	faces534Inactive = Join[Flatten[cells534Inactive, 1]];	
	edges534Inactive = Join[ Flatten[Map[getEdgesFromFace, faces534Inactive], 1]];
	(* vertices534 = Join[ Flatten[faces534, 1]]; *)


	Print["The counts of active and inactive cells for all components"];
	Print[Map[Length, componentCellsActive]];
	Print[Map[Length, componentCellsInactive]];
	(* Print[Map[Length, componentFaces]]; *)
	(* Print[Map[Length, componentEdges]]; *)
	(* Print[Map[Length, componentVertices]]; *)

	Print["The counts of active and inactive 534 cells, faces, edges"];
	Print[Map[Length, {cells534Active, cells534Inactive}]];
	Print[Map[Length, {faces534Active, faces534Inactive}]];
	Print[Map[Length, {edges534Active, edges534Inactive}]];

	projectTo3d[v_] := getPoincareExploded[v, HNorm[v]];

	(* splitCompEdges = {}; *)

(* 	For[compIndex=1, compIndex<=Length[componentEdges],compIndex++,
		edges = componentEdges[[compIndex]];
		splitEdges = {};
		For[eid = 1, eid <= Length[edges], eid++,
			edge = edges[[eid]];
			norm = HNorm[edge[[1]]];
			splitParts = splitEdge[edge, splitEdgeParts];
			splitParts = Map[HNormalize[#, norm]&, splitParts, {2}];
			splitEdges = Join[splitEdges, splitParts];
		];
		AppendTo[splitCompEdges, splitEdges]
	]; *)

	(* split 534 edges *)

	edges534Active = getSplitEdges[edges534Active];
	edges534Inactive = getSplitEdges[edges534Inactive];

	surfaces = {};
	cells534Centers = Map[getCellCenter, Join[cells534Active, cells534Inactive]];
	For[compIndex=1, compIndex<=Length[componentCellsActiveAndInactive],compIndex++,
		If[(displayCompIndex == 0 || displayCompIndex == compIndex),
			cellsTuple = componentCellsActiveAndInactive[[compIndex]];
			For[activeFlag = 2, activeFlag <= 2, activeFlag ++,
				opacity = 1;
				cells = cellsTuple[[activeFlag]];
				Print["Component {opacity, Length[cells]}"];
				Print[{opacity, Length[cells]}];

				For[cid = 1, cid <= Length[cells], cid++,
					sampleCell = cells[[cid]];
					sharing534Centers = Length[Intersection[{getCellCenter[sampleCell]}, cells534Centers, SameTest -> ApproxSamePoint ]] > 0;
					cellHeight = getCellHeight[sampleCell];

					triangles = Flatten[Map[splitFaceToTriangles, sampleCell], 1];
					For[triangleIndex = 1, triangleIndex <= Length[triangles], triangleIndex++,
						triangle = triangles[[triangleIndex]];
						norm = HNorm[triangle[[1]]];
						planarPoints[l1_, l2_] := 
							l2 (l1 triangle[[1]] + (1 - l1) triangle[[2]]) +
							(1 - l2) (l1 triangle[[1]] + (1 - l1) triangle[[3]]);
						normalizedPoints[l1_, l2_] := HNormalize[planarPoints[l1, l2], norm];
						poincarePoints[l1_, l2_] := getPoincareExploded[normalizedPoints[l1, l2], norm];
						color = Lighter[colors[[compIndex]], If[sharing534Centers, 0, 0.1]];
						surface = ParametricPlot3D[
							poincarePoints[l1, l2], {l1, 0, 1}, {l2, 0, 1},
							PlotStyle -> {color, Opacity[opacity]}, Mesh -> None
						];
						AppendTo[surfaces, surface];
					];

					If[Norm[getCellCenter[sampleCell][[Range[2,4]]] ] < epsilon, 
						range = Norm[poincarePoints[0, 0]] * rangeFactor;
					];
				];			
			];			
		];
	];

	If[show534,
		For[activeFlag = 2, activeFlag <= 2, activeFlag ++,
			cells = cells534ActiveAndInactive[[activeFlag]];
			opacity = faceOpacity;
			Print["534 {opacity, Length[cells]}"];
			Print[{opacity, Length[cells]}];

			For[cid = 1, cid <= Length[cells], cid++,
				sampleCell = cells[[cid]];
				cellHeight = getCellHeight[sampleCell];
				triangles = Flatten[Map[splitFaceToTriangles, sampleCell], 1];
				cellHeight = getCellHeight[sampleCell];

				For[triangleIndex = 1, triangleIndex <= Length[triangles], triangleIndex++,
					triangle = triangles[[triangleIndex]];
					norm = HNorm[triangle[[1]]];
					planarPoints[l1_, l2_] := 
						l2 (l1 triangle[[1]] + (1 - l1) triangle[[2]]) +
						(1 - l2) (l1 triangle[[1]] + (1 - l1) triangle[[3]]);
					normalizedPoints[l1_, l2_] := HNormalize[planarPoints[l1, l2], norm];
					poincarePoints[l1_, l2_] := getPoincareExploded[normalizedPoints[l1, l2], norm];
					color = color534;
					surface = ParametricPlot3D[
						poincarePoints[l1, l2], {l1, 0, 1}, {l2, 0, 1},
						PlotStyle -> {color, Opacity[opacity]}, Mesh -> None
					];
					AppendTo[surfaces, surface];
				];

				If[Norm[getCellCenter[sampleCell][[Range[2,4]]] ] < epsilon, 
					range = Norm[poincarePoints[0, 0]] * rangeFactor;
				];
			];
		];
	];

	Print["Populated surfaces. Count: " <> IntegerString[Length[surfaces]]];

	Print["Range of plot: " <> ToString[range]];
	thickness = range * 0.008;
	thicknessDecayExponent = 1.0;

	graphElements = {};

	edgesTuple = {edges534Active, edges534Inactive};
	For[activeFlag = 2, activeFlag <= 2, activeFlag++,
		edges = edgesTuple[[activeFlag]];
		opacity = 1 - faceOpacity;
		Print["534 edge {opacity, Length[edges]}"];
		Print[{opacity, Length[edges]}];
		AppendTo[graphElements,
			{
				Opacity[opacity], FaceForm[color534], 
				Table[
		            Cylinder[
		                Map[projectTo3d, edges[[edgeIndex]]], 
		                thickness * 5/(1 + (edges[[edgeIndex]][[1]][[1]]/HNorm[edges[[edgeIndex]][[1]]])^thicknessDecayExponent + (edges[[edgeIndex]][[2]][[1]]/HNorm[edges[[edgeIndex]][[2]]])^thicknessDecayExponent)
		            ], 
		            {edgeIndex, 1, Length[edges]}
		        ]
			}
		];
	];

	Print["Generated graph elements"];

	graph = Graphics3D[
		graphElements,
		ViewPoint -> {3, -3.5, 1}*30,
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
		Lighting -> lighting
	];

	If[Length[surfaces] > 1,
		graph = Apply[Show, Join[{graph}, surfaces]];
	];

	Print["Defined graph"];

	outputFrameFolder = FileNameJoin[{outputFolder, outputFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	outputFileName = "frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFrameFolder, outputFileName}];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[copyLast,
		lastOutputFileName = "frame_" <> IntegerString[frameIndex-1, 10, 4] <> ".png";
		lastOutputFileName = FileNameJoin[{outputFrameFolder, lastOutputFileName}];
		Print["Copying last file from " <> lastOutputFileName <> " to " <> outputFileName];
		If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
		CopyFile[lastOutputFileName, outputFileName],	

		If[exportToPov,
		  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
		  Print["Exporting to POV ray"];
		  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
		  Run[cropCommand];
		  ,
		  Export[outputFileName, graph];
		];
		Print["Saved an image to " <> outputFileName];
	];
];


exportAnimation[] := Module[
	{
		gifFileName,
		gifFolder,
		frameIndex,
		opacities,
		opacity,
		gifCommand
	},

	shapeWithParams = If[displayCompIndex == 0, "all_components", "component_" <> IntegerString[displayCompIndex]] <> "_transition_" <> IntegerString[heightThreshold];
	outputFrameFolderName = shapeWithParams <> "_frames";

	outputFrameFolder = FileNameJoin[{outputFolder, outputFrameFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	gifFileName = shapeWithParams;
	gifFolder = outputFolder;
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	frameIndex = 1;

	opacities1 = Table[
		1,
		{k, Ceiling[frameCount/3]}
	];

	opacities2 = Table[
		1 - k/(frameCount-1),
		{k, 0, frameCount - 1}
	];

	opacities = Join[opacities1, opacities2];

	Print["opacities"];
	Print[opacities];

	For[expId = 1, expId <= Length[opacities], expId++,
		opacity = opacities[[expId]];
		exportFrame[opacity, frameIndex++, outputFrameFolderName, False];
	];

	delay = If[frameCount >= 50, 3, Floor[500/frameCount]];

	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[delay] <> " -dispose 2 " <> outputFrameFolder <> "\\frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

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

(* exportFrame[0.0, 1, "test"]; *)
(* 
for ($num = 118 ; $num -le 134 ; $num++) { cp frame_0117.png "frame_0$num.png" }
 *)