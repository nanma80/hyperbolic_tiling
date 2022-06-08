#!/usr/bin/env wolframscript

(* models = {"hyperboloid", "klein", "poincare"}; *)
models = {"klein", "poincare"};
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
(* models = {"poincare"}; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2,
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 5; q = 3
];

cellThreshold = 50;
(* cellThreshold = 300; *)
(* cellThreshold = 1000; *)

If[Length[argv] >= 3,
	cellThreshold = ToExpression[argv[[3]]]
];

cellFirst = False;

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", " <> IntegerString[p] <> ", " <> IntegerString[q] <> "}"];
shape = "omni_ring_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> If[cellFirst, "_cell_first", ""];
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 0.8;

outputFolder = FileNameJoin[{imageFolder, "honeycombs", "omni_ring_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm]}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
exportToPov = True;
splitEdgeParts = 8;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

phi = (1 + Sqrt[5])/2;
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
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

Boost[boost_] := {{Cosh[boost], Sinh[boost], 0, 0}, {Sinh[boost], Cosh[boost], 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];
ApproxSamePointLoose[point1_, point2_] := Round[point1, 100 epsilon] == Round[point2, 100 epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Mean[N[set1]], Mean[N[set2]]];
sameCellCenter[cell1_, cell2_] := ApproxSamePoint[getCellCenter[cell1], getCellCenter[cell2]];
getCellCenter[cell_] := Total[Flatten[cell, 1]]/Length[Flatten[cell, 1]];
getCellHeight[cell_] := getCellCenter[cell][[1]];

getHyperboloid[v_] := {v[[2]], v[[3]], v[[4]]};
getKlein[v_] := {v[[2]], v[[3]], v[[4]]}/v[[1]];
getPoincare[v_] := {v[[2]], v[[3]], v[[4]]}/(1+v[[1]]);

colorMap := <|
	{5,3} -> Blue,
	{4,3} -> Green
|>;


filterCells[cells_, cellThreshold_] := Module[
	{cellCenterHeights, heightTally, tallyCounts, includeCount, heightThreshold, levelIndex},

	cellCenterHeights = Map[getCellCenter[#][[1]]&, cells];
	heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
	tallyCounts = Map[#[[2]]&, heightTally];

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


getCells[p_, q_, cellThreshold_] := Module[
	{},

	(* 4,3,4,3 *)
	mirrors = {
		{0, 0, 1, 0},
		{0, 0, -1, -1}/Sqrt[2],
		{0, 1, 0, 1}/Sqrt[2],
		{Sqrt[7]/2, -3/2, -1/2, 1/2}
	};

	phi = (1 + Sqrt[5])/2;

	If[{p, q} == {5, 3},
		mirrors = {
			{0, 0, 0, 2},
			{0, 1, 1/phi, -phi},
			{0, -2, 0, 0},
			{Sqrt[4 phi^4 + phi^2-3], phi, -2 phi^2, -1}
		}/2;
	];

	mirrorAngles = Table[Pi/ArcCos[HInner[mirrors[[i]], mirrors[[j]]]]//N, {i, 4}, {j, 4}];
	Map[Print, mirrorAngles];

	solutions = NSolve[ HNorm[{w,x,y,z}] == 1 && w > 0 &&
		HInner[{w, x, y, z}, mirrors[[1]]] == HInner[{w, x, y, z}, mirrors[[2]]] &&
		HInner[{w, x, y, z}, mirrors[[2]]] == HInner[{w, x, y, z}, mirrors[[3]]] &&
		HInner[{w, x, y, z}, mirrors[[3]]] == HInner[{w, x, y, z}, mirrors[[4]]], 
		{w, x, y, z}];
	Print[solutions];
	If[Length[solutions] < 1,
		Print["Cannot solve for initial vertex"];
		Exit[];
	];
	vertex1 = {w, x, y, z} /. solutions[[1]];
	Print[vertex1//N];

	cells = {};

	For[excludeMirrorIndex = 1, excludeMirrorIndex <= 4, excludeMirrorIndex++,
		orderedMirrors = mirrors[[Complement[Range[4], {excludeMirrorIndex}]]];

		faces = {};

		For[mirror1Index = 1, mirror1Index <= Length[orderedMirrors], mirror1Index++,
			mirror2Index = mirror1Index + 1;
			If[mirror2Index > Length[orderedMirrors], mirror2Index = 1];
			hexagon1 = {vertex1} // N;
			finished = False;
			For[i = 0, i < 20 && finished == False, i++,
			  beforeCount = Length[hexagon1];
			  newVertices1 = Map[HReflect[#, orderedMirrors[[mirror1Index]]] &, hexagon1];
			  newVertices2 = Map[HReflect[#, orderedMirrors[[mirror2Index]]] &, hexagon1];
			  hexagon1 = Union[hexagon1, newVertices1, newVertices2, SameTest -> ApproxSamePoint];
			  Print[{i, Length[hexagon1]}];
			  If[beforeCount == Length[hexagon1], finished = True]
			];

			hexagonCenter = Mean[hexagon1];
			proj = Orthogonalize[{{1, 2, 3, 4}, {1.23, 4.34, 2, 0}}];
			p2d = Map[(# - hexagonCenter) . Transpose[proj] &, hexagon1];
			args = Map[Arg[#[[1]] + #[[2]] I] &, p2d];
			hexagon1 = hexagon1[[Ordering[args]]];
			Print["Each type 1 face contains " <> IntegerString[Length[hexagon1]] <> " vertices"];
			
			AppendTo[faces, hexagon1];
		];
		
		cell1 = faces;
		countStart = -1;
		For[i=0, i < 100 && countStart != Length[cell1], i++,
			countStart = Length[cell1];
			
			For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
				cell1 = Union[cell1, Map[HReflect[#, orderedMirrors[[mirrorIndex]] ]&, cell1, {2}], SameTest -> sameCenter];
			];
		];
		Print["Each cell contains " <> IntegerString[Length[cell1]] <> " faces."];

		AppendTo[cells, cell1];
	];


	Print["only initial cells"];
	Print["deduping faces"];
	faces = Union[ Flatten[cells, 1], SameTest -> sameCenter];
	Print["deduping edges"];
	edges = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameCenter];
	Print["deduping vertices"];
	vertices = Union[ Flatten[faces,1], SameTest -> ApproxSamePoint];
	Print["Counts: cells, faces, edges, vertices"];
	Print[Map[Length, {cells, faces, edges, vertices}]];

	For[gIndex = 1, gIndex <= 6 && Length[cells] < cellThreshold, gIndex++,
		cells = Union[cells, Map[HReflect[#, mirrors[[4]] ]&, cells, {3}], SameTest -> sameCellCenter];

		countStart = -1;
		For[i = 0, i < 10 && countStart != Length[cells], i++,
			countStart = Length[cells];
			
			For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
				cells = Union[cells, Map[HReflect[#, mirrors[[mirrorIndex]] ]&, cells, {3}], SameTest -> sameCellCenter];
			];

			Print["gIndex = "  <> IntegerString[gIndex] <> "; i = "  <> IntegerString[i] <> ". Generated " <> IntegerString[Length[cells]] <> " cells so far"];

			If[Length[cells] > cellThreshold,
				cells = filterCells[cells, cellThreshold];
			];

			normalizedCellCenters = Map[HNormalize[getCellCenter[#],1]&, cells];
			normalizedCellCenterHeights = Map[#[[1]]&, normalizedCellCenters];
			normalizedHeightTally = Tally[normalizedCellCenterHeights, ApproxSamePoint]//Sort;
			Print["normalizedHeightTally"];
			Print[normalizedHeightTally];
		];
	];

	cells = filterCells[cells, cellThreshold];

	Print["Generated " <> IntegerString[Length[cells]] <> " cells in total."];
	cellCenterHeights = Map[getCellHeight, cells];
	heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
	tallyCounts = Map[#[[2]]&, heightTally];
	Print["heightTally without normalization"];
	Print[heightTally];

	normalizedCellCenters = Map[HNormalize[getCellCenter[#],1]&, cells];
	normalizedCellCenterHeights = Map[#[[1]]&, normalizedCellCenters];
	normalizedHeightTally = Tally[normalizedCellCenterHeights, ApproxSamePoint]//Sort;
	Print["normalizedHeightTally"];
	Print[normalizedHeightTally];

	cells
];

GetDataThroughCache[p, q] := Module[
	{},
	data = <||>;
	dataFileName = FileNameJoin[{dataFolder, shape <> "_" <> IntegerString[cellThreshold] <> ".wl"}];
	If[FileExistsQ[dataFileName], 
      Print["Reading from from " <> dataFileName];
      data = Get[dataFileName],
      (* file do not exist *)
      Print["Creating data for " <> dataFileName];

		cells = getCells[p, q, cellThreshold];
		cellCenterHeights = Map[getCellCenter[#][[1]]&, cells];
		heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
		tallyCounts = Map[#[[2]]&, heightTally];

		Print[heightTally];
		Print[tallyCounts];

		Print["deduping faces"];
		faces = Union[ Flatten[cells, 1], SameTest -> sameCenter];
		Print["deduping edges"];
		edges = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameCenter];
		Print["deduping vertices"];
		vertices = Union[ Flatten[faces,1], SameTest -> ApproxSamePoint];
		data = <|
			"cells" -> cells,
			"faces" -> faces,
			"edges" -> edges,
			"vertices" -> vertices
		|>;

    Save[dataFileName, data];
    Print["Saved data " <> dataFileName];
  ];
  data
];

data = GetDataThroughCache[p,q];
cells = data["cells"];
faces = data["faces"];
edges = data["edges"];
vertices = data["vertices"];
Print["Counts: cells, faces, edges, vertices"];
Print[Map[Length, {cells, faces, edges, vertices}]];

Print["unique cells"];
uniqueCells = Union[cells, SameTest -> sameCellCenter];
Print[Map[Length, {cells, uniqueCells}]];

normalizedCellCenters = Map[HNormalize[getCellCenter[#],1]&, cells];
normalizedCellCenterHeights = Map[#[[1]]&, normalizedCellCenters];
normalizedHeightTally = Tally[normalizedCellCenterHeights, ApproxSamePoint]//Sort;
Print["normalizedHeightTally"];
Print[normalizedHeightTally];


For[modelIndex = 1, modelIndex <= Length[models], modelIndex++,
	model = models[[modelIndex]];
	projectTo3d[v_] := If[
		model == "hyperboloid", 
			getHyperboloid[v],
		If[
			model == "klein",
				getKlein[v],
				getPoincare[v]
		]
	];

	If[model == "poincare",
		Print["There are " <> IntegerString[Length[edges]] <> " edges before split."];
		norm = HNorm[edges[[1]][[1]]];
		splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1];
		splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
		Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."];
		edges = splitEdges;
	];

	range = Max[Map[Norm[projectTo3d[#]]&, vertices]] * rangeFactor;
	thickness = range * 0.008;
	thicknessDecayExponent = 1;

	color = colorMap[{p,q}];

	graph = Graphics3D[
		{
			(* {
				Opacity[0.2], FaceForm[color], Map[Polygon, Map[projectTo3d, faces, {2}]]
			}, *)
			{
				Opacity[1], FaceForm[color], 
				Table[
		            Cylinder[
		                Map[projectTo3d, edges[[edgeIndex]]], 
		                thickness * 5/(1 + edges[[edgeIndex]][[1]][[1]]^thicknessDecayExponent + edges[[edgeIndex]][[2]][[1]]^thicknessDecayExponent)
		            ], 
		            {edgeIndex, 1, Length[edges]}
		        ]
			},
			{
				Opacity[1], FaceForm[Cyan], 
				Table[
		            Sphere[
		                projectTo3d[vertices[[vIndex]]],
		                thickness * 10/(1 + vertices[[vIndex]][[1]]^thicknessDecayExponent)
		            ], 
		            {vIndex, 1, Length[vertices]}
		        ]
			}
		},
		(* for real output *)
		ViewPoint -> {3, -4, 2}*30,
		(* for testing *)
		(* ViewPoint -> {0, 0, 1}*1000, ViewVertical -> {1, 0, 0}, *)
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
		Lighting -> lighting
	];

	outputFileName = shape <> "_" <> model <> "_" <> IntegerString[cellThreshold] <> ".png";

	outputFileName = FileNameJoin[{outputFolder, outputFileName}];
	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
		If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
		POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
		Run[cropCommand];
		,
		Export[outputFileName, graph]
	];
	Print["Saved the image to " <> outputFileName];
];

