#!/usr/bin/env wolframscript

(* models = {"hyperboloid", "klein", "poincare"}; *)
models = {"klein", "poincare"};
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
(* models = {"poincare"}; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 3, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]]; r = ToExpression[argv[[3]]],
	p = 5; q = 3; r = 4;
];

cellThreshold = 100;
(* cellThreshold = 300; *)
(* cellThreshold = 1000; *)

If[Length[argv] >= 4, 
	cellThreshold = ToExpression[argv[[4]]]
];

cellFirst = True;
(* cellFirst = False; *)

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "bitruncated_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm] <> If[cellFirst, "_cell_first", ""];
(* shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm] <> If[cellFirst, "_cell_first", ""]; *)
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 0.8;

outputFolder = FileNameJoin[{imageFolder, "honeycombs", "bitruncated_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm]}];
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
	{3,5,3} -> Red,
	{5,3,4} -> Blue,
	{4,3,5} -> Green,
	{5,3,5} -> Purple
|>;

getVertexFigureCoordinates[p_, q_, r_] := Module[
	{},

	vertexFigureCoordinates = {};

	If[q == 3 && r == 4,
		vertexFigureCoordinates = {
			{1, 0, 0},
			{0, 1, 0},
			{0, 0, 1},
			{0, 0, -1}
		};
	];

	If[q == 3 && r == 5,
		vertexFigureCoordinates = {
			{phi, 1, 0},
			{phi, -1, 0},
			{1, 0, phi},
			{1, 0, -phi}
		};
	];

	If[q == 5 && r == 3,
		vertexFigureCoordinates = {
			{phi, 1/phi, 0},
			{phi, -1/phi, 0},
			{1, 1, 1},
			{1, 1, -1}
		};
	];

	Map[#/Norm[#]&, vertexFigureCoordinates]
];

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

getCells[p_, q_, r_, cellThreshold_] := Module[
	{},
	piOverHpq = ArcCos[Sqrt[Cos[Pi/p]^2 + Cos[Pi/q]^2]]//Simplify;
	piOverHqr = ArcCos[Sqrt[Cos[Pi/q]^2 + Cos[Pi/r]^2]]//Simplify;

	(* chi: p0 and p2 *)
	coshchi = Cos[Pi/p] Cos[Pi/q] Cos[Pi/r]/Sin[piOverHpq]/Sin[piOverHqr]//Simplify;
	inner02 = coshchi;

	(* psi: p1 and p2 *)
	coshpsi = Sin[Pi/p] Cos[Pi/r]/Sin[piOverHpq]//Simplify;
	inner12 = coshpsi;

	(* phi: p0 and p1, 2phi: from p0 to p1 to another p0. edge length *)
	coshphi = Cos[Pi/p] Sin[Pi/r]/Sin[piOverHqr]//Simplify;
	inner01 = coshphi;
	cosh2phi = 2 coshphi^2 - 1//Simplify;
	sinh2phi = Sqrt[cosh2phi^2 - 1]//Simplify;

	normalizedVertexFigureCoordinates = getVertexFigureCoordinates[p, q, r];

	(* vertex first *)
	vertex1 = {1, 0, 0, 0};

	seedVertices = Map[{cosh2phi, sinh2phi #[[1]], sinh2phi #[[2]], sinh2phi #[[3]]}&, normalizedVertexFigureCoordinates];
	(* Print["seedVertices"]; *)
	mirrors = {
		seedVertices[[1]] - vertex1,
		seedVertices[[1]] - seedVertices[[2]],
		seedVertices[[2]] - seedVertices[[3]],
		seedVertices[[3]] - seedVertices[[4]]
	};

	(* normalized mirrors *)
	normalizedMirrors = Map[Re[HNormalize[#, I]] &, mirrors] // FullSimplify;
	solutions = Solve[HInner[{w, x, y, 0}, normalizedMirrors[[1]]] == 0 && 
		   HNorm[{w, x, y, 0}] == 1 && HInner[{w, x, y, 0}, normalizedMirrors[[2]]] == HInner[{w, x, y, 0}, normalizedMirrors[[3]]], {w, x, y}];
	Print[solutions];
	vertex1 = {w, x, y, 0} /. solutions[[2]] // FullSimplify;
	Print[vertex1//N];

	hexagon1 = {vertex1} // N;
	finished = False;
	For[i = 0, i < 20 && finished == False, i++,
	  beforeCount = Length[hexagon1];
	  newVertices1 = Map[HReflect[#, mirrors[[2]]] &, hexagon1];
	  newVertices2 = Map[HReflect[#, mirrors[[3]]] &, hexagon1];
	  hexagon1 = Union[hexagon1, newVertices1, newVertices2, SameTest -> ApproxSamePoint];
	  Print[{i, Length[hexagon1]}];
	  If[beforeCount == Length[hexagon1], finished = True]
	];

	hexagonCenter = Mean[hexagon1];
	proj = Orthogonalize[{{1, 2, 3, 4}, {1.23, 4.34, 2, 0}}];
	p2d = Map[(# - hexagonCenter) . Transpose[proj] &, hexagon1];
	args = Map[Arg[#[[1]] + #[[2]] I] &, p2d];
	hexagon1 = hexagon1[[Ordering[args]]];
	Print["Each type 1 face contains " <> IntegerString[Length[hexagon1]] <> " vertices. Expected value: " <> IntegerString[2 * q]];

	pentagon1 = {vertex1} // N;
	finished = False;
	For[i = 0, i < 20 && finished == False, i++, 
	  newVertex = FullSimplify[HDoubleReflect[pentagon1[[-1]], mirrors[[1]], mirrors[[2]]]];
	  If[!(ApproxSamePoint[newVertex, pentagon1[[1]]]), 
	  	AppendTo[pentagon1, newVertex], finished = True
	  ]
	];

	Print["Each type 2 face contains " <> IntegerString[Length[pentagon1]] <> " vertices. Expected value: " <> IntegerString[p]];

	cell1 = {hexagon1, pentagon1};
	countStart = -1;
	For[i=0, i < 100 && countStart != Length[cell1], i++,
		countStart = Length[cell1];
		
		For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
			cell1 = Union[cell1, Map[HReflect[#, mirrors[[mirrorIndex]] ]&, cell1, {2}], SameTest -> sameCenter];
		];
	];
	Print["Each cell contains " <> IntegerString[Length[cell1]] <> " faces."];
	
	Print["cell1 center"];
	cell1Center = getCellCenter[cell1];
	normalizedCell1Center = HNormalize[cell1Center,1];
	Print[cell1Center];
	Print[normalizedCell1Center];

	Print["flipping"];
	flipMirror = normalizedCell1Center - {1, 0, 0, 0};
	cell1Flip = Map[HReflect[#, flipMirror ]&, cell1, {2}];

	twist = ArrayFlatten[{{1, 0}, {0, RotationMatrix[Pi/q, normalizedCell1Center[[Range[2,4]]] ]}}];
	Print[twist//N];

	cell1Flip = Map[#.twist &, cell1Flip, {2}];

	Print[Length[cell1Flip]];
	Print[getCellCenter[cell1Flip]];

	cells = {cell1, cell1Flip};
	Print["only two cells"];
	Print["deduping faces"];
	faces = Union[ Flatten[cells, 1], SameTest -> sameCenter];
	Print["deduping edges"];
	edges = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameCenter];
	Print["deduping vertices"];
	vertices = Union[ Flatten[faces,1], SameTest -> ApproxSamePoint];
	Print["Counts: cells, faces, edges, vertices"];
	Print[Map[Length, {cells, faces, edges, vertices}]];



	cellThresholdStep1 = If[cellFirst, cellThreshold, cellThreshold];
	For[gIndex = 1, gIndex <= 6 && Length[cells] < cellThresholdStep1, gIndex++,
		cells = Union[cells, Map[HReflect[#, mirrors[[1]] ]&, cells, {3}], SameTest -> sameCellCenter];

		countStart = -1;
		For[i = 0, i < 10 && countStart != Length[cells], i++,
			countStart = Length[cells];
			
			For[mirrorIndex = 4, mirrorIndex >= 2, mirrorIndex--,
				cells = Union[cells, Map[HReflect[#, mirrors[[mirrorIndex]] ]&, cells, {3}], SameTest -> sameCellCenter];
			];

			Print["gIndex = "  <> IntegerString[gIndex] <> "; i = "  <> IntegerString[i] <> ". Generated " <> IntegerString[Length[cells]] <> " cells so far"];
			normalizedCellCenters = Map[HNormalize[getCellCenter[#],1]&, cells];
			normalizedCellCenterHeights = Map[#[[1]]&, normalizedCellCenters];
			normalizedHeightTally = Tally[normalizedCellCenterHeights, ApproxSamePoint]//Sort;
			Print["normalizedHeightTally"];
			Print[normalizedHeightTally];
		];
	];

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

	If[cellFirst,
		cells = filterCells[cells, cellThreshold];
		cellCenterHeights = Map[getCellHeight, cells];
		minHeight = Min[cellCenterHeights];
		Print["cellFirst"];
		Print[minHeight];
		sampleCell = Select[cells, getCellHeight[#] < minHeight + epsilon & ][[1]];

		sampleCellCenter = HNormalize[getCellCenter[sampleCell], 1];
		Print[sampleCellCenter];

		pentagons = Select[sampleCell, Mod[Length[#],5]==0&];
		zVector = HNormalize[Mean[pentagons[[1]]], 1];

		Print[zVector];
		zVectorToZaxis = ArrayFlatten[{{1, 0}, {0, RotationMatrix[{zVector[[Range[2,4]]], {0, 0, 1}}]}}];
		setup = Transpose[zVectorToZaxis];
		Map[Print, setup];

		cells = Map[#.setup &, cells, {3}];
	];

	cells
];

GetDataThroughCache[p, q, r] := Module[
	{},
	data = <||>;
	dataFileName = FileNameJoin[{dataFolder, shape <> "_" <> IntegerString[cellThreshold] <> ".wl"}];
	If[FileExistsQ[dataFileName], 
      Print["Reading from from " <> dataFileName];
      data = Get[dataFileName],
      (* file do not exist *)
      Print["Creating data for " <> dataFileName];

		cells = getCells[p, q, r, cellThreshold];
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

data = GetDataThroughCache[p,q,r];
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

	color = colorMap[{p,q,r}];

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
		(* ViewPoint -> {3, -4, 1}*30, *)
		(* for testing *)
		ViewPoint -> {0, 0, 1}*1000, ViewVertical -> {1, 0, 0},
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

(* 
cell centers height tally:
535 cell first, 300:
before normalization:
{{3.4270509831248432, 1}, {12.803695435310923, 12}, {27.975424859371593, 60}, {52.52379873561811, 120}, {67.69552815967873, 120}}
after normalization:
{{0.9999999999999999, 1}, {3.7360679774994017, 12}, {8.163118960624029, 60}, {15.326237921246358, 120}, {19.75328890437116, 120}}

535 vertex first, 300:
before normalization:
{{11.744678440936955, 20}, {36.293052317185825, 60}, {51.464781741247634, 120}, {76.0131556174965, 90}, {91.1848850415583, 120}}
after normlization:
{{3.4270509831248397, 20}, {10.590169943749803, 60}, {15.017220926873913, 120}, {22.180339887500928, 90}, {26.60739087061998, 120}}

bitruncated 535:
Why is the projection not having 5-fold symmetry?

Fixed flipped cell. After flipping, we need to twist. Then. with 300 cells:
without normalization:
{{1.7022169326902425, 1}, {5.83358421266784, 20}, {6.35959817298194, 12}, {13.895399318340015, 60}, {18.026766598317607, 60}, {25.562567743675686, 60}, {33.62438284934785, 60}, {37.75575012932545, 60}, {45.81756523499766, 12}}

with normalization:
{{1., 1}, {3.427050983124836, 20}, {3.736067977499712, 12}, {8.163118960624788, 60}, {10.59016994374946, 60}, {15.017220926873998, 60}, {19.753288904380668, 60}, {22.18033988749935, 60}, {26.916407864995026, 10}, {26.916407865003478, 2}}
Filter from 345 cells down to 153 cells. (1+20+12+60+60)

Heights 1 ~ 10.59 are correct. However, 15.01 should be 120 but we got only 60. Why?
Maybe there are multiple orbits of points contributing to 15.01. The symmetry operations we've done so far only reaches one orbit.

When we set the cell threshold to 500, applying one more mirror, we have:
{{1., 1}, {3.427050983124836, 20}, {3.736067977499712, 12}, {8.163118960624706, 60}, {10.59016994374946, 60}, {15.017220926873998, 78}, {15.326237921249557, 20}, {19.753288904380668, 60}, {22.18033988749935, 69}, {26.607390870624275, 20}, {26.916407864995026, 26}, {26.916407865003478, 6}, {29.343458848128716, 20}, ...



 *)