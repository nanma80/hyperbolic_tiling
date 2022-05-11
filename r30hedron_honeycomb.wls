#!/usr/bin/env wolframscript

(* models = {"hyperboloid", "klein", "poincare"}; *)
models = {"klein", "poincare"};
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
(* models = {"poincare"}; *)

argv = Rest @ $ScriptCommandLine;
p = 5; q = 3; r = 4;

cellThreshold = 100;
(* cellThreshold = 300; *)
(* cellThreshold = 1000; *)

If[Length[argv] >= 1, 
	cellThreshold = ToExpression[argv[[1]]]
];

cellFirst = True;
(* cellFirst = False; *)

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "r30hedron_honeycomb";
baseShape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm] <> If[cellFirst, "_cell_first", ""];
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 0.8;

outputFolder = FileNameJoin[{imageFolder, "honeycombs", shape }];
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
centralReflect[point_, mirror_] := - HReflect[point, mirror];
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

getMinInnerProduct[points_] := Module[
	{minInner, i, j, ijInner},
  minInner = 10^10;
  For[i = 1, i <= Length[points]-1, i++, 
      For[j = i + 1, j <= Length[points], j++, 
          ijInner = HInner[points[[i]], points[[j]]];
          If[ijInner < minInner,
              minInner = ijInner
          ]
      ];
  ];

  minInner
];


getAlternatedSubsets[points_] := Module[
	{
		subsets, minInner, allIndices, remainingIndices, index, neighborIndices
	},
	subsets = {{}, {}};

  minInner = getMinInnerProduct[points];

  If[Length[points] >= 1,
  	AppendTo[subsets[[1]], 1]
  ];

  allIndices = Range[Length[points]];
  remainingIndices = Complement[allIndices, subsets[[1]], subsets[[2]]];
  While[
  	Length[remainingIndices] > 0,
  	index = RandomChoice[remainingIndices];
  	neighborIndices = Select[allIndices, Abs[ HInner[ points[[#]], points[[index]] ] - minInner ] < epsilon &];

  	If[Length[Intersection[subsets[[1]], neighborIndices]] > 0 && Length[Intersection[subsets[[2]], neighborIndices]] == 0,
  		subsets[[2]] = Union[subsets[[2]], {index}];
  		subsets[[1]] = Union[subsets[[1]], neighborIndices];
  	];

  	If[Length[Intersection[subsets[[1]], neighborIndices]] == 0 && Length[Intersection[subsets[[2]], neighborIndices]] > 0,
  		subsets[[1]] = Union[subsets[[1]], {index}];
  		subsets[[2]] = Union[subsets[[2]], neighborIndices];
  	];

  	If[Length[Intersection[subsets[[1]], neighborIndices]] > 0 && Length[Intersection[subsets[[2]], neighborIndices]] > 0,
  		Print["Unexpected. Neighbor indices intersect with both subsets"];
  		Print[neighborIndices];
  		Print[subsets];
  	];

  	remainingIndices = Complement[allIndices, subsets[[1]], subsets[[2]]];
  ];

	subsets
];

getEdges[vertices_] := Module[
    {edges, inner, i, j},
    edges = {};
    minInner = getMinInnerProduct[vertices];
    For[i = 1, i <= Length[vertices]-1, i++, 
        For[j = i + 1, j <= Length[vertices], j++, 
            If[
                Abs[ HInner[vertices[[i]], vertices[[j]]] - minInner ] < 10 * epsilon,
                AppendTo[edges, {i, j}] 
            ];
        ];
    ];
    edges
];

getCells[] := Module[
	{},

	Print["testing"];

	baseDataFileName = FileNameJoin[{dataFolder, baseShape <> "_" <> IntegerString[cellThreshold] <> ".wl"}];
  baseData = Get[baseDataFileName];
  baseCells = baseData["cells"];
  Print[Length[baseCells]];

  cellCenters = Map[getCellCenter, baseCells];
	cellCenterHeights = Map[getCellHeight, baseCells];
	heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
	tallyCounts = Map[#[[2]]&, heightTally];
	Print["heightTally without normalization"];
	Print[heightTally];

	subsets = getAlternatedSubsets[cellCenters];
	Print[subsets];
	Print[Map[Length, subsets]];

	centerCellIndices = Select[Range[Length[cellCenters]], Norm[cellCenters[[#]][[ Range[2,4] ]]] < epsilon &];
	Print[centerCellIndices];
	centerCell = baseCells[[centerCellIndices[[1]]]];
	mainSubsetIndex = Select[Range[2], Length[Intersection[subsets[[#]], centerCellIndices ]] > 0&][[1]];
	mainSubset = subsets[[mainSubsetIndex]];
	otherSubset = subsets[[3 - mainSubsetIndex]];

	(* 
		For base honeycomb 534 cell threshold 100, the main subset has 91 = 1+30+60 cells. The other subset has 32=12+20 cells.
		For base honeycomb 534 cell threshold 300, the main subset has 163 = 1+30+60+72 cells. The other subset has 152=12+20+120 cells.
	*)
	Print["Selected a subset of cells from base cells. The subset has this many cells: " <> IntegerString[Length[mainSubset]] <> " the other subset has this many cells: " <> IntegerString[Length[otherSubset]]];

	mainSubsetCells = baseCells[[mainSubset]];
	otherSubsetCells = baseCells[[otherSubset]];

	cells = {};
	For[cid = 1, cid <=Length[mainSubsetCells], cid++,
		(* Print["cid: " <> IntegerString[cid]]; *)
		cell1 = mainSubsetCells[[cid]];
		center1 = getCellCenter[cell1];
		faces = cell1;
		faceCenters = Map[Mean, faces];
		newVertices = Map[centralReflect[center1, #]&, faceCenters];

		edges = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameCenter];
		(* Print[Length[edges]]; *)
		edgeCenterNorm = HNorm[Mean[edges[[1]]]];
		(* Map[Print[HNormalize[Mean[#], 1][[1]]]&, edges]; *)
		(* Print[edges[[1]]]; *)
		(* Map[Print[Mean[#][[1]]]&, edges]; *)

		newVertexPairIndices = getEdges[newVertices];
		newVerticesPairs = Map[newVertices[[#]]&, newVertexPairIndices];
		(* Print[Length[newVerticesPairs]]; *)
		newVerticesPairCenterNorm = HNorm[Mean[newVerticesPairs[[1]]]];
		
		(* Map[Print[HNormalize[Mean[#], 1][[1]]]&, newVerticesPairs]; *)

		(* Print[newVerticesPairs[[1]]]; *)

		newVerticesPairs = newVerticesPairs / newVerticesPairCenterNorm * edgeCenterNorm;
		(* Map[Print[Mean[#][[1]]]&, newVerticesPairs]; *)

		newFaces = {};
		For[eid = 1, eid <= Length[edges], eid++,
			edge = edges[[eid]];
			newVertexPairMatched = Select[newVerticesPairs, ApproxSamePointLoose[Mean[#], Mean[edge]]&];
			(* Print[Length[newVertexPairMatched]]; *)
			newVertexPair = newVertexPairMatched[[1]];		
			newFace = {edge[[1]], newVertexPair[[1]], edge[[2]], newVertexPair[[2]]};
			AppendTo[newFaces, newFace];
		];

		AppendTo[cells, newFaces];
	];

	Print["generated cells"];
	Print[Length[cells]];

(* 	Print["deduping faces"];
	faces = Union[ Flatten[cells, 1], SameTest -> sameCenter];
	Print["deduping edges"];
	edges = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameCenter];
	Print["deduping vertices"];
	vertices = Union[ Flatten[faces,1], SameTest -> ApproxSamePoint];
	Print["Counts: cells, faces, edges, vertices"];
	Print[Map[Length, {cells, faces, edges, vertices}]]; *)

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

		cells = getCells[];
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

