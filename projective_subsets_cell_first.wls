#!/usr/bin/env wolframscript

(* models = {"hyperboloid", "klein", "poincare"}; *)
(* models = {"klein", "poincare"}; *)
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
models = {"poincare"};

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 3, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]]; r = ToExpression[argv[[3]]],
	p = 3; q = 5; r = 3;
];

cellThreshold = 100;
(* cellThreshold = 300; *)
(* cellThreshold = 1000; *)

If[Length[argv] >= 4, 
	cellThreshold = ToExpression[argv[[4]]]
];

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm] <> "_cell_first";
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 0.8;

outputFolder = FileNameJoin[{imageFolder, "honeycombs"}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
exportToPov = True;
splitEdgeParts = 15;

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

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
centralReflect[point_, mirror_] := - HReflect[point, mirror];
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];
ApproxSamePointLoose[point1_, point2_] := Round[point1, 10 epsilon] == Round[point2, 10 epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
sameCellCenter[cell1_, cell2_] := sameCenter[Flatten[cell1, 1], Flatten[cell2, 1]];
getCellCenter[cell_] := Total[Flatten[cell, 1]]/Length[Flatten[cell, 1]];

getHyperboloid[v_] := {v[[2]], v[[3]], v[[4]]};
getKlein[v_] := {v[[2]], v[[3]], v[[4]]}/v[[1]];
getPoincare[v_] := {v[[2]], v[[3]], v[[4]]}/(1+v[[1]]);

colorMap := <|
	{3,5,3} -> Red,
	{5,3,4} -> Blue,
	{4,3,5} -> Green,
	{5,3,5} -> Purple
|>;


getEdges[vertices_, ignoreLessThan_] := Module[
    {edges, inner, i, j},
    edges = {};
    minInner = 10^10;
    For[i = 1, i <= Length[vertices]-1, i++, 
        For[j = i + 1, j <= Length[vertices], j++, 
            ijInner = HInner[vertices[[i]], vertices[[j]]];
            If[ijInner > ignoreLessThan + epsilon && ijInner < minInner,
                minInner = ijInner
            ]
        ];
    ];
    Print["minInner ="];
    Print[minInner];

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

getTriangles[vertices_] := Module[
    {edges, inner, i, j},
    edges = {};
    minInner = 10^10;
    For[i = 1, i <= Length[vertices]-1, i++, 
        For[j = i + 1, j <= Length[vertices], j++, 
            ijInner = HInner[vertices[[i]], vertices[[j]]];
            If[ijInner < minInner,
                minInner = ijInner
            ]
        ];
    ];
    Print["getTriangles: minInner ="];
    Print[minInner];

    For[i = 1, i <= Length[vertices]-2, i++, 
        For[j = i + 1, j <= Length[vertices]-1, j++, 
            For[k = j + 1, k <= Length[vertices], k++, 
                If[
                    Abs[HInner[vertices[[i]], vertices[[j]]] - minInner] < epsilon &&
                    Abs[HInner[vertices[[i]], vertices[[k]]] - minInner] < epsilon &&
                    Abs[HInner[vertices[[j]], vertices[[k]]] - minInner] < epsilon,
                    AppendTo[edges, {i, j, k}] 
                ];
            ];
        ];
    ];
    edges
];

getAllSubsetsThroughCache[vertices_, reflectionCenters_] := Module[
	{},
	allSubsets = {};
	subsetDataFileName = FileNameJoin[{dataFolder, shape <> "_" <> IntegerString[cellThreshold] <> "_subsets_cell_first.wl"}];
	If[FileExistsQ[subsetDataFileName], 
      Print["Reading from from " <> subsetDataFileName];
      allSubsets = Get[subsetDataFileName],
      (* file do not exist *)
      Print["Creating data for " <> subsetDataFileName];
      allSubsets = getAllSubsets[vertices, reflectionCenters];
      Save[subsetDataFileName, allSubsets];
      Print["Saved data to " <> subsetDataFileName];
  ];
  allSubsets
];

getAllSubsets[vertices_, reflectionCenters_] := Module[
	{maxHeight, discoveredIndices, allSubsets, remaining, seedIndex, subset1, k, verticesMap, subsetIndices},

	maxHeight = Max[Map[#[[1]]&, vertices]];
	
	allSubsets = {};

	verticesMap = <||>;
	For[vi = 1, vi <= Length[vertices],vi++,
		verticesMap[Round[ vertices[[vi]], epsilon] ] = vi
	];

	While[
		Length[Apply[Union, allSubsets]] < Length[vertices],
		discoveredIndices = Apply[Union, allSubsets];
		remaining = Sort[Complement[Range[Length[vertices]], discoveredIndices]];

		Print[];
		Print["remaining count"];
		Print[Length[remaining]];

		seedIndex = remaining[[1]];
		Print[seedIndex];
		subset1 = {vertices[[seedIndex]]};
		For[k = 1, k < 30, k++,
			beforeCount = Length[subset1];
			subset1 = Flatten[ Table[centralReflect[subset1[[i]], reflectionCenters[[j]] ], {i, Length[subset1]}, {j, Length[reflectionCenters]}], 1];
			subset1 = Select[subset1, #[[1]] <= maxHeight + epsilon &];
			subset1 = Union[subset1, SameTest -> ApproxSamePoint];
			subset1 = Intersection[vertices, subset1, SameTest -> ApproxSamePointLoose];
			Print[{k, Length[subset1]}];
			If[beforeCount == Length[subset1], Break[]];
		];

		subsetIndices = {};
		For[si = 1, si <= Length[subset1], si++,
			AppendTo[subsetIndices, verticesMap[Round[subset1[[si]], epsilon]]];
		];
		subsetIndices = Sort[subsetIndices];
		Print[{seedIndex, k, subset1//Length, subsetIndices//Length}];
		Print[subsetIndices];

		For[subIndex = 1, subIndex <= Length[allSubsets], subIndex++,
			If[Length[Intersection[ allSubsets[[subIndex]], subsetIndices ]] > 0,
				allSubsets[[subIndex]] = Sort[Union[allSubsets[[subIndex]], subsetIndices]];
				Print["Merge to an old component, #" <> IntegerString[subIndex]];
				Break[]
			]
		];

		If[subIndex > Length[allSubsets],
			Print["new component"];
			AppendTo[allSubsets, subsetIndices]
		];

		Print["Discovery progress: component count. before dedup, after dedup (should be same), final target"];
		Print[{Length[allSubsets], Total[Map[Length, allSubsets]], Length[Apply[Union, allSubsets]], Length[vertices]}];
	];

	allSubsets
];


dataFileName = FileNameJoin[{dataFolder, shape <> "_" <> IntegerString[cellThreshold] <> ".wl"}];
If[FileExistsQ[dataFileName], 
	data = Get[dataFileName],
	Print["Cannot find data file "<> dataFileName];
	Exit[]
];

cells = data["cells"];
faces = data["faces"];
edges = data["edges"];
vertices = data["vertices"];
Print["Counts: cells, faces, edges, vertices"];
Print[Map[Length, {cells, faces, edges, vertices}]];

normalizedCellCenters = Map[HNormalize[getCellCenter[#], HNorm[vertices[[1]]]]&, cells ];
reflectionCenters = Join[ vertices, normalizedCellCenters ];

allSubsets = getAllSubsetsThroughCache[normalizedCellCenters, reflectionCenters];

(* Print[allSubsets]; *)
Print["number of subsets:"];
Print[Length[allSubsets]];
Print["cells in all subsets:"];
Print[Map[Length, allSubsets]];

Print["All subsets before dedup, after dedup (should be same), all cell count"];
Print[{Total[Map[Length, allSubsets]], Length[Apply[Union, allSubsets]], Length[cells]}];

(* 
normalizedFaceCenters = Map[HNormalize[Mean[#], HNorm[vertices[[1]]]]&, faces ];
Print[Length[normalizedFaceCenters]];
facesSubsets = getAllSubsets[normalizedFaceCenters, reflectionCenters];
Print["number of face subsets:"];
Print[Length[facesSubsets]];
Print["faces in all subsets:"];
Print[Map[Length, facesSubsets]];

normalizedEdgeCenters = Map[HNormalize[Mean[#], HNorm[vertices[[1]]]]&, edges ];
Print[Length[normalizedEdgeCenters]];
edgesSubsets = getAllSubsets[normalizedEdgeCenters, reflectionCenters];
Print["number of face subsets:"];
Print[Length[edgesSubsets]];
Print["faces in all subsets:"];
Print[Map[Length, edgesSubsets]];
 *)