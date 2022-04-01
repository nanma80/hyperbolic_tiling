#!/usr/bin/env wolframscript

models = {"hyperboloid", "klein", "poincare"};
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
(* models = {"poincare"}; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 3, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]]; r = ToExpression[argv[[3]]],
	p = 5; q = 3; r = 4;
];

cellThreshold = 100;

Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"];
shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm];
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
imageSize = {4, 3}/3 * 720;
lighting = {{"Point", White, {10, -10, 10} * 5}};
rangeFactor = 0.8;

outputFolder = FileNameJoin[{imageFolder, "honeycombs"}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
exportToPov = True;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

phi = (1 + Sqrt[5])/2;
HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], Sin[t]},
  {0, -Sin[t], Cos[t]}
};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

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
	{5,3,5} -> Yellow
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

dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

getCells[p_, q_, r_] := Module[
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

	face1 = {vertex1}//N;
	finished = False;
	For[i = 0, i < 20 && finished == False, i++,
		newVertex = FullSimplify[HDoubleReflect[face1[[-1]], mirrors[[1]], mirrors[[2]]]];
		If[!(ApproxSamePoint[newVertex, face1[[1]]]),
			AppendTo[face1, newVertex],
			finished = True;
		]
	];

	Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

	cell1 = {face1};
	countStart = -1;
	For[i=0, i < 100 && countStart != Length[cell1], i++,
		countStart = Length[cell1];
		
		For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
			cell1 = Union[cell1, Map[HReflect[#, mirrors[[mirrorIndex]] ]&, cell1, {2}], SameTest -> sameCenter];
		];
	];
	Print["Each cell contains " <> IntegerString[Length[cell1]] <> " faces."];

	cells = {cell1};

	For[gIndex = 1, gIndex <= 6 && Length[cells] < cellThreshold, gIndex++,
		cells = Union[cells, Map[HReflect[#, mirrors[[1]] ]&, cells, {3}], SameTest -> sameCellCenter];

		countStart = -1;
		For[i = 0, i < 10 && countStart != Length[cells], i++,
			countStart = Length[cells];
			
			For[mirrorIndex = 4, mirrorIndex >= 2, mirrorIndex--,
				cells = Union[cells, Map[HReflect[#, mirrors[[mirrorIndex]] ]&, cells, {3}], SameTest -> sameCellCenter];
			];

			Print["gIndex = "  <> IntegerString[gIndex] <> "; i = "  <> IntegerString[i] <> ". Generated " <> IntegerString[Length[cells]] <> " cells so far"];
		];
	];

	Print["Generated " <> IntegerString[Length[cells]] <> " cells in total."];
	cells
];


cells = getCells[p, q, r];
cellCenterHeights = Map[getCellCenter[#][[1]]&, cells];
heightTally = Tally[cellCenterHeights, ApproxSamePoint]//Sort;
tallyCounts = Map[#[[2]]&, heightTally];

Print[heightTally];
Print[tallyCounts];

faces = Join[ Flatten[cells, 1]];
edges = Join[ Flatten[Map[getEdgesFromFace, faces], 1] ];
vertices = Join[ Flatten[faces,1]];

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

	range = Max[Map[Norm[projectTo3d[#]]&, vertices]] * rangeFactor;
	thickness = range * 0.008;

	color = colorMap[{p,q,r}];

	graph = Graphics3D[
		{
			(* {
				Opacity[0.2], FaceForm[color], Map[Polygon, Map[projectTo3d, faces, {2}]]
			}, *)
			{
				Opacity[1], FaceForm[color], Map[Cylinder[#, thickness]&, Map[projectTo3d, edges, {2}]]
			},
			{
				Opacity[1], FaceForm[color], Map[Sphere[#, thickness * 2]&, Map[projectTo3d, vertices]]
			}
		},
		ViewPoint -> {3, -4, 5}*30,
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

