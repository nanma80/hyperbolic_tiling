#!/usr/bin/env wolframscript

exportToPov = True;
highlightCenter = True;
topView = True;
plotPoints = 40;
maxRecursion = 2;
sphere = False;
doubleSheet = False;
If[sphere, doubleSheet = True];
mirror3Generation = 3;

splitEdgeParts = 8;

dataFolder = "data";
cosine = Cos[2 * Pi / 7];
ch2phi = (8/3 * Cos[Pi / 7]^2) - 1;
rsquare = (1 - cosine)/(ch2phi - 1);
x0square = rsquare + 1;
r = Sqrt[rsquare];
x0 = FullSimplify[Sqrt[x0square]];
epsilon = 0.0000000001;

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

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001;

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
getHyperboloid[v_] := {-v[[1]], v[[2]], v[[3]]};


mirrors = {
	{0, 0, 1},
	{0, -Sin[Pi/7], Cos[Pi/7]},
	{Sqrt[3]/2, Sqrt[3/4 + rsquare], 0}
};

allMirrors = mirrors[[{1, 2, 3}]];
For[j = 1, j <= mirror3Generation, j++,
	newMirrors3 = Map[HReflect[#, mirrors[[3]]]&, allMirrors];
	updatedMirrors = Union[allMirrors, newMirrors3, SameTest -> (Round[HNormalize[#1, 1], epsilon] == Round[HNormalize[#2, 1], epsilon]&)];
	Print[Map[Length, {allMirrors, newMirrors3, updatedMirrors}]];
	allMirrors = updatedMirrors//Simplify;

	For[i = 1, i <= 7, i++,
		newMirrors1 = Map[HReflect[#, mirrors[[1]]]&, allMirrors];
		newMirrors2 = Map[HReflect[#, mirrors[[2]]]&, allMirrors];
		updatedMirrors = Union[allMirrors, newMirrors1, newMirrors2, SameTest -> (Round[HNormalize[#1, 1], epsilon] == Round[HNormalize[#2, 1], epsilon]&)];
		Print[Map[Length, {allMirrors, newMirrors1, newMirrors2, updatedMirrors}]];
		allMirrors = updatedMirrors//Simplify;
	];
];

allMirrors = Union[allMirrors, (-1)*allMirrors, SameTest -> (Round[HNormalize[#1, 1], epsilon] == Round[HNormalize[#2, 1], epsilon]&)];
Print[Length[allMirrors]];
arbitraryVector = {1, E Pi, E / Pi};
allMirrors = Select[allMirrors, #.arbitraryVector >= 0 &];
Print[Length[allMirrors]];


plotRange = If[sphere, 1.1, 2];
opacity = 1;
zOffset = If[sphere, 0, If[doubleSheet, -plotRange - 0.9, 0]];

graph = ContourPlot3D[
	If[sphere, x^2 + y^2 + z^2, -x^2 - y^2 + z^2] == If[doubleSheet, 1, -1],
	{x, -plotRange, plotRange}, 
	{y, -plotRange, plotRange}, 
	{z, -plotRange + zOffset, plotRange + zOffset},
	Axes -> None, 
	BoundaryStyle -> Red, 
  PlotPoints -> plotPoints, 
  MaxRecursion -> maxRecursion,
	MeshFunctions -> {0 &}, 
	MeshShading -> { { Lighter[Blue] , Opacity[opacity] } },
	(* RegionBoundaryStyle -> None, *)
	RegionFunction -> Function[
                      {x, y, z},
                      Product[HInner[{z, x, y}, allMirrors[[k]]], {k, 1, Length[allMirrors]}] < 0
                  	]
];

outputFileName = If[sphere, "sphere_", If[doubleSheet, "double_", "single_"]] <> "fundamental_domain_" <> IntegerString[mirror3Generation] <> "_" <> IntegerString[plotPoints] <> "_" <> IntegerString[maxRecursion] <> ".png";
imageSize = {600, 450};

imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, "domains"}];
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
