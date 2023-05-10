#!/usr/bin/env wolframscript

nGenerations = 6;

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 6; q = 4; (* default *)
];

seedTheta = Pi/12;

shape = "hat";
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
exportToPov = True;
splitEdgeParts = 8;

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
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;

Boost[t_] :=
{
  {Cosh[t], Sinh[t], 0},
  {Sinh[t], Cosh[t], 0},
  {0,       0,       1}
};

generateFaces[p_, q_] := Module[
	{},
(* 
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
 *)
	(* new mirrors for hats tiling *)

	mirror1 = {0, 0, 1};
	mirror2old = {0, 1/2, Sqrt[3]/2};
	mirror2 = {0, -Sqrt[3]/2, 1/2};

	mirror3 = {Sqrt[rho ^ 2 + 1] - 1, rho, 0};
	mirror3 = mirror3 / HNorm[mirror3] * I;

	solutionRho = Solve[HInner[mirror3, mirror2old] == -Sqrt[2]/2, rho][[1]];
	mirror3 = mirror3/.solutionRho;
	(* Print[mirror3]; *)
	(* {1, Sqrt[2], 0} *)

	vertex1 = {Sqrt[r ^ 2 + 1], r Cos[seedTheta], r Sin[seedTheta]};

	solutionR = Solve[ HReflect[vertex1, mirror3] == vertex1, r][[1]];
	vertex1 = Simplify[vertex1/.solutionR] ;

	mirrors = { mirror1, mirror2, mirror3 }//N;

(* 	Print[mirrors];
	Print[vertex1];
 *)

	face1 = {vertex1}//N;
	finished = False;
	For[i = 0, i < 20 && finished == False, i++,
		For[mirrorIndex = 2, mirrorIndex >= 1, mirrorIndex--,
			newVertex = FullSimplify[HReflect[face1[[-1]], mirrors[[mirrorIndex]] ]];
			If[!(ApproxSamePoint[newVertex, face1[[1]]]),
				AppendTo[face1, newVertex],
				finished = True;
			];
		];
	];

	args = Map[Arg[#[[2]] + #[[3]] I] &, face1];
	face1 = face1[[Ordering[args]]];

	Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

	faces = {face1};
	countStart = -1;
	For[i = 1, i <= nGenerations && countStart != Length[faces], i++,
		Print["i = " <> IntegerString[i]];
		countStart = Length[faces];
		For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, mirrors[[mirrorIndex]] ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];

	Print["Making it more symmetric"];

	(* Making it symmetric *)
	centeredMirrors = Select[mirrors, #[[1]] == 0 &];

	countStart = -1;
	For[i = 0, i < 20 && countStart != Length[faces], i++,
		Print["i = " <> IntegerString[i]];
		countStart = Length[faces];
		For[mirrorIndex = Length[centeredMirrors], mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, centeredMirrors[[mirrorIndex]] ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];
	faces
];


plotShape3D[shape_, faces_, panning_] := Module[
	{allEdges, edges, allVertices, vertices},
	color = Blue;
	pannedFaces = Map[Boost[panning].# &, faces, {2}];
	displayFaces = Map[getPoincare, pannedFaces, {2}];

	Print["There are " <> IntegerString[Length[displayFaces]] <> " faces."];

	edges = Flatten[Map[getEdgesFromFace, faces], 1];
	Print["There are " <> IntegerString[Length[edges]] <> " edges before split."];

	norm = HNorm[edges[[1]][[1]]];

	splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1];
	splitEdges = Map[Re[HNormalize[#, norm]]&, splitEdges, {2}];
	displayEdges = Map[getPoincare, splitEdges, {2}];
	Print["There are " <> IntegerString[Length[displayEdges]] <> " edges after split."];

	vertices = Flatten[displayFaces,1];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	thickness = 0.006;
	lighting = {{"Point", White, {10, -10, 10}}};
	range = 0.8;

	graph = Graphics3D[
		{
(* 			{
				Opacity[0.8], FaceForm[color], Map[Polygon, displayFaces]
			}, *)
 			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[White], Map[Sphere[#, thickness * 2]&, vertices]
			}
		},
		ViewPoint -> {100, 0, 0},
		ViewVertical -> {0, 0, 1},
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
		Lighting -> lighting
	];
	graph
];

exportFrame[faces_, panning_, frameIndex_] := Module[
	{},

	graph = plotShape3D[shape, faces, panning];

	outputFolder = shape <> "_panning_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolder}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "panning_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFolder, outputFileName}];
	imageSize = {4, 3} * 720/4;


	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];
];


Print["Shape: " <> shape];
faces = generateFaces[p, q];
Print["There are " <> IntegerString[Length[faces]] <> " faces."];

exportFrame[faces, 0, 1];


(* dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
jsonDataFileName = FileNameJoin[{dataFolder, shape <> ".json"}];

If[FileExistsQ[dataFileName], DeleteFile[dataFileName]];
Save[dataFileName, faces];
Print["Saved faces to " <> dataFileName];

If[FileExistsQ[jsonDataFileName], DeleteFile[jsonDataFileName]];
Export[jsonDataFileName, faces, "RawJSON", Compact -> True];
Print["Saved faces to " <> jsonDataFileName];
 *)

