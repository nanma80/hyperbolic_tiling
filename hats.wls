#!/usr/bin/env wolframscript

faceHeightThreshold = 80;

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 6; q = 4; (* default *)
];

frameCount = 200;
(* frameCount = 5; *)

color = White;
diskColor = Blue;

seedTheta = Pi/12 * 1.0;
imageSize = {4, 3} * 720/3;
shape = "hat";
dataFolder = "data";
imageFolder = "images";
epsilon = 0.00000001;
exportToPov = True;
splitEdgeParts = 8;
nGenerations = 30;
shapeWithParams = shape <> "_" <> IntegerString[faceHeightThreshold] <> "_" <> IntegerString[frameCount];

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
Rotation[t_] :=
{
  {1, 0, 0},
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

getThickness[originalThickness_, v3d_] := originalThickness * ( 1 - Norm[v3d[[{2, 3} ]]] )^0.3;

Boost[t_] :=
{
  {Cosh[t], Sinh[t], 0},
  {Sinh[t], Cosh[t], 0},
  {0,       0,       1}
};

generateFaces[p_, q_] := Module[
	{},

	mirror1 = {0, 0, 1};
	mirror2 = {0, -Sqrt[3]/2, 1/2};

	mirror3 = {Sqrt[rho3 ^ 2 + 1] - 1, rho3, 0};
	mirror4 = {Sqrt[rho4 ^ 2 + 1] - 1, rho4/2, rho4/2 * Sqrt[3]};
	vertex1 = {Sqrt[r ^ 2 + 1], r Cos[seedTheta], r Sin[seedTheta]};

	solutions = NSolve[ 
		HReflect[vertex1, mirror3] == vertex1 
		&& HReflect[vertex1, mirror4] == vertex1 
		&& HInner[mirror3, mirror4] == 0 
		&& r > 0 && rho3 > 0 && rho4 > 0, {r, rho3, rho4}];

	solution = solutions[[1]];

	mirror3 = mirror3/.solution;
	mirror4 = mirror4/.solution;
	vertex1 = vertex1/.solution;

	mirrors = { mirror1, mirror2, mirror3, mirror4 }//N;

	rotations = {
		mirrors[[ {1, 2} ]],
		mirrors[[ {1, 3} ]],
		mirrors[[ {2, 4} ]]
	};

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

	(* construct hat *)
	origin = {1, 0, 0};
	mid34 = HNormalize[face1[[3]] + face1[[4]], 1];
	mid23 = HNormalize[face1[[2]] + face1[[3]], 1];
	mid45 = HNormalize[face1[[4]] + face1[[5]], 1];
	mid445 = HNormalize[face1[[4]] + mid45, 1];
	mid323 = HNormalize[face1[[3]] + mid23, 1];
	mid223 = HNormalize[face1[[2]] + mid23, 1];
	mid445323 = HNormalize[mid445 + mid323, 1];
	triangleCenter = HNormalize[origin + mid323 + mid445, 1];

	mid34Right = HDoubleReflect[ mid34, mirror4, mirror1];
	mid445323Right = HDoubleReflect[ mid445323, mirror4, mirror1];

	mid445o = HNormalize[mid445 + origin , 1];
	mid223o = HNormalize[mid223 + origin , 1];

	point45o = HNormalize[ ratio * mid45 + (1 - ratio) * origin, 1 ];
	solutionRatio = NSolve[ HInner[point45o, origin] == HInner[triangleCenter, origin] && ratio > 0 && ratio < 1, ratio ][[1]];
	point45o = point45o/.solutionRatio;
	point23o = HNormalize[ ratio * mid23 + (1 - ratio) * origin, 1 ]/.solutionRatio;

(* 	Print[ HInner[triangleCenter, origin] ];
	Print[ HInner[triangleCenter, mid323] ];
	Print[ HInner[triangleCenter, mid445] ];
 *)
	hat = {
		origin, mid223o, point23o, mid23, mid323, mid445323Right, mid34Right, mid34, mid445323, mid445, mid45, point45o, mid445o
	};	
	hat = Map[HNormalize[#, 1]&, hat];

	Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];
	Print["Each hat contains " <> IntegerString[Length[hat]] <> " vertices."];

	origin2 = HDoubleReflect[HReflect[origin, mirror4], mirror1, mirror2];
	Print[origin2];
	maxPanning = ArcCosh[ origin2[[1]] ];
	Print[maxPanning];

	faces = {hat};
	countStart = -1;
	For[i = 1, i <= nGenerations && countStart != Length[faces], i++,
		Print[DateString["ISODateTime"] <> " i = " <> IntegerString[i] <> " face count: " <> IntegerString[Length[faces]]];
		countStart = Length[faces];
		For[rid = Length[rotations], rid >= 1, rid--,
			newFaces = Map[HDoubleReflect[#, rotations[[rid]][[1]], rotations[[rid]][[2]]]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
			faces = Select[faces, Mean[#][[1]] < faceHeightThreshold &];
		];
	];

	Print["Making it more symmetric"];

	countStart = -1;
	For[i = 0, i < 20 && countStart != Length[faces], i++,
		Print[DateString["ISODateTime"] <> " i = " <> IntegerString[i] <> " face count: " <> IntegerString[Length[faces]]];
		countStart = Length[faces];
		For[rid = 1, rid >= 1, rid--,
			newFaces = Map[HDoubleReflect[#, rotations[[rid]][[1]], rotations[[rid]][[2]]]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
			faces = Select[faces, Mean[#][[1]] < faceHeightThreshold &];
		];
	];

	newFaces = Map[HDoubleReflect[#, rotations[[3]][[1]], rotations[[3]][[2]]]&, faces, {2}];
	newFaces = Map[HDoubleReflect[#, rotations[[1]][[1]], rotations[[1]][[2]]]&, newFaces, {2}];
	faces = Union[faces, newFaces, SameTest -> sameCenter];

	faces
];


plotShape3D[shape_, faces_, panning_] := Module[
	{allEdges, edges, allVertices, vertices},
	pannedFaces = Map[Rotation[Pi/3 * panning].Boost[panning * maxPanning].# &, faces, {2}];
	displayFaces = Map[getPoincare, pannedFaces, {2}];

	Print["There are " <> IntegerString[Length[displayFaces]] <> " faces."];

	edges = Flatten[Map[getEdgesFromFace, pannedFaces], 1];
	Print["There are " <> IntegerString[Length[edges]] <> " edges before split."];

	norm = HNorm[edges[[1]][[1]]];

	splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1];
	splitEdges = Map[Re[HNormalize[#, norm]]&, splitEdges, {2}];
	displayEdges = Map[getPoincare, splitEdges, {2}];
	Print["There are " <> IntegerString[Length[displayEdges]] <> " edges after split."];

	vertices = Flatten[displayFaces,1];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	thickness = 0.008;
	lighting = {{"Point", White, {10, -10, 10}}};
	range = 0.8;

	graph = Graphics3D[
		{
			{
				Opacity[1], FaceForm[diskColor], Cylinder[{{-0.001, 0, 0}, {0, 0, 0}}, 1 - 0.01]
			},
 			{
				Opacity[1], FaceForm[color], Map[Cylinder[#, getThickness[thickness, Mean[#]] ]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[color], Map[Sphere[#, getThickness[thickness, #] * 2 ]&, vertices]
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

	outputFolder = shapeWithParams <> "_panning_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolder}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "panning_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFolder, outputFileName}];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];
];

exportAnimation[] := Module[
	{
		outputFolder,
		gifFileName,
		gifFolder,
		maxExplosionFactor,
		explodeStep,
		explosionFactors,
		frameIndex,
		explosionFactor,
		gifCommand
	},

	outputFolder = shapeWithParams <> "_panning_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolder}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	gifFileName = shapeWithParams;
	gifFolder = FileNameJoin[{imageFolder, shape}];
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	Print["Shape: " <> shape];
	faces = generateFaces[p, q];
	Print["There are " <> IntegerString[Length[faces]] <> " faces."];

	For[frameIndex = 1, frameIndex <= frameCount, frameIndex++,
		Print[DateString["ISODateTime"] <> " frameIndex/frameCount: " <> IntegerString[frameIndex] <> " / " <> IntegerString[frameCount] ];
		exportFrame[faces, frameIndex/frameCount*1.0, frameIndex];
	];

	delay = If[frameCount >= 40, 3, Floor[300/frameCount]];

	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[delay] <> " -dispose 2 " <> outputFolder <> "\\panning_frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

	Print[gifCommand];
	Run[gifCommand];
	Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];
	
  mp4Command = "ffmpeg -i " <> gifFolder <> "\\" <> gifFileName <> ".gif" 
      (* <> " -movflags faststart -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" " *)
      <> " -movflags faststart -pix_fmt yuv420p "
      <> gifFolder <> "\\" <> gifFileName <> ".mp4";

  mp4FileName = gifFolder <> "\\" <> gifFileName <> ".mp4";
  If[FileExistsQ[mp4FileName], DeleteFile[mp4FileName]];
  Print[mp4Command];
  Run[mp4Command];
  Print["Saved the video to " <> gifFolder <> "\\" <> gifFileName <> ".mp4"];	
];

exportAnimation[];

(* For[fid = 0, fid <= 10, fid++,
	exportFrame[faces, fid/10, fid];
];
 *)

(* Print["Shape: " <> shape];
faces = generateFaces[p, q];
Print["There are " <> IntegerString[Length[faces]] <> " faces."];
 *)

(* exportFrame[faces, 0, 1]; *)
(* exportFrame[faces, 0.5, 2]; *)
(* exportFrame[faces, 1, 3]; *)

