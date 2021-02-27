#!/usr/bin/env wolframscript

(* nGenerations = 6; *)

shapes = {
	"h_7_3"
(* 	"h_3_7",
	"h_72_7",
	"h_7_72" *)
};


dataFolder = "data";

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

generateFaces[shape_] := Module[
	{},

	nGenerations = If[shape == "h_72_7", 5, If[shape == "h_3_7", 8, 8]];
	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

	configs = Get["configs.wl"];
	mirrors = configs[shape][[1]];
	vertex1 = configs[shape][[2]]//N;

	face1 = {vertex1};
	finished = False;
	For[i=0, i < 10 && finished == False, i++,
		newVertex = FullSimplify[HDoubleReflect[face1[[-1]], mirrors[[1]], mirrors[[2]]]];
		If[!(ApproxSamePoint[newVertex, face1[[1]]]),
			AppendTo[face1, newVertex],
			finished = True;
		]
	];

	Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

	faces = {face1};
	countStart = -1;
	For[i = 0, i <= nGenerations && countStart != Length[faces], i++,
		Print["i = " <> IntegerString[i]];
		countStart = Length[faces];
		For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, mirrors[[mirrorIndex]] ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];

	(* Making it symmetric *)
	centeredMirrors = Select[mirrors, #[[1]] == 0 &];

	countStart = -1;
	For[i = 0, i < 10 && countStart != Length[faces], i++,
		(* Print["i = " <> IntegerString[i]]; *)
		countStart = Length[faces];
		For[mirrorIndex = Length[centeredMirrors], mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, centeredMirrors[[mirrorIndex]] ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];
	faces
]

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
	shape = shapes[[shapeIndex]];
	Print[shape];

	faces = generateFaces[shape];
	Print["There are " <> IntegerString[Length[faces]] <> " faces."];
	
	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
	jsonDataFileName = FileNameJoin[{dataFolder, shape <> ".json"}];

	If[FileExistsQ[dataFileName], DeleteFile[dataFileName]];
	Save[dataFileName, faces];
	Print["Saved faces to " <> dataFileName];

	If[FileExistsQ[jsonDataFileName], DeleteFile[jsonDataFileName]];
	Export[jsonDataFileName, faces, "RawJSON", Compact -> True];
	Print["Saved faces to " <> jsonDataFileName];
]
