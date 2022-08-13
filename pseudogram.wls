#!/usr/bin/env wolframscript

exportToPov = True;

frameCount = 200;
sideCount = 13;

unit = 0.4;

epsilon = 0.00000001;
Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

imageSize = {4, 3}/3 * 720;
imageFolder = "images";
outputFolder = FileNameJoin[{imageFolder, "pseudogram"}];
axisColor = Darker[White, 0.5];
plotPoints = 20;
maxRecursion = 1;

If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

shapeWithParams = "pseudogram_" <> IntegerString[sideCount];

RoundEpsilon[v_] := Round[v, epsilon];
HInner[v_, u_] := 2 * v[[1]] * u[[1]] - Dot[v, u];
HNormSquare[v_] := HInner[v, v];
HNorm[v_] := Sqrt[HNormSquare[v]];
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], Sin[t]},
  {0, -Sin[t], Cos[t]}
};

Boost[t_] :=
{
  {Cosh[t], Sinh[t], 0},
  {Sinh[t], Cosh[t], 0},
  {0,       0,       1}
};

HReflect[point_, mirror_] := point - 2 * HInner[point, mirror]/HInner[mirror, mirror]*mirror;
HDoubleReflect[point_, mirror1_, mirror2_] := HReflect[HReflect[point, mirror1], mirror2];
ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001;
sameEdge[e1_, e2_] := sameCenter[e1, e2] && ((ApproxSamePoint[e1[[1]], e2[[1]]] && ApproxSamePoint[e1[[2]], e2[[2]]] ) || (ApproxSamePoint[e1[[1]], e2[[2]]] && ApproxSamePoint[e1[[2]], e2[[1]]] ));
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
getHyperboloid[v_] := {-v[[1]], v[[2]], v[[3]]};
getColor[center_] := Hue[ Mod[Arg[ center[[2]] + I center[[3]] ], 2 Pi]/2/Pi, (Tanh[center[[1]]/10 + 0.5 ] + 1)/2, 1];
(* getColor[center_] := Red; *)

scale[edge_, ratio_] := Module[
	{center},
	center = Mean[edge];
	Map[(# - center) * ratio + center&, edge]
];

scaleFace[face_, ratio_, stellationLevel_] := Module[
	{edges, edgeOrders},
	edges = Map[scale[#, ratio]&, getEdgesFromFace[face]];
	edgeOrders = Table[Mod[stellationLevel * k, 7]+1, {k, 0, 6}];
	Map[#[[1]]&, edges[[edgeOrders]]]
];

plotShape3D[progress_, frameIndex_] := Module[
	{displayFaces, allEdges, edges, allVertices, vertices},

	startDelay = 0.5;
	transitionDuration = 0.3;
	lighting = {{"Point", White, {50, -50, 50}}};
	range = 3;
	thickness = 0.01 * range;

	ks = Table[k, {k, -sideCount/2, sideCount/2, 1}];
	coefs = Table[k, {k, -sideCount/2, sideCount/2, 1}];
	absCoefs = Map[Abs[#+1/4]&, coefs];
	absCoefs = Map[#/Max[absCoefs] * startDelay &, absCoefs];
	(* Print[absCoefs//N]; *)

	vertices = Map[{Cosh[unit #], Sinh[unit #], 0}&, coefs];
	For[k = 1, k <= Length[vertices], k++,
		absCoef = absCoefs[[k]];
		oldVertex = vertices[[k]];
		If[Mod[k, 2] == 0,
			progressPerPoint = If[progress < absCoef, 
				0, 
				If[progress < absCoef + transitionDuration,
					(progress - absCoef)/transitionDuration,
					1
				]
			];
			
			newVertex = {Cos[progressPerPoint Pi] oldVertex[[1]], Cos[progressPerPoint Pi] oldVertex[[2]], Sin[progressPerPoint Pi] Norm[oldVertex[[{1,2}]]] - 0.001 Abs[oldVertex[[2]]] }
			,
			newVertex = oldVertex - {0, 0, 0.001 Abs[oldVertex[[2]] ]}
		];
		vertices[[k]] = newVertex;
	];
	edges = getEdgesFromFace[vertices];
	(* remove the closing edge *)
	edges = edges[[Range[Length[edges] - 1 ]]];

	elements = {};

	elements = Union[
		elements,
		{
			{FaceForm[axisColor], Cylinder[{{-2 range, 0, -1}, {2 range, 0, -1}}, thickness/2]},
			{FaceForm[axisColor], Cylinder[{{0, -2 range, -1}, {0, 2 range, -1}}, thickness/2]}
		}
	];

	elements = Union[
		elements,
		Table[
			{
				FaceForm[ If[Mod[k, 2] == 0, Red, Green] ],
				Sphere[vertices[[k]], thickness * 2]
			}, 
			{k, 1, Length[vertices]}
		]
	];

	elements = Union[
		elements,
		Table[
		{
				Opacity[1],
				FaceForm[Blue], 
				Cylinder[
					edges[[k]], 
					thickness
				]
		},
		{k, 1, Length[edges]}]
	];

	graph = Graphics3D[
		elements,
		ViewPoint -> {0, 0, 1000},
		ViewVertical -> {0, 1, 0},
		PlotRange -> {{range, -range}, {range, -range}, {range, -range}},
		Lighting -> lighting
	];

	background = ContourPlot3D[
		z == 2 range,
		{x, -2 range, 2 range}, 
		{y, -2 range, 2 range}, 
		{z, -2 range, 2 range},
		Axes -> None,
		BoundaryStyle -> None,
		MeshFunctions -> {0 &}, 
		MeshShading -> { { Lighter[Blue] , Opacity[0.3] } },
		RegionFunction -> Function[{x, y, z}, x^2 - y^2 > 1 * 1.1] 
	];

  graph = Show[graph, background];
	graph
];

exportFrame[progress_, frameIndex_, outputFolderName_] := Module[
	{outputFrameFolder},
	Print[DateString["ISODateTime"] <> " Generating frame for:"];
	Print[{frameIndex, progress}];
	(* Return[]; *)
	(* If[frameIndex <= 175, Return[]]; *)
	(* If[Mod[frameIndex, 2] != frameIndexMod, Return[]]; *)
	outputFrameFolder = FileNameJoin[{outputFolder, outputFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	outputFileName = "frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFrameFolder, outputFileName}];

	graph = plotShape3D[progress, frameIndex];


	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];

	Print["Saved image to " <> outputFileName];
];

exportAnimation[] := Module[
	{
		gifFileName,
		gifFolder,
		frameIndex,
		explosionFactor,
		gifCommand
	},

	outputFrameFolderName = shapeWithParams <> "_frames";

	outputFrameFolder = FileNameJoin[{outputFolder, outputFrameFolderName}];
	If[!DirectoryQ[outputFrameFolder], CreateDirectory[outputFrameFolder]];

	gifFileName = shapeWithParams;
	gifFolder = outputFolder;
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	progressIds = 
		Join[
			Table[0, {k, Ceiling[frameCount/6] }],
			Table[k, {k, 0, 1, 1/(frameCount - 1)}],
			Table[1, {k, Ceiling[frameCount/6] }]
		];
	Print[N[progressIds]];

	For[frameIndex = 1, frameIndex <= Length[progressIds], frameIndex++,
		progress = progressIds[[frameIndex]];
		Print["exportFrame"];
		Print[{progress, frameIndex, frameCount}];

		exportFrame[progress, frameIndex, outputFrameFolderName];
	];
	
	delay = If[frameCount >= 50, 3, Floor[500/frameCount]];

	(* gifCommand = "magick convert -duplicate 1,-2-1 -quiet +repage -alpha remove -loop 0 -delay " *)
	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " 
		<> IntegerString[delay] <> " -dispose 2 " <> outputFrameFolder <> "\\frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

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

(* exportFrame[0, 1, "test"]; *)
(* exportFrame[1/4, 2, "test"];
exportFrame[2/4, 3, "test"];
exportFrame[3/4, 4, "test"];
exportFrame[1, 5, "test"];
 *)