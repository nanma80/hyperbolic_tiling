#!/usr/bin/env wolframscript

exportToPov = True;

shapes = {
	"h_7_3"
	(* "h_3_7", *)
	(* "h_72_7", *)
	(* "h_7_72" *)
};

splitEdgeParts = 8;

dataFolder = "data";
imageFolder = "images";

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

Boost[t_] :=
{
  {Cosh[t], Sinh[t], 0},
  {Sinh[t], Cosh[t], 0},
  {0,       0,       1}
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
shrinkFace[face_, factor_] := Map[(factor) # + (1 - factor) Total[face]/Length[face] &, face];
shiftFace[face_, factor_] := Map[# + {Total[face][[1]], 0, 0} * factor &, face];
splitFace[face_, n_] := Module[
	{},
	edgesFace1 = Reverse[getEdgesFromFace[face]];
	splitEdges1 = Map[#[[2]]&, Map[splitEdge[#, n]&, edgesFace1], {2}];
	Catenate[splitEdges1]
];



plotShape3D[shape_, faces_, panning_] := Module[
	{allEdges, edges, allVertices, vertices},
	configs = Get["configs.wl"];
	color = configs[shape][[3]];
	norm = HNorm[faces[[1]][[1]]];

	pannedFaces = Map[Boost[panning].# &, faces, {2}];

	displayFaces = Map[getKlein, pannedFaces, {2}];

	Print["There are " <> IntegerString[Length[displayFaces]] <> " faces."];

	displayEdges = Flatten[Map[getEdgesFromFace, displayFaces], 1];
	Print["There are " <> IntegerString[Length[displayEdges]] <> " edges in total."];

	vertices = Flatten[displayFaces,1];
	Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."];

	thickness = 0.01;
	lighting = {{"Point", White, {10, -10, 10}}};
	range = 0.8;

	graph = Graphics3D[
		{
			{
				Opacity[0.8], FaceForm[color], Map[Polygon, displayFaces]
			},
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

exportFrame[shape_, panning_, frameIndex_] := Module[
	{},
	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
	Print["Reading data for " <> shape <> " from " <> dataFileName];

	faces = Get[dataFileName];

	graph = plotShape3D[shape, faces, panning];

	outputFolder = shape <> "_panning_klein_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolder}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "panning_klein_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFolder, outputFileName}];
	imageSize = {600, 450};


	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
	];

];

exportAnimation[shape_] := Module[
	{
		outputFolder,
		gifFileName,
		gifFolder,
		maxExplosionFactor,
		frameCount,
		explodeStep,
		explosionFactors,
		frameIndex,
		explosionFactor,
		gifCommand
	},

	outputFolder = shape <> "_panning_klein_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolder}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	gifFileName = shape <> "_panning_klein";
	gifFolder = FileNameJoin[{imageFolder, shape}];
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];


	cosine = Cos[2 * Pi / 7];
	ch2phi = (8/3 * Cos[Pi / 7]^2) - 1;
	rsquare = (1 - cosine)/(ch2phi - 1);
	x0square = rsquare + 1;
	r = Sqrt[rsquare];
	x0 = FullSimplify[Sqrt[x0square]];
	maxPanningParameters = ArcCoth[x0];
	
	(* maxPanningParameters = 0.5; *)
	frameCount = 30;
	step = maxPanningParameters/frameCount;

	epsilon = 0;
	panningParameters = Table[k, {k, epsilon, maxPanningParameters + epsilon, step}];
	For[ frameIndex = 1, frameIndex <= Length[panningParameters] , frameIndex++,
		Print["Frame #" <> IntegerString[frameIndex] <> "/" <> IntegerString[Length[panningParameters]]];
		panningParameter = panningParameters[[frameIndex]];
		Print[panningParameter];
		exportFrame[shape, panningParameter, frameIndex];
	];

	gifCommand = "magick convert -duplicate 1,-2-1 -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[Floor[150/frameCount]] <> " -dispose 2 " <> outputFolder <> "\\panning_klein_frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

	Print[gifCommand];
	Run[gifCommand];
	Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];
]


For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
	shape = shapes[[shapeIndex]];
	Print[shape];
	exportAnimation[shape];
]


(* exportFrame["h_7_3", maxPanningParameters, 1]; *)
