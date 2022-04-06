#!/usr/bin/env wolframscript

(* Objective: visualize 2{5, 3, 4}[5{4, 3, 5}] *)

(* models = {"hyperboloid", "klein", "poincare"}; *)
(* models = {"klein", "poincare"}; *)
(* models = {"hyperboloid"}; *)
(* models = {"klein"}; *)
models = {"poincare"};


cellThreshold = 100;
cellThreshold = 300;
(* cellThreshold = 1000; *)

plotTwoEdges = False;
(* plotTwoEdges = True; *)

highlightCycle = False;
(* highlightCycle = True; *)
(* cycleLength = 6; *)

(* Print["{" <> IntegerString[p] <> ", " <> IntegerString[q] <> ", "<> IntegerString[r] <> "}"]; *)
(* shape = "hyper_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm] <> "_" <> ToString[r, InputForm]; *)
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


dataFileName = FileNameJoin[{dataFolder, "hyper_5_3_4" <> "_" <> IntegerString[1000] <> ".wl"}];
If[FileExistsQ[dataFileName], 
	data = Get[dataFileName],
	Print["Cannot find data file "<> dataFileName];
	Exit[]
];

cells534 = data["cells"];
faces534 = data["faces"];
edges534 = data["edges"];
vertices534 = data["vertices"];
Print["534 Counts: cells, faces, edges, vertices"];
Print[Map[Length, {cells534, faces534, edges534, vertices534}]];

maxHeight534 = Max[Map[#[[1]]&, vertices534]];
Print[maxHeight534];


dataFileName = FileNameJoin[{dataFolder, "hyper_4_3_5" <> "_" <> IntegerString[cellThreshold] <> ".wl"}];
If[FileExistsQ[dataFileName], 
	data = Get[dataFileName],
	Print["Cannot find data file "<> dataFileName];
	Exit[]
];

cells435 = data["cells"];
faces435 = data["faces"];
edges435 = data["edges"];
vertices435 = data["vertices"];
Print["435 Counts: cells, faces, edges, vertices"];
Print[Map[Length, {cells435, faces435, edges435, vertices435}]];

maxHeight435 = Max[Map[#[[1]]&, vertices435]];
Print[maxHeight435];

only435 = Complement[vertices435, vertices534, SameTest -> ApproxSamePoint];
Print["only in 435"];
Print[Length[only435]];

Print[Min[Map[#[[1]]&, only435]]];
Print[Max[Map[#[[1]]&, only435]]];



intVertices = Intersection[vertices435, vertices534, SameTest -> ApproxSamePoint];
Print[Length[intVertices]];



Exit[];

normalizedCellCenters = Map[HNormalize[getCellCenter[#], HNorm[vertices[[1]]]]&, cells ];
reflectionCenters = Join[ vertices, normalizedCellCenters ];

allSubsets = getAllSubsetsThroughCache[vertices, reflectionCenters];

(* Print[allSubsets]; *)
Print["number of components:"];
Print[Length[allSubsets]];
Print["vertices in all components:"];
Print[Map[Length, allSubsets]];

Print["All components before dedup, after dedup (should be same), all vertex count"];
Print[{Total[Map[Length, allSubsets]], Length[Apply[Union, allSubsets]], Length[vertices]}];

component = vertices[[allSubsets[[componentIndex]]]];

edgeIndices = getEdges[component, -1];
(* no triangles *)
(* triangles = getTriangles[subset]; *)
Print["edge count"];
Print[edgeIndices//Length];
inner1 = Apply[HInner, component[[ edgeIndices[[1]] ]]];
Print[inner1];

edgeIndices2 = getEdges[component, inner1];

Print["secondary edge count"];
Print[edgeIndices2//Length];
inner2 = Apply[HInner, component[[ edgeIndices2[[1]] ]]];
Print[inner2];

g12 = Graph[Join[edgeIndices, edgeIndices2]];
cycle1 = FindCycle[g12, cycleLength];
Print["Found " <> IntegerString[Length[cycle1]] <> " cycles"];

If[Length[cycle1] > 0,
	Print["processing cycle"];
	Print[cycle1[[1]]];
	cycle1Indices = Map[#[[1]] &, cycle1[[1]]];
	cyclePoints = component[[cycle1Indices]];
	cycleEdges = getEdgesFromFace[cyclePoints],
	Print["no cycle detected"]
];

(* Print["edge inner"]; *)
(* Print[] *)

componentEdges = Map[component[[#]]&, edgeIndices];
componentEdges2 = Map[component[[#]]&, edgeIndices2];

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
		Print["There are " <> IntegerString[Length[componentEdges]] <> " component edges before split."];
		norm = HNorm[componentEdges[[1]][[1]]];
		splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, componentEdges], 1];
		splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
		Print["There are " <> IntegerString[Length[splitEdges]] <> " component edges after split."];
		componentEdges = splitEdges;

		Print["There are " <> IntegerString[Length[componentEdges2]] <> " secondary component edges before split."];
		norm = HNorm[componentEdges2[[1]][[1]]];
		splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, componentEdges2], 1];
		splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
		Print["There are " <> IntegerString[Length[splitEdges]] <> " secondary component edges after split."];
		componentEdges2 = splitEdges;

		If[Length[cycle1] > 0,
			Print["There are " <> IntegerString[Length[cycleEdges]] <> " cycle edges before split."];
			norm = HNorm[cycleEdges[[1]][[1]]];
			splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, cycleEdges], 1];
			splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}];
			Print["There are " <> IntegerString[Length[splitEdges]] <> " cycle edges after split."];
			cycleEdges = splitEdges;
		];
	];

	range = Max[Map[Norm[projectTo3d[#]]&, vertices]] * rangeFactor;
	thickness = range * 0.008;
	thicknessDecayExponent = 1;

	color = colorMap[{p,q,r}];

	graphElements = {
			{
				Opacity[1], FaceForm[color], 
				Table[
		            Cylinder[
		                Map[projectTo3d, componentEdges[[edgeIndex]]], 
		                thickness * 5/(1 + componentEdges[[edgeIndex]][[1]][[1]]^thicknessDecayExponent + componentEdges[[edgeIndex]][[2]][[1]]^thicknessDecayExponent)
		            ], 
		            {edgeIndex, 1, Length[componentEdges]}
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
		};

	altColor = Purple;
	If[color == Purple, altColor = Red];

	If[plotTwoEdges,
		AppendTo[graphElements, {
				Opacity[1], FaceForm[altColor], 
				Table[
		            Cylinder[
		                Map[projectTo3d, componentEdges2[[edgeIndex]]], 
		                thickness * 5/(1 + componentEdges2[[edgeIndex]][[1]][[1]]^thicknessDecayExponent + componentEdges2[[edgeIndex]][[2]][[1]]^thicknessDecayExponent)
		            ], 
		            {edgeIndex, 1, Length[componentEdges2]}
		        ]
			}]
	];

	If[highlightCycle && Length[cycle1] > 0,
		AppendTo[graphElements, {
				Opacity[1], FaceForm[Green], 
				Table[
		            Cylinder[
		                Map[projectTo3d, cycleEdges[[edgeIndex]]], 
		                thickness * 9/(1 + cycleEdges[[edgeIndex]][[1]][[1]]^(thicknessDecayExponent/2) + cycleEdges[[edgeIndex]][[2]][[1]]^(thicknessDecayExponent/2))
		            ], 
		            {edgeIndex, 1, Length[cycleEdges]}
		        ]
			}];
		AppendTo[graphElements, {
				Opacity[1], FaceForm[Darker[Green]], 
				Table[
		            Sphere[
		                projectTo3d[cyclePoints[[vIndex]]],
		                thickness * 18/(1 + cyclePoints[[vIndex]][[1]]^(thicknessDecayExponent/2))
		            ], 
		            {vIndex, 1, Length[cyclePoints]}
		        ]
			}]
	];

	graph = Graphics3D[
		graphElements,
		ViewPoint -> {3, -4, 1}*30,
		PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
		Lighting -> lighting
	];

	outputFileName = shape <> "_" <> model <> "_component_" <> IntegerString[componentIndex] <> "_" <> If[plotTwoEdges, "two_edges_", ""] <> If[highlightCycle, "cycle_", ""] <>  IntegerString[cellThreshold] <> ".png";

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


(* doesn't work well for 57-cell *)

(* 
100:
{13, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74}
1:
{1, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93}
2:
{2, 21, 28, 29, 36, 66, 71, 76, 124, 125, 143, 144, 148, 149, 159, 160, 161, 165, 178, 180, 219, 228, 229, 257, 259, 268, 283, 294, 295, 310, 312, 326, 356, 357, 361, 366, 376, 378, 423, 424, 426, 428, 444, 466, 467, 482, 483, 485, 486, 505, 540, 542, 543, 544, 545, 546, 618, 619, 620, 628, 642, 643, 648, 649, 657, 659, 681, 688, 694, 695, 698, 734, 751, 753}

300:
{333, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212}

1:
{1, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1614, 1615, 1616, 1617, 1618, 1619, 1620, 1621, 1622, 1623, 1624, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638, 1639, 1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649, 1650, 1651, 1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1669, 1670, 1671, 1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 1680, 1681, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 1690, 1691, 1692, 1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1701, 1702, 1703, 1704, 1705, 1706, 1707, 1708, 1709, 1710, 1711, 1712, 1713, 1714, 1715, 1716, 1717, 1718, 1719, 1720, 1721, 1722, 1723, 1724, 1725, 1726, 1727, 1728, 1729, 1730, 1731, 1732, 1733, 1734, 1735, 1736, 1737, 1738, 1739, 1740, 1741, 1742, 1743, 1744, 1745, 1746, 1747, 1748, 1749, 1750, 1751, 1752, 1753, 1754, 1755, 1756, 1757, 1758, 1759, 1760, 1761, 1762, 1763, 1764, 1765, 1766, 1767, 1768, 1769, 1770, 1771, 1772, 1773, 1774, 1775, 1776, 1777, 1778, 1779, 1780, 1781, 1782, 1783, 1784, 1785, 1786, 1787, 1788, 1789, 1790, 1791, 1792, 1793}
2:
{2, 21, 26, 27, 37, 79, 80, 81, 112, 113, 116, 117, 120, 121, 126, 128, 141, 186, 187, 195, 219, 230, 231, 257, 259, 268, 290, 291, 301, 306, 316, 318, 343, 354, 356, 371, 373, 376, 443, 445, 448, 449, 454, 472, 473, 484, 485, 488, 490, 516, 558, 580, 581, 600, 604, 611, 612, 613, 614, 647, 649, 650, 687, 688, 718, 720, 722, 723, 761, 763, 764, 765, 768, 769, 782, 783, 788, 789, 798, 800, 856, 857, 859, 861, 874, 875, 1014, 1015, 1018, 1057, 1059, 1073, 1085, 1094, 1095, 1111, 1113, 1122, 1182, 1183, 1188, 1189, 1198, 1199, 1205, 1206, 1210, 1211, 1242, 1243, 1260, 1261, 1274, 1275, 1282, 1283, 1343, 1345, 1353, 1357, 1366, 1367, 1390, 1391, 1396, 1397, 1400, 1401, 1406, 1407, 1416, 1470, 1471, 1476, 1500, 1501, 1514, 1515, 1522, 1523, 1584, 1585, 1591, 1593, 1606, 1607, 1836, 1837, 1844, 1845, 1846, 1848, 1860, 1861, 1862, 1863, 1870, 1873, 1914, 1916, 1917, 1970, 1971, 1973, 1992, 1993, 2010, 2011, 2022, 2024, 2044, 2045, 2055, 2056, 2074, 2076, 2128, 2129, 2136, 2137, 2143, 2145, 2162, 2164, 2170, 2172, 2178, 2180, 2243, 2244, 2267, 2269, 2270, 2272, 2275, 2277, 2278, 2279, 2302, 2304, 2356, 2357, 2360, 2361, 2386, 2388, 2399, 2401, 2426, 2427, 2430, 2432}


subset 1:
nearest neighbors:
12 vertices of an icosahedron {3,5}

Round[subset1[[Range[13]]]/(subset1[[2]][[4]]), epsilon]   
{
	{0.240495, 0., 0., 0.}, 
	{1.91726, -1.61803, 0., 1.}, 
	{1.91726, -1.61803, 0., -1.},
	{1.91726, -1., -1.61803, 0.}, 
	{1.91726, -1., 1.61803, 0.}, 
	{1.91726, 0., -1., -1.61803},
	{1.91726, 0., -1., 1.61803}, 
	{1.91726, 0., 1., -1.61803}, 
	{1.91726, 0., 1., 1.61803},
	{1.91726, 1., -1.61803, 0.}, 
	{1.91726, 1., 1.61803, 0.}, 
	{1.91726, 1.61803, 0., 1.},
	{1.91726, 1.61803, 0., -1.}
}

*)
