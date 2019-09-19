HyperbolicInner[v1_,v2_]:=v1[[1]] v2[[1]]-v1[[2]] v2[[2]]-v1[[3]] v2[[3]];
rotation7=(1	0	0
0	Cos[(2\[Pi])/7]	Sin[(2\[Pi])/7]
0	-Sin[(2\[Pi])/7]	Cos[(2\[Pi])/7]

);
ExtendByRotation[vertices_]:=Map[rotation7.#&,vertices];
ExtendByRotation[vertices_]:=Union[Flatten[Table[Map[MatrixPower[rotation7,k].#&,vertices],{k,0,6}],1]];
cosh2phi=8/3 Cos[\[Pi]/7]^2-1;
cosine=Cos[(2\[Pi])/7];
radiusSquare=(1-cosine)/(cosh2phi-1);
radius=Sqrt[radiusSquare]//FullSimplify;
x0=Sqrt[radiusSquare+1]//FullSimplify;
b=1+2 cosine;
a=Sqrt[radiusSquare+b^2];
NextPointInHeptagon[vertices_]:=b(vertices[[-1]]-vertices[[-2]])+vertices[[-3]];
(*NextFourPointsInHeptagon[vertices_]:=Module[{},For[i=0,i<4,i++,AppendTo[vertices,NextPointInHeptagon[vertices]]]];*)
CompleteHeptagon[vertices_]:=Module[{newVertices=vertices},For[i=0,i<4,i++,AppendTo[newVertices,NextPointInHeptagon[newVertices]]];TrigReduce[newVertices]];
CompleteHeptagonByIndex[vertices_,i_,j_,k_]:=CompleteHeptagon[{vertices[[i]],vertices[[j]],vertices[[k]]}];
TargetInner[vs_]:=Min[Table[N[Re[HyperbolicInner[vs[[1]],vs[[k]]]]],{k,2,Length[vs]}]];
v0s=Table[{x0,Cos[(2 k \[Pi])/7],Sin[(2 k \[Pi])/7]},{k,0,6}];
v10={a,b,0};
vertices=Union[v0s,{v10}];
vertices=Union[TrigReduce[ExtendByRotation[vertices]]];
vertices=Union[vertices,CompleteHeptagonByIndex[vertices,14,1,7]];
vertices=Union[TrigReduce[ExtendByRotation[vertices]]];
vertices=Union[vertices,CompleteHeptagonByIndex[vertices,26,14,34]];
vertices=Union[TrigReduce[ExtendByRotation[vertices]]];
(*vertices=Union[vertices,CompleteHeptagonByIndex[vertices,34,26,55]];
vertices=Union[TrigReduce[ExtendByRotation[vertices]]];*)
extended=vertices;
Print["Vertex count = ",Length[extended]]
Vertex count = 63
edgeIndices={};
target=TargetInner[v0s];
For[i=1,i<=Length[extended],i++,For[j=i+1,j<=Length[extended],j++,If[Abs[HyperbolicInner[N[extended[[i]]],N[extended[[j]]]]-target]<1*^-10,AppendTo[edgeIndices,{i,j}]]]];
Print["Edge count = ", Length[edgeIndices]]
edgeIndices
Edge count = 77
{{1,6},{1,7},{1,14},{2,3},{2,4},{2,8},{3,5},{3,9},{4,6},{4,12},{5,7},{5,13},{6,10},{7,11},{8,37},{8,44},{9,36},{9,45},{10,41},{10,49},{11,42},{11,47},{12,38},{12,43},{13,39},{13,46},{14,40},{14,48},{15,37},{15,45},{16,38},{16,44},{17,36},{17,46},{18,41},{18,43},{19,39},{19,47},{20,40},{20,49},{21,42},{21,48},{22,30},{22,50},{23,29},{23,51},{24,32},{24,52},{25,31},{25,53},{26,34},{26,55},{27,33},{27,54},{28,35},{28,56},{29,57},{30,58},{31,59},{32,60},{33,62},{34,61},{35,63},{36,58},{37,57},{38,59},{39,60},{40,63},{41,61},{42,62},{43,53},{44,51},{45,50},{46,52},{47,54},{48,56},{49,55}}
edgeThickness = 0.2;
vertexThickness = 2 edgeThickness;
edgesGraphics=Map[Cylinder[{extended[[#[[1]]]],extended[[#[[2]]]]},edgeThickness]&,edgeIndices];
verticesGraphics=Map[Sphere[#,vertexThickness]&,extended];
Print["Vertex count = ",Length[extended], "; edge count = ", Length[edgeIndices]]
graph=Graphics3D[Flatten[{verticesGraphics,edgesGraphics}]]
Vertex count = 63; edge count = 77
