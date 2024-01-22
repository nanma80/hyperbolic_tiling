#!/usr/bin/env wolframscript

(* 
parabola projection:

Start from the hyperboloid model
t^2 - x^2 - y^2 = 1
with projection center t=-1, x=y=0,
projection plane: x=t.

Sometimes, we call t as z. In this file, let's use t.

Projection:
(t, x, y) -> (x/(t-x+1), x/(t-x+1), y/(t-x+1))

if we tilt the coordinate by 45 deg in x-t plane, it becomes
(Sqrt[2] x/(t-x+1), y/(t-x+1))

if we change coordinate such that
t = Cosh[u]
x = Sinh[u] Cos[v]
y = Sinh[u] Sin[v]

Then we get:

*)

t[u_, v_] := Cosh[u];
x[u_, v_] := Sinh[u] Cos[v];
y[u_, v_] := Sinh[u] Sin[v];

px[u_, v_] := Sqrt[2] x[u,v]/(t[u,v] - x[u,v] + 1);
py[u_, v_] := y[u,v]/(t[u,v] - x[u,v] + 1);

(* px[u_, v_] := Sqrt[2] Sinh[u] Cos[v]/(Cosh[u] - Sinh[u] Cos[v] + 1 ); *)
(* py[u_, v_] := Sinh[u] Sin[v]/(Cosh[u] - Sinh[u] Cos[v] + 1 ); *)

Print[Sqrt[2] px[u,v] - py[u,v]^2 + 1//FullSimplify];

(* always positive *)
(* 
        2 (1 + Cosh[u])
-------------------------------
                              2
(1 + Cosh[u] - Cos[v] Sinh[u]) 
*)

solution = Solve[Sqrt[2] x/(t - x + 1) == px && y/(t - x + 1) == py && t^2-x^2-y^2==1, {x, y, t} ];
Print[solution];

(* 

                                2                                                        2     2
                 Sqrt[2] px + px          -2 py - Sqrt[2] px py       1 + Sqrt[2] px + px  + py
Out[5]= {{x -> --------------------, y -> ---------------------, t -> --------------------------}}
                                  2                           2                             2
               1 + Sqrt[2] px - py        -1 - Sqrt[2] px + py           1 + Sqrt[2] px - py
 *)

Solve[Sqrt[2] x/(z - x + 1) == px && y/(z - x + 1) == py && z^2-x^2-y^2==1, {x, y, z} ]