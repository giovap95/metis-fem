// Gmsh project created on Fri May 29 11:42:08 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {-0, -0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("Dirichlet") = {1, 2};
//+
Physical Curve("rock") = {1};
//+
