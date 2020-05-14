// Gmsh project created on Sat May 09 00:28:38 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {0, 10, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {2, 3};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Circle(5) = {5, 5, 0, 2, 0, 2*Pi};
//+
Curve Loop(1) = {1, -2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Surface(1) = {1, 2};
//+
Curve Loop(3) = {1, -2, 3, 4};
//+
Curve Loop(4) = {5};
//+
Plane Surface(1) = {3, 4};
//+
Physical Curve("Neumann", 6) = {4};
//+
Physical Curve("Dirichlet", 7) = {2};
//+
Physical Surface("matrix", 8) = {1};
//+
Characteristic Length {5} = 0.2;
//+

//+
Physical Curve("Neumann") += {4};
