// Gmsh project created on Sat May 09 00:28:38 2020
SetFactory("OpenCASCADE");

h = 2;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {10, 0, 0, h};
//+
Point(3) = {10, 10, 0,h};
//+
Point(4) = {0, 10, 0, h};
//+
Line(1) = {1,2};
//+
Line(2) = {2,3};
//+
Line(3) = {3,4};
//+
Line(4) = {4,1};
//+
Circle(5) = {5, 5, 0, 2, 0, 2*Pi};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};

Plane Surface(1) = {1,2};

Physical Surface("ALU 6061" , 100) = {1};

Physical Line("Neumann",10) = {2};
Physical Line("Dirichlet",20) = {4};
