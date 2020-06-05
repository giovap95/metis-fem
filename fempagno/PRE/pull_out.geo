// Gmsh project created on Sat May 09 00:28:38 2020
SetFactory("Built-In");

// Defining outer box
Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};

Point(5) = {10, 4, 0, 0.1};
Point(6) = {10, 6, 0, 0.1};
Point(7) = {0, 6, 0, 1.0};
Point(8) = {0, 4, 0, 1.0};

// Matrix
Line(1) = {1,2};
Line(2) = {2,5};
Line(3) = {5,8};
Line(4) = {8,1};

Line(9) = {7,6};
Line(10) = {6,3};
Line(11) = {3,4};
Line(12) = {4,7};

// Internal bar
Line(6) = {6,5}; // right face
Line(8) = {8,7}; // left face


// Closing lines
Curve Loop(1) = {1,2,3,4}; // matrix lower loop
Curve Loop(2) = {8,9,6,3}; // Inner loop
Curve Loop(3) = {9,10,11,12}; // matrix upper loop

//Defining Surfaces
Plane Surface(1) = {1}; // Matrix lower Surface
Plane Surface(3) = {3}; // matrix upper surface
Plane Surface(2) = -{2};    // Bar Surface

Physical Surface("matrix", 100) = {1};
Physical Surface("fiber", 200) = {2};
Physical Surface("matrix", 100) += {3};

// Boundary conditions
Physical Line("Neumann", 10) = {6};
Physical Line("Dirichlet", 20) = {4};
Physical Line("Dirichlet", 20) += {8};
Physical Line("Dirichlet", 20) += {12};
