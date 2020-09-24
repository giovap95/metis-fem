L = 3;

// Points
Point(1) = {0, 0, 0, L};
Point(2) = {200, 0, 0, L};
Point(3) = {200, 200, 0, L};
Point(4) = {0, 200, 0, L};
Point(5) = {25, 0, 0, L};
Point(6) = {0, 25, 0, L};


// Lines
Line(1) = {5, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 6};
Circle(5) = {6, 1, 5};

// Surfaces
Line Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// Physical groups
Physical Line("DirichletLeft") = {4};
Physical Line("DirichletBottom") = {1};
Physical Line("Neumann") = {2};
Physical Surface("AISI 316") = {1};
