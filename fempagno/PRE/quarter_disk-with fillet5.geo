L = 5;
rf = 2;
rh = 25;

// Points
Point(1) = {0, 0, 0, L};
Point(2) = {200, 0, 0, L};
Point(3) = {200, 200, 0, L};
Point(4) = {0, 200, 0, L};

Point(8) = {rf, rh + 2*rf, 0, L};
Point(6) = {0, rh + 2*rf, 0, L};
Point(7) = {rf, rh + rf, 0, L};

Point(10) = {rh + 2*rf, rf, 0, L};
Point(5) = {rh + 2*rf, 0, 0, L};
Point(9) = {rh + rf, rf, 0, L};

Point(11) = {rf, rf, 0, L};

// Lines
Line(1) = {5, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 6};
Circle(5) = {6, 8,  7};
Circle(6) = {7, 11, 9};
Circle(7) = {9, 10, 5};

// Surfaces
Line Loop(1) = {1, 2, 3, 4, 5,  6, 7};
Plane Surface(1) = {1};

// Physical groups
Physical Line("DirichletLeft") = {4};
Physical Line("DirichletBottom") = {1};
Physical Line("Neumann") = {2};
Physical Surface("AISI 316") = {1};
