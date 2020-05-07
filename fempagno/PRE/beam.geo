LX = 5; LY = 0.5; Nx = 40; Ny = 5;

Point(1) = {0, 0, 0, 1};
Point(2) = {LX, 0, 0, 1};
Point(3) = {LX, LY, 0, 1};
Point(4) = {0, LY, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1,3} = Nx+1;
Transfinite Line {2,4} = Ny+1;
Transfinite Surface{1};
Recombine Surface{1};

Physical Line("Fixed") = {4};
Physical Line("Forcing") = {2};
Physical Surface("Beam") = {1};

Mesh(2);

