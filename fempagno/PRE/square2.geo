h = 1;                         // Characteristic length of a mesh element
Point(1) = {0, 0, 0, h};       // Point construction
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Line(11) = {1,2};              //Lines
Line(12) = {2,3};
Line(13) = {3,4};
Line(14) = {4,1};
Curve Loop(1) = {11,12,13,14}; // A Boundary
Plane Surface(18) = {1};       // A Surface
Physical Surface(42) = {18};   // Setting a label to the Surface
Physical Curve(98) = {11}; // set physical tag 98 to Curve (=line) 1