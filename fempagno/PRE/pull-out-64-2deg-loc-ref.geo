// pull-out problem
Point(1) = {0,0,0,0.006325};
Point(2) = {10,0,0,10};

Line(1) = {1,2};

Physical Point("start") = {1};
Physical Point("end") = {2};

Physical Line("AISI 316") = {1};
