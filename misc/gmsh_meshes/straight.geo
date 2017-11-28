cl__1 = 1;
Point(1) = {-0.2, 0, 0, 1};
Point(2) = {0.2, 0, 0, 1};
Point(3) = {0.2, 0.6, 0, 1};
Point(4) = {-0.2, 0.6, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(12) = {1, 2, 3, 4};
Plane Surface(12) = {12};

Physical Surface("copper") = {12};
