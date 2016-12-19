cl__1 = 1;
Point(1) = {-0.2, 0, 0, 1};
Point(2) = {0.2, -0, 0, 1};
Point(3) = {0.2, 0.4, 0, 1};
Point(4) = {0.2, 0.6, 0, 1};
Point(5) = {-0.2, 0.6, 0, 1};
Point(6) = {-0.2, 0.2, 0, 1};

Line(1) = {2, 1};
Line(2) = {6, 1};
Line(3) = {2, 3};
Line(4) = {6, 3};
Line(5) = {5, 6};
Line(6) = {5, 4};
Line(7) = {4, 3};

Line Loop(12) = {2, -1, 3, -4};
Plane Surface(12) = {12};

Line Loop(13) = {5, 4, -7, -6};
Plane Surface(13) = {13};

Physical Surface("vacuum") = {13};
Physical Surface("copper") = {12};
