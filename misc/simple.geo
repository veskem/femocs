

Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 0, 0, 1.0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};

Plane Surface(5) = {5};

// physical id's

Physical Line(22) = {2};
Physical Line(33) = {4};

Physical Surface(44) = {5};


