
// -------------------------------------------------------
// Vacuum part

ms_far = 2.0;
ms_interm = 0.5;
ms_tip = 0.1;

box_x = 10;
box_h = 18;
tip_x = 1;
tip_h = 6;

// Box corners
Point(1) = {-box_x, 0, 0, ms_far};
Point(2) = {-box_x, box_h, 0, ms_far};
Point(3) = {box_x, box_h, 0, ms_far};
Point(4) = {box_x, 0, 0, ms_far};

Point(5) = {tip_x, 0, 0, ms_interm};

Point(6) = {tip_x, tip_h-tip_x, 0, ms_tip};
Point(7) = {0, tip_h-tip_x, 0, ms_tip};
Point(8) = {0, tip_h, 0, ms_tip};
Point(9) = {-tip_x, tip_h-tip_x, 0, ms_tip};

Point(10) = {-tip_x, 0, 0, ms_interm};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};

Circle(6) = {6, 7, 8};
Circle(7) = {8, 7, 9};

Line(8) = {9, 10};
Line(9) = {10, 1};

Line Loop(11) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Plane Surface(11) = {11};

// --------------------------------------------------------
// Copper part

bulk_h = 4;

Point(11) = {-box_x, -bulk_h, 0, ms_far};
Point(12) = {box_x, -bulk_h, 0, ms_far};

Line(10) = {1, 11};
Line(11) = {11, 12};
Line(12) = {12, 4};

Line Loop(12) = {-10, -11, -12, -4, -5, -6, -7, -8, -9};
Plane Surface(12) = {12};

// -------------------------------------------------------


// Physical id's

vacuum_top = 1;
bulk_bottom = 2;
copper_surface = 3;

Physical Line(vacuum_top) = {2};
//Physical Line(copper_surface) = {4, 5, 6, 7, 8, 9};
Physical Line(bulk_bottom) = {11};

vacuum_domain = 10;
copper_domain = 20;

Physical Surface(vacuum_domain) = {11};
Physical Surface(copper_domain) = {12};


// -----------------------------------------------------
// Meshing parameters

Mesh.RecombineAll = 1;

