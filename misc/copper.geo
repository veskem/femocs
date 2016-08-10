
ms_far = 2.0;
ms_interm = 0.5;
ms_tip = 0.1;

box_x = 10;
box_h = 14;
tip_x = 1;
tip_h = 4;

// Box corners
Point(1) = {-box_x, 0, 0, ms_far};
Point(2) = {box_x, 0, 0, ms_far};
Point(3) = {-box_x, box_h, 0, ms_far};
Point(4) = {box_x, box_h, 0, ms_far};

// emitter base points
Point(5) = {-tip_x, 0, 0, ms_interm};
Point(6) = {tip_x, 0, 0, ms_interm};

// Tip points
Point(7) = {-tip_x, tip_h-tip_x, 0, ms_tip};
Point(8) = {0, tip_h-tip_x, 0, ms_interm};
Point(9) = {tip_x, tip_h-tip_x, 0, ms_tip};
Point(10) = {0, tip_h, 0, ms_tip};

// Lines
Circle(1) = {9, 8, 10};
Circle(2) = {10, 8, 7};

Line(3) = {7, 5};
Line(4) = {5, 1};
Line(8) = {2, 6};
Line(9) = {6, 9};

// ----------------------------------------------
// Copper part

bulk_h = 4;

Point(11) = {-box_x, -bulk_h, 0, ms_far};
Point(12) = {box_x, -bulk_h, 0, ms_far};

Line(10) = {1, 11};
Line(11) = {11, 12};
Line(12) = {12, 2};

Line Loop(12) = {1, 2, 3, 4, 10, 11, 12, 8, 9};
Plane Surface(12) = {12};

// ----------------------------------------------
// Physical id-s

vacuum_top = 1;
copper_surface = 2;
bulk_bottom = 3;

Physical Line(copper_surface) = {1, 2, 3, 4, 9};
Physical Line(bulk_bottom) = {11};

vacuum_domain = 10;
copper_domain = 20;

Physical Surface(copper_domain) = {12};


