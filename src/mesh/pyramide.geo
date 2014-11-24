h=2;

Point(1) = {1,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {0,0,1,h};
Point(4) = {-1,0,0,h};
Point(5) = {0,-1,0,h};
Point(6) = {0,0,-1,h};

Line(1) = {1,2};
Line(2) = {3,1};
Line(3) = {3,2};
Line(4) = {2,4};
Line(5) = {3,4};
Line(6) = {3,5};
Line(7) = {1,5};
Line(8) = {4,5};
Line(9) = {6,1};
Line(10) = {6,2};
Line(11) = {6,4};
Line(12) = {6,5};

Line Loop(1) = {1,-3,2};
Line Loop(2) = {3,4,-5};
Line Loop(3) = {4,-11,10};
Line Loop(4) = {10,-1,-9};
Line Loop(5) = {2,7,-6};
Line Loop(6) = {5,8,-6};
Line Loop(7) = {8,-12,11};
Line Loop(8) = {9,7,-12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

Physical Surface(1) = {1,2,3,4,5,6,7,8};

Surface Loop(1) = {1,2,3,4,5,6,7,8};
Volume(1) = {1};
Physical Volume(1) = {1};

