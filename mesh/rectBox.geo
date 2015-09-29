h=0.1;

Point(1) = {-0.5,-0.4,-0.6,h};
Point(2) = {-0.5,-0.4,0.6,h};
Point(3) = {-0.5,0.4,-0.6,h};
Point(4) = {-0.5,0.4,0.6,h};
Point(5) = {0.5,-0.4,-0.6,h};
Point(6) = {0.5,-0.4,0.6,h};
Point(7) = {0.5,0.4,-0.6,h};
Point(8) = {0.5,0.4,0.6,h};

Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,5};
Line(4) = {2,4};
Line(5) = {2,6};
Line(6) = {3,4};
Line(7) = {3,7};
Line(8) = {4,8};
Line(9) = {5,6};
Line(10) = {5,7};
Line(11) = {6,8};
Line(12) = {7,8};

Line Loop(1) = {3,9,-5,-1};
Line Loop(2) = {10,12,-11,-9};
Line Loop(3) = {-7,6,8,-12};
Line Loop(4) = {-2,1,4,-6};
Line Loop(5) = {3,10,-7,-2};
Line Loop(6) = {5,11,-8,-4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Surface(1) = {1,2,3,4,5,6};
Physical Volume(1) = {1};


