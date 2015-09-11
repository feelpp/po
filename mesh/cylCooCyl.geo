h=0.1;

Point(1) = {0.,0.,0.,h};
Point(2) = {0.,0.,5.,h};
Point(3) = {0.,3.*Pi/4.,0.,h};
Point(4) = {0.,3.*Pi/4.,5.,h};
Point(5) = {1.,0.,0.,h};
Point(6) = {1.,0.,5.,h};
Point(7) = {1.,3.*Pi/4.,0.,h};
Point(8) = {1.,3.*Pi/4.,5.,h};

Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,5};
Line(4) = {5,1};
Line(5) = {3,4};
Line(6) = {4,8};
Line(7) = {8,7};
Line(8) = {7,3};
Line(9) = {1,3};
Line(10) = {2,4};
Line(11) = {6,8};
Line(12) = {5,7};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};
Line Loop(3) = {1,10,-5,-9};
Plane Surface(3) = {3};
Line Loop(4) = {-3,11,7,-12};
Plane Surface(4) = {4};
Line Loop(5) = {10,6,-11,-2};
Plane Surface(5) = {5};
Line Loop(6) = {9,-8,-12,4};
Plane Surface(6) = {6};

Physical Surface(1) = {1,2,3,4,5,6};

Surface Loop(1) = {1,2,3,4,5,6};

Volume(1) = {1};

Physical Volume(1) = {1};