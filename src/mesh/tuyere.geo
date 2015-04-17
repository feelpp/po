h=0.1;

l1=2;
l2=1;
l3=2;

r1=1;
r3=0.5;

// premier troncon
Point(1) = {0,0,0,h};
Point(2) = {-r1,0,0,h};
Point(3) = {0,r1,0,h};
Point(4) = {r1,0,0,h};
Point(5) = {0,-r1,0,h};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1) = {1,2,3,4};
Translate {0,0,l1}{ Duplicata { Line{1,2,3,4};}}

Line(9) = {2,6};
Line(10) = {3,8};
Line(11) = {4,13};
Line(12) = {5,18};

Line Loop(2) = {1,10,-5,-9};
Line Loop(3) = {2,11,-6,-10};
Line Loop(4) = {3,12,-7,-11};
Line Loop(5) = {4,9,-8,-12};

Dilate {{0,0,l1+l2},r3/r1}{ Translate {0,0,l2}{Duplicata { Line{5,6,7,8};}}}

Line(17) = {6,19};
Line(18) = {8,21};
Line(19) = {13,26};
Line(20) = {18,31};

Line Loop(6) = {5,18,-13,-17};
Line Loop(7) = {6,19,-14,-18};
Line Loop(8) = {7,20,-15,-19};
Line Loop(9) = {8,17,-16,-20};

Translate {0,0,l3}{ Duplicata { Line{13,14,15,16};}}

Line(25) = {19,32};
Line(26) = {21,34};
Line(27) = {26,39};
Line(28) = {31,44};

Line Loop(10) = {13,26,-21,-25};
Line Loop(11) = {14,27,-22,-26};
Line Loop(12) = {15,28,-23,-27};
Line Loop(13) = {16,25,-24,-28};

Line Loop(14) = {21,22,23,24};

Plane Surface(15) = {1};
Ruled Surface(16) = {2};
Ruled Surface(17) = {3};
Ruled Surface(18) = {4};
Ruled Surface(19) = {5};
Ruled Surface(20) = {6};
Ruled Surface(21) = {7};
Ruled Surface(22) = {8};
Ruled Surface(23) = {9};
Ruled Surface(24) = {10};
Ruled Surface(25) = {11};
Ruled Surface(26) = {12};
Ruled Surface(27) = {13};
Plane Surface(28) = {14};

Physical Surface("inflow") = {15};
Physical Surface("wall1") = {16,17,18,19};
Physical Surface("wall2") = {20,21,22,23};
Physical Surface("wall3") = {24,25,26,27};
Physical Surface("outflow") = {28};

Surface Loop(1) = {15,16,17,18,19,20,21,22,23,24,25,26,27,28};
Volume(1) = {1};
Physical Volume("volume") = {1};

