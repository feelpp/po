h=0.1;

Point(1) = {-Pi,-Pi,-Pi,h};
Point(2) = {-Pi,-Pi,Pi,h};
Point(3) = {-Pi,Pi,-Pi,h};
Point(4) = {-Pi,Pi,Pi,h};
Point(5) = {Pi,-Pi,-Pi,h};
Point(6) = {Pi,-Pi,Pi,h};
Point(7) = {Pi,Pi,-Pi,h};
Point(8) = {Pi,Pi,Pi,h};

Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {4,3};
Line(4) = {3,1};
Line(5) = {2,6};
Line(6) = {6,5};
Line(7) = {5,1};
Line(8) = {3,7};
Line(9) = {7,5};
Line(10) = {8,4};
Line(11) = {8,6};
Line(12) = {8,7};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {1,5,6,7};
Line Loop(3) = {-4,8,9,7};
Line Loop(4) = {8,-12,10,3};
Line Loop(5) = {12,9,-6,-11};
Line Loop(6) = {10,-2,5,-11};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};

Volume(1) = {1};

Physical Surface(1) = {1,2,3,4,5,6};
Physical Volume(1) = {1};
