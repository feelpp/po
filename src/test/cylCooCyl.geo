cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {0.5, 0, 0, 1};
Point(3) = {0, 0, 5, 1};
Point(4) = {0.5, 0, 5, 1};
Point(5) = {0, 6.283185307179586, 0, 1};
Point(6) = {0.5, 6.283185307179586, 5, 1};
Point(7) = {0, 6.283185307179586, 5, 1};
Point(8) = {0.5, 6.283185307179586, 0, 1};
Line(1) = {5, 1};
Line(2) = {1, 2};
Line(3) = {2, 8};
Line(4) = {8, 5};
Line(5) = {7, 3};
Line(6) = {3, 4};
Line(7) = {4, 6};
Line(8) = {6, 7};
Line(9) = {4, 2};
Line(10) = {1, 3};
Line(11) = {6, 8};
Line(12) = {5, 7};
Line Loop(14) = {5, -10, -1, 12};
Plane Surface(14) = {14};
Line Loop(16) = {7, 11, -3, -9};
Plane Surface(16) = {16};
Line Loop(18) = {6, 9, -2, 10};
Plane Surface(18) = {18};
Line Loop(20) = {2, 3, 4, 1};
Plane Surface(20) = {20};
Line Loop(22) = {4, 12, -8, 11};
Plane Surface(22) = {22};
Line Loop(24) = {8, 5, 6, 7};
Plane Surface(24) = {24};
Surface Loop(26) = {16, 24, 22, 20, 18, 14};
Volume(26) = {26};
Physical Surface(27) = {24};
Physical Surface(28) = {18};
Physical Surface(29) = {20};
Physical Surface(30) = {22};
Physical Surface(31) = {14};
Physical Surface(32) = {16};
Physical Volume(33) = {26};
