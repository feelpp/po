h=0.1;

Point(1) = {0.,0.,0.,h};
Point(2) = {1.,0.,0.,h};
Point(3) = {-Sqrt(2)/2,Sqrt(2)/2,0,h};
Circle(1) = {2,1,3};
Line(2) = {1,2};
Line(3) = {3,1};
Line Loop(4) = {2, 1, 3};
Plane Surface(5) = {4};
Extrude {0,0,5} {Surface{5}; }
Physical Surface(6) = {5,13,17,21,22};
Physical Volume(7) = {1};