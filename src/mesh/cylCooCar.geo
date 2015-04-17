h=0.1;

Point(1) = {0.,0.,0.,h};
Point(2) = {1.,0.,0.,h};
Point(3) = {-Sqrt(2)/2,Sqrt(2)/2,0,h};

Circle(1) = {2,1,3};
Line(2) = {1,2};
Line(3) = {3,1};

Line Loop(1) = {2, 1, 3};

Plane Surface(1) = {1};

Extrude {0,0,5} {Surface{1}; }
Physical Surface(1) = {1,11,15,19,20};
Physical Volume(1) = {1};