h=0.5;

Lt=18;
r=6;

angle=35;
rad=Pi*angle/180;
ac=16;
Lc=ac*Sin(rad);
bc=ac*Cos(rad);
L=Lt+Lc;

Point(1) = {0,0,-L,h};
Point(2) = {0,0,-Lc,h};
Point(3) = {0,0,Lc,h};
Point(4) = {0,0,L,h};

Point(5) = {r,0,-L,h};
Point(6) = {r,0,-Lc,h};
Point(7) = {r+bc,0,0,h};
Point(8) = {r,0,Lc,h};
Point(9) = {r,0,L,h};

Point(10) = {0,r,-L,h};
Point(11) = {0,r,-Lc,h};
Point(12) = {0,r+bc,0,h};
Point(13) = {0,r,Lc,h};
Point(14) = {0,r,L,h};

Point(15) = {-r,0,-L,h};
Point(16) = {-r,0,-Lc,h};
Point(17) = {-r-bc,0,0,h};
Point(18) = {-r,0,Lc,h};
Point(19) = {-r,0,L,h};

Point(20) = {0,-r,-L,h};
Point(21) = {0,-r,-Lc,h};
Point(22) = {0,-r-bc,0,h};
Point(23) = {0,-r,Lc,h};
Point(24) = {0,-r,L,h};

// center lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

i = 0;
For(1:4)

Line(4+i*11) = {5+i*5,1};
Line(5+i*11) = {6+i*5,2};
Line(6+i*11) = {8+i*5,3};
Line(7+i*11) = {9+i*5,4};

Line(8+i*11) = {5+i*5,6+i*5};
Circle(9+i*11) = {6+i*5,7+i*5,8+i*5};
Line(10+i*11) = {8+i*5,9+i*5};

Circle(11+i*11) = {5+i*5,1,5+((i+1)*5)%20};
Circle(12+i*11) = {6+i*5,2,6+((i+1)*5)%20};
Circle(13+i*11) = {8+i*5,3,8+((i+1)*5)%20};
Circle(14+i*11) = {9+i*5,4,9+((i+1)*5)%20};


Line Loop(1+i*6) = {1,-(5+i*11),-(8+i*11),4+i*11};
Line Loop(2+i*6) = {2,-(6+i*11),-(9+i*11),5+i*11};
Line Loop(3+i*6) = {3,-(7+i*11),-(10+i*11),6+i*11};

Line Loop(4+i*6) = {8+i*11,12+i*11,-(8+((i+1)*11)%44),-(11+i*11)};
Line Loop(5+i*6) = {9+i*11,13+i*11,-(9+((i+1)*11)%44),-(12+i*11)};
Line Loop(6+i*6) = {10+i*11,14+i*11,-(10+((i+1)*11)%44),-(13+i*11)};

Plane Surface(1+i*6) = {1+i*6};
Plane Surface(2+i*6) = {2+i*6};
Plane Surface(3+i*6) = {3+i*6};

i = i + 1;
EndFor

Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(10) = {10};
Ruled Surface(11) = {11};
Ruled Surface(12) = {12};
Ruled Surface(16) = {16};
Ruled Surface(17) = {17};
Ruled Surface(18) = {18};
Ruled Surface(22) = {22};
Ruled Surface(23) = {23};
Ruled Surface(24) = {24};

Line Loop(25) = {11,22,33,44};
Line Loop(26) = {12,23,34,45};
Line Loop(27) = {13,24,35,46};
Line Loop(28) = {14,25,36,47};

Plane Surface(25) = {25};
Plane Surface(26) = {26};
Plane Surface(27) = {27};
Plane Surface(28) = {28};

Surface Loop(1) = {4,5,6,10,11,12,16,17,18,22,23,24,25,28};
Volume(1) = {1};

Physical Surface("inflow") = {25};
Physical Surface("outflow") = {28};
Physical Surface("wall") = {4,5,6,10,11,12,16,17,18,22,23,24};
Physical Volume("domain") = {1};

