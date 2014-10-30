h1=0.01;

// sphere interieure
Point(1) = {0,0,0,h1};
Point(2) = {0.1,0,0,h1};
Point(3) = {0,0.1,0,h1};
Point(4) = {0,0,0.1,h1};
Point(5) = {-0.1,0,0,h1};
Point(6) = {0,-0.1,0,h1};
Point(7) = {0,0,-0.1,h1};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,2};
Circle(5) = {2,1,7};
Circle(6) = {7,1,5};
Circle(7) = {5,1,4};
Circle(8) = {4,1,2};
Circle(9) = {6,1,7};
Circle(10) = {7,1,3};
Circle(11) = {3,1,4};
Circle(12) ={4,1,6};

Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11}; 
Line Loop(3) = {3,-12,-7}; 
Line Loop(4) = {4,-8,12}; 
Line Loop(5) = {5,10,-1}; 
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9}; 
Line Loop(8) = {-4,9,-5};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

Surface Loop(1) = {1,2,3,4,5,6,7,8};

Volume(1) = {1};
Physical Volume("int") = {1};

// sphere exterieure
h=0.1;

Point(8) = {0,0,0,h};
Point(9) = {1,0,0,h};
Point(10) = {0,1,0,h};
Point(11) = {0,0,1,h};
Point(12) = {-1,0,0,h};
Point(13) = {0,-1,0,h};
Point(14) = {0,0,-1,h};

Circle(13) = {9,8,10};
Circle(14) = {10,8,12};
Circle(15) = {12,8,13};
Circle(16) = {13,8,9};
Circle(17) = {9,8,14};
Circle(18) = {14,8,12};
Circle(19) = {12,8,11};
Circle(20) = {11,8,9};
Circle(21) = {13,8,14};
Circle(22) = {14,8,10};
Circle(23) = {10,8,11};
Circle(24) ={11,8,13};

Line Loop(9) = {13,23,20};
Line Loop(10) = {14,19,-23}; 
Line Loop(11) = {15,-24,-19}; 
Line Loop(12) = {16,-20,24}; 
Line Loop(13) = {17,22,-13}; 
Line Loop(14) = {-14,-22,18};
Line Loop(15) = {-15,-18,-21}; 
Line Loop(16) = {-16,21,-17};

Ruled Surface(9) = {9};
Ruled Surface(10) = {10};
Ruled Surface(11) = {11};
Ruled Surface(12) = {12};
Ruled Surface(13) = {13};
Ruled Surface(14) = {14};
Ruled Surface(15) = {15};
Ruled Surface(16) = {16};

Surface Loop(2) = {9,10,11,12,13,14,15,16};

Volume(2) = {1,2};
Physical Volume("ext") = {2};
