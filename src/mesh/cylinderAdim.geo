h = 0.1;

Point(1) = {0.5, 0, 0, h};
Point(2) = {0.5, 0, 1, h};
Point(3) = {0.5, 0, 2, h};
Point(4) = {0.5, 0, 3, h};
Point(5) = {0.5, 0, 4, h};
Point(6) = {0.5, 0, 5, h};


Point(7) = {0, 0, 0, h};
Point(8) = {0, 0, 1, h};
Point(9) = {0, 0, 2, h};
Point(10) = {0, 0, 3, h};
Point(11) = {0, 0, 4, h};
Point(12) = {0, 0, 5, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,12};
Line(7) = {12,11};
Line(8) = {11,10};
Line(9) = {10,9};
Line(10) = {9,8};
Line(11) = {8,7};
Line(12) = {7,1};
Line(13) = {2,8};
Line(14) = {3,9};
Line(15) = {4,10};
Line(16) = {5,11};

Line Loop(1) = {1,13,11,12};
Line Loop(2) = {2,14,10,-13};
Line Loop(3) = {3,15,9,-14};
Line Loop(4) = {4,16,8,-15};
Line Loop(5) = {5,6,7,-16};
Plane Surface(6) = {1};
Plane Surface(7) = {2};
Plane Surface(8) = {3};
Plane Surface(9) = {4};
Plane Surface(10) = {5};

out1[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
  Surface{6,7,8,9,10};
};
out2[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
   Surface{out1[0], out1[5], out1[10], out1[15], out1[20]};
};
out3[] = Extrude{{0,0,1},{0,0,0},Pi/2} {
  Surface{out2[0], out2[5], out2[10], out2[15], out2[20]};
};
out4[] = Extrude{{0,0,1},{0,0,0},Pi/2} {
  Surface{out3[0], out3[5], out3[10], out3[15], out3[20]};
};

Physical Surface("inflow") = {out1[4], out2[4], out3[4], out4[4]};
Physical Surface("outflow") = {out1[23], out2[23], out3[23], out4[23]};
Physical Surface("cut1") = {6,7,8,9,10};
Physical Surface("cut2") = {out1[0], out1[5], out1[10], out1[15], out1[20]};
Physical Surface("cut3") = {out2[0], out2[5], out2[10], out2[15], out2[20]};
Physical Surface("cut4") = {out3[0], out3[5], out3[10], out3[15], out3[20]};
Physical Surface("wall") = {out1[2], out1[7], out1[12], out1[17], out1[22], out2[2], out2[7], out2[12], out2[17], out2[22], out3[2], out3[7], out3[12], out3[17], out3[22], out4[2], out4[7], out4[12], out4[17], out4[22]};

Physical Volume("fluid") = {out1[1], out1[6], out1[11], out1[16], out1[21], out2[1], out2[6], out2[11], out2[16], out2[21], out3[1], out3[6], out3[11], out3[16], out3[21], out4[1], out4[6], out4[11], out4[16], out4[21] };
