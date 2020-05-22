// Gmsh project created on Mon May  4 22:06:06 2020
//+
Point(1) = {-0, 0, 0, 1.0};
//+
Point(2) = {0, 0.1, 0, 1.0};
//+
Point(3) = {0.1, 0.1, 0, 1.0};
//+
Point(4) = {0.1, -0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.75} {
  Surface{1}; 
}
//+
Extrude {0, 0, -0.75} {
  Surface{1}; 
}

//+
Physical Volume("vol1") = {1};
//+
Physical Volume("vol2") = {2};
//+
Physical Surface("wall1") = {13, 35};
//+
Physical Surface("wall2") = {17, 39};
//+
Physical Surface("wall3") = {21, 43};
//+
Physical Surface("wall4") = {25, 47};
//+
Physical Surface("outlet") = {26};
//+
Physical Surface("inlet") = {48};
//+
Transfinite Surface {26};
