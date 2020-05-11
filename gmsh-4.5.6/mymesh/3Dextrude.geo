// Gmsh project created on Sun May  3 18:22:37 2020
//+
Point(1) = {-0, -0, 0, 1.0};
//+
Point(2) = {1, -0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 4} {
  Surface{1}; 
}
//+
Physical Surface("outlet") = {26};
//+
Physical Surface("inlet") = {1};
//+
Physical Surface("wall1") = {13};
//+
Physical Surface("wall2") = {25};
//+
Physical Surface("wall3") = {17};
//+
Physical Surface("wall4") = {21};
//+
Characteristic Length {3, 4, 1, 2, 6, 10, 14, 5} = 0.1;
