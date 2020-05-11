// Gmsh project created on Sat May  2 20:57:59 2020
//+
Point(1) = {-0, 0, 0, 1.0};
//+
Point(2) = {-0, 1, 0, 1.0};
//+
Point(3) = {2, 0, 0, 1.0};
//+
Point(4) = {2, 1, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet") = {4};
//+
Physical Curve("outlet") = {2};
//+
Physical Curve("top_wall") = {3};
//+
Physical Curve("bottom_wall") = {1};
//+
Characteristic Length {2, 4, 3, 1} = 0.1;
//+
Characteristic Length {2, 1, 4, 3} = 0.5;
