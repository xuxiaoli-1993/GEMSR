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
Characteristic Length {2, 4, 3, 1} = 0.1;
//+
Characteristic Length {2, 1, 4, 3} = 0.5;

//+
Physical Curve("p1") = {4, 2};
//+
Physical Curve("p2") = {3, 1};
