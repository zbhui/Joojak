ri=12;
ro=15;
li=20;
lo=200;
Point(1)={0,0,0,li};
Point(2)={ri,0,0,li};
Point(3)={-ri,0,0,li};
Point(4)={ro,0,0,lo};
Point(5)={-ro,0,0,lo};
Circle(1)={2,1,3};
Circle(2)={3,1,2};
Circle(3)={4,1,5};
Circle(4)={5,1,4};
Line Loop(5) = {4, 3};
Line Loop(6) = {2, 1};
Plane Surface(7) = {5, 6};
Physical Line(8) = {4, 3};
Physical Line(9) = {2, 1};
Physical Surface(10) = {7};
