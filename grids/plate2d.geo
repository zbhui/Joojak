scale = 0.2;
Point(1)={ -0.5,0,0,0.05*scale};
Point(2)={ 0,0,0,0.01*scale};
Point(3)={ 1,0,0,0.05*scale};
Point(4)={ 1,1,0,scale};
Point(5)={ -0.5,1,0,scale};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,1};

Line Loop(6)={1,2,3,4,5};

Plane Surface(7) = {6};

Physical Line(8) = {1};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Line(11) = {4};
Physical Line(12) = {5};

Physical Surface(13) = {7};
