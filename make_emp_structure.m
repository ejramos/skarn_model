function Emp_coords = make_emp_structure
%file: make_emp_structure.m

%author: Evan J. Ramos
%date:   08 Jan 2018

%Description: This script generates a structure containing the x- and y-
%coordinates for each feature in the A to A' transect of the Empire
%Mountain skarn in the Mineral King pendant. Subfield names are based off
%features from the Empire Mountain transect and simplified transect of
%Ramos et al., 2018 in prep.

%Main variable name:      Emp_coords
%Subfields of Emp_coords: i, ii, iii
%Subfield of i:           Kmrt
%Subfields of ii:         Kma, Kmrt
%Subfields of iii:        L, P, R
%Subfields of L and R:    Farewell, m1, Trma, cs1, Jmdt, cs2, m2, Kmrt

%Each subfield contains a 1 x 4 vector which contains the the x and y
%coordinates of the feature (all rectangular features), with the form:
%[xmin xmax ymin ymax]. All units are in meters.

%NOTE 1: SHIFT is added to the x-coordinates to account for the shift in 
%the origin and other crustal material to the left.

%NOTE 2: subfields I and II contain the overall geometry of subfields i and
%ii, respectively

%% Origin shift value

SHIFT = 2000; %[m] shift of x-values

%% Subfield i

Emp_coords.I = [0 444 0 1122];
Emp_coords.i.Kmrt = [0 386 93 1122];

%x-coord shift
Emp_coords.I(1:2) = Emp_coords.I(1:2) + SHIFT;
Emp_coords.i.Kmrt(1:2) = Emp_coords.i.Kmrt(1:2) + SHIFT;

%% Subfield ii

Emp_coords.II = [444 1496 0 1356];
Emp_coords.ii.Kma = [771 935 1076 1356];
Emp_coords.ii.Kmrt = [1122 1402 164 1356];

%x-coord shift
Emp_coords.II(1:2) = Emp_coords.II(1:2) + SHIFT;
Emp_coords.ii.Kma(1:2) = Emp_coords.ii.Kma(1:2) + SHIFT;
Emp_coords.ii.Kmrt(1:2) = Emp_coords.ii.Kmrt(1:2) + SHIFT;

%% Subfield iii

%NOTE: iii is non-rectangular and thus has no prescribed coordinates

%L
Emp_coords.iii.L.Farewell = [1496 2057 0 748];
Emp_coords.iii.L.m1       = [1496 2057 748 877];
Emp_coords.iii.L.Trma     = [1496 2057 877 1041];
Emp_coords.iii.L.cs1      = [1496 2057 1041 1134];
Emp_coords.iii.L.Jmdt     = [1496 2057 1134 1403];
Emp_coords.iii.L.cs2      = [1496 2057 1403 1543];
Emp_coords.iii.L.m2       = [1496 2057 1543 1660];
Emp_coords.iii.L.Kmrt     = [1496 2057 1660 1800];

%x-coord shift
Emp_coords.iii.L.Farewell(1:2) = Emp_coords.iii.L.Farewell(1:2) + SHIFT;
Emp_coords.iii.L.m1(1:2) = Emp_coords.iii.L.m1(1:2) + SHIFT;
Emp_coords.iii.L.Trma(1:2) = Emp_coords.iii.L.Trma(1:2) + SHIFT;
Emp_coords.iii.L.cs1(1:2) = Emp_coords.iii.L.cs1(1:2) + SHIFT;
Emp_coords.iii.L.Jmdt(1:2) = Emp_coords.iii.L.Jmdt(1:2) + SHIFT;
Emp_coords.iii.L.cs2(1:2) = Emp_coords.iii.L.cs2(1:2) + SHIFT;
Emp_coords.iii.L.m2(1:2) = Emp_coords.iii.L.m2(1:2) + SHIFT;
Emp_coords.iii.L.Kmrt(1:2) = Emp_coords.iii.L.Kmrt(1:2) + SHIFT;

%P
Emp_coords.iii.P = [2057 2805 0 1963];

%x-coord shift
Emp_coords.iii.P(1:2) = Emp_coords.iii.P(1:2) + SHIFT;

%R
Emp_coords.iii.R.Farewell = [2805 3880 0 982];
Emp_coords.iii.R.m1       = [2805 3880 982 1111];
Emp_coords.iii.R.Trma     = [2805 3880 1111 1275];
Emp_coords.iii.R.cs1      = [2805 3880 1275 1368];
Emp_coords.iii.R.Jmdt     = [2805 3880 1368 1637];
Emp_coords.iii.R.cs2      = [2805 3880 1637 1777];
Emp_coords.iii.R.m2       = [2805 3880 1777 1894];
Emp_coords.iii.R.Kmrt     = [2805 3880 1894 2034];

%x-coord shift
Emp_coords.iii.R.Farewell(1:2) = Emp_coords.iii.R.Farewell(1:2) + SHIFT;
Emp_coords.iii.R.m1(1:2) = Emp_coords.iii.R.m1(1:2) + SHIFT;
Emp_coords.iii.R.Trma(1:2) = Emp_coords.iii.R.Trma(1:2) + SHIFT;
Emp_coords.iii.R.cs1(1:2) = Emp_coords.iii.R.cs1(1:2) + SHIFT;
Emp_coords.iii.R.Jmdt(1:2) = Emp_coords.iii.R.Jmdt(1:2) + SHIFT;
Emp_coords.iii.R.cs2(1:2) = Emp_coords.iii.R.cs2(1:2) + SHIFT;
Emp_coords.iii.R.m2(1:2) = Emp_coords.iii.R.m2(1:2) + SHIFT;
Emp_coords.iii.R.Kmrt(1:2) = Emp_coords.iii.R.Kmrt(1:2) + SHIFT;
end