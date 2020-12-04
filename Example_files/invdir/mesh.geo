//2D mesh script for ResIPy (run the following in gmsh to generate a triangular mesh with topograpghy)
cl=0.77;//define characteristic length
//Define surface points
Point(1) = {-7.69,0.00,2968.10,cl};//topography point
Point(2) = {0.00,0.00,2968.10,cl};//electrode
Point(3) = {1.54,0.00,2966.60,cl};//electrode
Point(4) = {3.08,0.00,2964.90,cl};//electrode
Point(5) = {4.61,0.00,2964.10,cl};//electrode
Point(6) = {6.15,0.00,2963.30,cl};//electrode
Point(7) = {7.69,0.00,2962.60,cl};//electrode
Point(8) = {9.23,0.00,2961.50,cl};//electrode
Point(9) = {10.76,0.00,2960.70,cl};//electrode
Point(10) = {12.30,0.00,2960.00,cl};//electrode
Point(11) = {13.84,0.00,2959.30,cl};//electrode
Point(12) = {15.38,0.00,2958.70,cl};//electrode
Point(13) = {16.91,0.00,2958.10,cl};//electrode
Point(14) = {18.45,0.00,2957.30,cl};//electrode
Point(15) = {19.99,0.00,2956.60,cl};//electrode
Point(16) = {21.52,0.00,2955.80,cl};//electrode
Point(17) = {23.06,0.00,2955.00,cl};//electrode
Point(18) = {24.60,0.00,2954.30,cl};//electrode
Point(19) = {26.14,0.00,2953.50,cl};//electrode
Point(20) = {27.68,0.00,2952.80,cl};//electrode
Point(21) = {29.21,0.00,2951.90,cl};//electrode
Point(22) = {30.75,0.00,2951.10,cl};//electrode
Point(23) = {32.29,0.00,2950.20,cl};//electrode
Point(24) = {33.83,0.00,2949.50,cl};//electrode
Point(25) = {35.36,0.00,2948.70,cl};//electrode
Point(26) = {36.90,0.00,2948.10,cl};//electrode
Point(27) = {38.44,0.00,2947.30,cl};//electrode
Point(28) = {39.98,0.00,2946.70,cl};//electrode
Point(29) = {41.51,0.00,2946.10,cl};//electrode
Point(30) = {43.05,0.00,2945.50,cl};//electrode
Point(31) = {44.59,0.00,2945.10,cl};//electrode
Point(32) = {46.13,0.00,2944.50,cl};//electrode
Point(33) = {47.66,0.00,2944.10,cl};//electrode
Point(34) = {49.20,0.00,2943.60,cl};//electrode
Point(35) = {50.74,0.00,2943.40,cl};//electrode
Point(36) = {52.28,0.00,2943.00,cl};//electrode
Point(37) = {53.81,0.00,2942.50,cl};//electrode
Point(38) = {55.35,0.00,2942.10,cl};//electrode
Point(39) = {56.89,0.00,2941.60,cl};//electrode
Point(40) = {58.43,0.00,2941.00,cl};//electrode
Point(41) = {59.96,0.00,2940.60,cl};//electrode
Point(42) = {61.50,0.00,2940.10,cl};//electrode
Point(43) = {63.04,0.00,2939.50,cl};//electrode
Point(44) = {64.58,0.00,2939.10,cl};//electrode
Point(45) = {66.11,0.00,2939.20,cl};//electrode
Point(46) = {67.65,0.00,2939.20,cl};//electrode
Point(47) = {69.19,0.00,2939.10,cl};//electrode
Point(48) = {70.73,0.00,2938.90,cl};//electrode
Point(49) = {72.26,0.00,2938.60,cl};//electrode
Point(50) = {73.80,0.00,2938.30,cl};//electrode
Point(51) = {75.34,0.00,2938.50,cl};//electrode
Point(52) = {76.88,0.00,2938.70,cl};//electrode
Point(53) = {78.41,0.00,2938.60,cl};//electrode
Point(54) = {79.95,0.00,2938.80,cl};//electrode
Point(55) = {81.49,0.00,2938.90,cl};//electrode
Point(56) = {83.03,0.00,2939.10,cl};//electrode
Point(57) = {84.56,0.00,2939.10,cl};//electrode
Point(58) = {92.25,0.00,2939.10,cl};//topography point
//construct lines between each surface point
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,19};
Line(19) = {19,20};
Line(20) = {20,21};
Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,25};
Line(25) = {25,26};
Line(26) = {26,27};
Line(27) = {27,28};
Line(28) = {28,29};
Line(29) = {29,30};
Line(30) = {30,31};
Line(31) = {31,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,35};
Line(35) = {35,36};
Line(36) = {36,37};
Line(37) = {37,38};
Line(38) = {38,39};
Line(39) = {39,40};
Line(40) = {40,41};
Line(41) = {41,42};
Line(42) = {42,43};
Line(43) = {43,44};
Line(44) = {44,45};
Line(45) = {45,46};
Line(46) = {46,47};
Line(47) = {47,48};
Line(48) = {48,49};
Line(49) = {49,50};
Line(50) = {50,51};
Line(51) = {51,52};
Line(52) = {52,53};
Line(53) = {53,54};
Line(54) = {54,55};
Line(55) = {55,56};
Line(56) = {56,57};
Line(57) = {57,58};
//add points below surface to make a fine mesh region
Point(59) = {-7.69,0.00,2960.92,cl*2.00};//base of smoothed mesh region
Point(60) = {-2.14,0.00,2960.92,cl*2.00};//base of smoothed mesh region
Point(61) = {3.42,0.00,2957.55,cl*2.00};//base of smoothed mesh region
Point(62) = {8.97,0.00,2954.51,cl*2.00};//base of smoothed mesh region
Point(63) = {14.52,0.00,2951.86,cl*2.00};//base of smoothed mesh region
Point(64) = {20.07,0.00,2949.38,cl*2.00};//base of smoothed mesh region
Point(65) = {25.63,0.00,2946.59,cl*2.00};//base of smoothed mesh region
Point(66) = {31.18,0.00,2943.67,cl*2.00};//base of smoothed mesh region
Point(67) = {36.73,0.00,2940.99,cl*2.00};//base of smoothed mesh region
Point(68) = {42.28,0.00,2938.62,cl*2.00};//base of smoothed mesh region
Point(69) = {47.83,0.00,2936.87,cl*2.00};//base of smoothed mesh region
Point(70) = {53.39,0.00,2935.46,cl*2.00};//base of smoothed mesh region
Point(71) = {58.94,0.00,2933.69,cl*2.00};//base of smoothed mesh region
Point(72) = {64.49,0.00,2931.95,cl*2.00};//base of smoothed mesh region
Point(73) = {70.04,0.00,2931.81,cl*2.00};//base of smoothed mesh region
Point(74) = {75.60,0.00,2931.36,cl*2.00};//base of smoothed mesh region
Point(75) = {81.15,0.00,2931.70,cl*2.00};//base of smoothed mesh region
Point(76) = {86.70,0.00,2931.92,cl*2.00};//base of smoothed mesh region
Point(77) = {92.25,0.00,2931.92,cl*2.00};//base of smoothed mesh region
//make lines between base of fine mesh region points
Line(58) = {59,60};
Line(59) = {60,61};
Line(60) = {61,62};
Line(61) = {62,63};
Line(62) = {63,64};
Line(63) = {64,65};
Line(64) = {65,66};
Line(65) = {66,67};
Line(66) = {67,68};
Line(67) = {68,69};
Line(68) = {69,70};
Line(69) = {70,71};
Line(70) = {71,72};
Line(71) = {72,73};
Line(72) = {73,74};
Line(73) = {74,75};
Line(74) = {75,76};
Line(75) = {76,77};

//Adding boundaries
//end of boundaries.
//Add lines at leftmost side of fine mesh region.
Line(76) = {1,59};
//Add lines at rightmost side of fine mesh region.
Line(77) = {58,77};
//compile lines into a line loop for a mesh surface/region.
Line Loop(1) = {76, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, -77, -57, -56, -55, -54, -53, -52, -51, -50, -49, -48, -47, -46, -45, -44, -43, -42, -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1};

//Background region (Neumann boundary) points
cln=38.44;//characteristic length for background region
Point(78) = {-430.51,0.00,2968.10,cln};//far left upper point
Point(79) = {-430.51,0.00,-260.87,cln};//far left lower point
Point(80) = {515.07,0.00,2939.10,cln};//far right upper point
Point(81) = {515.07,0.00,-260.87,cln};//far right lower point
//make lines encompassing all the background points - counter clock wise fashion
Line(78) = {1,78};
Line(79) = {78,79};
Line(80) = {79,81};
Line(81) = {81,80};
Line(82) = {80,58};
//Add line loops and plane surfaces for the Neumann region
Line Loop(2) = {78, 79, 80, 81, 82, 77, -75, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65, -64, -63, -62, -61, -60, -59, -58, -76};
Plane Surface(1) = {1, 2};//Coarse mesh region surface

//Adding polygons
//end of polygons.
Plane Surface(2) = {1};//Fine mesh region surface

//Make a physical surface
Physical Surface(1) = {2, 1};

//End gmsh script
