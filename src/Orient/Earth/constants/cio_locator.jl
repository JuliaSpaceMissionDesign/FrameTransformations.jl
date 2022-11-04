
const COEFFS_CIO2006_SP = SVector(
    94.00, 
    3808.35, 
    -122.68,
    -72574.11, 
    27.98,
    15.62,
    )

const COEFFS_CIO2006_S = [
    # j = 0 
    [
        # 1-10
        IAUSeries(-2640.73,  0.39, SVector(0, 0,  0,  0,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(  -63.53,  0.02, SVector(0, 0,  0,  0,  2, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(  -11.75, -0.01, SVector(0, 0,  2, -2,  3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(  -11.21, -0.01, SVector(0, 0,  2, -2,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    4.57,  0.00, SVector(0, 0,  2, -2,  2, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -2.02,  0.00, SVector(0, 0,  2,  0,  3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -1.98,  0.00, SVector(0, 0,  2,  0,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    1.72,  0.00, SVector(0, 0,  0,  0,  3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    1.41,  0.01, SVector(0, 1,  0,  0,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    1.26,  0.01, SVector(0, 1,  0,  0, -1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),

        # 11-20
        IAUSeries(    0.63, 0.00, SVector(1, 0,  0,  0, -1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    0.63, 0.00, SVector(1, 0,  0,  0,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.46, 0.00, SVector(0, 1,  2, -2,  3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.45, 0.00, SVector(0, 1,  2, -2,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.36, 0.00, SVector(0, 0,  4, -4,  4, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    0.24, 0.12, SVector(0, 0,  1, -1,  1, 0, -8,  12, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.32, 0.00, SVector(0, 0,  2,  0,  0, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.28, 0.00, SVector(0, 0,  2,  0,  2, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.27, 0.00, SVector(1, 0,  2,  0,  3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.26, 0.00, SVector(1, 0,  2,  0,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        
        # 21-30
        IAUSeries(    0.21,  0.00, SVector(0, 0,  2, -2,  0, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.19,  0.00, SVector(0, 1, -2,  2, -3, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.18,  0.00, SVector(0, 1, -2,  2, -1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    0.10, -0.05, SVector(0, 0,  0,  0,  0, 0,  8, -13, 0, 0, 0, 0, 0, -1)),
        IAUSeries(   -0.15,  0.00, SVector(0, 0,  0,  2,  0, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    0.14,  0.00, SVector(2, 0, -2,  0, -1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(    0.14,  0.00, SVector(0, 1,  2, -2,  2, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.14,  0.00, SVector(1, 0,  0, -2,  1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.14,  0.00, SVector(1, 0,  0, -2, -1, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        IAUSeries(   -0.13,  0.00, SVector(0, 0,  4, -2,  4, 0,  0,   0, 0, 0, 0, 0, 0,  0)),
        
        # 31-40
        IAUSeries(    0.11,  0.00, SVector(0, 0,  2, -2,  4, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        IAUSeries(   -0.11,  0.00, SVector(1, 0, -2,  0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        IAUSeries(   -0.11,  0.00, SVector(1, 0, -2,  0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ],

    # j = 1
    [
        # 1-3
        IAUSeries(   -0.07,   3.57, SVector(0,  0,  0,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(    1.73,  -0.03, SVector(0,  0,  0,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(    0.00,   0.48, SVector(0,  0,  2, -2,  3,  0, 0,  0,  0, 0, 0, 0, 0, 0))
    ],

    # j = 2 
    [
        # 1-10
        IAUSeries(743.52, -0.17, SVector(0,  0,  0,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( 56.91,  0.06, SVector(0,  0,  2, -2,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  9.84, -0.01, SVector(0,  0,  2,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -8.85,  0.01, SVector(0,  0,  0,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -6.38, -0.05, SVector(0,  1,  0,  0,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -3.07,  0.00, SVector(1,  0,  0,  0,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  2.23,  0.00, SVector(0,  1,  2, -2,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  1.67,  0.00, SVector(0,  0,  2,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  1.30,  0.00, SVector(1,  0,  2,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  0.93,  0.00, SVector(0,  1, -2,  2, -2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),

        # 11-20
        IAUSeries(  0.68, 0.00, SVector(1,  0,  0, -2,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.55, 0.00, SVector(0,  0,  2, -2,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  0.53, 0.00, SVector(1,  0, -2,  0, -2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.27, 0.00, SVector(0,  0,  0,  2,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.27, 0.00, SVector(1,  0,  0,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.26, 0.00, SVector(1,  0, -2, -2, -2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.25, 0.00, SVector(1,  0,  0,  0, -1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  0.22, 0.00, SVector(1,  0,  2,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.21, 0.00, SVector(2,  0,  0, -2,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  0.20, 0.00, SVector(2,  0, -2,  0, -1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),

        # 21-25
        IAUSeries(  0.17,   0.00, SVector(0,  0,  2,  2,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(  0.13,   0.00, SVector(2,  0,  2,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.13,   0.00, SVector(2,  0,  0,  0,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.12,   0.00, SVector(1,  0,  2, -2,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries( -0.11,   0.00, SVector(0,  0,  2,  0,  0,  0, 0,  0,  0, 0, 0, 0, 0, 0))
    ],

    # j = 3
    [
        # 1-4
        IAUSeries(    0.30, -23.42, SVector(0,  0,  0,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(   -0.03,  -1.46, SVector(0,  0,  2, -2,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(   -0.01,  -0.25, SVector(0,  0,  2,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0)),
        IAUSeries(    0.00,   0.23, SVector(0,  0,  0,  0,  2,  0, 0,  0,  0, 0, 0, 0, 0, 0))
    ],

    # j = 4
    [
        # 1-1
        IAUSeries(   -0.26,  -0.01, SVector(0,  0,  0,  0,  1,  0, 0,  0,  0, 0, 0, 0, 0,0))
    ], 

    # j = 5
    [
        IAUSeries(   -0.0,  -0.00, SVector(0,  0,  0,  0,  0,  0, 0,  0,  0, 0, 0, 0, 0,0))
    ],
];
