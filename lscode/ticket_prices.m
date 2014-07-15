% Example: estimation of ticket prices
Y = [ 7.00 ;  11.55 ;  15.65 ;  18.35 ];
Phi = [ 1,   82.842 ;
        1,  147.480 ;
        1,  229.408 ;
        1,  266.341  ];
thetaLS = pinv(Phi)*Y
priceToVicenza = thetaLS(1) + thetaLS(2)*199.138
