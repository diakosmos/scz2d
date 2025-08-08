function ret0 = testf(in0)
     global x
     ret0 = 0;
     junk = subf(2);
     ret0 = x;

function ret1 = subf(in1)
     global x
     x = 2;
     ret1 = 3;
