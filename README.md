# pi-cruncher
Approaches to finding a large number of digits of pi
c++ methods use gmp

I was able to find 700,000,000 digits of pi in 20 minutes with the chudnovsky method, but ran out of memory on my device.
This is likely a limitation of the big-num methods used, and may be improved with a more memory-efficient,
less general-use big-num class.

The first method tried is shown in pi-atan.py which uses an algorithm that breaks pi into arctangents, and computes these to 
a specified level of precision.
This method was repeated in c++ using gmp big-num class.

The second method uses the chudnovsky algorithm and also depends on gmp.
