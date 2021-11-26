#include nothing lol
import time
from datetime import timedelta
import math

# calculate arctan of n / d
# arctan = x + x^3 / 3 + x^5 / 5 ...
def arctan(n, d, num_digits):
    print(num_digits)
    term = ((10**num_digits) * n)// d # decimal in form of integer
    term = (10**num_digits)
    term *= n
    print(f'd: {d}')
    term //= d
    print(term)
    total = term
    d2 = 0
    while(term != 0):
        d2 += 1
        term = term * (-n*n) // (d*d)
        total += term // (2*d2 + 1)
        #print(term)
    print(f'arctan of {n} / {d} took {d2} iterations')
    return total

# calculate pi up to n digits of accuracy
# method 1: pi = 4(4(arctan(1/5)) - arctan(1/239))
def calculatePiMethod1(num_digits):
    extra_digits = 2*int(math.log10(num_digits)) + 1
    pi = 4 * (4 * arctan(1,5,num_digits+extra_digits) - arctan(1,239,num_digits+extra_digits))
    pi = pi // (10**extra_digits)
    return pi

# calculate pi up to n digits of accuracy
# method 2: pi = 4 * (12arctan(1/49) + 32arctan(1/57) - 5arctan(1/239) + 12arctan(1/110443))
def calculatePiMethod2(num_digits):
    extra_digits = 2*int(math.log10(num_digits))
    digits = num_digits + extra_digits
    a49 = arctan(1, 49, digits)
    a57 = arctan(1, 57, digits)
    a239 = arctan(1, 239, digits)
    a110443 = arctan(1, 110443, digits)
    pi = 4*(12*a49 + 32*a57 - 5*a239 + 12*a110443)
    pi = pi // (10**extra_digits)
    return pi

# calculate pi up to n digits of accuracy
# method 3: pi = 4 * (183 arctan(1/239) + 32 arctan(1/1023) - 68 arctan(1/5832)
#                + 12 arctan(1 / 110443) - 12 arctan(1 / 4841182) - 100 arctan(1 / 6826318))
def calculatePiMethod3(num_digits):
    extra_digits = 2*int(math.log10(num_digits))
    print(extra_digits)
    digits = num_digits + extra_digits
    a239 = arctan(1, 239, digits)
    print(a239 * 183)
    a1023 = arctan(1, 1023, digits)
    print(a1023 * 32)
    a5832 = arctan(1, 5832, digits)
    print(a5832 * -68)
    a110443 = arctan(1, 110443, digits)
    print(a110443 * 12)
    a4841182 = arctan(1, 4841182, digits)
    print(a4841182 * -12)
    a6826318 = arctan(1, 6826318, digits)
    print(a6826318 * -100)
    pi = 4*(183*a239 + 32*a1023 - 68*a5832 + 12*a110443 - 12*a4841182 - 100*a6826318)
    pi = pi // (10**extra_digits)
    return pi

if __name__ == '__main__':
    num_digits = 30
    start = time.time()
    pi = calculatePiMethod3(num_digits)
    s1 = str(pi)
    print(s1)
    end = time.time()
    time_str = time.strftime("%Hh%Mm%Ss", time.gmtime(end-start))
    print(f'method 3: {num_digits} of pi calculate in {time_str}')
    with open(f'{num_digits}.txt', 'wt') as outfile:
        outfile.write(f'{s1[0]}.{s1[1:]}')
    