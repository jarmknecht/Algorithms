import random
import math


def prime_test(N, k):

    # Has a time complexity of O(n^6)
    # Has a space complexity of O(n) since it calls mod_exp and is_carmichael

    prime = None  # Bool to keep track if prime
    carmichael = None  # Bool to keep track if carmichael number
    numsUsed = []  # Keeps track of numbers used to ensure they are unique
    i = 1  # Starting counter for while loop

    while i <= k and i is not 0:  # Loop that does k random trials the not 0 ensures that we don't loop negative
        a = random.randint(1, N - 1)  # Creates a random number that is assigned to a for Fermat's Test
        if a not in numsUsed:  # checks to make sure the number hasn't been used
            i += 1
            numsUsed.append(a)
            mod = mod_exp(a, N - 1, N)  # saves the value from modular exponentiation

            if not mod == 1:  # Failed the Fermat test so it is composite
                prime = False
                break
            else:
                prime = True
        else:
            i -= 1

    numsUsed.clear()  # Clears the list so it can be used if need to check carmichael numbers
    i = 1
    if prime:  # May be prime or carmichael number so need to check k times again
        while i <= k and i is not 0:  # Loop that does k random trials the not 0 ensures that we don't loop negative
            a = random.randint(1, N - 1)
            if a not in numsUsed:  # checks to make sure the number hasn't been used so it is unique
                i += 1
                numsUsed.append(a)
                carmichael = is_carmichael(N, a)
                if carmichael:  # number is a carmichael number so break loop
                    prime = False
                    break
                else:
                    prime = True
            else:
                i -= 1

    if prime:
        return 'prime'
    elif not prime and carmichael:
        return 'carmichael'
    else:
        return 'composite'


def mod_exp(x, y, N):

    # Has a time O(n^3) n is the number of bits in x, y or N whichever is biggest
    # Space O(n) cause each recursive call of n is stored on the stack and then deleted when returned with each call

    if y == 0:  # If the exponent is 0 from flooring it stop
        return 1
    z = mod_exp(x, math.floor(y / 2), N)  # Recurse until you get y = 0
    if y % 2 == 0:  # if y is even
        return pow(z, 2) % N
    else:  # if y is odd
        return x * pow(z, 2) % N


def probability(k):

    # Has space complexity of this function is O(n) where n is the bit size of k
    # cause as k gets larger more space is needed to store the float
    # Has time complexity of O(k)

    return 1 - (1 / math.pow(2, k))  # calculates the probability that the number is prime


def is_carmichael(N, a):
    # The time complexity here is O(n^4) since it is O(n) and calls mod_exp which is O(n^3)
    # n is the bit size of y which is the number N - 1
    # Space complexity here is O(n) since it calls mod_exp
    y = N - 1
    while True:
        if y == 1:   # This breaks if y gets to one since 1 is odd
            break
        mod = mod_exp(a, y, N)
        if not mod == 1:  # if mod_exp is not one check if its prime or carmichael
            if mod == N - 1:  # This means it is prime
                return False
            else:  # This means it is a Carmichael number
                return True
        if not (y % 2) == 0:  # Checks to see if the exponent is divisible by two if not stop
            break
        y = y / 2  # Divide exponent by 2 which is the same as sqrt the whole equation
    return False
