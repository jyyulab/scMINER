#################################################################
## Implementation of Statistics::Distribution::uprob function  ##
## of perl in python.					       ##
## Author: Alireza Khatamian				       ##
## Email: akhatami@stjude.org				       ##
#################################################################
import math

def uprob(n):
        p = 0
        if abs(n) < 1.9:
                p = ( 1 + abs(n) * ( 0.049867347 + abs(n) * (0.0211410061 + abs(n) * 0.0032776263 + abs(n) * (0.0000380036 + abs(n) * (0.0000488906 + abs(n) * 0.000005383))))) ** (-16) / 2
        elif abs(n) <= 100:
                for i in range(18,0,-1):
                        p = i / (abs(n) + p)
                p = math.exp(-0.5 * abs(n) * abs(n)) / math.sqrt(2 * math.pi) / (abs(n) + p)
        if n < 0:
                p = 1 - p
        return p
