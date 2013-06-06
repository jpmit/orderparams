# testw.py
# get w value for fcc

import numpy as np

wig4s = np.array([ 0.104298,   # 4 -4  0
                  -0.164909,   # 4 -3 -1
                   0.186989,   # 4 -2 -2
                   0.156447,   # 3 -3  0
                  -0.0623298,  # 3 -2 -1
                  -0.0819482,  # 2 -2  0
                   0.141351,   # 2 -1 -1
                  -0.0670485,  # 1 -1  0
                   0.134097,   # 0  0  0
                   ])
wig4perms = np.array([6,12,6,6,12,6,6,6,1])
wig4is = [[4,-4, 0],
          [4,-3,-1],
          [4,-2,-2],
          [3,-3, 0],
          [3,-2,-1],
          [2,-2, 0],
          [2,-1,-1],
          [1,-1, 0],
          [0, 0, 0]]

# q4 for fcc
qlm = np.array([-0.0737554,
                0.0,
                0.0,
                0.0,
               -0.123416,
                0.0,
                0.0,
                0.0,
               -0.0737554
                ])

w4 = 0.0
for i in range(len(wig4s)):
    add = wig4s[i]*wig4perms[i]*qlm[wig4is[i][0]+4]*qlm[wig4is[i][1]+4]*qlm[wig4is[i][2]+4]
    print add
    w4 = w4 + add


denom = (sum(qlm**2))**1.5
print "top %.8f" %w4
print "bottom %.8f" %denom

print w4/denom # this gives -0.15932, which is correct!
