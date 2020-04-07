#----------------------------------------------
# LBG algorithm performance test
#
# SP5WWP 07/04/2020
#----------------------------------------------

from Q1 import q1
from Q2 import q2
from Q3 import q3

import lbg
import time

for size in range(1, 16+1):
    t1=time.perf_counter()

    #Q1 codebook
    f=open("q1.txt", "w+")

    cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q1, size)

    for i, c in enumerate(cb):
        f.write("%s\n" % (c))
        
    f.close()

    #Q2 codebook
    f=open("q2.txt", "w+")

    cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q2, size)

    for i, c in enumerate(cb):
        f.write("%s\n" % (c))
        
    f.close()

    #Q3 codebook
    f=open("q3.txt", "w+")

    cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q3, size)

    for i, c in enumerate(cb):
        f.write("%s\n" % (c))
        
    f.close()

    t2=time.perf_counter()

    print("%d\t%d\t%d\t%d" % (size, size, size, int(t2-t1)))
