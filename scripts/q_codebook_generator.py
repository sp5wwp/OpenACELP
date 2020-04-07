from Q1 import q1
from Q2 import q2
from Q3 import q3

import lbg
import time

t1=time.perf_counter()

#Q1 codebook
f=open("q1.txt", "w+")

cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q1, 16)

for i, c in enumerate(cb):
    f.write("%s\n" % (c))
    
f.close()

#Q2 codebook
f=open("q2.txt", "w+")

cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q2, 16)

for i, c in enumerate(cb):
    f.write("%s\n" % (c))
    
f.close()

#Q3 codebook
f=open("q3.txt", "w+")

cb, cb_abs_w, cb_rel_w = lbg.generate_codebook(q3, 16)

for i, c in enumerate(cb):
    f.write("%s\n" % (c))
    
f.close()

t2=time.perf_counter()

print("Time = %ds\n" % (int(t2-t1)))
