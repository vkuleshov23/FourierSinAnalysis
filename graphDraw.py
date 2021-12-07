import matplotlib.pyplot as p
import sys
import numpy as np

width = 4

for i in range(1, len(sys.argv)):
	words = open(sys.argv[i]).read().split()
	y = [float(v) for k, v in enumerate(words)]
	p.plot(y, linewidth = width)
	width = (width - 3)
	p.title(sys.argv[i])
p.grid()
p.show()