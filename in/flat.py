#!/usr/bin/python
import numpy as np
width = 200
Nwidth = 30

x = np.linspace(-width, width, Nwidth)
y = np.copy(x)
X,Y = np.meshgrid(x,y)
Z = np.copy(X) * 0

X = np.reshape(X,-1)
Y = np.reshape(Y,-1)
Z = np.reshape(Z,-1)

with open("flat.ckx","w") as f:
	f.write("%d\nMedium properties=type:I:1:pos:R:3\n"%len(X))
	for i in range(len(X)):
		f.write("2 %e %e %e\n"%(X[i],Y[i],Z[i]))
