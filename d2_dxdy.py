"""
Script to calculate the 2nd mixed derivative of data grids or functions
using the 5 points differentiaton formula 
Author: Aldo Arriola
email : aarriolac@uni.pe
"""
from __future__ import print_function

import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import pandas as pd

import time

def local_submatrix(A, r, c, SF=2):
	'''
	A main matrix
	(r, c) index of the center
	SF smooth factor
	(2*SF+1, 2*SF+1) shape of matrix
	'''
	return A[r-SF:r+SF+1 , c-SF: c+SF+1] # 
	
def dcentral5pts(y, h):
	''' 5 points formula'''
	#sh = abs(x[1]-x[0])
	return (-y[4]+8*y[3] - 8*y[1]+y[0])/(12*h)
	

def _2ndMixedDerivative(X, Y, Z, SF=2):
	rows, cols = X.shape
	dx = []
	dy = []
	dz = []
	for i in np.arange(SF, rows-SF):
		for j in np.arange(SF, cols-SF):
			Xs = local_submatrix(X, i, j)
			Ys = local_submatrix(Y, i, j)
			Zs = local_submatrix(Z, i, j)
			dy_ = []
			hy = abs(Ys[0][0] - Ys[1][0])
			hx = abs(Xs[0][0] - Xs[0][1])
			for idx in [0,1,3,4]: # index 2 is the center
				#print (Zs[:,idx])
				dy_.append(dcentral5pts(Zs[:,idx],hy))
			dy_ = dy_[0:2] + [ 0 ] +dy_[2:] 
			Fxy = dcentral5pts(dy_,hx)
			dx.append(X[i][j])
			dy.append(Y[i][j])
			dz.append(Fxy)
	
	return np.array(dx), np.array(dy), np.array(dz)
					
def close_floats(n1, n2, atol=0.01):
	return abs(n1-n2) <= atol

# 
def look_into(values, array, atol=0.01, extended=False):
	'''
	search for a value = [x, y] into a list = [[x1,y1,z1], [x2,y2,z2], ....]
	values is a tuple (x_i, y_i)'''
	for row in array:
		#print(row, values)
		if close_floats(values[0], row[0], atol=atol) and close_floats(values[1], row[1], atol=atol):
			return row[2]
			
	else :
		if extended:
			for row in array:
				if close_floats(values[0], row[0], atol=atol) and close_floats(values[0], row[1],atol=atol):
					return row[2]
		else:
			return False


		

def mesh_M(Ha, Hb, M, extended=False):
	'''creates 3 grids X, Y and Z(X,Y) from lists Ha, Hb, M'''
	
	Ha_grid, Hb_grid = np.meshgrid(np.unique(Ha), np.unique(Hb))
	Mmesh = np.zeros(Ha_grid.shape, dtype="float")*np.nan
	
	i, j = 0, 0
	tol=min(abs(Ha[0]-Ha[1]), abs(Hb[1]-Hb[0]))*0.1
	
	for row in Mmesh:
		j = 0
		for col in row:
			aux = look_into([Ha_grid[i][j], Hb_grid[i][j]] , np.transpose([Ha,Hb,M]), atol=tol, extended=extended)
			if aux:
				Mmesh[i][j] = aux
			j += 1
		i += 1
	
	return Ha_grid, Hb_grid, Mmesh


def function_F(x, y):
	# example function
	return np.power(x, 2)*np.power(y,3)


def function_d2Fdxdy(x,y):
	return 6*x*np.power(y,2)


def saveontxt(data):
	'''data is a list of 3D points'''
	f = open("2nd_deriv.txt", "w")
	
	for row in data:
		f.write(str(row[0])+","+str(row[1])+ ","+str(row[1])+"\n")
	f.close()
	
	
if __name__ == "__main__":
	savedata = False
	t0 = time.time()
	valuelim = 5
	npts = 51
	x = np.linspace(-valuelim, valuelim,npts)
	y = np.linspace(-valuelim*4/5, valuelim*4/5,npts)

	x_grid, y_grid = np.meshgrid(x, y)
	z_grid = function_F(x_grid, y_grid)
	
	fig = plt.figure("function")
	ax = plt.axes(projection="3d")
	ax.set_title("f(x, y)")
	ax.plot_surface(x_grid,y_grid,z_grid)#, c=z_grid)
	
	fig = plt.figure("analytical deriv function")
	ax2 = plt.axes(projection="3d")
	ax2.set_title("analytical d2f(x, y)")
	d2z_grid = function_d2Fdxdy(x_grid, y_grid)
	ax2.plot_surface(x_grid,y_grid,d2z_grid)#, c=z_grid)
	
	print('Calculating Second mixed derivative...', end=" " )
	d_x, d_y, d_z = _2ndMixedDerivative(x_grid, y_grid, z_grid)
	print("> Done.")
	#print("lenssss", len(datP), len(datP[0]))
	
	dxgrid, dygrid, dzgrid = mesh_M(d_x, d_y, d_z)
	
	fig = plt.figure("Numerical 2nd ord mix deriv")
	
	ax4 = plt.axes(projection='3d')
	ax4.set_title(" d2xy(x, y)")
	ax4.set_ylabel("x")
	ax4.set_xlabel("y")
	ax4.plot_surface(dxgrid,dygrid, dzgrid, cmap="RdGy")
	tf = time.time()
	print("TIME ELAPSED", (tf-t0)/60, "mins")

	if savedata :
		print("Saving data on txt...", end=" ")
		saveontxt(datP)
		print("> Done.")
	plt.show()

	
