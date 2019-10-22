import numpy as np
import sys
sys.path.append('/home/STG_04/WORK/SATUTILITIES/')
import satmaptools
from scipy.io.idl import readsav
from scipy.ndimage import convolve
from scipy.special import binom
import matplotlib.pyplot as plt
import h5py

import glob
#gas='NO2'
#gas='SO2'

def makeh5andeps(filenameidlfull,gas,datasetname,cbmax,cbinter):
	pathin='/home/D3_WORK1/claudia/INECC/Matrices_mapas/originales/%s/%s_idl/' %(gas,gas)
	pathout='/home/D3_WORK1/claudia/INECC/Matrices_mapas/smoothed/%s' %(gas)
	cblabel='%s  [molec./cm-2]' % (gas)
	cbmin=0.0

	filenameidl=filenameidlfull.split('/')[-1]
	filenameout=filenameidl.replace('.idl','.h5')

	filenameouteps=filenameidl.replace('.idl','.eps')
	dataorg=readsav(pathin+filenameidl)
	variables= dataorg.keys()


	def fill99(matrixorg):
		matrix=1.0*matrixorg
		nx=len(matrix[:,0])
		ny=len(matrix[0,:])
		for ix in  range(nx):
			for jy in  range(ny):
				if matrixorg[ix,jy]<-98.0:
					imin=ix-1
					if imin < 0:
						imin=0
					jmin=jy-1
					if jmin < 0:
						jmin=0

					imax=ix+2
					if imax > nx:
						imax=nx
					jmax=jy+2
					if jmax > ny:
						jmax=ny
					submatrix=matrixorg[imin:imax,jmin:jmax]
					index=np.where(submatrix[:,:] > -98.0)
					#print index
					matrix[ix,jy]=np.average(submatrix[index])
				



		return matrix


	def smoothed(matrixin,nsmooth):
		kernel=np.ones((nsmooth,nsmooth),dtype=float)
		kernel=kernel/np.sum(kernel[:,:])
		matrixout=convolve(matrixin,kernel,mode='nearest')
		return matrixout


	def smoothed_binom(matrixin,nsmooth):
		kernel=np.zeros((nsmooth,nsmooth),dtype=float)
		for i in range(nsmooth):
			for j in range(nsmooth):
				kernel[i,j]=binom(nsmooth,i)*binom(nsmooth,j)	
		kernel=kernel/np.sum(kernel[:,:])
		matrixout=convolve(matrixin,kernel,mode='nearest')
		return np.array(matrixout)
	










	fh5=h5py.File(pathout+'/'+filenameout,'w')
	for variable in variables:
		print variable
		fh5.create_dataset(variable,data=dataorg[variable])
	lats=dataorg['latnew']
	lons=dataorg['lonnew']
	#matrixorg=dataorg['no2matrix']
	#matrixorg=dataorg['so2matrix']
	matrixorg=dataorg[datasetname]
	print matrixorg
	matrixnew=fill99(matrixorg)
	matrixsmooth3=smoothed(matrixnew,3)
	matrixsmoothbinom5=smoothed_binom(matrixnew,5)

	fh5.create_dataset('smoothed_matrix',data=matrixsmoothbinom5)

	'''	
	plt.subplot(2,2,1)
	plt.contourf(lons,lats,matrixorg,100)
	plt.colorbar()
	plt.subplot(2,2,2)
	plt.contourf(lons,lats,matrixnew,100)
	plt.colorbar()
	plt.subplot(2,2,3)
	plt.contourf(lons,lats,matrixsmooth3,100)
	plt.colorbar()
	plt.subplot(2,2,4)
	plt.contourf(lons,lats,matrixsmoothbinom5,100)
	plt.colorbar()

	plt.show()
	'''	
	satmaptools.map_plot(lons,lats,matrixsmoothbinom5,cblabel=cblabel,cbmin=cbmin,cbmax=cbmax,cbinter=cbinter, alpha=0.75,nparalel=2,nmeridian=2,show=False)

	plt.savefig(filenameouteps)

	fh5.close()


gas='NO2'
cbmax=5.0E15
cbinter=1.0E14
datasetname='no2matrix'
filenameidl='no2matrix_TROP_RM_2004_2018_f.idl'
#makeh5andeps(filenameidl,gas,datasetname,cbmax,cbinter)

gas='SO2'
cbmax=10.0E16
cbinter=10.0E14
datasetname='so2matrix'
filenameidl='so2_PBL_matrix_TROP_RM_2004_2018_f.idl'
makeh5andeps(filenameidl,gas,datasetname,cbmax,cbinter)

gas='HCHO'
cbmax=3.0E16
cbinter=2.0E14
datasetname='hchomatrix_aver'
filenameidl='hchomatrix_TROP_RM_2004_2018_f.idl'
#makeh5andeps(filenameidl,gas,datasetname,cbmax,cbinter)
