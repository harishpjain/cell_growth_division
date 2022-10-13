
import numpy as np
import math
import os
import sys
import glob
import re
print('importing paraview')
from paraview import simple

#The below functions: calc_eccentricity and calc_circularity are from Dennis' work
def shape_prop(path, rank):
	"""
		calculate circularity for phase field simulations from vtu files
		circularity is defined according to S. Aland paper (see https://tu-dresden.de/mn/math/wir/ressourcen/dateien/forschung/publikationen/pdf2012/benchmark_computations_of_diffuse_interface_models.pdf?lang=de)
		here circularity is 1 for a circle, smaller than 1 for ellipses
		many other define circularity as the inverse of this one (so it is not bounded anymore btw 0 and 1)
		t_min, t_max and jump are extracted so that b.c. do not cause any problem (might be problem dependent)
	"""
 
	def ResetSession():
		pxm = simple.servermanager.ProxyManager()
		pxm.UnRegisterProxies()
		del pxm
		simple.Disconnect()
		simple.Connect()
  
	pvd_file = path + '/phase_p' + str(rank) + '_.pvd'
	reader_pvd = simple.OpenDataFile(pvd_file)
	#save in a list all the timestep
	t_step = reader_pvd.TimestepValues
	print('step is ', t_step)
	t = t_step#*np.arange(t_min,t_max,jump)

	#==================================================================================================
 
	perimeter_temp = []
	a_temp = []
	b_temp = []
	times_temp = []
	perimeter_arr = np.zeros(len(t))
 
 
	simple.Delete(reader_pvd)
	del reader_pvd 
	
	#ResetSession()
	


	#calculating the eccentricity for different time step and then take the average value (it takes some time)
	for indt, time in enumerate(t):
		#ResetSession()
		filename = path + "data/phase_p" + str(rank) + "_" + '{:06.3f}'.format(time) + ".vtu"
		#read the vtu file
		reader=simple.XMLUnstructuredGridReader(FileName=filename)
		reader.UpdatePipeline()
  
		#=====================calculate the area of the domain where psi is bigger than zero================	
		#apply a clip to isolate only the values where phi is bigger than -0.7 (play with this value?)
		Clip1 = simple.Clip(Input = reader, ClipType="Scalar")
		Clip1.Scalars = ['POINTS', 'phi']
		Clip1.Value = -0.7

		#===================================================================================================
		#extract the value of the x coords with the calculator and take its range
		#Calculator1 = simple.Calculator()
		#Calculator1.Function = 'coordsX' #coorsdX'
		#bb = Calculator1.PointData
		#qq = bb.GetArray('Result')
		#xRange = qq.GetRange()


		#extract the value of the y coords with the calculator and take its range
		#Calculator2 = simple.Calculator()
		#Calculator2.Function = 'coordsY'
		#bb1 = Calculator2.PointData
		#qq1 = bb1.GetArray('Result')
		#yRange = qq1.GetRange()

		#a = abs(yRange[1]-yRange[0])/2.0; b = abs(xRange[1] - xRange[0])/2.0
		#a_temp.append(a)
		#b_temp.append(b)
		#===================================================================================================

		#set everything to one with the calculator
		#Calculator = simple.Calculator()
		#Calculator.Function = '1' 

		#integrate and take the value of interest
		#IntegrateVariable = simple.IntegrateVariables()
		#integral_data = simple.servermanager.Fetch(IntegrateVariable)
		#integral_point = integral_data.GetPointData()
		#area = integral_point.GetArray('Result').GetTuple(0)[0] #area (integral of 1 on domain where psi > 0)
		#===================================================================================================

		#===========================calculate the perimeter of the contour with value zero==================
		Contour = simple.Contour(Input = reader, PointMergeMethod="Uniform Binning")
		Contour.Isosurfaces = [-0.7]
		Contour.ContourBy = ['POINTS', 'phi']

		IntegrateVariable1 = simple.IntegrateVariables()
		integral_data1 = simple.servermanager.Fetch(IntegrateVariable1)
		integral_cell = integral_data1.GetCellData()
		try:
			perimeter_arr[indt] = integral_cell.GetArray('Length').GetTuple(0)[0] #perimeter of contour 1
		except:
			perimeter_arr[indt] = 0.0
		#===================================================================================================
			
		times_temp.append(time)
		simple.Delete(reader)
		simple.Delete(Clip1)
		simple.Delete(Contour)
		simple.Delete(IntegrateVariable1)
		#simple.Delete(integral_cell)
		#simple.Delete(integral_data1)
  
		del reader
		del Clip1
		del Contour
		del IntegrateVariable1
		del integral_cell
		del integral_data1

	return t,  perimeter_arr

def shape_prop2(path, rank):
	"""
		calculate circularity for phase field simulations from vtu files
		circularity is defined according to S. Aland paper (see https://tu-dresden.de/mn/math/wir/ressourcen/dateien/forschung/publikationen/pdf2012/benchmark_computations_of_diffuse_interface_models.pdf?lang=de)
		here circularity is 1 for a circle, smaller than 1 for ellipses
		many other define circularity as the inverse of this one (so it is not bounded anymore btw 0 and 1)
		t_min, t_max and jump are extracted so that b.c. do not cause any problem (might be problem dependent)
	"""
 
	def ResetSession():
		pxm = simple.servermanager.ProxyManager()
		pxm.UnRegisterProxies()
		del pxm
		simple.Disconnect()
		simple.Connect()
  
	pvd_file = path + '/phase_p' + str(rank) + '_.pvd'
	reader_pvd = simple.OpenDataFile(pvd_file)
	#save in a list all the timestep
	t_step = reader_pvd.TimestepValues
	print('step is ', t_step)
	t = t_step#*np.arange(t_min,t_max,jump)

	#==================================================================================================
 
	perimeter_arr = np.zeros(len(t))
 
 
	simple.Delete(reader_pvd)
	del reader_pvd 
	
	#ResetSession()
	filename = path + "data/phase_p" + str(rank) + "_" + '{:06.3f}'.format(0.0) + ".vtu"
	#read the vtu file
	reader=simple.XMLUnstructuredGridReader(FileName=filename)
	#=====================calculate the area of the domain where psi is bigger than zero================	
	#apply a clip to isolate only the values where phi is bigger than -0.7 (play with this value?)
	Clip1 = simple.Clip(Input = reader, ClipType="Scalar")
	Clip1.Scalars = ['POINTS', 'phi']
	Clip1.Value = -0.7

	#===========================calculate the perimeter of the contour with value zero==================
	Contour = simple.Contour(Input = reader, PointMergeMethod="Uniform Binning")
	Contour.Isosurfaces = [-0.7]
	Contour.ContourBy = ['POINTS', 'phi']

	IntegrateVariable1 = simple.IntegrateVariables()
	integral_data1 = simple.servermanager.Fetch(IntegrateVariable1)
	integral_cell = integral_data1.GetCellData()

	#calculating the eccentricity for different time step and then take the average value (it takes some time)
	for indt, time in enumerate(t):
		filename = path + "data/phase_p" + str(rank) + "_" + '{:06.3f}'.format(time) + ".vtu"
		reader.FileName(filename)
		reader.UpdatePipeline()
		try:
			perimeter_arr[indt] = integral_cell.GetArray('Length').GetTuple(0)[0] #perimeter of contour 1
		except:
			perimeter_arr[indt] = 0.0
		#===================================================================================================
  
	return t,  perimeter_arr


if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'phase_p*.pvd'
    out_dir = sys.argv[1] + 'shape'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
#num_cells = sys.argv[2]

ranks = []
for filename in glob.glob(file_pattern):
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))
ranks = np.array(ranks)

sorted_index = np.argsort(ranks)
def sort_ranks(elem):
    return elem.iloc[0]['rank']
ranks = ranks[sorted_index]

for ind_r, rank in enumerate(ranks):
	print(ind_r)
	ti, pe = shape_prop(sys.argv[1], rank)
	print(pe.shape)
	np.save(out_dir + '/perimeter_' + str(rank) + '.npy', pe)
	if rank == 0:
		np.save(out_dir + '/times.npy', ti)
