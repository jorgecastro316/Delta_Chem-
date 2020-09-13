import re
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import subprocess as sub
from scipy.signal import find_peaks
from scipy import optimize
import time 
from scipy.optimize import curve_fit
from pylab import *
from prettytable import PrettyTable
from itertools import combinations
import sys


def scaling_transformation(vector,scale_factor):
	S = np.matrix([[scale_factor,0,0,0],[0,scale_factor,0,0],[0,0,scale_factor,0],[0,0,0,1]])
	vector = np.append(vector,1)
	scale = np.dot(vector,S)
	scaled_vector = []
	scaled_vector.append(scale.item(0))
	scaled_vector.append(scale.item(1))
	scaled_vector.append(scale.item(2))
	return scaled_vector
if __name__ == "__scaling_transformation__":
    scaling_transformation()
    
def translation_transformation(vector,direction):
	T = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[direction[0],direction[1],direction[2],1]])
	vector = np.append(vector,1)
	translate = np.dot(vector,T)
	translated_vector = []
	translated_vector.append(translate.item(0))
	translated_vector.append(translate.item(1))
	translated_vector.append(translate.item(2))
	return translated_vector
if __name__ == "__translation_transformation__":
    translation_transformation()

def rotation_matrix(a,b):
	c = np.dot(a,b)
	if abs(c + 1) <= 0.1:
		a = a*(-1)
	v = np.cross(a,b)
	s = np.linalg.norm(v)
	
	vx = np.matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
	I = np.matrix([[1,0,0],[0,1,0],[0,0,1]])
	R = I + vx + (np.dot(vx,vx)*(1/(1+c)))

	v = np.cross(a,b)
	s = np.linalg.norm(v)
	c = np.dot(a,b)
	vx = np.matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
	I = np.matrix([[1,0,0],[0,1,0],[0,0,1]])
	R = I + vx + (np.dot(vx,vx)*(1/(1+c)))
	return R
if __name__ == "__rotation_matrix__":
    rotation_matrix()

def reflection_matrix(normal_vector):
	px = normal_vector[0]
	py = normal_vector[1]
	pz = normal_vector[2]
	RM = np.matrix([[(-(px**2) + (pz**2) + (py**2)),(-2*px*py),(-2*px*pz)],
		[(-2*py*px),(-(py**2) + (px**2) + (pz**2)),(-2*py*pz)],
		[(-2*pz*px),(-2*pz*py),(-(pz**2) + (py**2) + (px**2))]])
	return(RM)

if __name__ == "__reflection_matrix__":
    rotation_matrix()
    
def Volume_of_Tetrahedron(v0,v1,v2,v3):
	#Compute the volume of a tetrahedron given the 4 vertices
	tetrahedron = []
	tetrahedron.append(v0)
	tetrahedron.append(v1)
	tetrahedron.append(v2)
	tetrahedron.append(v3)
	tet = np.asarray(tetrahedron)
	tet = np.insert(tet,3,1,axis=1)
	matrix_tet = np.matrix(tet)
	tetrehedron_volume = (np.linalg.det(matrix_tet))/6
	return(tetrehedron_volume)
if __name__ == "__Volume_of_Tetrahedron_":
    Volume_of_Tetrahedron()

def Baricenter(vertices):
	#Calculates the baricenter of a set of points
	lenghts=[]
	number_of_points = len(vertices)
	x = 0
	y = 0
	z = 0
	baricenter = []
	for vertex in vertices:
		x = x + vertex[0]
		y = y + vertex[1]
		z = z + vertex[2]

	x = x/number_of_points
	y = y/number_of_points
	z = z/number_of_points

	baricenter.append(x)
	baricenter.append(y)
	baricenter.append(z)

	baricenter = np.asarray(baricenter)

	#Compute the distance between the baricenter and each point

	for vertex in vertices:
		vertex = np.asarray(vertex)
		lenghts.append(np.linalg.norm(vertex - baricenter))

	#Obtain the maximun distance from the baricenter to the set of points

	max_distance = max(lenghts)

	return(baricenter,max_distance)
if __name__ == "__Baricenter_":
    Baricenter()

def Slerp(P0,P1,t):
	 #Compute the cosine of the angle between the two vectors.
	dot = np.dot(P0,P1)
	# Compute the angle
	Omega_0 = np.arccos(dot) #angle between input vectors
	sin_omega_0 = np.sin(Omega_0) #angle between v0 and result
	#print(sin_omega_0)
	#Create an array of the angle times the values of the parameter t
	Omega = Omega_0*t
	sin_Omega = np.sin(Omega)

	s0 = np.cos(Omega) - dot*sin_Omega/sin_omega_0
	s1 = sin_Omega/sin_omega_0 
	path = (s0[:,np.newaxis] * P0[np.newaxis,:]) + (s1[:,np.newaxis] * P1[np.newaxis,:])
	x = []
	y = []
	z = []
	for i in path:
		x.append(i[0])
		y.append(i[1])
		z.append(i[2])
	return(x,y,z)
if __name__ == "__Slerp_":
    Slerp()

def PSD_Solve(structure_file,psd_time,psd_grid_resolution,psd_test_particle_size):
	process = sub.Popen(("psd",structure_file,str(psd_grid_resolution),str(psd_test_particle_size)))
	time.sleep(int(psd_time))
	process.kill()

	
if __name__ == "__PSD_Solve__":
    PSD_Solve()

def PSD_analysis(pore_diameter_file,peak_finder_height,peak_finder_prominence):
	Pore_Diameter = []
	Probability_Density = []
	proba_pore = {}
	#proba_pore is a dictionary whose keys are the probability density 
	#and the values are the pore diameters
	pores = open(pore_diameter_file)
	for line in pores:
		pore = re.compile(r'\s+')
		pore = pore.split(line)
		pore.pop(-1)
		pore_diameter = float(pore[0])
		probability_density = float(pore[1])
		proba_pore[probability_density] = pore_diameter
		Pore_Diameter.append(pore_diameter)
		Probability_Density.append(probability_density)
	x = Pore_Diameter
	y = Probability_Density
	Pore_Diameter = np.array(Pore_Diameter)
	Probability_Density = np.array(Probability_Density)
	peaks,properties = find_peaks(Probability_Density,prominence = float(peak_finder_prominence),height = float(peak_finder_height))
	number_of_porosities = len(peaks)
	print("Number of Porosities: ",number_of_porosities)
	A = properties['peak_heights']
	mu = []
	for i in range(number_of_porosities):
		mu.append(proba_pore[A[i]])

	def gauss(x,mu,sigma,A):
		return A*exp(-(x-mu)**2/2/sigma**2)
	
	def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
		return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

	def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3):
		return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

	def tetramodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3,mu4,sigma4,A4):
		return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)+gauss(x,mu4,sigma4,A4)

	expected = []
	for j in range(number_of_porosities):
		expected.append(mu[j])
		expected.append(0.2)
		expected.append(A[j])


	if number_of_porosities == 1:
		params,cov = curve_fit(gauss,x,y,expected)
		model = 1

	elif number_of_porosities == 2:
		params,cov = curve_fit(bimodal,x,y,expected)
		model = 2

	elif number_of_porosities == 3:
		params,cov = curve_fit(trimodal,x,y,expected)
		model = 3

	elif number_of_porosities == 4:
		params,cov = curve_fit(tetramodal,x,y,expected)
		model = 4

	sigma=sqrt(diag(cov))
	table = PrettyTable()
	table.field_names = [" ", "Parameters","Sigma"]
	j,k,l = 0,1,2
	pores = []
	dsv = []
	for m in range(1,int((len(params)/3)+1)):
		table.add_row(["Pore radius "+str(m),params[j],sigma[j]])
		table.add_row(["sigma"+str(m),params[k],sigma[k]])
		table.add_row(["Height Peak "+str(m),params[l],sigma[l]])
		pores.append(params[j])
		dsv.append(sigma[j])
		j+=3
		k+=3
		l+=3

	print(table)

	if number_of_porosities == 1:
		plt.plot(x,gauss(x,*params),color='k',lw=2,label='Gaussian Fit')

	elif number_of_porosities == 2:
		plt.plot(x,bimodal(x,*params),color='k',lw=2,label='Gaussian Fit')

	elif number_of_porosities == 3:
		plt.plot(x,trimodal(x,*params),color='k',lw=2,label='Gaussian Fit')

	elif number_of_porosities == 4:
		plt.plot(x,tetramodal(x,*params),color='k',lw=2,label='Gaussian Fit')


	plt.bar(Pore_Diameter,Probability_Density,0.025,color = 'm',alpha = 0.25,edgecolor = 'k',linewidth = 0.5)
	plt.xlabel('Pore Radius')
	plt.ylabel('Probability Density')
	plt.show()
	pores_dat = {}
	#pores = [6,10,12.5]
	for i in range(len(pores)):
		pores_dat[pores[i]*2] = 2*dsv[i]

	return(pores_dat)
	
if __name__ == "__PSD_analysis__":
    PSD_analysis()


def grouping (pore_diameter,points):
	pore_neighbourhood = []
	control_list = []
	points = {tuple(i) : None for i in points}
	points_list = list(points.keys())
	for point in points_list:
		neighbours = []
		neighbours.append(point)
		neigh_list = [point for point in points_list]
		neigh_list.remove(point)
		for neigh in neigh_list:
			dist_mid_neigh = np.linalg.norm(np.array(neigh) -np.array(point))
			if dist_mid_neigh <= (pore_diameter):
				neighbours.append(neigh)
				points_list.remove(neigh)
	
		control_list.append(len(neighbours))
		pore_neighbourhood.append(neighbours)

	#Get the baricenters
	baricenters_neighbourhood = []
	for neighbours in pore_neighbourhood:
		baricenter,max_distance = Baricenter(neighbours)
		baricenters_neighbourhood.append(tuple(baricenter))
			
	#Check if there only one point per neighbourhood
	counter = 0
	for number_of_neighbours in control_list:
		if number_of_neighbours == 1:
			counter += 1

	if counter == len(control_list):
		return baricenters_neighbourhood

	else:
		return grouping(pore_diameter,baricenters_neighbourhood)
if __name__ == "__grouping__":
    grouping()

def Polyhedron_Alignment(objective_polyhedron,reference_polyhedron,radius_c,c_vertex,reference_vertices,sphr_vdw_corr):
	#Choose the biggest regular face of the hull
	c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces_pr,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = reference_polyhedron
	faces_polyhedron,polyhedron_vertices,num_faces_po,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas = objective_polyhedron
	irreducible_volumes_points = []
	faces_list = list(c_hull_v_area.keys())

	if polytope_type == 'Cantellation':
		num_small_reg_faces = abs(num_faces_pr[0] - num_faces_pr[1])
		if len(polyhedron_vertices) == num_small_reg_faces:
			interpretation = 'VERTICES'
		else:
			interpretation = 'FACES'
					
	if polytope_type == 'Truncation':
		if len(polyhedron_vertices) < num_vertices_polytope:
			interpretation = 'VERTICES'
		elif len(polyhedron_vertices) >= num_vertices_polytope:
			interpretation = 'EDGES'

	#Create a vector normal to the surface of the hull's max polygon
	vertices_max_polygon = {}
	for edge in max_polygon:
		e1 = reference_vertices[edge[0]]
		e2 = reference_vertices[edge[1]]
		vertices_max_polygon[e1] = None
		vertices_max_polygon[e2] = None	
	vertices_max_polygon = list(vertices_max_polygon.keys())
	ia = 0
	ib = 1
	ic = 2
	norm = 0.0
	while norm == 0.0:
		vp1 = np.array(vertices_max_polygon[ia])
		vp2 = np.array(vertices_max_polygon[ib])
		vp3 = np.array(vertices_max_polygon[ic])
		vp2_1 = np.subtract(vp1,vp2)
		vp2_3 = np.subtract(vp3,vp2)
		p_vector = np.cross(vp2_1,vp2_3)
		norm = float(np.linalg.norm(p_vector))
		ib += 1
		ic += 1

	#Center the vector at the c_vertex by translation
	c_direction = np.subtract(np.array(c_vertex),np.array((0,0,0)))
	p_vector = translation_transformation(p_vector,c_direction)
	bar_max_polygon,usls = Baricenter(vertices_max_polygon)

	#Create a vector normal to the surface of the hull's min regular polygon
	vertices_min_reg_polygon = {}
	for edge in min_reg_polygon:
		e1 = reference_vertices[edge[0]]
		e2 = reference_vertices[edge[1]]
		vertices_min_reg_polygon[e1] = None
		vertices_min_reg_polygon[e2] = None
	vertices_min_reg_polygon = list(vertices_min_reg_polygon.keys())

	mia = 0
	mib = 1
	mic = 2
	mnorm = 0.0
	while mnorm == 0.0:
		vp1 = np.array(vertices_min_reg_polygon[mia])
		vp2 = np.array(vertices_min_reg_polygon[mib])
		vp3 = np.array(vertices_min_reg_polygon[mic])
		vp2_1 = np.subtract(vp1,vp2)
		vp2_3 = np.subtract(vp3,vp2)
		mp_vector = np.cross(vp2_1,vp2_3)
		mnorm = float(np.linalg.norm(mp_vector))
		mib += 1
		mic += 1

	#Center the vector at the c_vertex by translation
	mp_vector = translation_transformation(mp_vector,c_direction)
	bar_min_reg_polygon,usls = Baricenter(vertices_min_reg_polygon)

	if polytope_type == 'Regular_b':
		polytope_dir_vector = (np.subtract(mp_vector,c_vertex))*-1

	else:
		polytope_dir_vector = np.subtract(p_vector,c_vertex)

	polytope_dir_vector = (polytope_dir_vector/np.linalg.norm(polytope_dir_vector))

	#Hasta aquí todo bien

	#Analysis of the Fit Polyhedron
	polyhedron_i_baricenters = {}
	for face in faces_polyhedron:
		face_verts = [polyhedron_vertices[i] for i in face]
		face_baricenter, usls = Baricenter(face_verts)
		polyhedron_i_baricenters[tuple(face)] = face_baricenter

	area_diff = {}
	for area in areas_faces:
		area_d = abs(area - areas[0])
		area_diff[area] = area_d

	min_area_diff = min(list(area_diff.values()))

	for area in area_diff:
		if area_diff[area] == min_area_diff:
			fit_area = area
	bar_polygon = bar_max_polygon

	if polytope_type == 'Regular_b':
		fit_area = small_face_area
		bar_polygon = bar_min_reg_polygon
		

	if polytope_type == 'Cantellation':
		fit_area = big_face_area
		bar_polygon = bar_max_polygon

	fit_face_list = areas_faces[fit_area] #List containing the reference faces selected by their area
	# bar_polygon ---- > Baricenter of polytope's reference face. 
	#Select the face in the polyhedron which baricenter is the closest to the baricenter of the polytope's reference face
	face_b_dist = {}
	for face in fit_face_list:
		face_bar = np.array(polyhedron_i_baricenters[face])
		dist_bars = np.linalg.norm(face_bar - bar_polygon)
		face_b_dist[dist_bars] = face

	min_b_dist = min(list(face_b_dist.keys()))
	fit_face = face_b_dist[min_b_dist]#Polyhedron Reference face
	bar_f = polyhedron_i_baricenters[fit_face]
	fit_face_verts = {}

	fit_face_reg_b= fit_face
	fit_face_ct = fit_face
	for point in fit_face:
 		fit_face_verts[tuple(polyhedron_vertices[point])] = None

	fit_face_verts = list(fit_face_verts.keys())

	#Create a vector perpendicular to the face
	#Choose any 3 points/position vectors
	fa = 0
	fb = 1
	fc = 2
	normf = 0.0
	while normf == 0.0:
		vf1 = np.array(fit_face_verts[fa])
		vf2 = np.array(fit_face_verts[fb])
		vf3 = np.array(fit_face_verts[fc])
		vf2_1 = np.subtract(vf1,vf2)
		vf2_3 = np.subtract(vf3,vf2)
		f_vector = np.cross(vf2_1,vf2_3)
		normf = float(np.linalg.norm(f_vector))
		fb += 1
		fc += 1
	#Center the vector at the c_vertex by translation
	c_direction = np.subtract(np.array(c_vertex),np.array((0,0,0)))
	f_vector = translation_transformation(f_vector,c_direction)
	polyhedron_dir_vector = np.subtract(f_vector,c_vertex)
	polyhedron_dir_vector = polyhedron_dir_vector/np.linalg.norm(polyhedron_dir_vector)
########################################################################FIRST_ROTATION################################################################################
	R = rotation_matrix(polyhedron_dir_vector,polytope_dir_vector)
	#Hasta aquí todo bien
	#R = np.array([[1,0,0],[0,1,0],[0,0,1]])

	#translation to the origin of coordinates
	origin_verts = []
	o_direction = np.subtract(np.array((0,0,0)),np.array(c_vertex))
	for vertex in polyhedron_vertices:
		o_vert = translation_transformation(vertex,o_direction)
		origin_verts.append(o_vert)

	#Rotation to align vectors
	aligned_vertices = []
	for vertex in origin_verts:
		align = np.dot(R,np.transpose(vertex))
		aligned_vertex = []
		aligned_vertex.append(align.item(0))
		aligned_vertex.append(align.item(1))
		aligned_vertex.append(align.item(2))
		aligned_vertices.append(aligned_vertex)
	
	#Geometric correction ----> the points corresponding to the polyhedron must be treated now as spheres with wdw radii
	radius_d = radius_c - float(sphr_vdw_corr)
	#radius_d = radius_c
	scale_factor = radius_d/radius_c
	scaled_vertices = []
	for vertex in aligned_vertices:
		scaled_vertex = scaling_transformation(vertex,scale_factor)
		scaled_vertices.append(scaled_vertex)
###############################################################PILAS
	
	#Back translation
	transformed_vertices = []
	for vertex in scaled_vertices:
		bt_vertex = translation_transformation(vertex,o_direction*(-1))
		transformed_vertices.append(bt_vertex)
 ####################################################END OF THE FIRST ROTATION###########################################################
	#After the first transformation redefine the position of the fit face

	#Create a new fit list with the polygons of the polyhedron which area is similar to the biggest regular polygon of the polytope
	area_diff_reg = {}
	for area in areas_faces:
		area_d = abs(area - areas[1])
		area_diff_reg[area] = area_d

	min_area_diff_reg = min(list(area_diff_reg.values()))


	for area in area_diff_reg:
		if area_diff_reg[area] == min_area_diff_reg:
			fit_area_reg = area

	fit_face_reg_list = areas_faces[fit_area_reg]
	polyhedron_i_baricenters_2 = {}

	for face in faces_polyhedron:
		face_verts = []
		for vert in face:
			face_verts.append(transformed_vertices[vert])
		face_baricenter, usls = Baricenter(face_verts)
		polyhedron_i_baricenters_2[tuple(face)] = face_baricenter

	face_b_dist_2 = {}
	for face in fit_face_list:
		face_bar = np.array(polyhedron_i_baricenters_2[face])
		dist_bars = np.linalg.norm(face_bar - bar_max_polygon)
		face_b_dist_2[dist_bars] = face

	min_b_dist_2 = min(list(face_b_dist_2.keys()))
	fit_face_2 = face_b_dist_2[min_b_dist_2]
	bar_f_2 = polyhedron_i_baricenters_2[fit_face_2]
	fit_face_verts_2 = {}

	for point in fit_face_2:
 		fit_face_verts_2[tuple(transformed_vertices[point])] = None

	fit_face_verts_2 = list(fit_face_verts_2.keys())

	poly_pair_dist = {}
	for pair in max_polygon:
		mp1 = reference_vertices[pair[0]]
		mp2 = reference_vertices[pair[1]]
		pair_verts = (mp1,mp2)
		dist_mp = np.linalg.norm(np.array(mp1) - np.array(mp2))
		poly_pair_dist[tuple(pair_verts)] = dist_mp

	max_dist_mp = max(list(poly_pair_dist.values()))

	large_edges_list = []
	for pair in poly_pair_dist:
		if abs(poly_pair_dist[pair] - max_dist_mp) < 0.01:
			large_edges_list.append(pair)

	max_poly_edge = large_edges_list[0]
	mxp1 = np.array(max_poly_edge[0])
	mxp2 = np.array(max_poly_edge[1])
	mid_mp = np.array([0.5*(mxp1[0] + mxp2[0]),0.5*(mxp1[1] + mxp2[1]),0.5*(mxp1[2] + mxp2[2])])
	mp_vector = np.subtract(mid_mp,bar_polygon)
	c_direction = np.subtract(np.array(c_vertex),np.array((0,0,0)))
	mp_vector = translation_transformation(mp_vector,c_direction)
	mp_dir_vector = np.subtract(mp_vector,c_vertex)
	mp_dir_vector = mp_dir_vector/np.linalg.norm(mp_dir_vector)
	#############################################Hasta aqui todo bien##################################################
	#The last weird theory
	#Choose the edge that is associated with the face representign the biggest regular polyhedron
	fit_baricenter_2 = polyhedron_i_baricenters_2[fit_face_2]
	fit_comb = list(combinations(fit_face_2,2))
	fit_pair_dist = {}
	if polytope_type == 'Truncation':
		for pair in fit_comb:
			for fit_face in fit_face_reg_list:
				if pair[0] in fit_face and pair[1] in fit_face:
					p1 = np.array(transformed_vertices[pair[0]])
					p2 = np.array(transformed_vertices[pair[1]])
					fit_pair = pair
					fit_pair_verts = (p1,p2)
					break
	else:
		if len(fit_face_2) > len(max_polygon): 
			for pair in fit_comb:
				for fit_face in fit_face_list:
					if pair[0] in fit_face and pair[1] in fit_face: 
						p1 = np.array(transformed_vertices[pair[0]])
						p2 = np.array(transformed_vertices[pair[1]])
						fit_pair = pair
						fit_pair_verts = (p1,p2)
						break
		else:
			for pair in fit_comb:
				for fit_face in fit_face_reg_list:
					if pair[0] in fit_face and pair[1] in fit_face: 
						p1 = np.array(transformed_vertices[pair[0]])
						p2 = np.array(transformed_vertices[pair[1]])
						fit_pair = pair
						fit_pair_verts = (p1,p2)
						break

	fe1 = fit_pair_verts[0]
	fe2 = fit_pair_verts[1]
	mid_fe = np.array([0.5*(fe1[0] + fe2[0]),0.5*(fe1[1] + fe2[1]),0.5*(fe1[2] + fe2[2])])

	fe_vector = np.subtract(mid_fe,fit_baricenter_2)
	c_direction = np.subtract(np.array(c_vertex),np.array((0,0,0)))
	fe_vector = translation_transformation(fe_vector,c_direction)
	fe_dir_vector = np.subtract(fe_vector,c_vertex)
	fe_dir_vector = fe_dir_vector/np.linalg.norm(fe_dir_vector)
	Rr = rotation_matrix(fe_dir_vector,mp_dir_vector)
	Rf = Rr

	#Look for the biggest faces of the polytope
	poly_list_r = []
	for face in c_hull_i_area:
		if len(face) == len(max_polygon):
			poly_list_r.append(face)

	surrounding_faces = {}
	max_polygon_face = {}
	for edge in max_polygon:
		max_polygon_face[edge[0]] = None
		max_polygon_face[edge[1]] = None
		for face in poly_list_r:
			if edge[0] in face and edge[1] in face:
				surrounding_faces[face] = None
	max_polygon_face = list(max_polygon_face.keys())
	surrounding_faces = list(surrounding_faces.keys())
	#The surrounding faces are the faces with the same number of vertices connected to the fit_face
	#Eliminate the fit face from the list and leave only surrounding faces
	for face in surrounding_faces:
		count = 0
		for point in face:
			if point in max_polygon_face:
				count += 1
		if count == len(face):
			surrounding_faces.remove(face)

	if polytope_type == 'Cantellation':
		print('Interpretation selected: ',interpretation)
		#Select the small regular faces surroundign the max polygon
		if interpretation == 'VERTICES':
			#Find all the small regular faces
			poly_list_ct = [] 
			poly_list_not_ct =[] #List of all the faces that are not the min_reg_polygon of the polytope
			for face in c_hull_i_area:
				if len(face) == len(min_reg_polygon) and abs(c_hull_i_area[face] - areas[1]) < 0.1:
					poly_list_ct.append(face)
				else:
					poly_list_not_ct.append(face)

			#Select the small regular faces that are connected with the reference largest face
			ct_surrounding_faces = {}
			for face in poly_list_ct:
				for point in face:
					for pair in max_polygon:
						if point in pair:
							ct_surrounding_faces[face] = None
			
			#Select the that face are connected with the reference largest face and that are not the small regular faces
			ct_surrounding_other_faces = {}
			for face in poly_list_not_ct:
				for point in face:
					for pair in max_polygon:
						if point in pair:
							ct_surrounding_other_faces[face] = None

					
			ct_surrounding_faces = list(ct_surrounding_faces.keys())
			ct_surrounding_other_faces = list(ct_surrounding_other_faces.keys())
			small_face_list = areas_faces[small_face_area]
			#Find the baricenter of each surrounding small regular faces
			#Compute a vector normal to each of these faces
			ct_surround_bar = {}
			ct_surround_vec = {}
			for face in ct_surrounding_faces:
				face_verts = []
				for point in face:
					face_verts.append(reference_vertices[point])
				face_bar,usls = Baricenter(face_verts)
				f1 = face_verts[0]
				f2 = face_verts[1]
				f1_b = np.subtract(f1,face_bar)
				f2_b = np.subtract(f2,face_bar)
				vector_plane = np.cross(f1_b,f2_b)
				vector_plane = vector_plane/np.linalg.norm(vector_plane)
				ct_surround_bar[tuple(face_bar)] = face
				ct_surround_vec[face] = vector_plane
			
			ct_surround_vec = list(ct_surround_vec.values())#List containing the vectors perpendicular to the small regular faces conected to the largest face in the polytope
			#Find the baricenter of each surrounding NOT small regular faces
			#Compute a vector normal to each of these faces
			ct_surround_not_bar = {}
			ct_surround_not_vec = {}
			for face in ct_surrounding_other_faces:
				face_verts = []
				for point in face:
					face_verts.append(reference_vertices[point])
				face_bar,usls = Baricenter(face_verts)
				f1 = face_verts[0]
				f2 = face_verts[1]
				f1_b = np.subtract(f1,face_bar)
				f2_b = np.subtract(f2,face_bar)
				vector_plane = np.cross(f1_b,f2_b)
				vector_plane = vector_plane/np.linalg.norm(vector_plane)
				ct_surround_not_bar[tuple(face_bar)] = face
				ct_surround_not_vec[face] = vector_plane
			
			ct_surround_not_vec = list(ct_surround_not_vec.values())#List containing the vectors perpendicular to the NOT small regular faces conected to the largest face in the polytope
			fit_comb_ct = list(combinations(fit_face_ct,2))
			faces_wo_ct = [list(face) for face in faces_polyhedron]#Polyhedron
			faces_wo_ct.remove(list(fit_face_ct))

			#create the edges (pairs of vertices) of the polyhedron's reference face 
			poly_h_big_face_edges = {}
			for pair in fit_comb_ct:
				for face in faces_wo_ct:
					if pair[0] in face and pair[1] in face:
						poly_h_big_face_edges[pair] = None
						break

			poly_h_big_face_edges = list(poly_h_big_face_edges.keys())
			
			#Identify the indices of the polyhedron's reference face
			poly_h_big_face_inds = {}
			for pair in poly_h_big_face_edges:
				poly_h_big_face_inds[pair[0]]=None
				poly_h_big_face_inds[pair[1]]=None

			poly_h_big_face_inds = tuple(poly_h_big_face_inds.keys())

			for face in small_face_list:
				f_count = 0
				for point in poly_h_big_face_inds:
					if point in face:
						f_count +=1
				if f_count == len(poly_h_big_face_inds):
					small_face_list.remove(face)
					break

			#Select the edges that connect with the small faces
			cool_pairs = {}
			for pair in poly_h_big_face_edges:
				for face in small_face_list:
					if pair[0] in face and pair[1] in face:
						cool_pairs[pair] = None

			cool_pairs = list(cool_pairs.keys())
			cool_verts = {}
			cool_vectors = {}
			for pair in cool_pairs:
				pv1 = transformed_vertices[pair[0]]
				pv2 = transformed_vertices[pair[1]]
				vec = np.subtract(pv1,pv2)
				vec = vec/np.linalg.norm(vec)
				cool_vectors[pair] = vec
				cool_verts[tuple(pv1)] = None
				cool_verts[tuple(pv2)] = None

			bar_small_polyh_face,usls = Baricenter(list(cool_verts.keys()))
			#dx.scatter(bar_small_polyh_face[0],bar_small_polyh_face[1],bar_small_polyh_face[2],color='r')
			#dx.scatter(bar_polygon[0],bar_polygon[1],bar_polygon[2],color='k')

			target_vectors = []
			#for pair in cool_pairs:
			for cool_vertex in cool_verts:
				vector = np.subtract(np.array(cool_vertex),bar_small_polyh_face)
				#dx.quiver(mid_pm[0],mid_pm[1],mid_pm[2],vector[0],vector[1],vector[2])
				vector = vector/np.linalg.norm(vector)
				target_vectors.append(vector)

			cool_vectors_list = list(cool_vectors.values())
			
			if len(max_polygon)%2 == 0 and len(max_polygon) != 4:
				ct_surrounding_edges = {}	
				for pair in max_polygon:
					for face in poly_list_ct:
						if pair[0] in face and pair[1] in face:
							ct_surrounding_edges[pair] = None

				ct_surrounding_edges = list(ct_surrounding_edges.keys())
				ct = ct_surrounding_edges[0]
				ct_1 = reference_vertices[ct[0]]
				ct_2 = reference_vertices[ct[1]]
				mid_ct = np.array([0.5*(ct_1[0] + ct_2[0]),0.5*(ct_1[1] + ct_2[1]),0.5*(ct_1[2] + ct_2[2])])
				#dx.scatter(mid_ct[0],mid_ct[1],mid_ct[2],color = 'y')
				ct_dir_vector = np.subtract(mid_ct,bar_polygon)
				ct_dir_vector = ct_dir_vector/np.linalg.norm(ct_dir_vector)
				#dx.quiver(mid_ct[0],mid_ct[1],mid_ct[2],ct_dir_vector[0],ct_dir_vector[1],ct_dir_vector[2])
			else:
				ct = max_polygon[0]
				ct_1 = reference_vertices[ct[0]]
				ct_2 = reference_vertices[ct[1]]
				ct_dir_vector = np.subtract(ct_1,bar_polygon)
				ct_dir_vector = ct_dir_vector/np.linalg.norm(ct_dir_vector)

	######################################Todo bien todo correcto
			fit_count = 0
			for fit_vector in cool_vectors_list:
				for s_vector in ct_surround_not_vec:
					test = np.dot(fit_vector,s_vector)
					if  -0.1 <= test <= 0.1:
						fit_count += 1

			if fit_count != len(cool_vectors_list) or fit_count != 2*len(cool_vectors_list):
				for t_vector in target_vectors:
					Rt = rotation_matrix(t_vector,ct_dir_vector)
					t_aligned_vectors = []
					for c_vector in cool_vectors_list:
						align = np.dot(Rt,np.transpose(np.array(c_vector)))
						aligned_vector = []
						aligned_vector.append(align.item(0))
						aligned_vector.append(align.item(1))
						aligned_vector.append(align.item(2))
						t_aligned_vectors.append(aligned_vector)
					
					fit_count_2 = 0
					for fit_vector in t_aligned_vectors:
						for s_vector in ct_surround_not_vec:
							test = np.dot(fit_vector,s_vector)
							if  -0.1 <= test <= 0.1:
								fit_count_2 += 1

					#Multiply times 2 because this vectors can be perpendicular to the vector defining two parallel planes.
					if fit_count_2 >= 1:
						Rf = Rt
						break

		if interpretation == 'FACES':
			#Find all the small regular faces
			poly_list_ct = []
			for face in c_hull_i_area:
				if len(face) == len(min_reg_polygon) and abs(c_hull_i_area[face] - areas[1]) < 0.1:
					poly_list_ct.append(face)

			#Select the small regular faces that are connected with the reference largest face
			ct_surrounding_faces = {}
			for face in poly_list_ct:
				for point in face:
					for pair in max_polygon:
						if point in pair:
							ct_surrounding_faces[face] = None
					
			ct_surrounding_faces = list(ct_surrounding_faces.keys())#Polytope
			small_face_list = areas_faces[small_face_area]#Polyhedron

			#Find the baricenter of each surrounding small regular faces
			#Compute a vector normal to each of these faces
			ct_surround_bar = {}
			ct_surround_vec = {}
			for face in ct_surrounding_faces:
				face_verts = []
				for point in face:
					face_verts.append(reference_vertices[point])
				face_bar,usls = Baricenter(face_verts)
				f1 = face_verts[0]
				f2 = face_verts[1]
				f1_b = np.subtract(f1,face_bar)
				f2_b = np.subtract(f2,face_bar)
				vector_plane = np.cross(f1_b,f2_b)
				vector_plane = vector_plane/np.linalg.norm(vector_plane)
				ct_surround_bar[tuple(face_bar)] = face
				ct_surround_vec[face] = vector_plane

			ct_surround_vec = list(ct_surround_vec.values())#List containing the vectors perpendicular to the small regular faces conected to the largest face in the polytope
			fit_comb_ct = list(combinations(fit_face_ct,2))
			faces_wo_ct = [list(face) for face in faces_polyhedron]#Polyhedron
			faces_wo_ct.remove(list(fit_face_ct))

			#create the edges (oairs of vertices) of the polyhedron's reference face 
			poly_h_big_face_edges = {}
			for pair in fit_comb_ct:
				for face in faces_wo_ct:
					if pair[0] in face and pair[1] in face:
						poly_h_big_face_edges[pair] = None
						break

			poly_h_big_face_edges = list(poly_h_big_face_edges.keys())
			
			#Identify the indices of the polyhedron's reference face
			poly_h_big_face_inds = {}
			for pair in poly_h_big_face_edges:
				poly_h_big_face_inds[pair[0]]=None
				poly_h_big_face_inds[pair[1]]=None

			poly_h_big_face_inds = tuple(poly_h_big_face_inds.keys())

			for face in small_face_list:
				f_count = 0
				for point in poly_h_big_face_inds:
					if point in face:
						f_count +=1
				if f_count == len(poly_h_big_face_inds):
					small_face_list.remove(face)
					break

			#Select the edges that conect with the small faces
			cool_pairs = {}
			for pair in poly_h_big_face_edges:
				for face in small_face_list:
					if pair[0] in face and pair[1] in face:
						cool_pairs[pair] = None

			cool_pairs = list(cool_pairs.keys())
			#Create a set of direction vectors parallel to the edges of the large faee conecting with small faces in the polyhedron:

			cool_verts = {}
			cool_vectors = {}
			for pair in cool_pairs:
				pv1 = transformed_vertices[pair[0]]
				pv2 = transformed_vertices[pair[1]]
				vec = np.subtract(pv1,pv2)
				vec = vec/np.linalg.norm(vec)
				cool_vectors[pair] = vec
				cool_verts[tuple(pv1)] = None
				cool_verts[tuple(pv2)] = None

			bar_small_polyh_face,usls = Baricenter(list(cool_verts.keys()))
			#dx.scatter(bar_small_polyh_face[0],bar_small_polyh_face[1],bar_small_polyh_face[2],color='r')
			#dx.scatter(bar_polygon[0],bar_polygon[1],bar_polygon[2],color='k')

			target_vectors = []
			for pair in cool_pairs:
				pm1 = np.array(transformed_vertices[pair[0]])
				pm2 = np.array(transformed_vertices[pair[1]])
				mid_pm = np.array((0.5*(pm1[0] + pm2[0]),0.5*(pm1[1] + pm2[1]),0.5*(pm1[2] + pm2[2])))
				#dx.scatter(mid_pm[0],mid_pm[1],mid_pm[2],color = 'm')
				vector = np.subtract(mid_pm,bar_small_polyh_face)
				#dx.quiver(mid_pm[0],mid_pm[1],mid_pm[2],vector[0],vector[1],vector[2])
				vector = vector/np.linalg.norm(vector)
				target_vectors.append(vector)

			cool_vectors_list = list(cool_vectors.values())
			
			for pair in max_polygon:
				pmp1 = reference_vertices[pair[0]]
				pmp2 = reference_vertices[pair[1]]
				#dx.scatter(pmp1[0],pmp1[1],pmp1[2],color = 'k')
				#dx.scatter(pmp2[0],pmp2[1],pmp2[2],color = 'k')

			
			if len(max_polygon)%2 == 0 and len(max_polygon) != 4:
				ct_surrounding_edges = {}	
				for pair in max_polygon:
					for face in poly_list_ct:
						if pair[0] in face and pair[1] in face:
							ct_surrounding_edges[pair] = None

				ct_surrounding_edges = list(ct_surrounding_edges.keys())
				ct = ct_surrounding_edges[0]
				ct_1 = reference_vertices[ct[0]]
				ct_2 = reference_vertices[ct[1]]
				mid_ct = np.array([0.5*(ct_1[0] + ct_2[0]),0.5*(ct_1[1] + ct_2[1]),0.5*(ct_1[2] + ct_2[2])])
				#dx.scatter(mid_ct[0],mid_ct[1],mid_ct[2],color = 'y')
				ct_dir_vector = np.subtract(mid_ct,bar_polygon)
				ct_dir_vector = ct_dir_vector/np.linalg.norm(ct_dir_vector)
				#dx.quiver(mid_ct[0],mid_ct[1],mid_ct[2],ct_dir_vector[0],ct_dir_vector[1],ct_dir_vector[2])
			else:
				ct = max_polygon[0]
				ct_1 = reference_vertices[ct[0]]
				ct_2 = reference_vertices[ct[1]]
				ct_dir_vector = np.subtract(ct_1,bar_polygon)
				ct_dir_vector = ct_dir_vector/np.linalg.norm(ct_dir_vector)			
				#dx.quiver(ct_1[0],ct_1[1],ct_1[2],ct_dir_vector[0],ct_dir_vector[1],ct_dir_vector[2])
	######################################Todo bien todo correcto
			fit_count = 0
			for fit_vector in cool_vectors_list:
				for s_vector in ct_surround_vec:
					test = np.dot(fit_vector,s_vector)
					if  -0.1 <= test <= 0.1:
						fit_count += 1

			if fit_count != len(cool_vectors_list) or fit_count != 2*len(cool_vectors_list):
				for t_vector in target_vectors:
					Rt = rotation_matrix(t_vector,ct_dir_vector)
					t_aligned_vectors = []
					for c_vector in cool_vectors_list:
						align = np.dot(Rt,np.transpose(np.array(c_vector)))
						aligned_vector = []
						aligned_vector.append(align.item(0))
						aligned_vector.append(align.item(1))
						aligned_vector.append(align.item(2))
						t_aligned_vectors.append(aligned_vector)
					
					fit_count_2 = 0
					for fit_vector in t_aligned_vectors:
						for s_vector in ct_surround_vec:
							test = np.dot(fit_vector,s_vector)
							if  -0.1 <= test <= 0.1:
								fit_count_2 += 1

					#Multiply times 2 because this vectors can be perpendicular to the vector defining two parallel planes.
					if fit_count_2 >= 1:
						Rf = Rt
						break

	if polytope_type == 'Regular_a' and len(surrounding_faces) != 0: 
	#Define a vector perpendicular for each of the surounding faces and their baricenters
		surround_bar = {}
		surround_vec = {}
		for face in surrounding_faces:
			face_verts = []
			for point in face:
				face_verts.append(reference_vertices[point])
			face_bar,usls = Baricenter(face_verts)
			f1 = face_verts[0]
			f2 = face_verts[1]
			f1_b = np.subtract(f1,face_bar)
			f2_b = np.subtract(f2,face_bar)
			vector_plane = np.cross(f1_b,f2_b)
			vector_plane = vector_plane/np.linalg.norm(vector_plane)
			surround_bar[tuple(face_bar)] = face
			surround_vec[face] = vector_plane

		#for s_vector in list(surround_vec.values()):
		#	cx.quiver(c_vertex[0],c_vertex[1],c_vertex[2],s_vector[0],s_vector[1],s_vector[2],length = 10.0)
		
		#Look for the faces in fit face list that are conected to the fit face 2 (Fit polyhedron)
		fit_list_r = [face for face in fit_face_list]
		fit_list_r.remove(fit_face_2)
		fit_rounding_faces = []
		cool_pairs = []
		for pair in fit_comb:
			for fit_face in fit_list_r:
				if pair[0] in fit_face and pair[1] in fit_face: 
					fit_rounding_faces.append(fit_face)
					cool_pairs.append(pair)


		#Create a set of direction vectors:
		target_vectors = []
		for pair in max_polygon:
			pm1 = np.array(reference_vertices[pair[0]])
			pm2 = np.array(reference_vertices[pair[1]])
			mid_pm = np.array((0.5*(pm1[0] + pm2[0]),0.5*(pm1[1] + pm2[1]),0.5*(pm1[2] + pm2[2])))
			#cx.scatter(mid_pm[0],mid_pm[1],mid_pm[2],color = 'm')
			vector = np.subtract(mid_pm,bar_max_polygon)
			vector = vector/np.linalg.norm(vector)
			target_vectors.append(vector)

		cool_vectors = {}
		for pair in cool_pairs:
			pv1 = transformed_vertices[pair[0]]
			pv2 = transformed_vertices[pair[1]]
			vec = np.subtract(pv1,pv2)
			vec = vec/np.linalg.norm(vec)
			cool_vectors[pair] = vec

		cool_vectors_list = list(cool_vectors.values())
		fit_count = 0
		for fit_vector in cool_vectors_list:
			for s_vector in list(surround_vec.values()):
				test = np.dot(fit_vector,s_vector)
				if  -0.1 <= test <= 0.1:
					fit_count += 1

		if fit_count != len(cool_vectors_list):
			for t_vector in target_vectors:
				Rt = rotation_matrix(t_vector,mp_dir_vector)
				t_aligned_vectors = []
				for c_vector in cool_vectors_list:
					align = np.dot(Rt,np.transpose(np.array(c_vector)))
					aligned_vector = []
					aligned_vector.append(align.item(0))
					aligned_vector.append(align.item(1))
					aligned_vector.append(align.item(2))
					t_aligned_vectors.append(aligned_vector)
					fit_count_2 = 0
				for fit_vector in t_aligned_vectors:
					for s_vector in list(surround_vec.values()):
						test = np.dot(fit_vector,s_vector)
						if  -0.1 <= test <= 0.1:
							fit_count_2 += 1

				#Multiply time 2 because this vectors can be perpendicular to the vector defining two parallel planes
				if fit_count_2 == len(cool_vectors_list)*2:
					Rf = Rt
					break

	#Look for the smallest faces of the polytope
	small_reg_faces_list = []
	for face in c_hull_i_area:
		if len(face) == len(min_reg_polygon):
			small_reg_faces_list.append(face)
#########################################################################################################################################################################################################
	min_surrounding_faces = {}
	min_polygon_face = {}
	for edge in min_reg_polygon:
		for face in c_hull_i_area:
			if edge[0] in face and edge[1] in face:
				min_surrounding_faces[face] = None
				min_polygon_face[edge[0]] = None
				min_polygon_face[edge[1]] = None
	min_polygon_face = list(min_polygon_face.keys())
	min_surrounding_faces = list(min_surrounding_faces.keys())

	#The surrounding faces in this second case are the faces with the same number of vertices of the smallest polygons connected to the fit_face
	#Eliminate the fit face from th list. only surrounding faces
	for face in min_surrounding_faces:
		count = 0
		for point in face:
			if point in min_polygon_face:
				count += 1
		if count == len(face):
			min_surrounding_faces.remove(face)
	
	if polytope_type == 'Regular_b' and len(min_surrounding_faces) != 0: 
	#Define a vector perpendicular for each of the surounding faces and their baricenters################################################################################
		min_surround_bar = {}
		min_surround_vec = {}
		for face in min_surrounding_faces:
			face_verts = []
			for point in face:
				face_verts.append(reference_vertices[point])
			face_bar,usls = Baricenter(face_verts)
			f1 = face_verts[0]
			f2 = face_verts[1]
			f1_b = np.subtract(f1,face_bar)
			f2_b = np.subtract(f2,face_bar)
			vector_plane = np.cross(f1_b,f2_b)
			vector_plane = vector_plane/np.linalg.norm(vector_plane)
			min_surround_bar[tuple(face_bar)] = face
			min_surround_vec[face] = vector_plane

		#for s_vector in list(surround_vec.values()):
		#	cx.quiver(c_vertex[0],c_vertex[1],c_vertex[2],s_vector[0],s_vector[1],s_vector[2],length = 10.0)
		fit_comb_reg_b = list(combinations(fit_face_reg_b,2))
		faces_wo_b = [face for face in faces_polyhedron]
		faces_wo_b.remove(list(fit_face_reg_b))

		poly_h_small_face_edges = {}
		for pair in fit_comb_reg_b:
			for face in faces_wo_b:
				if pair[0] in face and pair[1] in face:
					poly_h_small_face_edges[pair] = None
					break

		cool_pairs = list(poly_h_small_face_edges.keys())
		#Create a set of direction vectors:

		cool_verts = {}
		cool_vectors = {}
		for pair in cool_pairs:
			pv1 = transformed_vertices[pair[0]]
			pv2 = transformed_vertices[pair[1]]
			vec = np.subtract(pv1,pv2)
			vec = vec/np.linalg.norm(vec)
			cool_vectors[pair] = vec
			cool_verts[tuple(pv1)] = None
			cool_verts[tuple(pv2)] = None

		bar_small_polyh_face,usls = Baricenter(list(cool_verts.keys()))
		#cx.scatter(bar_small_polyh_face[0],bar_small_polyh_face[1],bar_small_polyh_face[2],color='r')
		#cx.scatter(bar_polygon[0],bar_polygon[1],bar_polygon[2],color='k')

		target_vectors = []
		for pair in cool_pairs:
			pm1 = np.array(transformed_vertices[pair[0]])
			pm2 = np.array(transformed_vertices[pair[1]])
			mid_pm = np.array((0.5*(pm1[0] + pm2[0]),0.5*(pm1[1] + pm2[1]),0.5*(pm1[2] + pm2[2])))
			#cx.scatter(mid_pm[0],mid_pm[1],mid_pm[2],color = 'm')
			vector = np.subtract(mid_pm,bar_small_polyh_face)
			#cx.quiver(mid_pm[0],mid_pm[1],mid_pm[2],vector[0],vector[1],vector[2])
			vector = vector/np.linalg.norm(vector)
			target_vectors.append(vector)

		cool_vectors_list = list(cool_vectors.values())
		reg_b = min_reg_polygon[0]
		r_b_1 = reference_vertices[reg_b[0]]
		r_b_2 = reference_vertices[reg_b[1]]
		mid_reg_b = np.array([0.5*(r_b_1[0] + r_b_2[0]),0.5*(r_b_1[1] + r_b_2[1]),0.5*(r_b_1[2] + r_b_2[2])])
		#cx.scatter(mid_reg_b[0],mid_reg_b[1],mid_reg_b[2],color = 'g')
		reg_b_dir_vector = np.subtract(mid_reg_b,bar_polygon)
		reg_b_dir_vector = reg_b_dir_vector/np.linalg.norm(reg_b_dir_vector)
		#cx.quiver(mid_reg_b[0],mid_reg_b[1],mid_reg_b[2],reg_b_dir_vector[0],reg_b_dir_vector[1],reg_b_dir_vector[2])
		fit_count = 0
		for fit_vector in cool_vectors_list:
			for s_vector in list(min_surround_vec.values()):
				test = np.dot(fit_vector,s_vector)
				if  -0.1 <= test <= 0.1:
					fit_count += 1

		if fit_count != len(cool_vectors_list) or fit_count != 2*len(cool_vectors_list):
			for t_vector in target_vectors:
				Rt = rotation_matrix(t_vector,reg_b_dir_vector)
				t_aligned_vectors = []
				for c_vector in cool_vectors_list:
					align = np.dot(Rt,np.transpose(np.array(c_vector)))
					aligned_vector = []
					aligned_vector.append(align.item(0))
					aligned_vector.append(align.item(1))
					aligned_vector.append(align.item(2))
					t_aligned_vectors.append(aligned_vector)
				
				fit_count_2 = 0
				for fit_vector in t_aligned_vectors:
					for s_vector in list(min_surround_vec.values()):
						test = np.dot(fit_vector,s_vector)
						if  -0.1 <= test <= 0.1:
							fit_count_2 += 1

				#Multiply times 2 because this vectors can be perpendicular to the vector defining two parallel planes.
				if fit_count_2 >= 1:
					Rf = Rt
					break

	##################################################################SECOND ROTATION###################################################################
	#translation to the origin of coordinates
	f_origin_verts = []
	f_o_direction = np.subtract(np.array((0,0,0)),np.array(c_vertex))
	for vertex in transformed_vertices:
		f_o_vert = translation_transformation(vertex,f_o_direction)
		f_origin_verts.append(f_o_vert)

	#Rotation to align vectors
	f_aligned_vertices = []
	for vertex in f_origin_verts:
		align = np.dot(Rf,np.transpose(vertex))
		aligned_vertex = []
		aligned_vertex.append(align.item(0))
		aligned_vertex.append(align.item(1))
		aligned_vertex.append(align.item(2))
		f_aligned_vertices.append(aligned_vertex)

	#Back translation
	f_transformed_vertices = []
	for vertex in f_aligned_vertices:
		bt_vertex = translation_transformation(vertex,f_o_direction*(-1))
		f_transformed_vertices.append(bt_vertex)

	#Final translation to center both reference and objective polyhedrons in the same baricenter

	ref_bar,usls = Baricenter(reference_vertices)
	obj_bar,usls = Baricenter(f_transformed_vertices)
	ft_dir_vector = np.subtract(ref_bar,obj_bar)

	ft_transformed_vertices = []
	for vertex in f_transformed_vertices:
		ft_vertex = translation_transformation(vertex,ft_dir_vector)
		ft_transformed_vertices.append(ft_vertex)


	faces_list=[]
	for f in range(num_faces_po):
		polygon = []
		for j in faces_polyhedron[f]:
			#polygon.append(transformed_vertices[j])#After first rotation
			polygon.append(ft_transformed_vertices[j])#After second rotation
			#polygon.append(polyhedron_vertices[j])#Without alignment
		polygon = np.array(polygon)
		faces_list.append(polygon)

	#solidp = Poly3DCollection(faces_list)
	#solidp.set_edgecolor('k')
	#solidp.set_facecolor('r')
	#solidp.set_alpha(0.15)
	##dx.add_collection3d(solidp)
	#plt.show()
	return(ft_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction)

if __name__ == "__Polyhedron_Alignment__":
    Polyhedron_Alignment()

def Polygonalization(convex_hull,vertices,plane_tolerance):
	triangles = convex_hull.Export_Triangles()
	Xp = []
	Yp = []
	Zp = []
	for vertex in vertices:
		Xp.append(vertex[0])
		Yp.append(vertex[1])
		Zp.append(vertex[2])

	#ax.plot_trisurf(Xp, Yp, triangles, Zp,edgecolors='k',linewidth=0.5,cmap = 'autumn',alpha=0)
	
	triangles2 = convex_hull.Export_Triangles()
	ts = {}
	for triangle in triangles:
		ts[tuple(triangle)] = None
	triangles = list(ts.keys())
	in_plane_indices_dict = {}
	c_hull_v_area = {}
	c_hull_i_area = {}
	q = 0
	hull_polygons = {}
	hull_poly_ratio = {}
	hull_polygons_area = {}
	for triangle in triangles:
		tris = [tri for tri in triangles]
		tris.remove(triangle)
		tetrahedron = []
		in_plane_triangles = []
		#Get all the vertices that are coplanar with those of the given triangle
		in_plane_vertices = []
		in_plane_indices = []
		ver_l = [vertex for vertex in vertices]
		in_plane_triangles.append(triangle)
		for i in triangle:
			tetrahedron.append(vertices[i])
			in_plane_vertices.append(vertices[i])
			in_plane_indices.append(i)
			ver_l.remove(vertices[i])
		#Create tetrahedrons with the remaining vertices in the hull	
		for vertex in ver_l:
			tetrahedron.append(vertex)
			tet = np.asarray(tetrahedron)
			tet = np.insert(tet,3,1,axis = 1)
			matrix_tet = np.matrix(tet)
			tetrahedron_volume = (np.linalg.det(matrix_tet))/6
			if abs(tetrahedron_volume) <= float(plane_tolerance):
				in_plane_vertices.append(vertex)
				in_plane_indices.append(vertices.index(vertex))
			
			tetrahedron.remove(vertex)
		in_plane_indices_dict[tuple(in_plane_indices)] = in_plane_vertices

		#The plane for a triangle has been created
		#Find all the triangles that belongs to the same plane

		for trgl in tris:
			if trgl[0] in in_plane_indices and trgl[1] in in_plane_indices and trgl[2] in in_plane_indices:
				in_plane_triangles.append(trgl)
				triangles.remove(trgl)
		#Convex hull of the -----------> POLYGON <------------- 

		conv_hull_verts = {}
		conv_hull_inds = {}
		edges = {}
		area_of_plane = 0
		for triangle in in_plane_triangles:
			#Compute the area of the plane 
			p1 = np.array(vertices[triangle[0]])
			p2 = np.array(vertices[triangle[1]])
			p3 = np.array(vertices[triangle[2]])
			p1p2 = np.subtract(p1,p2)
			p1p3 = np.subtract(p1,p3)
			cross_p = np.cross(p1p2,p1p3)
			area_of_triangle = np.linalg.norm(cross_p)/2.0
			area_of_plane += area_of_triangle
				
			#Generate the edges of the triangle 
			tri_comb = list(combinations(triangle,2))
			for edge in tri_comb:
				tri_edge = 0
				for triangle in in_plane_triangles:
					if edge[0] in triangle and edge[1] in triangle:
						tri_edge += 1
				if tri_edge == 1:
					edges[edge] = None
					ae1 = vertices[edge[0]]
					ae2 = vertices[edge[1]]
					xe = [ae1[0],ae2[0]]
					ye = [ae1[1],ae2[1]]
					ze = [ae1[2],ae2[2]]
					#ax.plot(xe,ye,ze,linewidth = 1,color = 'r')
			#An edge belongs to the convex hull iff it is part of only one triangle

		#Eliminate vertices internal to the edges
		for edge in list(edges.keys()):
			conv_hull_inds[edge[0]] = None
			conv_hull_inds[edge[1]] = None
			conv_hull_verts[tuple(vertices[edge[0]])] = None
			conv_hull_verts[tuple(vertices[edge[1]])] = None

		conv_hull_verts = list(conv_hull_verts.keys())
		conv_hull_inds = list(conv_hull_inds.keys())
		plane_ind_combs = list(combinations(in_plane_indices,3))
		to_eliminate = {}
		for triplet in plane_ind_combs:
			t1 = np.array(vertices[triplet[0]])
			t2 = np.array(vertices[triplet[1]])
			t3 = np.array(vertices[triplet[2]])
			ta = np.subtract(t1,t2)
			tb = np.subtract(t2,t3)
			tc = np.subtract(t3,t1)
			t_cross = np.cross(ta,tb)
			area_t = np.linalg.norm(t_cross)/2.0
			if abs(area_t) <= 0.25:
				t_dict = {}
				t_dict[np.linalg.norm(ta)] = (triplet[0],triplet[1])
				t_dict[np.linalg.norm(tb)] = (triplet[1],triplet[2])
				t_dict[np.linalg.norm(tc)] = (triplet[2],triplet[0])
				max_n = max(list(t_dict.keys()))
				good_t = t_dict[max_n]
				for i in triplet:
					if i not in good_t:
						to_eliminate[i] = None
		to_eliminate = list(to_eliminate.keys())
		actual_inds = []
		actual_verts = []

		for ind in conv_hull_inds:
			if ind not in to_eliminate:
				actual_inds.append(ind)
				actual_verts.append(tuple(vertices[ind]))

		actual_inds = tuple(actual_inds)
		actual_verts = tuple(actual_verts)
		c_hull_v_area[actual_verts] = area_of_plane########################################### conv_hull_verts
		c_hull_i_area[actual_inds] = area_of_plane
		
		edges = list(edges.keys())
	
		actual_edges_v = {}
		#The last weird theory
		#Align the plane perpendicular to z-axis
		al1 = np.array(actual_verts[0])
		al2 = np.array(actual_verts[1])
		al3 = np.array(actual_verts[2])
		aa = np.subtract(al1,al2)
		ab = np.subtract(al1,al3)
		normal_v = np.cross(aa,ab)
		normal_v = normal_v/np.linalg.norm(normal_v)
		R = rotation_matrix(normal_v,np.array((0,0,1)))
		aligned_plane_verts = []
		for vertex in actual_verts:
			align = np.dot(R,np.transpose(vertex))
			aligned_vertex = []
			aligned_vertex.append(align.item(0))
			aligned_vertex.append(align.item(1))
			aligned_plane_verts.append(aligned_vertex)

		triangulation_polygon = Delaunay(np.array(aligned_plane_verts))################################
		triangles_polygon = list(triangulation_polygon.simplices)
		for triangle in triangles_polygon:
			tri_comb = list(combinations(triangle,2))
			for edge in tri_comb:
				tri_edge = 0
				for triangle in triangles_polygon:
					if edge[0] in triangle and edge[1] in triangle:
						tri_edge += 1
				if tri_edge == 1:
					actual_edge_vert1 = tuple(actual_verts[edge[0]])
					actual_edge_vert2 = tuple(actual_verts[edge[1]])
					actual_edge_verts = (actual_edge_vert1,actual_edge_vert2)
					actual_edges_v[actual_edge_verts] = None
		# actual edges_v ----> pairs of verts of the polytope
		actual_edges_i = {}
		for edge in actual_edges_v: 
			i1 = vertices.index(edge[0])
			i2 = vertices.index(edge[1])
			ult_edge = (i1,i2)
			actual_edges_i[ult_edge] = None		
	
		edges_lenght = {}
		for edge in actual_edges_v:
			e1 = np.array(edge[0])
			e2 = np.array(edge[1])
			d_e = np.linalg.norm(e1 - e2)
			edges_lenght[edge] = d_e		

		max_edge = max(list(edges_lenght.values()))
		min_edge = min(list(edges_lenght.values()))

		proportion = min_edge/max_edge

		actual_edges_i = tuple(list(actual_edges_i.keys()))

		hull_polygons_area[actual_edges_i] = area_of_plane

		if proportion >= 0.9:
			hull_polygons[actual_edges_i] = 'regular'

		if 0.45 <= proportion < 0.9:
			hull_polygons[actual_edges_i] = 'semiregular'

		if proportion < 0.45:
			hull_polygons[actual_edges_i] = 'irregular'

		hull_poly_ratio[actual_edges_i] = proportion
		q += 1

	#Number of vertices in pore polytope
	planes_inds = list(c_hull_i_area.keys())
	no_rep_inds = {}
	for plane in planes_inds:
		for ind in plane:
			no_rep_inds[ind] = None
	num_vertices_polytope = len(list(no_rep_inds.keys()))


	hull_polygons_num_edges = {}
	for polygon in hull_polygons:
		points = len(polygon)
		hull_polygons_num_edges[polygon] = points

	#Organize the polygons according to the number of vertices
	hull_types_of_polygons = {}
	for points in list(hull_polygons_num_edges.values()):
		hull_types_of_polygons[points] = None
	
	hull_types_of_polygons = list(hull_types_of_polygons.keys())

	types_dict = {}
	for polygon_type in hull_types_of_polygons:
		list_t_polygons = []
		for polygon in hull_polygons_num_edges:
			if hull_polygons_num_edges[polygon] == polygon_type:
				list_t_polygons.append(polygon)
		types_dict[polygon_type] = list_t_polygons

	#Organize the polygons according to the areas
	hull_types_of_areas = {}
	for area in list(hull_polygons_area.values()):
		hull_types_of_areas[area] = None
	hull_types_of_areas = list(hull_types_of_areas.keys())

	area_dict = {}
	for area in hull_types_of_areas:
		list_a_polygons = []
		for polygon in hull_polygons_area:
			if hull_polygons_area[polygon] == area:
				list_a_polygons.append(polygon)
		area_dict[area] = list_a_polygons
	
	reg_min_areas = {}
	for polygon in hull_polygons:
		if hull_polygons[polygon] != 'irregular':
			area = hull_polygons_area[polygon]
			reg_min_areas[area] = polygon

	min_reg_area = min(list(reg_min_areas.keys()))
	min_reg_polygon = reg_min_areas[min_reg_area]

	max_poly_area = max(hull_types_of_areas)
	max_polygon = area_dict[max_poly_area]
	max_polygon = max_polygon[0]
	re_count_a = 0
	re_count_b = 0
	for polygon in hull_polygons:
		if hull_polygons[polygon] == 'regular':
			re_count_a += 1
		#Include semiregular if the area is representative
		if hull_polygons[polygon] == 'semiregular':
			a = hull_polygons_area[polygon]
			if a >= max_poly_area*(3/4):
				re_count_b += 1
				re_count_a += 1
	if re_count_a == len(hull_polygons):
		if re_count_b == 0:
			polytope_type = 'Regular_a'
		else:
			polytope_type = 'Regular_b'
			
		sugestions = hull_types_of_polygons
		num_faces = len(hull_polygons)
		areas = [max_poly_area,max_poly_area]
		print('Regular')
		return(c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio)

	else:
		#Compute the biggest irregular polyhedron
		irre_areas = {}
		for polygon in hull_polygons:
			if hull_polygons[polygon] == 'irregular':
				area = hull_polygons_area[polygon]
				irre_areas[area] = polygon

		max_irre_area = max(list(irre_areas.keys()))
		max_irre_polygon = irre_areas[max_irre_area]

		#Compute the biggest regular polyhedron 
		reg_areas = {}
		for polygon in hull_polygons:
			if hull_polygons[polygon] == 'regular':
				area = hull_polygons_area[polygon]
				reg_areas[area] = polygon

		if len(reg_areas) != 0:		
			max_reg_area = max(list(reg_areas.keys()))
			max_reg_polygon = reg_areas[max_reg_area]
			#Compute the  smallest regular polyhedron
			reg_min_areas = {}
			for polygon in hull_polygons:
				if hull_polygons[polygon] != 'irregular':
					area = hull_polygons_area[polygon]
					reg_min_areas[area] = polygon

			min_reg_area = min(list(reg_min_areas.keys()))
			min_reg_polygon = reg_min_areas[min_reg_area]

		else:
			semireg_areas = {}
			for polygon in hull_polygons:
				if hull_polygons[polygon] == 'semiregular':
					area = hull_polygons_area[polygon]
					semireg_areas[area] = polygon
			#Compute the  smallest semiregular polyhedron
			semireg_min_areas = {}
			for polygon in hull_polygons:
				if hull_polygons[polygon] != 'irregular':
					area = hull_polygons_area[polygon]
					semireg_min_areas[area] = polygon

			min_semireg_area = min(list(reg_min_areas.keys()))
			max_semireg_area = max(list(reg_min_areas.keys()))
			min_semireg_polygon = reg_min_areas[min_semireg_area]
			max_semireg_polygon = reg_min_areas[max_semireg_area]
		
		cantellation = 'off'
		truncation = 'off'
		semiregular = 'off'
		
		if len(reg_areas) != 0:
			if abs(min_reg_area - max_reg_area) < 0.01 and len(min_reg_polygon) == len(max_reg_polygon):
				truncation = 'on'
			else:
				cantellation = 'on'

		else:
			semiregular = 'on'

		if cantellation == 'on':
			sugestions = []
			num_faces = 0
			num_faces_2 = 0

			if len(max_irre_polygon) > 4:
				sugestions.append(len(max_reg_polygon))
				sugestions.append(len(max_irre_polygon))
			elif len(min_reg_polygon) > 4:
				sugestions.append(len(max_reg_polygon)*2)
				sugestions.append(len(min_reg_polygon)/2)

			for polygon in hull_polygons:
				#Account for the small regular polygons
				if len(polygon) == len(min_reg_polygon) and abs(hull_polygons_area[polygon] - min_reg_area) < 0.01 and len(min_reg_polygon) > 4:
					num_faces_2 += 1
				#Account for the biggest irregular faces
				if len(polygon) == len(max_irre_polygon) and abs(hull_polygons_area[polygon] - max_irre_area) < 0.01  and  max_irre_area > max_reg_area:
					num_faces += 1
					num_faces_2 +=1
				#Account for the biggest regular faces
				if len(polygon) == len(max_reg_polygon) and abs(hull_polygons_area[polygon] - max_reg_area) < 0.01:
					num_faces += 1
					num_faces_2 +=1	

			num_faces = [num_faces,num_faces_2]
			areas = [max_poly_area,min_reg_area]
			#i.e if the biggest irregular doesn't contribute
			polytope_type = 'Cantellation'
			print('Cantellation')
			return(c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio)

		if truncation == 'on':
			sugestions = [len(max_reg_polygon),len(max_irre_polygon)]
			num_faces = 0
			for polygon in hull_polygons:
				if hull_polygons[polygon] == 'regular':
					num_faces += 1
				if hull_polygons[polygon] == 'irregular':
					if len(polygon) == len(max_irre_polygon):
						num_faces += 1

			if int(len(max_irre_polygon)/2)not in sugestions:
				sugestions.append(int(len(max_irre_polygon)/2))
	
			polytope_type = 'Truncation'
			areas = [max_poly_area,max_reg_area]
			print('Truncation')
			return(c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio)

		if semiregular == 'on':
			polytope_type = 'Semiregular'
			num_faces = len(hull_polygons)
			sugestions = [len(max_irre_polygon),len(max_semireg_polygon),len(min_semireg_polygon)]
			areas = [max_poly_area,max_semireg_area]
			print('Semiregular')
			return(c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,max_semireg_polygon,polytope_type,hull_polygons,hull_poly_ratio)

if __name__ == "__Polygonalization__":
    Polygonalization()

def General_Polygonalization(convex_hull,vertices):
	X = []
	Y = []
	Z = []
	for vert in vertices:
		X.append(vert[0])
		Y.append(vert[1])
		Z.append(vert[2])
	triangles = convex_hull.Export_Triangles()
	triangles2 = convex_hull.Export_Triangles()
	ts = {}
	for triangle in triangles:
		ts[tuple(triangle)] = None

	triangles = list(ts.keys())
	in_plane_indices_dict = {}
	c_hull_v_area = {}
	c_hull_i_area = {}
	hull_polygons = {}
	hull_polygons_area = {}

	q = 0
	for triangle in triangles:
		tris = [tri for tri in triangles]
		tris.remove(triangle)
		tetrahedron = []
		in_plane_triangles = []
		#Get all the vertices that are coplanar with those of the given triangle
		in_plane_vertices = []
		in_plane_indices = []
		ver_l = [vertex for vertex in vertices]
		in_plane_triangles.append(triangle)
		for i in triangle:
			tetrahedron.append(vertices[i])
			in_plane_vertices.append(vertices[i])
			in_plane_indices.append(i)
			ver_l.remove(vertices[i])
			
		for vertex in ver_l:
			tetrahedron.append(vertex)
			tet = np.asarray(tetrahedron)
			tet = np.insert(tet,3,1,axis = 1)
			matrix_tet = np.matrix(tet)
			tetrahedron_volume = (np.linalg.det(matrix_tet))/6
			if abs(tetrahedron_volume) <= 0.001:
				in_plane_vertices.append(vertex)
				in_plane_indices.append(vertices.index(vertex))
			
			tetrahedron.remove(vertex)
			
		in_plane_indices_dict[tuple(in_plane_indices)] = in_plane_vertices


		#The plane for a triangle has been created
		#Find all the triangles that belongs to the same plane

		for trgl in tris:
			if trgl[0] in in_plane_indices and trgl[1] in in_plane_indices and trgl[2] in in_plane_indices:
				in_plane_triangles.append(trgl)
				triangles.remove(trgl)
		#Convex hull of the -----------> POLYGON <------------- 
		conv_hull_verts = {}
		conv_hull_inds = {}
		edges = {}
		area_of_plane = 0
		for triangle in in_plane_triangles:
			#Compute the area of the plane 
			p1 = np.array(vertices[triangle[0]])
			p2 = np.array(vertices[triangle[1]])
			p3 = np.array(vertices[triangle[2]])
			p1p2 = np.subtract(p1,p2)
			p1p3 = np.subtract(p1,p3)
			cross_p = np.cross(p1p2,p1p3)
			area_of_triangle = np.linalg.norm(cross_p)/2.0
			area_of_plane += area_of_triangle
				
			#Generate the edges of the triangle 
			tri_comb = list(combinations(triangle,2))
			for edge in tri_comb:
				tri_edge = 0
				for triangle in in_plane_triangles:
					if edge[0] in triangle and edge[1] in triangle:
						tri_edge += 1
				if tri_edge == 1:
					edges[edge] = None
					ae1 = vertices[edge[0]]
					ae2 = vertices[edge[1]]
					xe = [ae1[0],ae2[0]]
					ye = [ae1[1],ae2[1]]
					ze = [ae1[2],ae2[2]]

		
		#An edge belongs to the convex hull iff it belong only to one triangle

		for edge in list(edges.keys()):
			conv_hull_inds[edge[0]] = None
			conv_hull_inds[edge[1]] = None
			conv_hull_verts[tuple(vertices[edge[0]])] = None
			conv_hull_verts[tuple(vertices[edge[1]])] = None

		conv_hull_verts = list(conv_hull_verts.keys())
		conv_hull_inds = list(conv_hull_inds.keys())
		plane_ind_combs = list(combinations(in_plane_indices,3))
		to_eliminate = {}
		for triplet in plane_ind_combs:
			t1 = np.array(vertices[triplet[0]])
			t2 = np.array(vertices[triplet[1]])
			t3 = np.array(vertices[triplet[2]])
			ta = np.subtract(t1,t2)
			tb = np.subtract(t2,t3)
			tc = np.subtract(t3,t1)
			t_cross = np.cross(ta,tb)
			area_t = np.linalg.norm(t_cross)/2.0
			if abs(area_t) <= 0.015:
				t_dict = {}
				t_dict[np.linalg.norm(ta)] = (triplet[0],triplet[1])
				t_dict[np.linalg.norm(tb)] = (triplet[1],triplet[2])
				t_dict[np.linalg.norm(tc)] = (triplet[2],triplet[0])
				max_n = max(list(t_dict.keys()))
				good_t = t_dict[max_n]
				for i in triplet:
					if i not in good_t:
						to_eliminate[i] = None
		to_eliminate = list(to_eliminate.keys())
		actual_inds = []
		actual_verts = []

		for ind in conv_hull_inds:
			if ind not in to_eliminate:
				actual_inds.append(ind)
				actual_verts.append(tuple(vertices[ind]))

		actual_inds = tuple(actual_inds)
		actual_verts = tuple(actual_verts)
		c_hull_v_area[actual_verts] = area_of_plane
		c_hull_i_area[actual_inds] = area_of_plane
		edges = list(edges.keys())
	
		actual_edges_v = {}
		#The last weird theory
		#Align the plane perpendicular to z-axis
		al1 = np.array(actual_verts[0])
		al2 = np.array(actual_verts[1])
		al3 = np.array(actual_verts[2])
		aa = np.subtract(al1,al2)
		ab = np.subtract(al1,al3)
		normal_v = np.cross(aa,ab)
		normal_v = normal_v/np.linalg.norm(normal_v)
		R = rotation_matrix(normal_v,np.array((0,0,1)))
		aligned_plane_verts = []
		for vertex in actual_verts:
			align = np.dot(R,np.transpose(vertex))
			aligned_vertex = []
			aligned_vertex.append(align.item(0))
			aligned_vertex.append(align.item(1))
			aligned_plane_verts.append(aligned_vertex)

		triangulation_polygon = Delaunay(np.array(aligned_plane_verts))################################
		triangles_polygon = list(triangulation_polygon.simplices)
		for triangle in triangles_polygon:
			tri_comb = list(combinations(triangle,2))
			for edge in tri_comb:
				tri_edge = 0
				for triangle in triangles_polygon:
					if edge[0] in triangle and edge[1] in triangle:
						tri_edge += 1
				if tri_edge == 1:
					actual_edge_vert1 = tuple(actual_verts[edge[0]])
					actual_edge_vert2 = tuple(actual_verts[edge[1]])
					actual_edge_verts = (actual_edge_vert1,actual_edge_vert2)
					actual_edges_v[actual_edge_verts] = None
		# actual edges_v ----> pairs of verts of the polytope
		actual_edges_i = {}
		for edge in actual_edges_v: 
			i1 = vertices.index(edge[0])
			i2 = vertices.index(edge[1])
			ult_edge = (i1,i2)
			actual_edges_i[ult_edge] = None		
	
		edges_lenght = {}
		for edge in actual_edges_v:
			e1 = np.array(edge[0])
			e2 = np.array(edge[1])
			d_e = np.linalg.norm(e1 - e2)
			edges_lenght[edge] = d_e		

		max_edge = max(list(edges_lenght.values()))
		min_edge = min(list(edges_lenght.values()))

		proportion = min_edge/max_edge

		actual_edges_i = tuple(list(actual_edges_i.keys()))

		hull_polygons_area[actual_edges_i] = area_of_plane
		hull_polygons[actual_edges_i] = None
		q += 1

	polygonalization_verts = []
	for polygon in hull_polygons:
		polygon_verts = []
		for edge in polygon:
			p1 = vertices[edge[0]]
			p2 = vertices[edge[1]]
			pair = (p1,p2)
			polygon_verts.append(pair)
		polygonalization_verts.append(polygon_verts)

	polygonalization_inds = list(hull_polygons.keys())

	hull_polygons_num_edges = {}
	for polygon in hull_polygons:
		points = len(polygon)
		hull_polygons_num_edges[polygon] = points

	#Organize the polygons according to the number of vertices
	hull_types_of_polygons = {}
	for points in list(hull_polygons_num_edges.values()):
		hull_types_of_polygons[points] = None
	
	hull_types_of_polygons = list(hull_types_of_polygons.keys())

	types_dict = {}
	for polygon_type in hull_types_of_polygons:
		list_t_polygons = []
		for polygon in hull_polygons_num_edges:
			if hull_polygons_num_edges[polygon] == polygon_type:
				list_t_polygons.append(polygon)
		types_dict[polygon_type] = list_t_polygons

	#Organize the polygons according to the areas
	hull_types_of_areas = {}
	for area in list(hull_polygons_area.values()):
		hull_types_of_areas[area] = None
	hull_types_of_areas = list(hull_types_of_areas.keys())

	area_dict = {}
	for area in hull_types_of_areas:
		list_a_polygons = []
		for polygon in hull_polygons_area:
			if hull_polygons_area[polygon] == area:
				list_a_polygons.append(polygon)
		area_dict[area] = list_a_polygons

	return(polygonalization_verts,polygonalization_inds,types_dict,area_dict)
if __name__ == "__General_Polygonalization__":
	General_Polygonalization()

def Polyhedron_Dimensional_Parameters(triangles,vertices,o_vertex):
	#To Calculate the volume of a polytope whose convex hull
	#is that obtained from chull.py. Simplices are passed as 
	#'triangles', among a lists of coordinates of each vertex of the convex hull 
	#and the center op the polytope
	poly_vol = []
	Polyhedron_Total_Volume = 0
	for j in range(len(triangles)):
		tetrahedron = [o_vertex]
		for i in triangles[j]:
			tetrahedron.append(vertices[i])
		tet = np.asarray(tetrahedron)
		tet = np.insert(tet,3,1,axis = 1)
		matrix_tet = np.matrix(tet)
		tetrahedron_volume = (np.linalg.det(matrix_tet))/6
		poly_vol.append(tetrahedron_volume)

	for k in poly_vol:
		Polyhedron_Total_Volume += abs(k)


	Polyhedron_Total_Area = 0
	for triangle in triangles:
		p1 = np.array(vertices[triangle[0]])
		p2 = np.array(vertices[triangle[1]])
		p3 = np.array(vertices[triangle[2]])
		p1p2 = np.subtract(p1,p2)
		p1p3 = np.subtract(p1,p3)
		cross_p = np.cross(p1p2,p1p3)
		area_of_triangle = np.linalg.norm(cross_p)/2.0
		Polyhedron_Total_Area += area_of_triangle
	return(Polyhedron_Total_Volume,Polyhedron_Total_Area)
if __name__ == "__Polyhedron_Dimensional_Parameters__":
    Polyhedron_Dimensional_Parameters()

def Ghost_Protocol(parameters_file):
    sub.Popen(("cp",parameters_file,"parachutes.txt"))
    ghost_in_the_shell = sub.Popen(("gedit","parachutes.txt"))
    ghost_in_the_shell.wait()
    return ("parachutes.txt")
if __name__ == "__ Ghost_Protocol__":
	Ghost_Protocol() 
