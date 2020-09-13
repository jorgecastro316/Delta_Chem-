import re
from scipy.spatial import Voronoi
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import subprocess as sub
from scipy.signal import find_peaks
from scipy import optimize
import time 
from scipy.optimize import curve_fit
from pylab import *
from prettytable import PrettyTable
from itertools import combinations
import statistics 
import sys
import numpy as np
from math import *
from random import random  
from deltachem.Convex_Hull_I import Vector,Hull,Face,Edge,Vertex  
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import sympy
from deltachem.irreducible_volume import Irreducible_Volume
#from celluloid import Camera
from itertools import combinations
from prettytable import PrettyTable
from deltachem.DC_library import *
from deltachem.DC_library import PSD_Solve
from deltachem.DC_library import Polygonalization
from deltachem.DC_polyhedrons import *
import os

class Delta_Chem():
	def __init__(self,basename,dat,plots,parameters_file):
		self.basename = basename
		self.dat = dat
		self.plots = plots
		self.parameters_file = parameters_file
		self.Structure = None
		self.Atomic_Species = None
		self.Atomic_Species_Grouped = None
		self.Convex_Hulls = {}
		self.Polytopes = []
		self.Centers = []
		self.Pore_radii = []
		self.Volumes = []
		self.Areas = []
		self.Contributing_Atoms = []
		self.Fittest_Polyhedrons = []
		self.Irreducible_Volumes = []
		self.Generator_Triangles = []
		self.Porosity_Volumes = []
		self.Porosity_Areas = []
		self.Porosity_num_faces = []


	def Leganto(self,structure_file):
		colors = ['r','k','y','m']
		size = [15.2,17.00,12.00,13.90]
		#Open the structure file .xyz
		strct = open(structure_file)
		strct = strct.readlines()
		strct = strct[2:]
		elements = []
		coords = []
		coordinates = []
		structure = {}
		#Separate the chemical symbols and the coordinates of each atom in
		#the structure
		for line in strct:
			symbol = re.match('([a-zA-Z]*)(.*)',line,re.M|re.I)
			if symbol:
				elements.append(symbol.group(1))
				coords.append(symbol.group(2))
		#Separate each component of the coordinates to create a position vector
		#for each atom
		for coord in coords:
			point = re.compile(r'\s+')
			point = point.split(coord)
			point.pop(0)
			px = float(point[0])
			py = float(point[1])
			pz = float(point[2])
			point = (px,py,pz)
			coordinates.append(point)

		#Create a lisit of the elements present in the structure
		list_of_elements = {}
		for i in range(len(elements)):
			list_of_elements[elements[i]] = coordinates[i]
			structure[coordinates[i]] = elements[i]

		elements_set = list_of_elements
		list_of_elements = list(list_of_elements.keys())
		self.Atomic_Species = list_of_elements
		#Gruop the atoms in the structure clasified by their symbol
		grouped_atoms = {}
		for atom in list_of_elements:
			group = []
			for i in structure:
				if structure[i] == atom:
					group.append(i)
			grouped_atoms[atom] = group
		self.Atomic_Species_Grouped = grouped_atoms
		atoms_list = list(grouped_atoms.values())
		lst = []
		for elem in atoms_list:
			X=[]
			Y=[]
			Z=[]
			for single_atom in elem:
				X.append(single_atom[0])
				Y.append(single_atom[1])
				Z.append(single_atom[2])
			lst_elem =[X,Y,Z]
			lst.append(lst_elem)

		#Open the file with de VdW radii data
		#and create a dictionary whose keys are the elements and 
		#the values are the vdw radii 
		radii_list = open('radii_list.dat')
		VdW_radii = {}
		for line in radii_list:
			elmt = re.match('([a-zA-Z]*)',line,re.M|re.I)
			VdW_radius = re.search('([0-9].[0-9]*)',line,re.M|re.I)
			if elmt and VdW_radius:
				VdW_radii[elmt.group(1)] = float(VdW_radius.group(1))

		#Create a list enumerating the atoms in the structure
		indices = []
		for i in range(len(coordinates)):
			indices.append(i)

		#Compute all the possible combinations between pairs of atoms
		pairs_combination = list(combinations(indices,2))
		#Calculate the distances between the atoms of each pair, the create a dict whose keys are
		#the distances and the values are the pairs
		distances = {}
		for pair in pairs_combination:
			elem_pair = []
			#P1 and P2 are the coordinates of each atom corresponding with the indices pair
			P1 = coordinates[pair[0]]
			P2 = coordinates[pair[1]]
			#elem_pair contains the element that corresponds with each atom of the pair
			elem_pair.append(structure[P1])
			elem_pair.append(structure[P2])
			#Obtains de VdW radius of each atom in the pair
			vdwrad_1 = VdW_radii[elem_pair[0]]
			vdwrad_2 = VdW_radii[elem_pair[1]]
			#Sum both radii and subtract this quantity from the distance between the points
			#rad_correct = vdwrad_1 + vdwrad_2
			rad_correct = 0
			p1 = np.array(P1)
			p2 = np.array(P2)
			distance1_2 = np.linalg.norm(p1-p2)
			distance1_2 = distance1_2 - rad_correct
			distances[pair] = distance1_2 
			elem_pair = []
			self.Structure = structure
		return(coordinates,distances,lst,structure,VdW_radii,grouped_atoms)

	if __name__ == "__Leganto__":
	    Leganto()

	def Kerno(self,basename,dat,plots,parameters_file):
		path_tst = os.path.dirname(__file__)
		xyz = basename + '.xyz'
		parameters = open("./" + parameters_file,'r')
		parameters = parameters.readlines()
		for line in parameters:
			match_time = re.match(r'PSD_solve_time = (.*)',line,re.M|re.I)
			if match_time:
				execution_time = match_time.group(1)
				break
		for line in parameters:
			match_resolution = re.match(r'PSD_grid_res = (.*)',line,re.M|re.I)
			if match_resolution:
				grid_resolution = match_resolution.group(1)
				break
		for line in parameters:
			match_size = re.match(r'PSD_test_particle = (.*)',line,re.M|re.I)
			if match_size:
				test_particle_size = match_size.group(1)
				break
		for line in parameters:
			match_height = re.match(r'peak_finder_height = (.*)',line,re.M|re.I)
			if match_height:
				height = match_size.group(1)
				break				
		for line in parameters:
			match_prominence = re.match(r'peak_finder_prominence = (.*)',line,re.M|re.I)
			if match_prominence:
				prominence = match_prominence.group(1)
				break
		for line in parameters:
			match_overlap = re.match(r'overlap_sphr_extra = (.*)',line,re.M|re.I)
			if match_overlap:
				overlap_sphr_extra = match_overlap.group(1)
				break	

		PSD_Solve(xyz,execution_time,grid_resolution,test_particle_size)
		pores_diameters_dsv = PSD_analysis(dat,height,prominence)
		coordinates,distances,lst,structure,VdW_radii,grouped_atoms = self.Leganto(xyz)
		elements_list = list(grouped_atoms.keys())  
		pores_diameters = list(pores_diameters_dsv.keys())

		self.Pore_radii = [radius for radius in pores_diameters]
		if os.path.isfile(basename + '_kerno_prev_bar.txt') == True:
			print('Kerno files found: Avoiding centers computation')
			kerno_bars_file = basename + '_kerno_prev_bar.txt'
			kerno_prev_bars = open(kerno_bars_file,'r')
			kerno_prev_bars = kerno_prev_bars.readlines()
			kerno_bar_coordinates = []
			
			for line in kerno_prev_bars:
				kerno_bar_coords = re.match('(.*)',line,re.M|re.I)
				kerno_bar_coords = kerno_bar_coords.group()
				kerno_bar_coordinates.append(kerno_bar_coords)
		
		#Separate each component of the coordinates to create a position vector
			centers_of_porosities = []
			for coord in kerno_bar_coordinates:
				centers_per_pore = []
				point = re.compile(r'\s+')
				point = point.split(coord)
				point.pop(0)
				px = float(point[0])
				py = float(point[1])
				pz = float(point[2])
				point = (px,py,pz)
				centers_per_pore.append(point)
				centers_of_porosities.append(centers_per_pore)
		else:
			xc = []
			yc = []
			zc = []
			for coord in coordinates:
				xc.append(coord[0])
				yc.append(coord[1])
				zc.append(coord[2])

			maxx = max(xc)
			maxy = max(yc)
			maxz = max(zc)

			#################################################
			vertices = coordinates
			position_vectors = [] 
			for vertex in vertices:
				v0 = np.around(vertex[0],3)
				v1 = np.around(vertex[1],3)
				v2 = np.around(vertex[2],3)
				position_vectors.append(Vector(v0,v1,v2))
			convex_hull = Hull(position_vectors)
			triangles = (convex_hull.Export_Triangles())
			X = []
			Y = []
			Z = []
			for j in vertices:
				X.append(j[0])
				Y.append(j[1])
				Z.append(j[2])


			voro = Voronoi(coordinates)
			#################################################
			
			print('Original candidates',len(voro.vertices))

			# filter 1: The size of the system. Reject the Voronoi vertices that are too far from the system
			F1_voro_verts = []
			for vertex in voro.vertices:
				if abs(vertex[0]) < maxx and abs(vertex[1]) < maxy and abs(vertex[2]) < maxz:
					v1 = np.around(float(vertex[0]),5)
					v2 = np.around(float(vertex[1]),5)
					v3 = np.around(float(vertex[2]),5)
					F1_voro_verts.append((v1,v2,v3))
			print('Candidates after F1: ',len(F1_voro_verts))
			################################################################################################
			# filter 2: Select only the candidates inside the convex hull
			F2_voro_verts = []
			maxc = []
			epsilon = 1e-5
			for F1_Candidate in F1_voro_verts:
				candidate_count = 0
				tetrahedron = []
				for triangle in triangles:
					tetrahedron.append(coordinates[triangle[0]])
					tetrahedron.append(coordinates[triangle[1]])
					tetrahedron.append(coordinates[triangle[2]])	
					tetrahedron.append(F1_Candidate)
					tet = np.asarray(tetrahedron)
					tet = np.insert(tet,3,1,axis = 1)
					matrix_tet = np.matrix(tet)
					tetrahedron_volume = (np.linalg.det(matrix_tet))/6
					if tetrahedron_volume < epsilon:
						candidate_count += 1
					#print(tetrahedron_volume)
					tetrahedron = []
					maxc.append(candidate_count)

				if candidate_count == len(triangles):
					F2_voro_verts.append(F1_Candidate)

			if len(F2_voro_verts) < len(F1_voro_verts)/4:
				F2_voro_verts = F1_voro_verts

			print('Candidates after F2: ',len(F2_voro_verts))

			########################################################################################################################################		
			#For each candidate compute the closest atom distance

			F2_dist_dict = {}
			F2_vdw_dict = {}
			dist_wo_vdw = {}
			for F2_Candidate in F2_voro_verts:
				distances = {}
				distances_wo_vdw = {}
				for atom_position in coordinates:
					dist_cand_atom = np.linalg.norm(np.array(atom_position) - np.array(F2_Candidate))
					vdw_corr = VdW_radii[structure[atom_position]]
					d_wo_vdw =dist_cand_atom
					dist_cand_atom = abs(dist_cand_atom - vdw_corr)
					distances[dist_cand_atom] = vdw_corr
					distances_wo_vdw[d_wo_vdw] = None
				min_dist = min(list(distances.keys()))
				min_dist_wo = min(list(distances_wo_vdw.keys()))
				F2_dist_dict[F2_Candidate] = min_dist
				F2_vdw_dict[F2_Candidate] = distances[min_dist]
				dist_wo_vdw[F2_Candidate] = min_dist_wo

			#########################################################################################################################################
			#Filter 3: 
			F3_all_candidates = {}
			#Ascendant order
			pores_diameters.sort()
			for F2_Candidate in F2_dist_dict:
				low_bound_a = pores_diameters[0]/2.0 - F2_vdw_dict[F2_Candidate]/2
				up_bound_a = pores_diameters[-1]/2.0 + F2_vdw_dict[F2_Candidate]
				if low_bound_a <= F2_dist_dict[F2_Candidate] <= up_bound_a:
					F3_all_candidates[F2_Candidate] = None 
						
			print('Candidates after F3: ',len(F3_all_candidates))
			centers_of_porosities = []

			if len(pores_diameters) == 1:
				centers_per_pore = []
				candidates = list(F3_all_candidates.keys())
				distances = {}
				for candidate in candidates:
					cad = dist_wo_vdw[candidate]
					distances[candidate] = cad
				LES = max(list(distances.values()))
				for candidate in candidates:
					if distances[candidate] == LES:
						centers_per_pore.append(candidate)
						break
				print(print('Candidates after F4: ', len(centers_per_pore)))
				centers_of_porosities.append(centers_per_pore)

			#If there are more than one type of porosity
			else:
				#Prior to the next filter, classify the candidates accordin to the porosit they probably belong to
				#Descendant order
				pores_diameters.sort(reverse = True)
				F3_all_candidates = list(F3_all_candidates.keys())
				F3_classified = {}
				for F3_candidate in F3_all_candidates:
					cad = F2_dist_dict[F3_candidate]
					diffs = {}
					for pore_diameter in pores_diameters:
						diff = abs(cad - pore_diameter/2)
						diffs[pore_diameter] = diff
					diffs_l = list(diffs.values())
					min_diff = min(diffs_l)
					for pore_diameter in pores_diameters:
						if diffs[pore_diameter] == min_diff:
							F3_classified[F3_candidate] = pore_diameter
							break

				#Order the candidates according to the assigned pore diameter
				F3_per_pore = []
				F3_all_candidates = {}
				#Ascendant order
				pores_diameters.sort()
				for pore_diameter in pores_diameters:
					candidates_per_pore = []
					for candidate in F3_classified:
						if F3_classified[candidate] == pore_diameter:
							candidates_per_pore.append(candidate)
							F3_all_candidates[candidate] = pore_diameter
					F3_per_pore.append(candidates_per_pore)

			###################################################################################################################################################
				# Filter 4: Organize the candidates by porosity and eliminate those that don't belong to a determined cavity
				F4_all_candidates_list = [candidate for candidate in reversed(list(F3_all_candidates.keys()))]
				F4_all_candidates_pds = [pds for pds in reversed(list(F3_all_candidates.values()))]
				l = 0
				for F3_candidate in F4_all_candidates_list:
					F4_all_others_list = [candidate for candidate in F4_all_candidates_list]
					F4_all_others_list.remove(F3_candidate)
					coo = 0
					for F3_other_candidate in F4_all_others_list:
						if F3_classified[F3_other_candidate] != F3_classified[F3_candidate]:
							dist_bar_otb = np.linalg.norm(np.array(F3_other_candidate) - np.array(F3_candidate))
							if dist_bar_otb <= dist_wo_vdw[F3_candidate]:
								coo += 1
								rmv_indx = F4_all_candidates_list.index(F3_other_candidate)
								F4_all_candidates_list.remove(F3_other_candidate)
								F4_all_candidates_pds.remove(F4_all_candidates_pds[rmv_indx])

				F4_candidates = []
				F4_all_candidates = dict(zip(F4_all_candidates_list,F4_all_candidates_pds))
				for pore_diameter in pores_diameters:
					F4_candidates_per_pore = []
					for F4_candidate in F4_all_candidates:
						if abs(F4_all_candidates[F4_candidate] - pore_diameter)<0.1:
							F4_candidates_per_pore.append(F4_candidate)
					F4_candidates.append(F4_candidates_per_pore)

				print('Candidates after F4: ', len(F4_all_candidates))
				#########################################################################################################################################################
				# Filter 5: Grouping
				F5_candidates = []
				j = 0
				lenf5 = 0
				for F4_candidates_per_pore in F4_candidates:
					pore_diameter = pores_diameters[j]
					pore_radius = pore_diameter/2.0
					F5_candidates_per_pore = grouping(pore_diameter/2,F4_candidates_per_pore)
					F5_candidates_per_pore = {tuple(candidate): None for candidate in F5_candidates_per_pore}
					F5_candidates_per_pore = list(F5_candidates_per_pore.keys())
					lenf5 += len(F5_candidates_per_pore)
					F5_candidates.append(F5_candidates_per_pore)
					j += 1

				print('Candidates after F5: ', lenf5)

				###################################################################################################################################################################
				#Filter 6: Largest Overlaping Spherical Shells. 
				j = 0
				lenf6 = 0
				for F5_candidates_per_pore in F5_candidates:
					pore_diameter = pores_diameters[j]
					pore_radius = pore_diameter/2.0
					can_dict_1 = {}
					can_dict_2 = {}
					can_dict_3 = {}
					can_dict_4 = {}
					can_dict_5 = {}
					temp_list = []
					centers_per_pore = []
					atoms_per_pore = []
				#Obtain the closest atom distances
					for F5_candidate in F5_candidates_per_pore:
						distances = {}
						overlaping_atoms_number = 0
						overlaping_atoms = []
						for atom_position in coordinates:
							dist_atom_F5 = np.linalg.norm(np.array(F5_candidate) - np.array(atom_position))
							vdw_corr = VdW_radii[structure[atom_position]]
							dist_atom_F5 = dist_atom_F5 - vdw_corr
							distances[atom_position] = dist_atom_F5

						closest_atom_distance = min(list(distances.values()))
						diff = abs(closest_atom_distance - (pore_diameter/2))
						if diff < 1.0: 
							for atom_position in coordinates:
								if  closest_atom_distance <= distances[atom_position] <= closest_atom_distance + 3.6:
									overlaping_atoms_number += 1	
									overlaping_atoms.append(atom_position)

						can_dict_1[F5_candidate] = overlaping_atoms_number
						can_dict_2[F5_candidate] = closest_atom_distance
						can_dict_3[F5_candidate] = diff

					max_surround_count = max(list(can_dict_1.values()))
					for F5_candidate in can_dict_1:
						if can_dict_1[F5_candidate] == max_surround_count:
							temp_list.append(F5_candidate)

					for center in temp_list:
						can_dict_4[can_dict_3[center]] = center

					min_diff = min(list(can_dict_4.keys()))
					for diff in can_dict_4:
						if diff == min_diff:
							centers_per_pore.append(can_dict_4[diff])
							break
					centers_of_porosities.append(centers_per_pore)
					lenf6 += len(centers_per_pore)
					j += 1

				print('Candidates after F6: ', lenf6)
			with open(basename + '_kerno_prev_bar.txt', 'w') as fp:
				for centers_per_pore in centers_of_porosities:
					for center in centers_per_pore:
						fp.write('{0:12.6f} {1:12.6f} {2:12.6f}\n'.format(center[0],center[1],center[2]))
		#########################################################################################################################################
		#Select the contributing atoms
		#Select the contributing atoms
		contributing_atoms = []
		pc = 0
		for centers_per_pore in centers_of_porosities:
			for center in centers_per_pore:
				self.Centers.append(center)
				atoms_per_pore = []
				distances = {}
				for atom_position in coordinates:
					distance = np.linalg.norm(np.array(atom_position) - np.array(center))
					vdw_corr = VdW_radii[structure[atom_position]]
					distace = distance - vdw_corr
					distances[atom_position] = distance
				cad = min(list(distances.values()))
				for atom_position in coordinates:
					if distances[atom_position] <= cad + float(overlap_sphr_extra):
						atoms_per_pore.append(atom_position)
				self.Contributing_Atoms.append(atoms_per_pore)
				contributing_atoms.append(atoms_per_pore)
		
		if plots == 'y':
			fig_b= plt.figure()
			bx = fig_b.gca(projection='3d')
			colors_P = ['g','m','b','y','k','r']
			bb = 0
			cb = 1
			Xc = []
			Yc = []
			Zc = []
			for j in coordinates:
				Xc.append(j[0])
				Yc.append(j[1])
				Zc.append(j[2])

			bx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in Structure')
			for center in self.Centers:
				bx.scatter(center[0],center[1],center[2],color = colors_P[bb-cb],alpha = 0.7,label = 'Center porosity '+ str(cb))
				cb += 1
			#bx.plot_trisurf(X, Y, triangles, Z,edgecolors='k',linewidth=1.5,cmap = 'winter',alpha=0.25)
			for atoms_per_pore in self.Contributing_Atoms:
				xa = []
				ya = []
				za = []
				for atom in atoms_per_pore:
					xa.append(atom[0])
					ya.append(atom[1])
					za.append(atom[2])
				bx.scatter(xa,ya,za,color = colors_P[bb],alpha = 0.7,label ='Atoms contributing to porosity ' + str(bb+1))
				bb += 1
			bx.legend()
			plt.title('Kerno')
			plt.show()

		return(coordinates,contributing_atoms,centers_of_porosities,lst,pores_diameters,structure,VdW_radii,elements_list)

	if __name__ == "__Kerno__":
	    Kerno()


	def Mozaiko(self,vertices,parameters_file):
		parameters = open(parameters_file,'r')
		parameters = parameters.readlines()
		for line in parameters:
			match_plane_tol = re.match(r'Tol_plane_tet = (.*)',line,re.M|re.I)
			if match_plane_tol:
				plane_tol = match_plane_tol.group(1)
				break

		position_vectors = [] 

		for vertex in vertices:
			position_vectors.append(Vector(vertex[0],vertex[1],vertex[2]))

		convex_hull = Hull(position_vectors)
		triangles = (convex_hull.Export_Triangles())

		X = []
		Y = []
		Z = []
		for j in vertices:
			X.append(j[0])
			Y.append(j[1])
			Z.append(j[2])

		#Get the faces
		self.Convex_Hulls[convex_hull] = [X,Y,Z]
		#Ignore the position_vectors that haven't been considered on the convex hull

		ver_left = []
		for i in convex_hull.Export_Vertices():
			ver_left.append(vertices[i])

		#######################################Analyze the surface of the polytope##########################################
		c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = Polygonalization(convex_hull,vertices,plane_tol)
		self.Porosity_num_faces.append(num_faces)
		Polytope = [c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio]
		self.Polytopes.append(Polytope)

		b_vertex,b_radius = Baricenter(ver_left)
		polytope_volume,polytope_area = Polyhedron_Dimensional_Parameters(triangles,ver_left,b_vertex)
		self.Porosity_Volumes.append(polytope_volume)
		self.Porosity_Areas.append(polytope_area)

		return(X,Y,Z,convex_hull,Polytope,ver_left,polytope_volume,polytope_area)
	if __name__ == "__Mozaiko__":
	    Mozaiko()

	def Fragmento(self,atoms_of_porosities,centers_of_porosities,atoms_in_structure,pores_diameters,structure,Vdw_radii,elements_list,parameters_file,plots):
		ls = [i for i in range(len(centers_of_porosities))]
		comb_ls = list(combinations(ls,2))
		coordinates = list(structure.keys())
		irreducible_volumes_points = []
		c = 0
		for atoms_per_pore in atoms_of_porosities:
			X,Y,Z,convex_hull,Polytope,ver_left,polytope_volume,polytope_area = self.Mozaiko(atoms_per_pore,parameters_file)
			triangles = (convex_hull.Export_Triangles())
			c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = Polytope
			if plots == 'y':
				fig_cg = plt.figure(1)
				#############################################################################################################################################
				cgx = fig_cg.gca(projection='3d')
				xcg = []
				ycg = []
				zcg = []
				for atom in atoms_per_pore:
					xcg.append(atom[0])
					ycg.append(atom[1])
					zcg.append(atom[2])
					#cgx.scatter(xcg,ycg,zcg,color='r',alpha = 1,linewidths = 2)
				cgx.plot_trisurf(X, Y, triangles, Z,edgecolors='k',linewidth=1.5,cmap = 'winter',alpha=0.25)
				plt.axis('off')
				plt.title('Mozaiko: Delaunay Triangulation stage' + '\n' + self.basename + ' Porosity # ' + str(c+1))
	            ##############################################################################################################################################
				fig_d = plt.figure(2)
				dx = fig_d.gca(projection='3d')
				d = 0
				for polygon in hull_polygons:
					xe = []
					ye = []
					ze = []
					for edge in polygon:
						p1 = atoms_per_pore[edge[0]]
						p2 = atoms_per_pore[edge[1]]
						xe = [p1[0],p2[0]]
						ye = [p1[1],p2[1]]
						ze = [p1[2],p2[2]]
						dx.plot(xe,ye,ze,color = 'k',linewidth = 1.5) 
				d += 1
				plt.axis('off')
				plt.title('Mozaiko: Polygonalization stage' + '\n' + self.basename + ' Porosity # ' + str(c+1))
				plt.show()
	            
			########################################Fitting of the pore cavity########################################################

			print("////////////////////Computing minumum enclosing sphere # " + str(c+1) + "/////////////////////////")
			b_vertex,b_radius = Baricenter(ver_left)
			print("////////////////////Calculating volume of Polytope # " + str(c+1) + "/////////////////////////////")
			print("/////////////////////////////////Fitting with regular polyhedrons/////////////////////////////////")
			polyhedrons = {"Tetrahedron":Tetrahedron,
						"Cube":Cube,
						"Octahedron":Octahedron,
						"Dodecahedron":	Dodecahedron,
				    	"Icosahedron": Icosahedron,
				    	"Truncated_Tetrahedron": Truncated_Tetrahedron,
				    	"Truncated_Cube": Truncated_Cube,
				    	"Truncated_Octahedron": Truncated_Octahedron,
				    	"Truncated_Dodecahedron": Truncated_Dodecahedron,
				    	"Truncated_Icosahedron": Truncated_Icosahedron,
				    	"Cuboctahedron": Cuboctahedron,
				    	"Icosidodecahedron": Icosidodecahedron,
				    	"Small_Rhombicuboctahedron": Small_Rhombicuboctahedron,
				    	"Small_Rhombicosidodecahedron": Small_Rhombicosidodecahedron,
				    	"Great_Rhombicuboctahedron" :Great_Rhombicuboctahedron,
				    	"Great_Rhombicosidodecahedron": Great_Rhombicosidodecahedron,
				    	}

			volume_difference_p = {}
			volume_difference_b = {}
			area_difference_p = {}
			area_difference_b = {}
			polyhedrons_filters_123 = {}
			polyhedrons_area_diff = {}
			polyhedrons_vol_diff = {}
			c1_filter = {}
			c2_filter = {}
			r_filter = {}
			t_filter = {}
			s_filter = {}
			########################################Polyhedron fitting####################################
			for polyhedron in list(polyhedrons.keys()):
				#p_faces,p_vertices,p_num_faces,p_Volume_of_polyhedron,p_Area_of_polyhedron,p_max_area_face = polyhedrons[polyhedron](p_radius,p_vertex)
				b_faces,b_vertices,b_num_faces,b_Volume_of_polyhedron,b_Area_of_polyhedron,b_max_area_face,b_min_area_face,b_areas_faces,b_faces_areas = polyhedrons[polyhedron](b_radius,b_vertex)
				#vol_diff_p = abs (polytope_volume - p_Volume_of_polyhedron)
				vol_diff_b = abs (polytope_volume - b_Volume_of_polyhedron)
				#area_diff_p = abs(polytope_area - p_Area_of_polyhedron)
				area_diff_b = abs(polytope_area - b_Area_of_polyhedron)
				polyhedrons_area_diff[polyhedron] = area_diff_b 
				polyhedrons_vol_diff[polyhedron] = vol_diff_b
				type_faces_polyhedron = {}
				for face in b_faces:
					type_faces_polyhedron[len(face)] = None
				type_faces_polyhedron = list(type_faces_polyhedron.keys())


				if polytope_type == 'Cantellation':
					num_small_reg_faces = abs(num_faces[0] - num_faces[1])
					#Filter 1: Number of faces
					if b_num_faces == num_faces[0]:
						c1_filter[polyhedron] = abs(b_max_area_face - areas[0])

					if b_num_faces == num_faces[1]:
						s_polyt_num_faces = abs(num_faces[0] - num_faces[1])
						s_polyh_num_faces = len(b_areas_faces[b_min_area_face])
						if s_polyt_num_faces == s_polyh_num_faces:
							c2_filter[polyhedron] = abs(b_min_area_face - areas[1])
							
				if polytope_type == 'Regular_a' or polytope_type == 'Regular_b':
					
					#Filter 2: Number of faces
					if b_num_faces == num_faces:
						face_count = 0
						if len(type_faces_polyhedron) == len(sugestions):
							for face in type_faces_polyhedron:
								if face in sugestions:
									face_count += 1

						if face_count == len(sugestions):
							r_filter[polyhedron] = abs(b_max_area_face - areas[0])


				if polytope_type == 'Truncation':

					if b_num_faces == num_faces:				
						if sugestions[0] in type_faces_polyhedron:
							sug_count = 0
							for face in type_faces_polyhedron:
								if face in sugestions:
									sug_count += 1
							if sug_count > 1:
								t_filter[polyhedron] = abs(b_max_area_face - areas[0])

				if polytope_type == 'Semiregular':
					if b_num_faces == num_faces:				
						if sugestions[0] in type_faces_polyhedron:
							sug_count = 0
							for face in type_faces_polyhedron:
								if face in sugestions:
									sug_count += 1
							if sug_count > 1:
								s_filter[polyhedron] = abs(b_max_area_face - areas[0])

			if len(c1_filter) != 0:
				min_C1 = min(list(c1_filter.values()))
				for poly in c1_filter:
					if c1_filter[poly] == min_C1:
						polyhedrons_filters_123[poly] = min_C1

			if len(c2_filter) != 0:
				min_C2 = min(list(c2_filter.values()))
				for poly in c2_filter:
					if c2_filter[poly] == min_C2:
						polyhedrons_filters_123[poly] = min_C2

			if len(r_filter) != 0:
				min_r = min(list(r_filter.values()))
				for poly in r_filter:
					if r_filter[poly] == min_r:
						polyhedrons_filters_123[poly] = min_r

			if len(t_filter) != 0:
				for poly in t_filter:
					polyhedrons_filters_123[poly] = t_filter[poly]

			if len(s_filter) != 0:
				min_s = min(list(s_filter.values()))
				for poly in s_filter:
					if s_filter[poly] == min_s:
						polyhedrons_filters_123[poly] = min_s
			
			list_differences = list(polyhedrons_filters_123.values())
			min_diff = min(list_differences)
			fit_polyhedrons = list(polyhedrons_filters_123.keys())

			min_diff_v = min(list(polyhedrons_vol_diff.values()))
			min_diff_a = min(list(polyhedrons_area_diff.values()))
			radius_c = b_radius
			c_vertex = b_vertex

			parameters = open(parameters_file,'r')
			parameters = parameters.readlines()
			for line in parameters:
				match_correction = re.match(r'pore_sphr_vdw_corr = (.*)',line,re.M|re.I)
				if match_correction:
					sphr_vdw_corr= match_correction.group(1)
					break			

			if sphr_vdw_corr == 'VdW':
				radii = {}
				for atom in atoms_per_pore:
					symbol = structure[atom]
					vdw_r = Vdw_radii[symbol]
					radii[vdw_r] = None
				radii = list(radii.keys())
				sphr_vdw_corr = max(radii)

			if polytope_type == 'Regular_a' or polytope_type == 'Regular_b' :
				for polyhedron in polyhedrons_filters_123:
					if polyhedrons_filters_123[polyhedron] == min_diff:
						fit_polyhedron = polyhedron
				faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas = polyhedrons[fit_polyhedron](radius_c,c_vertex)
				Polyhedron = [faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas]
				diff_v = abs(Volume_of_polyhedron - polytope_volume)
				diff_a = abs(Area_of_polyhedron - polytope_area)
				volume_fit_percentage = 100 - (diff_v/polytope_volume)*100
				area_fit_percentage = 100 - (diff_a/polytope_area)*100
				table = PrettyTable()
				table.field_names = ["Fittest Polyhedron # " +str(c+1) , "Volume_of_polyhedron", "Area_of_polyhedron", 'Volume_of_polytope # '+str(c+1), 'Area_of_polytope # '+str(c+1), "Volume fit %", "Area fit %"]
				table.add_row([fit_polyhedron, Volume_of_polyhedron , Area_of_polyhedron, polytope_volume, polytope_area, volume_fit_percentage, area_fit_percentage])
				print(table)
				f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction = Polyhedron_Alignment(Polyhedron,Polytope,radius_c,c_vertex,atoms_per_pore,sphr_vdw_corr)
				fittest_polyhedron = [f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas]
				self.Fittest_Polyhedrons.append(fittest_polyhedron)


			elif polytope_type != 'Regular_a' and len(fit_polyhedrons) == 1 or polytope_type != 'Regular_b' and len(fit_polyhedrons) == 1 :
				fit_polyhedron = fit_polyhedrons[0]
				faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas = polyhedrons[fit_polyhedron](radius_c,c_vertex)
				Polyhedron = [faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas]
				diff_v = abs(Volume_of_polyhedron - polytope_volume)
				diff_a = abs(Area_of_polyhedron - polytope_area)
				volume_fit_percentage = 100 - (diff_v/polytope_volume)*100
				area_fit_percentage = 100 - (diff_a/polytope_area)*100
				table = PrettyTable()
				table.field_names = ["Fittest Polyhedron # " +str(c+1) , "Volume_of_polyhedron", "Area_of_polyhedron", 'Volume_of_polytope # '+str(c+1), 'Area_of_polytope # '+str(c+1), "Volume fit %", "Area fit %"]
				table.add_row([fit_polyhedron, Volume_of_polyhedron , Area_of_polyhedron, polytope_volume, polytope_area, volume_fit_percentage, area_fit_percentage])
				print(table)
				f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction = Polyhedron_Alignment(Polyhedron,Polytope,radius_c,c_vertex,atoms_per_pore,sphr_vdw_corr)
				fittest_polyhedron = [f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas]
				self.Fittest_Polyhedrons.append(fittest_polyhedron)

			elif polytope_type != 'Regular_a' and len(fit_polyhedrons) != 1 or polytope_type != 'Regular_b' and len(fit_polyhedrons) == 1 :
				print('More than one polyhedron can be used as fittest polyhedron ') 

				if polytope_type == 'Cantellation':
					print('The shape of the porosity suggests that Morfemo may interprete the '+ '\n' + 
					'small regular faces in the pore polytope as any, FACES or VERTICES ' + '\n' + 
					'of the fittest polyhedron.')
				if polytope_type == 'Truncation':
					print('The shape of the porosity suggests that Morfemo may interprete the '+ '\n' + 
					'small regular faces in the pore polytope as any, EDGES or VERTICES ' + '\n' + 
					'of the fittest polyhedron.')


				for i in range(len(fit_polyhedrons)):
					interpretation = None
					fit_polyhedron = fit_polyhedrons[i]
					faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas = polyhedrons[fit_polyhedron](radius_c,c_vertex)
					Polyhedron = [faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas]
					 #Interpretation of the shape
					if polytope_type == 'Cantellation':
						if len(polyhedron_vertices) == num_small_reg_faces:
							interpretation = 'VERTICES'
						else:
							interpretation = 'FACES'
					
					if polytope_type == 'Truncation':
						if len(polyhedron_vertices) < num_vertices_polytope:
							interpretation = 'VERTICES'
						elif len(polyhedron_vertices) >= num_vertices_polytope:
							interpretation = 'EDGES'


					print('Interpretation as ' + interpretation)
					diff_v = abs(Volume_of_polyhedron - polytope_volume)
					diff_a = abs(Area_of_polyhedron - polytope_area)
					volume_fit_percentage = 100 - (diff_v/polytope_volume)*100
					area_fit_percentage = 100 - (diff_a/polytope_area)*100
					table = PrettyTable()
					print('Number of Faces: ' + str(num_faces))
					table.field_names = ["Fittest Polyhedron # " +str(c+1) , "Volume_of_polyhedron", "Area_of_polyhedron", 'Volume_of_polytope # '+str(c+1), 'Area_of_polytope # '+str(c+1), "Volume fit %", "Area fit %"]
					table.add_row([fit_polyhedron, Volume_of_polyhedron , Area_of_polyhedron, polytope_volume, polytope_area, volume_fit_percentage, area_fit_percentage])
					print(table)

				print('Please introduce the name of one of the polyhedra above: ')
				user_selection = str(input())
				fit_polyhedron = user_selection
				faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas = polyhedrons[user_selection](radius_c,c_vertex)
				Polyhedron = [faces_polyhedron,polyhedron_vertices,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,areas_faces,faces_areas]
				###########################################Alignment#################################################
				f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction = Polyhedron_Alignment(Polyhedron,Polytope,radius_c,c_vertex,atoms_per_pore,sphr_vdw_corr)
				fittest_polyhedron = [f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas]
				self.Fittest_Polyhedrons.append(fittest_polyhedron)

			#####################################################################################################


			faces_list_o = []
			for face_i in faces_list:
				face_o = []
				for vert_i in face_i:
					vert_vec = np.array(vert_i)
					vert_vec = translation_transformation(vert_vec,o_direction)
					face_o.append(vert_vec)
				faces_list_o.append(np.array(face_o))

			fig_e = plt.figure()
			ex = fig_e.gca(projection='3d')
			solid = Poly3DCollection(faces_list_o)
			solid.set_edgecolor('b')
			solid.set_facecolor('k')
			solid.set_alpha(0.15)
			ex.add_collection3d(solid)

			e = 0
			for polygon in hull_polygons:
				xeg = []
				yeg = []
				zeg = []
				for edge in polygon:
					p1 = atoms_per_pore[edge[0]]
					p2 = atoms_per_pore[edge[1]]
					p1 = translation_transformation(np.array(p1),o_direction)
					p2 = translation_transformation(np.array(p2),o_direction)
					xeg = [p1[0],p2[0]]
					yeg = [p1[1],p2[1]]
					zeg = [p1[2],p2[2]]
					ex.plot(xeg,yeg,zeg,color = 'k',linewidth = 1.5) 
				e += 1
			if plots == 'y':
				plt.title('Morfemo: ' + self.basename + ' Porosity # ' + str(c+1))
				plt.axis('off')
				plt.show()
			
			elemental_triangles = {"Tetrahedron":[2,3,3,'c',2,radius_c,c_vertex],
								"Cube":[2,4,3,'c',2,radius_c,c_vertex],
								"Octahedron":[4,2,3,'c',0,radius_c,c_vertex],
								"Dodecahedron":[5,3,2,'c',1,radius_c,c_vertex],
				    			"Icosahedron":[2,3,5,'c',2,radius_c,c_vertex],
				    			"Truncated_Tetrahedron": [2,3,3,'s',2,radius_c,c_vertex],
				    			"Truncated_Cube": [2,4,3,'s',1,radius_c,c_vertex],
				    			"Truncated_Octahedron":[2,4,3,'s',2,radius_c,c_vertex],
				    			"Truncated_Dodecahedron": [2,5,3,'s',1,radius_c,c_vertex],
				    			"Truncated_Icosahedron": [2,5,3,'s',2,radius_c,c_vertex],
				    			"Cuboctahedron": [2,4,3,'c',0,radius_c,c_vertex],
				    			"Icosidodecahedron": [2,5,3,'c',0,radius_c,c_vertex],
				    			"Small_Rhombicuboctahedron": [2,4,3,'s',0,radius_c,c_vertex],
				    			"Small_Rhombicosidodecahedron": [2,5,3,'s',0,radius_c,c_vertex],
				    			"Great_Rhombicuboctahedron":[2,4,3,'i',radius_c,c_vertex],
				    			"Great_Rhombicosidodecahedron": [2,5,3,'i',radius_c,c_vertex],
				    			}

			triangle_parameters = elemental_triangles[fit_polyhedron]
			generator,relevant_points,triangle_type,perimeter_points = Irreducible_Volume(triangle_parameters,c_vertex)
			pore_radius = pores_diameters[c]/2.0

			#translation to the origin of coordinates
			#R = np.array([[1,0,0],[0,1,0],[0,0,1]])
			#Rf = np.array([[1,0,0],[0,1,0],[0,0,1]])
			Riv = np.dot(Rf,R)
			o_points = []
			
			generator = translation_transformation(generator,o_direction)

			for point in relevant_points:
				o_point = translation_transformation(point,o_direction)
				o_points.append(o_point)

			g_points = []
			for point in perimeter_points:
				g_point = translation_transformation(point,o_direction)
				g_points.append(g_point)


			#Rotation to align vectors
			gen = np.dot(Riv,np.transpose(generator))
			generator = np.array([gen.item(0),gen.item(1),gen.item(2)])
			generator = scaling_transformation(generator,scale_factor)
			aligned_points = []
			for point in o_points:
				align = np.dot(Riv,np.transpose(point))
				aligned_point = []
				aligned_point.append(align.item(0))
				aligned_point.append(align.item(1))
				aligned_point.append(align.item(2))
				aligned_points.append(aligned_point)

			aligned_points2 = []
			for point in g_points:
				align = np.dot(Riv,np.transpose(point))
				aligned_point = []
				aligned_point.append(align.item(0))
				aligned_point.append(align.item(1))
				aligned_point.append(align.item(2))
				aligned_points2.append(aligned_point)

			scaled_points = []
			for point in aligned_points:
				scaled_point = scaling_transformation(point,scale_factor)
				scaled_points.append(scaled_point)
#################################################################################
			scaled_points2 = []
			for point in aligned_points2:
				scaled_point = scaling_transformation(point,scale_factor)
				scaled_points2.append(scaled_point)

			#Back translation
			transformed_points = []
			generator_1 = generator#generator before back translation for the origin of corrdinates to the center of the porosity
			generator = translation_transformation(generator,o_direction*(-1))
			for point in scaled_points:
				t_point = translation_transformation(point,o_direction*(-1))
				transformed_points.append(t_point)
	##########################################################################3PILAAAAAAAAAAAAAAAAAAAAAs scaled aligned
			transformed_points2 = []
			for point in scaled_points2:
				t_point = translation_transformation(point,o_direction*(-1))
				transformed_points2.append(t_point)

			#Create the faces of the irreducible volume
###########################################################################PERSPECTIVE PLOTS#################################################################################################
			if triangle_type == 'i':
				P_proj = transformed_points[0]
				Q_proj = transformed_points[1]
				R_proj = transformed_points[2]
				p_proj = transformed_points[3]
				q_proj = transformed_points[4]
				r_proj = transformed_points[5]

			#Plane P
				x_plane_P = np.array([generator[0],q_proj[0],P_proj[0],r_proj[0]])
				y_plane_P = np.array([generator[1],q_proj[1],P_proj[1],r_proj[1]])
				z_plane_P = np.array([generator[2],q_proj[2],P_proj[2],r_proj[2]])
				verts_plane_P = []
				verts_plane_P.append(list(zip(x_plane_P,y_plane_P,z_plane_P)))

			#Plane Q
				x_plane_Q = np.array([generator[0],p_proj[0],Q_proj[0],r_proj[0]])
				y_plane_Q = np.array([generator[1],p_proj[1],Q_proj[1],r_proj[1]])
				z_plane_Q = np.array([generator[2],p_proj[2],Q_proj[2],r_proj[2]])
				verts_plane_Q = []
				verts_plane_Q.append(list(zip(x_plane_Q,y_plane_Q,z_plane_Q)))

			#Plane R
				x_plane_R = np.array([generator[0],q_proj[0],R_proj[0],p_proj[0]])
				y_plane_R = np.array([generator[1],q_proj[1],R_proj[1],p_proj[1]])
				z_plane_R = np.array([generator[2],q_proj[2],R_proj[2],p_proj[2]])
				verts_plane_R = []
				verts_plane_R.append(list(zip(x_plane_R,y_plane_R,z_plane_R)))

				plane_P = Poly3DCollection(verts_plane_P,facecolor='b',alpha=0.8)
				plane_Q = Poly3DCollection(verts_plane_Q,facecolor='y',alpha=0.8)
				plane_R = Poly3DCollection(verts_plane_R,facecolor='r',alpha=0.8)

				planes = [plane_P,plane_Q,plane_R]

			if triangle_type == 'c1':
				proj = transformed_points[0]
				E1_proj = transformed_points[1]
				E2_proj = transformed_points[2]

				#Plane of the Opposite Reference Point
				x_plane_E1 = np.array([generator[0],proj[0],E1_proj[0]])
				y_plane_E1 = np.array([generator[1],proj[1],E1_proj[1]])
				z_plane_E1 = np.array([generator[2],proj[2],E1_proj[2]])
				verts_plane_E1 = []
				verts_plane_E1.append(list(zip(x_plane_E1,y_plane_E1,z_plane_E1)))

				#Plane of the Secondary Reference Point
				x_plane_E2 = np.array([generator[0],proj[0],E2_proj[0]])
				y_plane_E2 = np.array([generator[1],proj[1],E2_proj[1]])
				z_plane_E2 = np.array([generator[2],proj[2],E2_proj[2]])
				verts_plane_E2 = []
				verts_plane_E2.append(list(zip(x_plane_E2,y_plane_E2,z_plane_E2)))
				plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='r',alpha=0.9)
				plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='y',alpha=0.9)
				planes = [plane_E1,plane_E2]

			if triangle_type == 'c2':
				projR = transformed_points[0]
				M_proj = transformed_points[1]

			#Plane Referece Point
				x_plane = np.array([generator[0],projR[0],M_proj[0]])
				y_plane = np.array([generator[1],projR[1],M_proj[1]])
				z_plane = np.array([generator[2],projR[2],M_proj[2]])
				verts_plane = []
				verts_plane.append(list(zip(x_plane,y_plane,z_plane)))
				plane = Poly3DCollection(verts_plane,facecolor='r',alpha=0.9)
				planes = [plane]
				

			if triangle_type == 's1':
				Ref_proj = transformed_points[0]
				E1_proj = transformed_points[1]
				E2_proj = transformed_points[2]
				e1_proj = transformed_points[3]
				e2_proj = transformed_points[4]
				
				#Plane Referece Point
				x_plane_Ref = np.array([generator[0],e2_proj[0],Ref_proj[0],e1_proj[0]])
				y_plane_Ref = np.array([generator[1],e2_proj[1],Ref_proj[1],e1_proj[1]])
				z_plane_Ref = np.array([generator[2],e2_proj[2],Ref_proj[2],e1_proj[2]])
				verts_plane_Ref = []
				verts_plane_Ref.append(list(zip(x_plane_Ref,y_plane_Ref,z_plane_Ref)))

				#Plane E1
				x_plane_E1 = np.array([generator[0],E1_proj[0],e2_proj[0]])
				y_plane_E1 = np.array([generator[1],E1_proj[1],e2_proj[1]])
				z_plane_E1 = np.array([generator[2],E1_proj[2],e2_proj[2]])
				verts_plane_E1 = []
				verts_plane_E1.append(list(zip(x_plane_E1,y_plane_E1,z_plane_E1)))

				#Plane E2
				x_plane_E2 = np.array([generator[0],e1_proj[0],E2_proj[0]])
				y_plane_E2 = np.array([generator[1],e1_proj[1],E2_proj[1]])
				z_plane_E2 = np.array([generator[2],e1_proj[2],E2_proj[2]])
				verts_plane_E2 = []
				verts_plane_E2.append(list(zip(x_plane_E2,y_plane_E2,z_plane_E2)))
				plane_Ref = Poly3DCollection(verts_plane_Ref,facecolor='b',alpha=0.9)
				plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='y',alpha=0.9)
				plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='r',alpha=0.9)
				planes = [plane_Ref,plane_E1,plane_E2]

			if triangle_type == 's2':
				proj = transformed_points[0]
				Op_proj = transformed_points[1]
				Sec_proj = transformed_points[2]
				projRef = transformed_points[3]

			#Plane Oppsite Reference Point
				x_plane_Op = np.array([generator[0],proj[0],Op_proj[0],projRef[0]])
				y_plane_Op = np.array([generator[1],proj[1],Op_proj[1],projRef[1]])
				z_plane_Op = np.array([generator[2],proj[2],Op_proj[2],projRef[2]])
				verts_plane_Op = []
				verts_plane_Op.append(list(zip(x_plane_Op,y_plane_Op,z_plane_Op)))

			#Plane Secondary Reference Point
				x_plane_Sec = np.array([generator[0],proj[0],Sec_proj[0]])
				y_plane_Sec = np.array([generator[1],proj[1],Sec_proj[1]])
				z_plane_Sec = np.array([generator[2],proj[2],Sec_proj[2]])
				verts_plane_Sec = []
				verts_plane_Sec.append(list(zip(x_plane_Sec,y_plane_Sec,z_plane_Sec)))
				plane_Op = Poly3DCollection(verts_plane_Op,facecolor='r',alpha=0.9)
				plane_Sec = Poly3DCollection(verts_plane_Sec,facecolor='y',alpha=0.9)
				planes = [plane_Op,plane_Sec]

#######################################################################ZERO PLOTS################################################################################################
			if triangle_type == 'i':
				#P_proj = transformed_points[0]
				#Q_proj = transformed_points[1]
				#R_proj = transformed_points[2]
				#p_proj = transformed_points[3]
				#q_proj = transformed_points[4]
				#r_proj = transformed_points[5]
				P_proj = scaled_points[0]
				Q_proj = scaled_points[1]
				R_proj = scaled_points[2]
				p_proj = scaled_points[3]
				q_proj = scaled_points[4]
				r_proj = scaled_points[5]
			#Plane P
				x_plane_P = np.array([generator_1[0],q_proj[0],P_proj[0],r_proj[0]])
				y_plane_P = np.array([generator_1[1],q_proj[1],P_proj[1],r_proj[1]])
				z_plane_P = np.array([generator_1[2],q_proj[2],P_proj[2],r_proj[2]])
				verts_plane_P = []
				verts_plane_P.append(list(zip(x_plane_P,y_plane_P,z_plane_P)))

			#Plane Q
				x_plane_Q = np.array([generator_1[0],p_proj[0],Q_proj[0],r_proj[0]])
				y_plane_Q = np.array([generator_1[1],p_proj[1],Q_proj[1],r_proj[1]])
				z_plane_Q = np.array([generator_1[2],p_proj[2],Q_proj[2],r_proj[2]])
				verts_plane_Q = []
				verts_plane_Q.append(list(zip(x_plane_Q,y_plane_Q,z_plane_Q)))

			#Plane R
				x_plane_R = np.array([generator_1[0],q_proj[0],R_proj[0],p_proj[0]])
				y_plane_R = np.array([generator_1[1],q_proj[1],R_proj[1],p_proj[1]])
				z_plane_R = np.array([generator_1[2],q_proj[2],R_proj[2],p_proj[2]])
				verts_plane_R = []
				verts_plane_R.append(list(zip(x_plane_R,y_plane_R,z_plane_R)))

				plane_P = Poly3DCollection(verts_plane_P,facecolor='b',alpha=0.8)
				plane_Q = Poly3DCollection(verts_plane_Q,facecolor='y',alpha=0.8)
				plane_R = Poly3DCollection(verts_plane_R,facecolor='r',alpha=0.8)

				z_planes = [plane_P,plane_Q,plane_R]

			if triangle_type == 'c1':
				#proj = transformed_points[0]
				#E1_proj = transformed_points[1]
				#E2_proj = transformed_points[2]
				proj = scaled_points[0]
				E1_proj = scaled_points[1]
				E2_proj = scaled_points[2]

				#Plane of the Opposite Reference Point
				x_plane_E1 = np.array([generator_1[0],proj[0],E1_proj[0]])
				y_plane_E1 = np.array([generator_1[1],proj[1],E1_proj[1]])
				z_plane_E1 = np.array([generator_1[2],proj[2],E1_proj[2]])
				verts_plane_E1 = []
				verts_plane_E1.append(list(zip(x_plane_E1,y_plane_E1,z_plane_E1)))

				#Plane of the Secondary Reference Point
				x_plane_E2 = np.array([generator_1[0],proj[0],E2_proj[0]])
				y_plane_E2 = np.array([generator_1[1],proj[1],E2_proj[1]])
				z_plane_E2 = np.array([generator_1[2],proj[2],E2_proj[2]])
				verts_plane_E2 = []
				verts_plane_E2.append(list(zip(x_plane_E2,y_plane_E2,z_plane_E2)))
				plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='r',alpha=0.9)
				plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='y',alpha=0.9)
				z_planes = [plane_E1,plane_E2]

			if triangle_type == 'c2':
				#projR = transformed_points[0]
				#M_proj = transformed_points[1]
				projR = scaled_points[0]
				M_proj = scaled_points[1]
			#Plane Referece Point
				x_plane = np.array([generator_1[0],projR[0],M_proj[0]])
				y_plane = np.array([generator_1[1],projR[1],M_proj[1]])
				z_plane = np.array([generator_1[2],projR[2],M_proj[2]])
				verts_plane = []
				verts_plane.append(list(zip(x_plane,y_plane,z_plane)))
				plane = Poly3DCollection(verts_plane,facecolor='r',alpha=0.9)
				z_planes = [plane]
				

			if triangle_type == 's1':
				#Ref_proj = transformed_points[0]
				#E1_proj = transformed_points[1]
				#E2_proj = transformed_points[2]
				#e1_proj = transformed_points[3]
				#e2_proj = transformed_points[4]
				Ref_proj = scaled_points[0]
				E1_proj = scaled_points[1]
				E2_proj = scaled_points[2]
				e1_proj = scaled_points[3]
				e2_proj = scaled_points[4]
				
				#Plane Referece Point
				x_plane_Ref = np.array([generator_1[0],e2_proj[0],Ref_proj[0],e1_proj[0]])
				y_plane_Ref = np.array([generator_1[1],e2_proj[1],Ref_proj[1],e1_proj[1]])
				z_plane_Ref = np.array([generator_1[2],e2_proj[2],Ref_proj[2],e1_proj[2]])
				verts_plane_Ref = []
				verts_plane_Ref.append(list(zip(x_plane_Ref,y_plane_Ref,z_plane_Ref)))

				#Plane E1
				x_plane_E1 = np.array([generator_1[0],E1_proj[0],e2_proj[0]])
				y_plane_E1 = np.array([generator_1[1],E1_proj[1],e2_proj[1]])
				z_plane_E1 = np.array([generator_1[2],E1_proj[2],e2_proj[2]])
				verts_plane_E1 = []
				verts_plane_E1.append(list(zip(x_plane_E1,y_plane_E1,z_plane_E1)))

				#Plane E2
				x_plane_E2 = np.array([generator_1[0],e1_proj[0],E2_proj[0]])
				y_plane_E2 = np.array([generator_1[1],e1_proj[1],E2_proj[1]])
				z_plane_E2 = np.array([generator_1[2],e1_proj[2],E2_proj[2]])
				verts_plane_E2 = []
				verts_plane_E2.append(list(zip(x_plane_E2,y_plane_E2,z_plane_E2)))
				plane_Ref = Poly3DCollection(verts_plane_Ref,facecolor='b',alpha=0.9)
				plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='y',alpha=0.9)
				plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='r',alpha=0.9)
				z_planes = [plane_Ref,plane_E1,plane_E2]

			if triangle_type == 's2':
				#proj = transformed_points[0]
				#Op_proj = transformed_points[1]
				#Sec_proj = transformed_points[2]
				#projRef = transformed_points[3]
				proj = scaled_points[0]
				Op_proj = scaled_points[1]
				Sec_proj = scaled_points[2]
				projRef = scaled_points[3]

			#Plane Oppsite Reference Point
				x_plane_Op = np.array([generator_1[0],proj[0],Op_proj[0],projRef[0]])
				y_plane_Op = np.array([generator_1[1],proj[1],Op_proj[1],projRef[1]])
				z_plane_Op = np.array([generator_1[2],proj[2],Op_proj[2],projRef[2]])
				verts_plane_Op = []
				verts_plane_Op.append(list(zip(x_plane_Op,y_plane_Op,z_plane_Op)))

			#Plane Secondary Reference Point
				x_plane_Sec = np.array([generator_1[0],proj[0],Sec_proj[0]])
				y_plane_Sec = np.array([generator_1[1],proj[1],Sec_proj[1]])
				z_plane_Sec = np.array([generator_1[2],proj[2],Sec_proj[2]])
				verts_plane_Sec = []
				verts_plane_Sec.append(list(zip(x_plane_Sec,y_plane_Sec,z_plane_Sec)))
				plane_Op = Poly3DCollection(verts_plane_Op,facecolor='r',alpha=0.9)
				plane_Sec = Poly3DCollection(verts_plane_Sec,facecolor='y',alpha=0.9)
				z_planes = [plane_Op,plane_Sec]
#################################################################################################################################################################################
			
			fig_f = plt.figure()
			fx = fig_f.gca(projection='3d')
			fx.set_xlim((radius_c,-radius_c))
			fx.set_ylim((radius_c,-radius_c))
			fx.set_zlim((radius_c,-radius_c))
			f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = fittest_polyhedron

			faces_list_o = []
			for face_i in faces_list:
				face_o = []
				for vert_i in face_i:
					vert_vec = np.array(vert_i)
					vert_vec = translation_transformation(vert_vec,o_direction)
					face_o.append(vert_vec)
				faces_list_o.append(np.array(face_o))
			
			solid = Poly3DCollection(faces_list_o)
			solid.set_edgecolor('k')
			solid.set_facecolor('k')
			solid.set_alpha(0.15)
			fx.add_collection3d(solid)
			e = 0

			for polygon in hull_polygons:
				xeg = []
				yeg = []
				zeg = []
				for edge in polygon:
					p1 = atoms_per_pore[edge[0]]
					p2 = atoms_per_pore[edge[1]]
					p1 = translation_transformation(np.array(p1),o_direction)
					p2 = translation_transformation(np.array(p2),o_direction)
					xeg = [p1[0],p2[0]]
					yeg = [p1[1],p2[1]]
					zeg = [p1[2],p2[2]]
					#fx.plot(xeg,yeg,zeg,color = 'k',linewidth = 0.5) 
			for plane in z_planes:
				fx.add_collection3d(plane)
#############################################################PERSPECTIVE IRREDUCIBLE PLOT##################################################################################################							
			ver = {}
			irreducible_vol = []

			for point in transformed_points2:
				ver[(point[0],point[1],point[2])] = None
				irreducible_vol.append(Vector(point[0],point[1],point[2]))

			ver[(c_vertex[0],c_vertex[1],c_vertex[2])] = None
			ver[(generator[0],generator[1],generator[2])] = None
			irreducible_vol.append(Vector(c_vertex[0],c_vertex[1],c_vertex[2]))
			irreducible_vol.append(Vector(generator[0],generator[1],generator[2]))

			ver = list(ver.keys())
			xIV = []
			yIV = []
			zIV = []
			for vert in ver:
				xIV.append(vert[0])
				yIV.append(vert[1])
				zIV.append(vert[2])

			irreducible_hull = Hull(irreducible_vol)
			triangles = irreducible_hull.Export_Triangles()
			irred_vol_verts,irred_vol_inds,irred_types_dict,irred_area_dict = General_Polygonalization(irreducible_hull,ver)
			self.Irreducible_Volumes.append(irred_vol_verts)
			self.Generator_Triangles.append(planes)
			irreducible_volumes_points.append(ver)
###############################################################ZEROS IRREDUCIBLE PLOT############################################################################################
			z_ver = {}
			z_irreducible_vol = []
			for point in scaled_points2:
				z_ver[(point[0],point[1],point[2])] = None
				z_irreducible_vol.append(Vector(point[0],point[1],point[2]))

			z_c_vertex = np.array((0,0,0))
			z_ver[(z_c_vertex[0],z_c_vertex[1],z_c_vertex[2])] = None
			z_ver[(generator_1[0],generator_1[1],generator_1[2])] = None
			z_irreducible_vol.append(Vector(z_c_vertex[0],z_c_vertex[1],z_c_vertex[2]))
			z_irreducible_vol.append(Vector(generator_1[0],generator_1[1],generator_1[2]))

			z_ver = list(z_ver.keys())
			xIV = []
			yIV = []
			zIV = []
			for vert in z_ver:
				xIV.append(vert[0])
				yIV.append(vert[1])
				zIV.append(vert[2])

			z_irreducible_hull = Hull(z_irreducible_vol)
			z_triangles = z_irreducible_hull.Export_Triangles()
			z_irred_vol_verts,z_irred_vol_inds,z_irred_types_dict,z_irred_area_dict = General_Polygonalization(z_irreducible_hull,z_ver)
			z_irred_verts_coords = {}
			for plane in z_irred_vol_verts:
				for edge in plane:
					p1 = edge[0]
					p2 = edge[1]
					z_irred_verts_coords[p1] = None
					z_irred_verts_coords[p2] = None
					xg = [p1[0],p2[0]]
					yg = [p1[1],p2[1]]
					zg = [p1[2],p2[2]]
					fx.plot(xg,yg,zg,color = 'k',linewidth = 1.5)
			
			z_irred_verts_coords = list(z_irred_verts_coords.keys())
			xir = []
			yir = []
			zir = []
			for vert_coords in z_irred_verts_coords:
				xir.append(vert_coords[0])
				yir.append(vert_coords[1])
				zir.append(vert_coords[2])

			fx.scatter(xir,yir,zir,color = 'b',alpha = 0.5, label = 'Critical Points')
			fx.scatter(generator_1[0],generator_1[1],generator_1[2],color = 'k',alpha=1,label = 'Generator vertex')

			if plots == 'y':
				plt.legend()
				plt.title('Fragmento: ' + self.basename + ' Porosity # ' + str(c+1))
				plt.axis('off')
				plt.show()

			print("////////////////////////Irreducible Volume Generated////////////////////////")
			c += 1

		return(irreducible_volumes_points)
	if __name__ == "__Fragmento__":
	    Fragmento()

	def find_critical_points(self):
		coordinates,atoms_of_porosities,centers_of_porosities,lst,pores_diameters,structure,VdW_radii,elements_list = self.Kerno(self.basename,self.dat,self.plots,self.parameters_file)
		centers_of_porosities = [center[0] for center in centers_of_porosities]
		irreducible_volume = self.Fragmento(atoms_of_porosities,centers_of_porosities,lst,pores_diameters,structure,VdW_radii,elements_list,self.parameters_file,self.plots)
		return(irreducible_volume)
	if __name__ == "__find_ad_sites__":
	    find_ad_sites()

	def kerno_only(self):
		coordinates,atoms_of_porosities,centers_of_porosities,lst,pores_diameters,structure,VdW_radii,elements_list = self.Kerno(self.basename,self.dat,self.plots,self.parameters_file)
		centers_of_porosities = [center[0] for center in centers_of_porosities]
		return(centers_of_porosities,atoms_of_porosities)

	def default(self):
		coordinates,atoms_of_porosities,centers_of_porosities,lst,pores_diameters,structure,VdW_radii,elements_list = self.Kerno(self.basename,self.dat,self.plots,self.parameters_file)
		centers_of_porosities = [center[0] for center in centers_of_porosities]
		convex_hulls = []
		Polytopes = []
		c = 0
		for atoms_per_pore in atoms_of_porosities:
			X,Y,Z,convex_hull,Polytope,ver_left,polytope_volume,polytope_area = self.Mozaiko(atoms_per_pore,self.parameters_file)
			convex_hulls.append(convex_hull)
			Polytopes.append(Polytope)
			if self.plots == 'y':
				fig_cg = plt.figure(1)
				#############################################################################################################################################
				cgx = fig_cg.gca(projection='3d')
				triangles = (convex_hull.Export_Triangles())
				cgx.plot_trisurf(X, Y, triangles, Z,edgecolors='k',linewidth=1.5,cmap = 'winter',alpha=0.25)
				plt.axis('off')
				plt.title('Mozaiko: Delaunay Triangulation stage' + '\n' + self.basename + ' Porosity # ' + str(c+1))
				#plt.show()
	            ##############################################################################################################################################
				c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = Polytope

				fig_d = plt.figure(2)
				dx = fig_d.gca(projection='3d')
				d = 0
				for polygon in hull_polygons:
					xe = []
					ye = []
					ze = []
					for edge in polygon:
						p1 = atoms_per_pore[edge[0]]
						p2 = atoms_per_pore[edge[1]]
						xe = [p1[0],p2[0]]
						ye = [p1[1],p2[1]]
						ze = [p1[2],p2[2]]
						dx.plot(xe,ye,ze,color = 'k',linewidth = 1.5) 
				d += 1
				plt.axis('off')
				plt.title('Mozaiko: Polygonalization stage' + '\n' + self.basename + ' Porosity # ' + str(c+1))
				plt.show()
			c += 1
		return(centers_of_porosities,atoms_of_porosities,convex_hulls,Polytopes)
