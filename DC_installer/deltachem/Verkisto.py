import re
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import subprocess as sub
import time 
from pylab import *
import statistics 
import sys
import numpy as np
from math import *
from itertools import combinations
import sympy
import pyfiglet
from prettytable import *
import subprocess as sub
import time
import os
from pylab import *
#Import required delta_chem modules
import deltachem.DC_library
from deltachem.DC_library import *
import deltachem.Delta_Chem
from deltachem.Delta_Chem import *
import deltachem.Convex_Hull_I 
from deltachem.Convex_Hull_I import Vector,Hull,Face,Edge,Vertex
import deltachem.irreducible_volume 
from deltachem.irreducible_volume import Irreducible_Volume
import deltachem.DC_polyhedrons
from deltachem.DC_polyhedrons import *
import deltachem.Wythoff
from deltachem.Wythoff import Wythoff

if len(sys.argv) < 3:
	print ('Usage: delta_chem <structure_file_basename> <psd_file_name.dat> <options>')
	sys.exit(1)

basename = sys.argv[1]
dat = sys.argv[2]
options = []
if len(sys.argv) >= 4:
	options = sys.argv[3:]

parameters = 'parameters.txt'
ghost_protocol = 'off'
path_tst = os.path.dirname(__file__)

if 'dc_par' in options:
	ghost_protocol = 'on'
	ghostbuster = open(path_tst + "/" + 'who_you_gonna_call.txt')
	ghostbuster = ghostbuster.readlines()
	ghost_title = pyfiglet.figlet_format('Ghost Protocol Activated',font =  "alphabet" )
	ghost_warning = pyfiglet.figlet_format('You have 30 seconds to introduce changes',font =  "digital" )
	print(ghost_title)
	for line in ghostbuster:
		print(line)
	
	parameters = Ghost_Protocol(parameters)
	time.sleep(30)

if 'no_plots' in options:
	plots = 'n'
else:
	plots = 'y'
#Initialize the class Delta Chem
if os.path.isfile('parachutes.txt') == True:
	parameters = 'parachutes.txt'
material = Delta_Chem(basename,dat,plots,parameters)

if 'kerno_only' in options:
	centers_of_porosities,atoms_of_porosities = material.kerno_only()
	#Get the atributes of the Delta Chem object
	structure = getattr(material,'Structure')
	coordinates = list(structure.keys())
	atoms_in_material = list(structure.values())
	atomic_species_list = getattr(material,'Atomic_Species')
	atomic_species_grouped = getattr(material,'Atomic_Species_Grouped')
	centers_of_porosities = getattr(material,'Centers')
	contributing_atoms = getattr(material,'Contributing_Atoms')
	pores_radii = getattr(material,'Pore_radii')

elif 'dc_approach' in options:
	#Execute Delta Chem Analysis
	critical_points = material.find_critical_points()
	structure = getattr(material,'Structure')
	coordinates = list(structure.keys())
	atoms_in_material = list(structure.values())
	atomic_species_list = getattr(material,'Atomic_Species')
	atomic_species_grouped = getattr(material,'Atomic_Species_Grouped')
	centers_of_porosities = getattr(material,'Centers')
	contributing_atoms = getattr(material,'Contributing_Atoms')
	pores_radii = getattr(material,'Pore_radii')
	polytopes = getattr(material,'Polytopes')
	convex_hulls = getattr(material,'Convex_Hulls')
	fittest_polyhedrons = getattr(material,'Fittest_Polyhedrons')
	irreducible_volumes = getattr(material,'Irreducible_Volumes')
	generator_triangles = getattr(material,'Generator_Triangles')
	porosity_volumes = getattr(material,'Porosity_Volumes')
	porosity_areas = getattr(material,'Porosity_Areas')
	porosity_num_faces = getattr(material,'Porosity_num_faces')

else:
	centers_of_porosities,atoms_of_porosities,convex_hull,Polytope = material.default()
	structure = getattr(material,'Structure')
	coordinates = list(structure.keys())
	atoms_in_material = list(structure.values())
	atomic_species_list = getattr(material,'Atomic_Species')
	atomic_species_grouped = getattr(material,'Atomic_Species_Grouped')
	centers_of_porosities = getattr(material,'Centers')
	contributing_atoms = getattr(material,'Contributing_Atoms')
	pores_radii = getattr(material,'Pore_radii')
	polytopes = getattr(material,'Polytopes')
	convex_hulls = getattr(material,'Convex_Hulls')
	porosity_volumes = getattr(material,'Porosity_Volumes')
	porosity_areas = getattr(material,'Porosity_Areas')
	porosity_num_faces = getattr(material,'Porosity_num_faces')


###########################################################Plots anf GIFs###############################################################
########################################################################################################################################
################################################################ gif ####################################################################

if 'gif' in options:
	fig_ag= plt.figure()
	agx = fig_ag.gca(projection ='3d')
	colors = ['m','k','y','b','g']
	colors_P = ['r','b','g']
	ag = 0
	for element in atomic_species_grouped:
		xe = []
		ye = []
		ze = []
		atoms_of_element = atomic_species_grouped[element]
		for atom in atoms_of_element:
			xe.append(atom[0])
			ye.append(atom[1])
			ze.append(atom[2])
		agx.scatter(xe,ye,ze,color = colors[ag],alpha = 1,label = atomic_species_list[ag])
		ag += 1

	for angle in range(0,90,2):
		agx.view_init(30,angle)
		filename = 'legando_stage'+ str(angle) + basename +'.png'
		plt.legend()
		plt.savefig(filename, dpi = 96)
		plt.gca()

	sub.call(('convert','-delay','10','legando_stage*.png',basename +'_legando.gif'))
	args = ('rm', '-rf', 'legando_stage*.png')
	sub.call('%s %s %s' % args, shell=True)

	
#########################################################################################################################################
########################################################Plot the results of Kerno########################################################
u = np.linspace(0,  2*np.pi, 25)
v = np.linspace(0, np.pi, 25)
alp = [0.75,0.75,0.75]
size = [15.2,17.00,12.00,13.90]
colors_P = ['k','m','b','r','y','g']
#Plot the centers of the porosities
if 'perspective' in options:
	fig_b = plt.figure()
	bx = fig_b.gca(projection='3d')
	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	bx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')

	b = 0
	for center in centers_of_porosities:
		pore_radius = pores_radii[b]/2
		bx.scatter(center[0],center[1],center[2],color='r',alpha = 1,label='Center of porosity '+str(b+1))
		xS = pore_radius * np.outer(np.cos(u), np.sin(v)) + center[0]
		yS = pore_radius * np.outer(np.sin(u), np.sin(v)) + center[1]
		zS = pore_radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
		bx.plot_surface(xS, yS, zS, cmap = 'flag',alpha = 0.25)
		b += 1

	bb = 0
	for atoms_per_pore in contributing_atoms:
		xa = []
		ya = []
		za = []
		for atom in atoms_per_pore:
			xa.append(atom[0])
			ya.append(atom[1])
			za.append(atom[2])
		bx.scatter(xa,ya,za,color=colors_P[bb],alpha = 0.7,label='Atoms contributing to porosity '+str(bb+1))
		bb += 1
	#plt.title('Kerno')
	plt.axis('off')
	plt.legend()
	plt.show()

################################################################ gif #######################################################################

if 'gif' in options:
	fig_bg = plt.figure()
	bgx = fig_bg.gca(projection='3d')

	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	bgx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in structure')

	bg = 0
	for center in centers_of_porosities:
		pore_radius = pores_radii[bg]/2
		bgx.scatter(center[0],center[1],center[2],color='r',alpha = 1,label='Center of porosity '+str(bg+1))
		xS = pore_radius * np.outer(np.cos(u), np.sin(v)) + center[0]
		yS = pore_radius * np.outer(np.sin(u), np.sin(v)) + center[1]
		zS = pore_radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
		bgx.plot_surface(xS, yS, zS, cmap = 'summer',alpha = 0.25)
		bg += 1

	bbg = 0
	for atoms_per_pore in contributing_atoms:
		xa = []
		ya = []
		za = []
		for atom in atoms_per_pore:
			xa.append(atom[0])
			ya.append(atom[1])
			za.append(atom[2])
		bgx.scatter(xa,ya,za,color=colors_P[bbg],alpha = 0.7,label='Atoms contributing to porosity '+str(bbg+1))
		bbg += 1
	for angle in range(90,180,2):
		bgx.view_init(30,angle)
		filename = 'kerno_stage'+ str(angle) + basename +'.png'
		plt.legend()
		plt.savefig(filename, dpi = 96)
		plt.gca()
	sub.call(('convert','-delay','10','kerno_stage*.png',basename +'_kerno.gif'))
	args = ('rm', '-rf', 'kerno_stage*.png')
	sub.call('%s %s %s' % args, shell=True)

############################################################################################################################################
#################################################### Plot the results of Mozaiko ###########################################################
########################################################## a) Triangulation ################################################################
if 'perspective' in options and 'kerno_only' not in options:
	fig_c = plt.figure()
	cx = fig_c.gca(projection='3d')

	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	cx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')
	colorsss = ['b','r','g']

	sss=0
	for atoms_per_pore in contributing_atoms:
		xcg = []
		ycg = []
		zcg = []
		for atom in atoms_per_pore:
			xcg.append(atom[0])
			ycg.append(atom[1])
			zcg.append(atom[2])
		cx.scatter(xcg,ycg,zcg,color=colorsss[sss],alpha = 1,linewidths = 2,label = 'Contributing atoms Pore # ' + str(sss+1))
		sss += 1

	for convex_hull in convex_hulls:
		X,Y,Z = convex_hulls[convex_hull]
		triangles = (convex_hull.Export_Triangles())
		cx.plot_trisurf(X, Y, triangles, Z,edgecolors='k',linewidth=0.75,cmap = 'flag',alpha=0.25)
	#plt.title('Mozaiko: Delaunay Triangulation stage' + '\n' + basename)
	plt.axis('off')
	plt.legend()
	plt.show()

################################################################ gif #######################################################################
if 'gif' in options and 'kerno_only' not in options:
	fig_cg = plt.figure()
	cgx = fig_cg.gca(projection='3d')

	for j in coordinates:
		cgx.scatter(j[0],j[1],j[2],color='k',alpha = 0.15)

	for atoms_per_pore in contributing_atoms:
		xcg = []
		ycg = []
		zcg = []
		for atom in atoms_per_pore:
			xcg.append(atom[0])
			ycg.append(atom[1])
			zcg.append(atom[2])
		cgx.scatter(xcg,ycg,zcg,color='r',alpha = 1,linewidths = 2)

	for convex_hull in convex_hulls:
		X,Y,Z = convex_hulls[convex_hull]
		triangles = (convex_hull.Export_Triangles())
		cgx.plot_trisurf(X, Y, triangles, Z,edgecolors='k',linewidth=1.5,cmap = 'autumn',alpha=0.25)
	for angle in range(180,270,2):
		cgx.view_init(30,angle)
		filename = 'mozaiko_tri_stage'+ str(angle) + basename +'.png'
		plt.savefig(filename, dpi = 96)
		plt.gca()

	sub.call(('convert','-delay','10','mozaiko_tri_stage*.png',basename +'_mozaiko_tri.gif'))
	args = ('rm', '-rf', 'mozaiko_tri_stage*.png')
	sub.call('%s %s %s' % args, shell=True)
############################################################ b) Polygonalization ################################################################
if 'perspective' in options and 'kerno_only' not in options:
	fig_d = plt.figure()
	dx = fig_d.gca(projection='3d')

	d = 0
	for polytope in polytopes:
		atoms = contributing_atoms[d]
		c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytope
		for polygon in hull_polygons:
			xe = []
			ye = []
			ze = []
			for edge in polygon:
				p1 = atoms[edge[0]]
				p2 = atoms[edge[1]]
				xe = [p1[0],p2[0]]
				ye = [p1[1],p2[1]]
				ze = [p1[2],p2[2]]
				dx.plot(xe,ye,ze,color = 'k',linewidth = 1.5) 
		d += 1

	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	dx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')
	plt.title('Mozaiko: Polygonalization stage' + '\n' + basename)
	plt.legend()
	plt.show()

################################################################ gif ###########################################################################

if 'gif' in options and 'kerno_only' not in options:
	fig_dg = plt.figure()
	dgx = fig_dg.gca(projection='3d')

	dg = 0
	for polytope in polytopes:
		atoms = contributing_atoms[dg]
		c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytope
		for polygon in hull_polygons:
			xdg = []
			ydg = []
			zdg = []
			for edge in polygon:
				p1 = atoms[edge[0]]
				p2 = atoms[edge[1]]
				xdg = [p1[0],p2[0]]
				ydg = [p1[1],p2[1]]
				zdg = [p1[2],p2[2]]
				dgx.plot(xdg,ydg,zdg,color = 'b',linewidth = 1.5) 
		dg += 1

	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	dgx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')

	for angle in range(270,360,2):
		dgx.view_init(30,angle)
		filename = 'mozaiko_poly_stage'+ str(angle) + basename +'.png'
		plt.savefig(filename, dpi = 96)
		plt.gca()

	sub.call(('convert','-delay','10','mozaiko_poly_stage*.png',basename +'_mozaiko_poly.gif'))
	args = ('rm', '-rf', 'mozaiko_poly_stage*.png')
	sub.call('%s %s %s' % args, shell=True)

################################################################################################################################################
############################################################ Results of Fragmento ##############################################################
########################################################## a) Fitting and Alignment ############################################################
if 'perspective' in options and 'dc_approach' in options:
	fig_e = plt.figure()
	ex = fig_e.gca(projection='3d')
	for polyhedron in fittest_polyhedrons:
		f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = polyhedron
		solid = Poly3DCollection(faces_list)
		solid.set_edgecolor('k')
		solid.set_facecolor('r')
		solid.set_alpha(0.15)
		ex.add_collection3d(solid)

	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	ex.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')

	e = 0
	for polytope in polytopes:
		atoms = contributing_atoms[e]
		c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytope
		for polygon in hull_polygons:
			xeg = []
			yeg = []
			zeg = []
			for edge in polygon:
				p1 = atoms[edge[0]]
				p2 = atoms[edge[1]]
				xeg = [p1[0],p2[0]]
				yeg = [p1[1],p2[1]]
				zeg = [p1[2],p2[2]]
				ex.plot(xeg,yeg,zeg,color = 'k',linewidth = 1.5) 
		e += 1
	plt.title('Morfemo: ' + basename)
	plt.legend()
	plt.show()
################################################################ gif ###########################################################################
if 'gif' in options and 'dc_approach' in option:
	fig_eg = plt.figure()
	egx = fig_eg.gca(projection='3d')
	for polyhedron in fittest_polyhedrons:
		f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = polyhedron
		solid = Poly3DCollection(faces_list)
		solid.set_edgecolor('k')
		solid.set_facecolor('r')
		solid.set_alpha(0.15)
		egx.add_collection3d(solid)
	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	egx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')

	eg = 0
	for polytope in polytopes:
		atoms = contributing_atoms[eg]
		c_hull_v_area,c_hull_i_area,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytope
		for polygon in hull_polygons:
			xeg = []
			yeg = []
			zeg = []
			for edge in polygon:
				p1 = atoms[edge[0]]
				p2 = atoms[edge[1]]
				xeg = [p1[0],p2[0]]
				yeg = [p1[1],p2[1]]
				zeg = [p1[2],p2[2]]
				egx.plot(xeg,yeg,zeg,color = 'k',linewidth = 1.5) 
		eg += 1
	for angle in range(0,90,2):
		egx.view_init(30,angle)
		filename = 'fragmento_fit_stage'+ str(angle) + basename +'.png'
		plt.savefig(filename, dpi = 96)
		plt.gca()

	sub.call(('convert','-delay','10','fragmento_fit_stage*.png',basename +'_fragmento_fit.gif'))
	args = ('rm', '-rf', 'fragmento_fit_stage*.png')
	sub.call('%s %s %s' % args, shell=True)

########################################################## b) Generator Triangles ##############################################################
if 'perspective' in options and 'dc_approach' in options:
	fig_f = plt.figure()
	fx = fig_f.gca(projection='3d')
	for polyhedron in fittest_polyhedrons:
		f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = polyhedron
		solid = Poly3DCollection(faces_list)
		solid.set_edgecolor('k')
		solid.set_facecolor('r')
		solid.set_alpha(0.15)
		fx.add_collection3d(solid)

	for generator_triangle in generator_triangles:
		for plane in generator_triangle:
			fx.add_collection3d(plane)
		
	Xc = []
	Yc = []
	Zc = []
	for j in coordinates:
		Xc.append(j[0])
		Yc.append(j[1])
		Zc.append(j[2])

	fx.scatter(Xc,Yc,Zc,color= 'k',alpha = 0.15,label = 'Atoms in ' + basename + ' structure')

	for irreducible_volume in irreducible_volumes:
		for plane in irreducible_volume:
			for edge in plane:
				p1 = edge[0]
				p2 = edge[1]
				xg = [p1[0],p2[0]]
				yg = [p1[1],p2[1]]
				zg = [p1[2],p2[2]]
				fx.plot(xg,yg,zg,color = 'k',linewidth = 1.5)
	xir = []
	yir = []
	zir = []
	for list_i in critical_points:
		for i in list_i:
			xir.append(i[0])
			yir.append(i[1])
			zir.append(i[2])
	fx.scatter(xir,yir,zir,color = 'g',alpha = 1, label = 'Critical Points')
	
	plt.title('Fragmento: ' + basename)
	plt.legend()
	plt.show()

################################################################ gif ###########################################################################
if 'gif' in options and 'dc_approach' in options:
	fig_fg = plt.figure()
	fgx = fig_fg.gca(projection='3d')
	for polyhedron in fittest_polyhedrons:
		f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = polyhedron
		solid = Poly3DCollection(faces_list)
		solid.set_edgecolor('k')
		solid.set_facecolor('r')
		solid.set_alpha(0.15)
		fgx.add_collection3d(solid)

	for generator_triangle in generator_triangles:
		for plane in generator_triangle:
			fgx.add_collection3d(plane)
		
	for j in coordinates:
		fgx.scatter(j[0],j[1],j[2],color='k',alpha = 0.15)

	for irreducible_volume in irreducible_volumes:
		for plane in irreducible_volume:
			for edge in plane:
				p1 = edge[0]
				p2 = edge[1]
				xg = [p1[0],p2[0]]
				yg = [p1[1],p2[1]]
				zg = [p1[2],p2[2]]
				fgx.plot(xg,yg,zg,color = 'k',linewidth = 1.5)

	for angle in range(90,180,2):
		fgx.view_init(30,angle)
		filename = 'fragmento_gen_stage'+ str(angle) + basename +'.png'
		plt.savefig(filename, dpi = 96)
		plt.gca()

	sub.call(('convert','-delay','10','fragmento_gen_stage*.png',basename +'_fragmento_gen.gif'))
	args = ('rm', '-rf', 'fragmento_gen_stage*.png')
	sub.call('%s %s %s' % args, shell=True)

########################################################################################################################################################
############################################################### Write the DC_approach output ###############################################################
if 'dc_approach' in options:
	irreducible_volume_points = []
	for list_i in critical_points:
		for i in list_i:
			irreducible_volume_points.append(i)
	parameters = open(path_tst + "/" + parameters,'r')
	parameters = parameters.readlines()
	for line in parameters:
		symbol_par = re.match(r'critical_point_element = (.*)',line,re.M|re.I)
		if symbol_par:
			critical_point_element = symbol_par.group(1)
			break
	for line in parameters:
		fit_symbol_par = re.match(r'fittest_element = (.*)',line,re.M|re.I)
		if fit_symbol_par:
			fittest_element = fit_symbol_par.group(1)
			break

	with open(path_tst + "/" + basename + '_delta_Chem.xyz', 'w') as fp:
		len_polyh = 0
		polyhedra_vertices = []
		for fittest_polyhedron in fittest_polyhedrons:
			f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = fittest_polyhedron
			len_plh = len(f_transformed_vertices)
			len_polyh += len_plh
			polyhedra_vertices.append(f_transformed_vertices)		

		fp.write(str(len(structure) + len(irreducible_volume_points) + len_polyh) +'\n')
		fp.write('XYZ file: irreducible volume points as ' + critical_point_element + ' atoms. ' + 'Fittest polyhedra vertices shown as ' + fittest_element + ' atoms' + '\n')

		for atom_position in coordinates:
			fp.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(structure[atom_position],atom_position[0],atom_position[1],atom_position[2]))

		for coord in irreducible_volume_points:
			fp.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(str(critical_point_element),coord[0], coord[1], coord[2]))

		for fittest_polyhedron in fittest_polyhedrons:
			f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = fittest_polyhedron
			for vertex in f_transformed_vertices:
				fp.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(str(fittest_element),vertex[0], vertex[1], vertex[2]))

#######################################################################################################################################################
############################################################### Write the general output ##############################################################

if len(sys.argv) >= 4:
	output_file = open("./" + basename + '_dc_info.out','w')
	dc_title = pyfiglet.figlet_format("Delta Chem",font = "isometric1")
	output_file.write(str(dc_title) + '\n')
	output_file.write('Delta Chem by J.Castro and T.Terencio.' + '\n')
	output_file.write('Yachay Tech, Urcuqui, Ecuador.' + '\n')
	output_file.write('\n')
	output_file.write('Material: ' + str(basename) + '\t' + 'Atoms in the unit cell: ' + str(len(structure)) + '\n')
	output_file.write('Number of atomic species: ' + str(len(atomic_species_list)) + '\t' + ','.join(atomic_species_list) + '\n')
	output_file.write('Kerno' + '\n')
	c_title = pyfiglet.figlet_format('Centers of Porosities')
	output_file.write(str(c_title) + '\n')
	cx = PrettyTable(['Center of Porosity','PSD_radius'])
	c = 0
	for center in centers_of_porosities:
		cx.add_row([center,pores_radii[c]])
		c += 1

	output_file.write(str(cx) + '\n')
	print(c_title)
	print(cx)
	ct_title = pyfiglet.figlet_format('Contributing Atoms')
	output_file.write(str(ct_title) + '\n')
	pa = 0
	for pore_atoms in contributing_atoms:
		ci_title = pyfiglet.figlet_format('Atoms of Porosity # ' + str(pa + 1),font = 'digital')
		output_file.write(str(ci_title) + '\n')
		for atom in pore_atoms:
			output_file.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(structure[atom],atom[0],atom[1],atom[2]))	
		pa += 1

	if 'triangles' in options and 'kerno_only' not in options:
		output_file.write('Selected option: triangles ' + '\n')
		ch = 1
		t_title = pyfiglet.figlet_format("Triangles")
		output_file.write(str(t_title) + '\n')
		print(t_title)
		for convex_hull in convex_hulls:
			#Create the faces with the coordinates of each vertex
			X,Y,Z = convex_hulls[convex_hull]
			verts = []
			for i in range(len(X)):
				coord = [X[i], Y[i], Z[i]]
				verts.append(coord)
			triangles = (convex_hull.Export_Triangles())
			faces = []

			#Create table for simplices
			x = PrettyTable(['A','B','C'])
			x.title = "Simplices pore # " + str(ch)
			for triangle in triangles:
				simplex = [verts[triangle[0]], verts[triangle[1]], verts[triangle[2]]]
				faces.append(simplex)
				x.add_row(triangle)

			#Create table for faces
			y = PrettyTable(['A','B','C'])
			y.title= "Faces from triangulation of pore # " + str(ch)
			for face in faces:
				y.add_row(face)

			#Write the tables in the output file
			output_file.write(str(x)+'\n')
			output_file.write(str(y)+'\n')

			#Display the output in the console
			print("Simplices pore # " + str(ch))
			print(triangles)
			print(y)
			ch += 1

	if 'polygons' in options and 'kerno_only' not in options: 
		output_file.write('Selected option: polygons ' + '\n')
		p = 0
		p_title = pyfiglet.figlet_format("Polygons")
		output_file.write(str(p_title) + '\n')
		print(p_title)
		for polytope in polytopes:
			atoms = contributing_atoms[p]
			c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytope
			c_hull_v = list(c_hull_v_area.keys())
			c_hull_i = list(c_hull_i_area.keys())
			c_hull_e = list(hull_polygons.keys())
			c_hull_a = list(c_hull_i_area.values())

			xp = PrettyTable(['Simplices','Edges','Areas'])
			xp.title = 'Polygonalization of pore # ' + str(p+1)
			for i in range(len(c_hull_e)):
				xp.add_row([c_hull_i[i],c_hull_e[i],c_hull_a[i]])
			output_file.write(str(xp)+'\n')
			e_title = pyfiglet.figlet_format('Edges Porosity # '+str(p+1),font = 'bubble')
			print(e_title)
			print(c_hull_e)
			f_title = pyfiglet.figlet_format('Faces Porosity # '+str(p+1),font = 'bubble')
			print(f_title)
			print(c_hull_v)
			s_title = pyfiglet.figlet_format('Simplices Porosity # '+str(p+1),font = 'bubble')
			print(s_title)
			print(c_hull_i)
			yp = PrettyTable(['Simplices','Faces','Areas'])
			yp.title = 'Faces Polygonalization of pore # ' + str(p+1)
			c_hull_v = [np.array(v) for v in c_hull_v]
			chv = 0
			for v in c_hull_v:
				yp.add_row([c_hull_i[chv],v,c_hull_a[chv]])
				chv += 1
			output_file.write(str(yp)+'\n')
			p += 1

	if 'dc_approach' in options:
		output_file.write('Selected option: irreducible ' + '\n')
		i_title = pyfiglet.figlet_format('Irreducible Volumes')
		output_file.write(str(i_title) + '\n')
		print(i_title)
		i = 0
		for ad_pore in critical_points:
			c_title = pyfiglet.figlet_format('Critical Points Porosity # '+str(i+1),font = "digital")
			output_file.write(str(c_title) + '\n')
			print(c_title)
			for coord in ad_pore:
				output_file.write('{0:12.6f} {1:12.6f} {2:12.6f}\n'.format(coord[0], coord[1], coord[2]))
				print('{0:12.6f} {1:12.6f} {2:12.6f}\n'.format(coord[0], coord[1], coord[2]))
			i += 1
		iv = 0	
		for irreducible_volume in irreducible_volumes:
			iv_title = pyfiglet.figlet_format('Irreducible Volume Planes Porosity # ' + str(iv + 1),font = 'bubble')
			output_file.write(str(iv_title) + '\n')
			xv = PrettyTable(['Irreducible_Volume_Planes'])
			xv.hrules = ALL
			irreducible_volume_planes = []
			for plane in irreducible_volume:
				points_in_plane = {}
				for edge in plane:
					p1 = tuple(edge[0])
					p2 = tuple(edge[1])
					points_in_plane[p1] = None
					points_in_plane[p2] = None
				points_in_plane = list(points_in_plane.keys())
				irreducible_volume_planes.append(np.array(points_in_plane))

			for plane in irreducible_volume_planes:
				xv.add_row([plane])
			output_file.write(str(xv)+'\n')
			iv += 1
	
	if 'porosity_par' in options and 'kerno_only' not in options:
		output_file.write('Selected option: porosity_par ' + '\n')
		po_p_title = pyfiglet.figlet_format('Geometric Analysis of Porosities')
		output_file.write(str(po_p_title) + '\n')
		print(po_p_title)
		for i in range(len(porosity_volumes)):
			table = PrettyTable()
			pp_title = pyfiglet.figlet_format('Analysis Porosity # ' + str(i + 1),font = 'digital')
			output_file.write(str(pp_title) + '\n')
			print(pp_title)
			table.field_names = ['Volume_Porosity # ' + str(i + 1), 'Area_int_surf_Porosity # ' + str(i + 1),'# of faces Porosity #' + str(i + 1)]
			table.add_row([porosity_volumes[i], porosity_areas[i], porosity_num_faces[i]])
			output_file.write(str(table) + '\n')
			print(table)
			table_2 = PrettyTable()
			table_2.field_names = ['Simplices','Uniformity Ratio','Clasification','Area']
			c_hull_v_area,c_hull_i_area,num_vertices_polytope,sugestions,num_faces,areas,max_polygon,min_reg_polygon,polytope_type,hull_polygons,hull_poly_ratio = polytopes[i]
			c_hull_i = list(c_hull_i_area.keys())
			c_hull_r = list(hull_poly_ratio.values())
			c_hull_a = list(c_hull_i_area.values())
			c_hull_e = list(hull_polygons.keys())
			c_hull_c = list(hull_polygons.values())
			for j in range(len(c_hull_i)):
				table_2.add_row([c_hull_i[j],c_hull_r[j],c_hull_c[j],c_hull_a[j]])
			
			output_file.write(str(table_2) + '\n')
			print(table_2)

	if 'fittest_par' in options and 'dc_approach' in options:
		output_file.write('Selected option: porosity_par ' + '\n')
		py_p_title = pyfiglet.figlet_format('Geometric Analysis of Fittest Polyhedrons')
		output_file.write(str(py_p_title) + '\n')
		print(py_p_title)
		fc = 0
		for fittest_polyhedron in fittest_polyhedrons:
			f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = fittest_polyhedron
			f_poly_title = pyfiglet.figlet_format('Fittest Polyhedron # ' + str(fc + 1) + ': ' + str(fit_polyhedron), font = 'digital')
			fp_table = PrettyTable(['Volume','Surface Area','# Faces','Biggest Face Area','Smallest Face Area',])
			fp_table.add_row([Volume_of_polyhedron,Area_of_polyhedron,num_faces,big_face_area,small_face_area])
			output_file.write(str(f_poly_title) + '\n')
			output_file.write(str(fp_table) + '\n')
			print(f_poly_title)
			print(fp_table)
			ver_title = pyfiglet.figlet_format('Vertices',font = 'digital')
			output_file.write(str(ver_title) + '\n')
			print(ver_title)
			vc = 1
			for vertex in f_transformed_vertices:
				output_file.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(str(vc),vertex[0],vertex[1],vertex[2]))
				print('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(str(vc),vertex[0],vertex[1],vertex[2]))
				vc += 1
			fc += 1
			ff_table = PrettyTable(['Face','Area'])
			for face in faces_polyhedron:
				f_area = faces_areas[tuple(face)]
				ff_table.add_row([face,f_area])
			output_file.write(str(ff_table) + '\n')
			print(ff_table)

	if 'Wythoff' in options and 'dc_approach' in options:
		output_file.write('Selected option: Wythoff' + '\n')
		py_p_title = pyfiglet.figlet_format('All the critical points of the Fittest Polyhedrons')
		output_file.write(str(py_p_title) + '\n')
		count_p = 0
		all_points = []
		all_gens = []
		for fittest_polyhedron in fittest_polyhedrons:
			cen = centers_of_porosities[count_p]
			f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = fittest_polyhedron
			triangles,planes_gen,generators,all_sites,points_gens = Wythoff(fit_polyhedron,radius_d,[cen[0],cen[1],cen[2]])
			Riv = np.dot(Rf,R)
			aligned_generators = []
			for generator in generators:
				generator = translation_transformation(generator,o_direction)
				align_g = np.dot(Riv,np.transpose(generator))
				aligned_point = []
				aligned_point.append(align_g.item(0))
				aligned_point.append(align_g.item(1))
				aligned_point.append(align_g.item(2))
				generator = aligned_point
				generator = translation_transformation(generator,o_direction*(-1))
				aligned_generators.append(generator)
			all_gens.append(aligned_generators)	

			aligned_sites = {}
			for site in all_sites:
				site = translation_transformation(site,o_direction)
				align_s = np.dot(Riv,np.transpose(site))
				aligned_point = []
				aligned_point.append(align_s.item(0))
				aligned_point.append(align_s.item(1))
				aligned_point.append(align_s.item(2))
				site = aligned_point
				site = translation_transformation(site,o_direction*(-1))
				aligned_sites[tuple(site)] = None
			aligned_sites = list(aligned_sites.keys())
			all_points.append(aligned_sites)
			count_p += 1
			fc = 0
			for pore_points in all_points:
				all_title = pyfiglet.figlet_format('Critical Points Porosity # ' + str(fc + 1), font = 'digital')
				output_file.write(str(all_title) + '\n')
				for site in pore_points:
					output_file.write('{0:6.2s} {1:12.6f} {2:12.6f} {3:12.6f}\n'.format(str(critical_point_element),site[0],site[1],site[2]))
				fc += 1

		if 'perspective' in options:
			fig_d = plt.figure()
			dx = fig_d.gca(projection='3d')
			colorsd = ['m','g','b','y','r','b']
			d = 0
			for pore_points in all_points:
				wx = []
				wy = []
				wz = []
				for site in pore_points:
					wx.append(site[0])
					wy.append(site[1])
					wz.append(site[2])
				dx.scatter(wx,wy,wz,color = colorsd[d], alpha = 0.5, label = 'Wythoff vertices Polyhedron # ' + str(d+1))
				d += 1
			
			for cen in centers_of_porosities:
				dx.scatter(cen[0],cen[1],cen[2],color = 'r', alpha = 1)
			for polyhedron in fittest_polyhedrons:
				f_transformed_vertices,faces_list,R,Rf,radius_d,scale_factor,o_direction,faces_polyhedron,num_faces,Volume_of_polyhedron,Area_of_polyhedron,big_face_area,small_face_area,fit_polyhedron,faces_areas = polyhedron
				solid = Poly3DCollection(faces_list)
				solid.set_edgecolor('k')
				solid.set_facecolor('r')
				solid.set_alpha(0.0)
				dx.add_collection3d(solid)
			
			cwx = []
			cwy = []
			cwz = []
			for j in coordinates:
				cwx.append(j[0])
				cwy.append(j[1])
				cwz.append(j[2])
			dx.scatter(cwx,cwy,cwz,color='k',alpha = 0.25)
			plt.legend()
			plt.show()
	output_file.close()

