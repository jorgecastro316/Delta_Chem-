import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import *
import sympy
from deltachem.general_rotation import General_Rotation
from deltachem.DC_library import *
from deltachem.DC_polyhedrons import *

def Irreducible_Volume(list_of_parameters,vertices):
	p = list_of_parameters[0]
	q = list_of_parameters[1]
	r = list_of_parameters[2]
	ap = float(sympy.pi/p)
	aq = float(sympy.pi/q)
	ar = float(sympy.pi/r)
	triangle_type = list_of_parameters[3]
	radius = list_of_parameters[-2]
	center = list_of_parameters[-1]
	#Order the vertices of the triangle according with the angle at that point 

	Reference_point = np.array(vertices[0])
	E1 = np.array(vertices[1])
	E2 = np.array(vertices[2])

	if triangle_type == 'c':
		corner = list_of_parameters[4]
		generator,points,triangle_type,perimeter_points,planes = Elemental_Triangle_Cornered(p,q,r,corner,radius,center,Reference_point,E1,E2)

	if triangle_type == 's':
		side = list_of_parameters[4]
		generator,points,triangle_type,perimeter_points,planes = Elemental_Triangle_Sided(p,q,r,side,radius,center,Reference_point,E1,E2)

	if triangle_type == 'i':
		generator,points,triangle_type,perimeter_points,planes = Elemental_Triangle_Incentered(p,q,r,radius,center,Reference_point,E1,E2)


	return(generator,points,triangle_type,perimeter_points,planes)
if __name__ == "__Irreducible_Volume__":
    Irreducible_Volume()

def Elemental_Triangle_Cornered(p,q,r,corner,radius,center,Reference_point,E1,E2):
	#The three numbers in Wythoff's symbol, p, q and r, 
	#represent the corners of the Schwarz triangle which angles are
	rad = radius
	radius = 1
	center_of_sphere = np.array([0,0,0])
	tp = sympy.pi / p
	tq = sympy.pi / q
	tr = sympy.pi / r
	cp = float(sympy.cos(tp))
	sp = float(sympy.sin(tp))
	cq = float(sympy.cos(tq))
	cr = float(sympy.cos(tr))

	if p == 2:
		tref = tp
		te1 = tq
		te2 = tr
		if corner == 0:
			generator = Reference_point
		elif corner == 1:
			generator = E1
		elif corner == 2:
			generator = E2
		
	elif q == 2:

		tref = tq
		te1 = tr
		te2 = tp
		if corner == 0:
			generator = E2
		elif corner == 1:
			generator = Reference_point
		elif corner == 2:
			generator = E1

	elif r == 2:
		tref = tr
		te1 = tp
		te2 = tq
		if corner == 0:
			generator = E1
		elif corner == 1:
			generator = E2
		elif corner == 2:
			generator = Reference_point

	#Angles of the arcs forming the edges of the triangle
	ang_arcRef_E1 = acos(np.dot(Reference_point,E1)/(np.linalg.norm(Reference_point)*np.linalg.norm(E1)))
	ang_arcRef_E2 = acos(np.dot(Reference_point,E2)/(np.linalg.norm(Reference_point)*np.linalg.norm(E2)))
	ang_arcE1_E2 = acos(np.dot(E1,E2)/(np.linalg.norm(E1)*np.linalg.norm(E2)))

	if generator[0] == Reference_point[0] and generator[1] == Reference_point[1] and generator[2] == Reference_point[2] :
		triangle_type = 'c1'
		#In this case the generator must be projected orthogonal to the opposite side
		arc_Ref_E1 = radius*ang_arcRef_E1
		arc_proj = math.atan(math.tan(arc_Ref_E1)*math.cos(float(te1)))
		ang_proj = arc_proj/radius
		projection = General_Rotation(center_of_sphere,np.cross(E1,E2),ang_proj,E1)

		#Duplicate the triangle:
		Refd = Reference_point
		E1d = E1
		E2d = General_Rotation(center_of_sphere,np.cross(Reference_point,E2),-ang_arcRef_E2,Reference_point)
		projectiond = General_Rotation(center_of_sphere,np.cross(E1d,E2d),ang_proj,E1d)

		#steps
		ang_g_proj = acos(np.dot(generator,projection)/(np.linalg.norm(generator)*np.linalg.norm(projection)))
		step_proj = General_Rotation(center_of_sphere,np.cross(generator,projection),2*ang_g_proj,generator)
		step_projd = General_Rotation(center_of_sphere,np.cross(generator,projectiond),2*ang_g_proj,generator)
		step_projm = General_Rotation(center_of_sphere,np.cross(generator,projectiond),-2*ang_g_proj,generator)

		#Projections on the plane
		proj = np.array([((step_proj[0]-generator[0])/2)+generator[0],((step_proj[1]-generator[1])/2)+generator[1],((step_proj[2]-generator[2])/2)+generator[2]])
		projd = np.array([((step_projd[0]-generator[0])/2)+generator[0],((step_projd[1]-generator[1])/2)+generator[1],((step_projd[2]-generator[2])/2)+generator[2]])
		projm = np.array([((step_projm[0]-generator[0])/2)+generator[0],((step_projm[1]-generator[1])/2)+generator[1],((step_projm[2]-generator[2])/2)+generator[2]])

		#Define the vertices of partitional faces

		p_g = np.subtract(proj,generator)
		pd_g = np.subtract(projd,generator)
		pm_g = np.subtract(projm,generator)
	
		normal_plane_E1 = np.cross(p_g,pd_g)
		normal_plane_E1 = normal_plane_E1 * 1/np.linalg.norm(normal_plane_E1) 
		normal_plane_E2 = np.cross(pm_g,p_g)
		normal_plane_E2 = normal_plane_E2 * 1/np.linalg.norm(normal_plane_E2)
	
		E1_proj = E1 - np.dot((E1 - generator),normal_plane_E1)*normal_plane_E1
		E2_proj = E2 - np.dot((E2 - generator),normal_plane_E2)*normal_plane_E2

		points = [generator,proj,E1_proj,E2_proj]
		#points = Scaling_Traslation_Point(points,center,rad)
		generator = points[0]
		perimeter_points = points
		points = points[1:]

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
		plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='k',alpha=0.5)
		plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='m',alpha=0.5)
		planes = [plane_E1,plane_E2]


	else:
		triangle_type = 'c2'
		if generator[0] == E1[0] and generator[1] == E1[1] and generator[2] == E1[2]:
			mobile_point = E2
		if generator[0] == E2[0] and generator[1] == E2[1] and generator[2] == E2[2]:
			mobile_point = E1
		#Duplicate:
		#Look for the larger edge
		lenRef_mob = np.linalg.norm(Reference_point - mobile_point)
		lenRef_gen = np.linalg.norm(Reference_point - generator)

		if lenRef_gen > lenRef_mob:
			rot = generator
			stat = mobile_point
		elif lenRef_gen < lenRef_mob:
			rot = mobile_point
			stat = generator
		else: 
			rot = mobile_point
			stat = generator

		#1)Perform a step of arc R_rot in the direction rot X Ref 
		ang_R_rot =  acos(np.dot(Reference_point,rot)/(np.linalg.norm(Reference_point)*np.linalg.norm(rot)))
		generatorp = General_Rotation(center_of_sphere,np.cross(Reference_point,rot),-2*ang_R_rot,rot)
		#2)Transport the rotative (rot) point over an arc Rstat on the direction generatorp X stat
		ang_R_stat =  acos(np.dot(Reference_point,stat)/(np.linalg.norm(Reference_point)*np.linalg.norm(stat)))
		Refd = General_Rotation(center_of_sphere,np.cross(generatorp,stat),ang_R_stat,stat)

		#steps
		ang_R_g =  acos(np.dot(Reference_point,generator)/(np.linalg.norm(Reference_point)*np.linalg.norm(generator)))
		step_Ref = General_Rotation(center_of_sphere,np.cross(generator,Reference_point),2*ang_R_g,generator)
		step_Refd = General_Rotation(center_of_sphere,np.cross(generator,Refd),2*ang_R_g,generator)


		#Projections on the plane
		projR = np.array([((step_Ref[0] - generator[0])/2) + generator[0],((step_Ref[1]-generator[1])/2) + generator[1],((step_Ref[2] - generator[2])/2) + generator[2]])
		projRd = np.array([((step_Refd[0] - generator[0])/2) + generator[0],((step_Refd[1]-generator[1])/2) + generator[1],((step_Refd[2] - generator[2])/2) + generator[2]])

		#Define the vertices of faces

		pr_g = np.subtract(projR,generator)
		prd_g = np.subtract(projRd,generator)

		normal_plane = np.cross(pr_g,prd_g)
		normal_plane = normal_plane * 1/np.linalg.norm(normal_plane) 

		M_proj = mobile_point - np.dot((mobile_point - generator),normal_plane)*normal_plane

		points = [generator,projR,M_proj]
		#points = Scaling_Traslation_Point(points,center,rad)
		generator = points[0]
		perimeter_points = points[1:]
		x_plane = np.array([generator[0],projR[0],M_proj[0]])
		y_plane = np.array([generator[1],projR[1],M_proj[1]])
		z_plane = np.array([generator[2],projR[2],M_proj[2]])
		verts_plane = []
		verts_plane.append(list(zip(x_plane,y_plane,z_plane)))
		plane = Poly3DCollection(verts_plane,facecolor='m',alpha=0.5,linewidth = 1.5)
		planes = [plane]
	#plt.show()
	return(generator,points,triangle_type,perimeter_points,planes)
if __name__ == "__Elemental_Triangle_Cornered__":
    Elemental_Triangle_Cornered()

def Elemental_Triangle_Incentered(p,q,r,radius,center,Reference_point,E1,E2):
	#The three numbers in Wythoff's symbol, p, q and r, 
	#represent the corners of the Schwarz triangle which angles are
	# pi/p, pi/q and pi/r radians respectively.
	triangle_type = 'i'
	rad = radius
	radius = 1
	center_of_sphere = np.array([0,0,0]) 
	tp = sympy.pi / p
	tq = sympy.pi / q
	tr = sympy.pi / r
	cp = float(sympy.cos(tp))
	sp = float(sympy.sin(tp))
	cq = float(sympy.cos(tq))
	cr = float(sympy.cos(tr))

	#here we go again with another weird theory 	
	#Calculate the great circles angles
	arcA_B = math.acos(np.dot(Reference_point,E2)/(np.linalg.norm(Reference_point)*np.linalg.norm(E2)))
	arcA_C = math.acos(np.dot(Reference_point,E1)/(np.linalg.norm(Reference_point)*np.linalg.norm(E1)))
	arcB_C = math.acos(np.dot(E2,E1)/(np.linalg.norm(E2)*np.linalg.norm(E1)))

	#Calculate the lenght of the arcs
	a = radius*arcB_C
	b = radius*arcA_C
	c = radius*arcA_B
	# s is the half sum of the arcs
	s = (a + b + c)/2
	ss = math.sin(s)
	ssa = math.sin(s - a)
	ssb = math.sin(s - b)
	ssc = math.sin(s - c)
	n = math.sqrt(ss*ssa*ssb*ssc)
	#compute the angular radius of the circle inscribed on the spherical triangle
	ang_rad = math.atan2(n,ss)

	#Calculate the difference between the lenght of the arc of the side opposite to the angle pi/2
	#and the angular radius. 
	arc_left = arcA_B - ang_rad
	#Get the angles 
	rotation_angle = -ang_rad/radius
	angle_left = arc_left/radius

	if p == 2:
		p_axis = np.cross(E2,E1)
		q_axis = np.cross(Reference_point,E2)
		r_axis = np.cross(E1,Reference_point)
		i_axis = np.cross(Reference_point,E2)
		p_angle = angle_left
		q_angle = -rotation_angle
		r_angle = rotation_angle
		i_angle = -rotation_angle
		p_origin = E2
		q_origin = Reference_point
		r_origin = Reference_point


	if q == 2: 
		p_axis = np.cross(E1,Reference_point)
		q_axis = np.cross(E2,E1)
		r_axis = np.cross(Reference_point,E2)
		i_axis = np.cross(E1,Reference_point)
		p_angle = rotation_angle
		q_angle = angle_left
		r_angle = -rotation_angle
		i_angle = rotation_angle
		q_origin = E2
		p_origin = Reference_point
		r_origin = Reference_point

	if  r == 2:
		p_axis = np.cross(Reference_point,E2)
		q_axis = np.cross(E1,Reference_point)
		r_axis = np.cross(E2,E1)
		i_axis = np.cross(Reference_point,E2)
		p_angle = -rotation_angle
		q_angle = rotation_angle
		r_angle = angle_left
		i_angle = -rotation_angle
		p_origin = Reference_point
		q_origin = Reference_point
		r_origin = E2


	#To obtain the incenter and the orthogonal projections over the edges of the triangle
	#the position vectors of the vertices will be rotated on angle proportional to the lenght
	#of the angular radius. 

	#proyection over the edge r

	Proj_on_r = General_Rotation(center_of_sphere,r_axis,r_angle,r_origin)

	#proyection over the edge q
	
	Proj_on_q = General_Rotation(center_of_sphere,q_axis,q_angle,q_origin)

	#proyection over the edge p

	Proj_on_p = General_Rotation(center_of_sphere,p_axis,p_angle,p_origin)

	#To obtain the incenter by rotating one of the proyections over an arc 
	#proportional to the angular radius

	if r == 2:
		Proj_to_i = Proj_on_q
	else:
		Proj_to_i = Proj_on_r

	spherical_incenter = General_Rotation(center_of_sphere,i_axis,i_angle,Proj_to_i)

	#########################Get the faces of the irreducible volume polytope#######################
	#One step at time
	#Use the edges of the triangle as reflection planes to move the generator vertex along the 
	#arcs joining it with its projections over the edges.

	if p == 2:
		P = Reference_point
		Q = E1
		R = E2
		step_p_angle = -2*rotation_angle 
		step_q_angle = -2*rotation_angle
		step_r_angle = 2*rotation_angle
		ref_step_p = spherical_incenter
		ref_step_q = spherical_incenter
		ref_step_r = P
		proj_step_p = Proj_on_p
		proj_step_q = Proj_on_q
		proj_step_r = Proj_on_q
	
	elif q == 2:
		Q = Reference_point
		R = E1
		P = E2
		step_p_angle = 2*rotation_angle 
		step_q_angle = -2*rotation_angle
		step_r_angle = 2*rotation_angle
		ref_step_p = Q
		ref_step_q = spherical_incenter
		ref_step_r = Q
		proj_step_p = Proj_on_r
		proj_step_q = Proj_on_q
		proj_step_r = Proj_on_p

	
	elif r == 2:
		R = Reference_point
		P = E1
		Q = E2
		step_p_angle = 2*rotation_angle 
		step_q_angle = 2*rotation_angle
		step_r_angle = -2*rotation_angle
		ref_step_p = R
		ref_step_q = R
		ref_step_r = spherical_incenter
		proj_step_p = Proj_on_q
		proj_step_q = Proj_on_p
		proj_step_r = Proj_on_r


	step_p = General_Rotation(center_of_sphere,np.cross(ref_step_p,proj_step_p),step_p_angle,spherical_incenter)
	step_q = General_Rotation(center_of_sphere,np.cross(ref_step_q,proj_step_q),step_q_angle,spherical_incenter)
	step_r = General_Rotation(center_of_sphere,np.cross(ref_step_r,proj_step_r),step_r_angle,spherical_incenter)

	#Define thre planes using the three step points and the incenter
	#Obtain the normal vectors for the planes
	p_i = np.subtract(step_p,spherical_incenter)
	q_i = np.subtract(step_q,spherical_incenter)
	r_i = np.subtract(step_r,spherical_incenter)

	normal_plane_P = np.cross(r_i,q_i)
	normal_plane_P = normal_plane_P * 1/np.linalg.norm(normal_plane_P) 
	normal_plane_Q = np.cross(r_i,p_i)
	normal_plane_Q = normal_plane_Q * 1/np.linalg.norm(normal_plane_Q)
	normal_plane_R = np.cross(q_i,p_i) 
	normal_plane_R = normal_plane_R * 1/np.linalg.norm(normal_plane_R)


	#Project the vertices of the triangle P,Q,R 
	#over the surface of the plane along the direction of the normal vector

	P_proj = P - np.dot((P - spherical_incenter),normal_plane_P)*normal_plane_P
	Q_proj = Q - np.dot((Q - spherical_incenter),normal_plane_Q)*normal_plane_Q
	R_proj = R - np.dot((R - spherical_incenter),normal_plane_R)*normal_plane_R
	p_proj = np.array([((step_p[0]-spherical_incenter[0])/2)+spherical_incenter[0],((step_p[1]-spherical_incenter[1])/2)+spherical_incenter[1],((step_p[2]-spherical_incenter[2])/2)+spherical_incenter[2]])
	q_proj = np.array([((step_q[0]-spherical_incenter[0])/2)+spherical_incenter[0],((step_q[1]-spherical_incenter[1])/2)+spherical_incenter[1],((step_q[2]-spherical_incenter[2])/2)+spherical_incenter[2]])
	r_proj = np.array([((step_r[0]-spherical_incenter[0])/2)+spherical_incenter[0],((step_r[1]-spherical_incenter[1])/2)+spherical_incenter[1],((step_r[2]-spherical_incenter[2])/2)+spherical_incenter[2]])

	points = [spherical_incenter,P_proj,Q_proj,R_proj,p_proj,q_proj,r_proj]
	#points = Scaling_Traslation_Point(points,center,rad)
	generator = np.array(points[0])
	perimeter_points = points[1:]
#####################################################Graphics#########################################3####33333
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

	plane_P = Poly3DCollection(verts_plane_P,facecolor='k',alpha=0.5)
	plane_Q = Poly3DCollection(verts_plane_Q,facecolor='b',alpha=0.5)
	plane_R = Poly3DCollection(verts_plane_R,facecolor='m',alpha=0.5)
	planes = [plane_P,plane_Q,plane_R]

	return(generator,points,triangle_type,perimeter_points,planes)

if __name__ == "__Elemental_Triangle_Incentered__":
    Elemental_Triangle_Incentered()

def Elemental_Triangle_Sided(p,q,r,side,radius,center,Reference_point,E1,E2):
	#The three numbers in Wythoff's symbol, p, q and r, 
	#represent the corners of the Schwarz triangle which angles are
	# pi/p, pi/q and pi/r radians respectively.
	rad = radius
	radius = 1
	center_of_sphere = np.array([0,0,0])
	tp = sympy.pi / p
	tq = sympy.pi / q
	tr = sympy.pi / r
	cp = float(sympy.cos(tp))
	sp = float(sympy.sin(tp))
	cq = float(sympy.cos(tq))
	cr = float(sympy.cos(tr))
	if p == 2:
		reference = 0
		tref = tp
		t1 = tq
		t2 = tr

	if q == 2:
		reference = 1
		tref = tq
		t1 = tr
		t2 = tp

	if r == 2:
		reference = 2
		tref = tr
		t1 = tp 
		t2 = tq
	#Angles of the arcs forming the edges pf the triangle
	ang_arcRef_E1 = acos(np.dot(Reference_point,E1)/(np.linalg.norm(Reference_point)*np.linalg.norm(E1)))
	ang_arcRef_E2 = acos(np.dot(Reference_point,E2)/(np.linalg.norm(Reference_point)*np.linalg.norm(E2)))
	ang_arcE1_E2 = acos(np.dot(E1,E2)/(np.linalg.norm(E1)*np.linalg.norm(E2)))
	
	if reference == side:
		triangle_type = 's1'
		c = radius*ang_arcRef_E1
		A = tref/2
		B = t1

		#Napier's analogies to obtain the generator vertex

		tnc2 = float(sympy.tan(c/2))

		ApB = ((A + B)/2)
		AmB = ((A - B)/2)
		tnapb = (math.cos(AmB)/math.cos(ApB))*tnc2
		tnamb = (math.sin(AmB)/math.sin(ApB))*tnc2

		apb = math.atan(tnapb)*2
		amb = math.atan(tnamb)*2

		a = (apb + amb)/2
		b = apb - a

		ang_arc_a =  math.radians(math.degrees(a/radius))
		generator = General_Rotation(center_of_sphere,np.cross(E1,E2),ang_arc_a,E1)
		#Duplicate the triangle
		E1d = E1
		E2d = E2
		ang_arcRef_g = math.acos(np.dot(Reference_point,generator)/(np.linalg.norm(Reference_point)*np.linalg.norm(generator)))
		Refd = General_Rotation(center_of_sphere,np.cross(Reference_point,generator),ang_arcRef_g,generator)

		#Obtain the projections of the generator over the remaining edges
		proj_r_C = A
		arc_proj_e2 = math.atan(math.tan(b)*math.cos(proj_r_C))
		arc_proj_e1 = math.atan(math.tan(b)*math.cos(proj_r_C))
		ang_proj_e2 = arc_proj_e2/radius
		ang_proj_e1 = arc_proj_e1/radius
		proj_on_e2 = General_Rotation(center_of_sphere,np.cross(Reference_point,E1),ang_proj_e2,Reference_point)
		proj_on_e1 = General_Rotation(center_of_sphere,np.cross(Reference_point,E2),ang_proj_e1,Reference_point)
		proj_on_e1d = General_Rotation(center_of_sphere,np.cross(Refd,E1d),ang_proj_e2,Refd)
		proj_on_e2d = General_Rotation(center_of_sphere,np.cross(Refd,E2d),ang_proj_e1,Refd)

		#Steps
		ang_step_e2 = math.acos(np.dot(generator,proj_on_e2)/(np.linalg.norm(generator)*np.linalg.norm(proj_on_e2)))
		step_e2 = General_Rotation(center_of_sphere,np.cross(generator,proj_on_e2),ang_step_e2,proj_on_e2)
		ang_step_e1 = math.acos(np.dot(generator,proj_on_e1)/(np.linalg.norm(generator)*np.linalg.norm(proj_on_e1)))
		step_e1 = General_Rotation(center_of_sphere,np.cross(generator,proj_on_e1),ang_step_e1,proj_on_e1)
		step_e2d = General_Rotation(center_of_sphere,np.cross(generator,proj_on_e2d),ang_step_e2,proj_on_e2d)
		step_e1d = General_Rotation(center_of_sphere,np.cross(generator,proj_on_e1d),ang_step_e1,proj_on_e1d)

		#Projections on the plane

		e1_proj = np.array([((step_e1[0]-generator[0])/2)+generator[0],((step_e1[1]-generator[1])/2)+generator[1],((step_e1[2]-generator[2])/2)+generator[2]])
		e2_proj = np.array([((step_e2[0]-generator[0])/2)+generator[0],((step_e2[1]-generator[1])/2)+generator[1],((step_e2[2]-generator[2])/2)+generator[2]])
		e1d_proj = np.array([((step_e1d[0]-generator[0])/2)+generator[0],((step_e1d[1]-generator[1])/2)+generator[1],((step_e1d[2]-generator[2])/2)+generator[2]])
		e2d_proj = np.array([((step_e2d[0]-generator[0])/2)+generator[0],((step_e2d[1]-generator[1])/2)+generator[1],((step_e2d[2]-generator[2])/2)+generator[2]])

		#Built the faces

		e1_g = np.subtract(e1_proj,generator)
		e2_g = np.subtract(e2_proj,generator)
		e1d_g = np.subtract(e1d_proj,generator)
		e2d_g = np.subtract(e2d_proj,generator)
		normal_plane_Ref = np.cross(e2_g,e1_g)
		normal_plane_Ref = normal_plane_Ref * 1/np.linalg.norm(normal_plane_Ref) 
		normal_plane_E1 = np.cross(e2_g,e1d_g)
		normal_plane_E1 = normal_plane_E1 * 1/np.linalg.norm(normal_plane_E1)
		normal_plane_E2 = np.cross(e1_g,e2d_g) 
		normal_plane_E2 = normal_plane_E2 * 1/np.linalg.norm(normal_plane_E2)
		Ref_proj = Reference_point - np.dot((Reference_point - generator),normal_plane_Ref)*normal_plane_Ref
		E1_proj = E1 - np.dot((E1 - generator),normal_plane_E1)*normal_plane_E1
		E2_proj = E2 - np.dot((E2 - generator),normal_plane_E2)*normal_plane_E2

		#Scale all the relevant points of the triangle
		points = [generator,Ref_proj,E1_proj,E2_proj,e1_proj,e2_proj]
		#đpoints = Scaling_Traslation_Point(points,center,rad)
		generator = points[0]
		perimeter_points = points
		points = points[1:]

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
		plane_Ref = Poly3DCollection(verts_plane_Ref,facecolor='k',alpha=0.5)
		plane_E1 = Poly3DCollection(verts_plane_E1,facecolor='b',alpha=0.5)
		plane_E2 = Poly3DCollection(verts_plane_E2,facecolor='m',alpha=0.5)
		planes = [plane_Ref,plane_E1,plane_E2]

	else:
		triangle_type = 's2'
		if side == 0:
			if q == 2:
				second_ref = E1
				opposed_ref = E2
				sec_ang = tr
				sign = 1
			if r == 2:
				second_ref = E2
				opposed_ref = E1
				sec_ang = tq
				sign = -1
			op_ang = tp
			A = op_ang/2

		if side == 1:
			if p == 2:
				second_ref = E2
				sec_ang = tr
				opposed_ref = E1
				sign = -1
			if r == 2:
				second_ref = E1
				opposed_ref = E2
				sec_ang = tp
				sign = 1
			op_ang = tq
			A = op_ang/2

		if side == 2:
			if p == 2:
				second_ref = E1
				opposed_ref = E2
				sec_ang = tq
				sign = 1 
			if q == 2:
				second_ref = E2
				opposed_ref = E1
				sec_ang = tp
				sign = -1
			op_ang = tr
			A = op_ang/2


		#Napier's analogies to obtain the generator vertex
		ang_arcRef_Op = math.acos(np.dot(Reference_point,opposed_ref)/(np.linalg.norm(Reference_point)*np.linalg.norm(opposed_ref)))

		#generator = General_Rotation(center_of_sphere,np.cross(Reference_point,second_ref),ang_arc_a,Reference_point)
		c = radius*ang_arcRef_Op
		B = tref
		tnc2 = float(sympy.tan(c/2))
		ApB = ((A + B)/2)
		AmB = ((A - B)/2)
		tnapb = (math.cos(AmB)/math.cos(ApB))*tnc2
		tnamb = (math.sin(AmB)/math.sin(ApB))*tnc2
		apb = math.atan(tnapb)*2
		amb = math.atan(tnamb)*2
		a = (apb + amb)/2
		b = apb - a
		ang_arc_a = math.radians(math.degrees(a/radius))
		ang_arc_b = math.radians(math.degrees(b/radius))
		generator = General_Rotation(center_of_sphere,np.cross(Reference_point,second_ref),ang_arc_a,Reference_point)

		#Duplicate the triangle
		E1d = second_ref
		Refd = Reference_point
		ang_arcRef_Op = math.acos(np.dot(Reference_point,opposed_ref)/(np.linalg.norm(Reference_point)*np.linalg.norm(opposed_ref)))
		E2d = General_Rotation(center_of_sphere,np.cross(Reference_point,opposed_ref),-ang_arcRef_Op,Reference_point)

		#Obtain the projections of the generator over the remaining edges
		proj_r_C = sec_ang
		ang_arc_Sec_g = math.acos(np.dot(second_ref,generator)/(np.linalg.norm(second_ref)*np.linalg.norm(generator)))
		arc_Sec_g = ang_arc_Sec_g/radius
		arc_proj_ref = math.atan(math.tan(arc_Sec_g)*math.cos(proj_r_C))
		ang_proj_ref = sign*arc_proj_ref/radius
		proj_on_ref = General_Rotation(center_of_sphere,np.cross(E1,E2),ang_proj_ref,second_ref)
		proj_on_refd = General_Rotation(center_of_sphere,np.cross(E1d,E2d),sign*ang_proj_ref,second_ref)

		#Steps
		ang_step_proj = math.acos(np.dot(generator,proj_on_ref)/(np.linalg.norm(generator)*np.linalg.norm(proj_on_ref)))
		step_proj = General_Rotation(center_of_sphere,np.cross(generator,proj_on_ref),2*ang_step_proj,generator)
		ang_step_projd = math.acos(np.dot(generator,proj_on_refd)/(np.linalg.norm(generator)*np.linalg.norm(proj_on_refd)))
		step_projd = General_Rotation(center_of_sphere,np.cross(generator,proj_on_refd),2*ang_step_projd,generator)
		ang_stepRef = math.acos(np.dot(Reference_point,generator)/(np.linalg.norm(Reference_point)*np.linalg.norm(generator)))
		step_Ref = General_Rotation(center_of_sphere,np.cross(generator,Reference_point),ang_stepRef,Reference_point)

		#Projections on the plane

		proj = np.array([((step_proj[0]-generator[0])/2)+generator[0],((step_proj[1]-generator[1])/2)+generator[1],((step_proj[2]-generator[2])/2)+generator[2]])
		projd = np.array([((step_projd[0]-generator[0])/2)+generator[0],((step_projd[1]-generator[1])/2)+generator[1],((step_projd[2]-generator[2])/2)+generator[2]])
		projRef = np.array([((step_Ref[0]-generator[0])/2)+generator[0],((step_Ref[1]-generator[1])/2)+generator[1],((step_Ref[2]-generator[2])/2)+generator[2]])

		#Define the vertices of partitional faces

		p_g = np.subtract(proj,generator)
		pd_g = np.subtract(projd,generator)
		pref_g = np.subtract(projRef,generator)
	
		normal_plane_Op = np.cross(p_g,pref_g)
		normal_plane_Op = normal_plane_Op * 1/np.linalg.norm(normal_plane_Op) 
		normal_plane_sec = np.cross(p_g,pd_g)
		normal_plane_sec = normal_plane_sec * 1/np.linalg.norm(normal_plane_sec)
	
		Op_proj = opposed_ref - np.dot((opposed_ref - generator),normal_plane_Op)*normal_plane_Op
		Sec_proj = second_ref - np.dot((second_ref - generator),normal_plane_sec)*normal_plane_sec

		#Scale all the relevant points of the triangle
		points = [generator,proj,Op_proj,Sec_proj,projRef]
		#points = Scaling_Traslation_Point(points,center,rad)
		generator = points[0]
		perimeter_points = points[1:]
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
		plane_Op = Poly3DCollection(verts_plane_Op,facecolor='m',alpha=0.5)
		plane_Sec = Poly3DCollection(verts_plane_Sec,facecolor='k',alpha=0.5)
		planes = [plane_Op,plane_Sec]

	return(generator,points,triangle_type,perimeter_points,planes)
if __name__ == "__Elemental_Triangle_Sided__":
    Elemental_Triangle_Sided()

def Scaling(vertices,radius):
	origin = np.array([0,0,0])
	distances = []  
	for i in vertices:
		distances.append(np.linalg.norm(i - origin))
	max_distance = max(distances)

	scale_factor = radius/max_distance
	#Built the scaling matrix 4x4
	S = np.matrix([[scale_factor,0,0,0],[0,scale_factor,0,0],[0,0,scale_factor,0],[0,0,0,1]])
	new_verts = []
	for j in vertices:
		j = np.append(j,1)
		new_verts.append(np.dot(j,S))

	transformed_verts = []
	vert = []

	for k in range(0,len(new_verts)):
		#new_verts[k]
		vert.append(new_verts[k].item(0))
		vert.append(new_verts[k].item(1))
		vert.append(new_verts[k].item(2))
		transformed_verts.append(vert)
		vert = []

	return(transformed_verts)

if __name__ == "__Scaling__":
    Scaling()
#Create the initial elemental triangle
def Wythoff(Polyhedron,radius,center):
	center = np.array(center)
	elemental_triangles = {
		"Tetrahedron":[2,3,3,'c',2,radius,center],
		"Cube":[4,2,3,'c',2,radius,center],
		"Octahedron":[4,2,3,'c',0,radius,center],
		"Dodecahedron":[2,3,5,'c',1,radius,center],
		"Icosahedron":[2,3,5,'c',2,radius,center],
		"Truncated_Tetrahedron": [2,3,3,'s',2,radius,center],
		"Truncated_Cube": [2,4,3,'s',1,radius,center],
		"Truncated_Octahedron":[2,4,3,'s',2,radius,center],
		"Truncated_Dodecahedron": [2,5,3,'s',1,radius,center],
		"Truncated_Icosahedron": [2,5,3,'s',2,radius,center],
		"Cuboctahedron": [2,4,3,'c',0,radius,center],
		"Icosidodecahedron": [2,5,3,'c',0,radius,center],
		"Small_Rhombicuboctahedron": [2,4,3,'s',0,center,center],
		"Small_Rhombicosidodecahedron": [2,5,3,'s',0,radius,center],
		"Great_Rhombicuboctahedron" :[2,4,3,'i',radius,center],
		"Great_Rhombicosidodecahedron": [2,5,3,'i',radius,center]}

	triangle_parameters = elemental_triangles[Polyhedron]
	p = triangle_parameters[0]
	q = triangle_parameters[1]
	r = triangle_parameters[2]
	pres = 7
	center_of_sphere = np.array([0,0,0])
	#fix the bug of the dodecahedron. Python works in mystherious ways
	if p == 5:
		p5 = p
		q5 = q
		r5 = r
		tp5 = sympy.pi / p5
		tq5 = sympy.pi / q5
		tr5 = sympy.pi / r5
		cp5 = float(sympy.cos(tp5))
		sp5 = float(sympy.sin(tp5))
		cq5 = float(sympy.cos(tq5))
		cr5 = float(sympy.cos(tr5))
		#Edges are labeled according to the vertex that each of them is opposed to
		lr5 = np.array([1,0,0])
		lq5 = np.array([-cp5,sp5,0])
		lp5 = np.array([-cq5,-(cr5 + cp5*cq5)/sp5])
		lpz5 = sqrt(-(cp5*cq5 + cr5)**2/sp5**2 - cq5**2 + 1)
		lp5 = np.append(lp5,lpz5)

		#Compute the coordinates of the vertices of the triangle
		P5 = np.round((np.cross(lr5,lq5))/(np.linalg.norm(np.cross(lr5,lq5))),pres)
		Q5 = np.round((np.cross(lp5,lr5))/(np.linalg.norm(np.cross(lp5,lr5))),pres)
		R5 = np.round((np.cross(lq5,lp5))/(np.linalg.norm(np.cross(lq5,lp5))),pres)
		p = r5
		r = p5


	tp = sympy.pi / p
	tq = sympy.pi / q
	tr = sympy.pi / r
	cp = float(sympy.cos(tp))
	sp = float(sympy.sin(tp))
	cq = float(sympy.cos(tq))
	cr = float(sympy.cos(tr))
	#Edges are labeled according to the vertex that each of them is opposed to
	lr = np.array([1,0,0])
	lq = np.array([-cp,sp,0])
	lp = np.array([-cq,-(cr + cp*cq)/sp])
	lpz = sqrt(-(cp*cq + cr)**2/sp**2 - cq**2 + 1)
	lp = np.append(lp,lpz)

	#Compute the coordinates of the vertices of the triangle
	P = np.round((np.cross(lr,lq))/(np.linalg.norm(np.cross(lr,lq))),pres)
	Q = np.round((np.cross(lp,lr))/(np.linalg.norm(np.cross(lp,lr))),pres)
	R = np.round((np.cross(lq,lp))/(np.linalg.norm(np.cross(lq,lp))),pres)

	#Calculate the surface area of the spherical triangle.
	#First calculate the spherical excess
	E = tp + tq + tr - sympy.pi
	A = radius*radius*E
	sphere_area = 4*sympy.pi*radius*radius
	tnt = int(sphere_area/A)

	if p == 2:
		Reference_point = P
		E1 = Q
		E2 = R
		tref = tp
		te1 = tq
		te2 = tr
			
	elif q == 2:
		Reference_point = Q
		E1 = R
		E2 = P
		tref = tq
		te1 = tr
		te2 = tp

	elif r == 2:
		Reference_point = R
		E1 = P
		E2 = Q
		tref = tr
		te1 = tp
		te2 = tq

	#Angles of the arcs forming the edges of the triangle
	ang_arcRef_E1 = acos(np.dot(Reference_point,E1)/(np.linalg.norm(Reference_point)*np.linalg.norm(E1)))
	ang_arcRef_E2 = acos(np.dot(Reference_point,E2)/(np.linalg.norm(Reference_point)*np.linalg.norm(E2)))
	ang_arcE1_E2 = acos(np.dot(E1,E2)/(np.linalg.norm(E1)*np.linalg.norm(E2)))

	##########################################################################################################
	#Set the data structure
	sphere_tiling = {}
	tiling_edges = {}
	eye = {}
	triangles = []
	edges_p_t = []
	edges = {}
	#Set the intial values
	#Number of triangles that will cover the sphere to generate the polyhedron
	nt = tnt
	initial_triangle = (tuple(Reference_point),tuple(E1),tuple(E2))
	initial_edges = [(tuple(E1),tuple(E2)),(tuple(Reference_point),tuple(E1)),(tuple(Reference_point),tuple(E2))]
	sphere_tiling[initial_triangle] = initial_edges
	edges_p_t.append(initial_edges)
	triangles.append(initial_triangle)
	for edge in initial_edges:
		tiling_edges[edge] = 'Available'
		edges[edge] = None 

	eye[initial_triangle] = {initial_edges[0]:'Available',initial_edges[1]:'Available',initial_edges[2]:'Available'}

	tri_count = 1
	for triangle in triangles:
		if tri_count < nt:
			Refo = triangle[0]
			E1o = triangle[1]
			E2o = triangle[2]
			sides = sphere_tiling[triangle]
			for side in sides:
				if tiling_edges[side] == 'Available':
					#Count in hw many triangles is the side included
					count_t = 0
					for ttri in sphere_tiling:
						if side[0] in ttri and side[1] in ttri:
							count_t += 1

					#Select the vertex to rotate
					if count_t > 1:
						tiling_edges[side] = 'Occupied'
						pass

					elif count_t == 1:
						for vertex in triangle:
							if vertex not in side:
								to_rotate = vertex
						
						if to_rotate == Refo:
							#Duplicate the triangle by rotatiog the reference using the hipotenuse as mirror plane
							E1d = E1o
							E2d = E2o
							E2c = np.subtract(center_of_sphere,E2d)
							E1c = np.subtract(center_of_sphere,E1d)
							norm_mirror = np.cross(E2c,E1c)
							norm_mirror = norm_mirror/np.linalg.norm(norm_mirror)
							#Obtain the reflection matrix
							RM = reflection_matrix(norm_mirror)
							Rd = np.dot(RM,np.transpose(Refo))
							Refd = np.array([Rd.item(0),Rd.item(1),Rd.item(2)])
							#Built New triangle
							new_triangle = (tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres)))
							new_edges = {(tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres))):'Occupied',(tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres))):None,(tuple(np.round(Refd,pres)),tuple(np.round(E2d,pres))):None}

						if to_rotate == E1o:	
							#Duplicate the triangle by rotating the reference using the hipotenuse as mirror plane
							Refd = Refo
							E2d = E2o
							ang = acos(np.dot(Refo,E1o)/(np.linalg.norm(Refo)*np.linalg.norm(E1o)))
							E1d = General_Rotation(center_of_sphere,np.cross(Refo,E1o),-2*ang,E1o)
							#Built New triangle
							new_triangle = (tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres)))
							new_edges = {(tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres))):None,(tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres))):'Occupied',(tuple(np.round(Refd,pres)),tuple(np.round(E2d,7))):None}


						if to_rotate == E2o:
						#Duplicate the triangle by rotatiog the reference using the hipotenuse as mirror plane
							Refd = Refo
							E1d = E1o
							ang = acos(np.dot(Refo,E2o)/(np.linalg.norm(Refo)*np.linalg.norm(E2o)))
							E2d = General_Rotation(center_of_sphere,np.cross(Refo,E2o),-2*ang,E2o)
						#Built New triangle
							new_triangle = (tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres)))
							new_edges = {(tuple(np.round(E1d,pres)),tuple(np.round(E2d,pres))):None,(tuple(np.round(Refd,pres)),tuple(np.round(E1d,pres))):None,(tuple(np.round(Refd,pres)),tuple(np.round(E2d,pres))):'Occupied'}

						
						if tri_count < nt:
						#Determine if the new formed triangle is repeated
							ctc = 0
							for t_triangle in triangles:
								if new_triangle[0] in t_triangle and new_triangle[1] in t_triangle and new_triangle[2] in t_triangle:
									ctc += 1

						#add new triangle to spherical tiling only if it is not already there.
							
							if ctc == 0:
								#Classify the new_edges
								for n_edge in new_edges:
									if new_edges[n_edge] != 'Occupied':
										cc = 0
										to_update = {}
										for triangle_b in triangles:
											if n_edge[0] in triangle_b and n_edge[1] in triangle_b:
												cc += 1
												to_update[triangle_b] = n_edge
										if cc == 0:
											new_edges[n_edge] = 'Available'
										else:
											new_edges[n_edge] = 'Occupied'
											for triangle_u in to_update:
												d_eye = eye[triangle_u]
												edge_u = to_update[triangle_u]
												d_eye.update({edge_u:'Occupied'})
												eye.update({triangle_u:d_eye})


								n_edges = list(new_edges.keys())
								sphere_tiling[new_triangle] = n_edges
								triangles.append(new_triangle)
								for edge in n_edges:
									edges[edge] = None
								#Actualize the ocupation of the edge in the current triangle
								tiling_edges.update({side:'Occupied'})###
								#########################################

								for edge in new_edges:
									tiling_edges.update({edge:new_edges[edge]})

								act_eye = {}
								for edge_t in sides:
									act_eye[edge_t] = tiling_edges[edge_t]
								eye.update({triangle:act_eye})
								eye[new_triangle] = new_edges
								tri_count += 1
						else:
							break
		else:
			break

	if len(triangles) < tnt:
		if p == 4 or q == 4 or r == 4:
			lt1 = triangles[34]
			lt2 = triangles[35]
			lt3 = triangles[36]
			lt4 = triangles[37]
			hole_triangles = [lt1,lt2,lt3,lt4]

		elif q == 5:
			lt1 = triangles[103]
			lt2 = triangles[105]
			lt3 = triangles[106]
			lt4 = triangles[102]
			hole_triangles = [lt1,lt2,lt3,lt4]

		elif r == 5:
			lt1 = triangles[103]
			lt2 = triangles[104]
			lt3 = triangles[106]
			lt4 = triangles[107]
			hole_triangles = [lt1,lt2,lt3,lt4]
		
		lazy_vertices = {}
		lazy_border = {}
		lazy_references = {}
		for hole_triangle in hole_triangles:
			h_t = [t for t in hole_triangles]
			h_t.remove(hole_triangle)
			for h_vert in hole_triangle:
				for h_tri in h_t:
					if h_vert in h_tri:
						lazy_vertices[h_vert] = None
						ind = h_tri.index(h_vert)
						if ind != 0:
							lazy_border[h_vert] = None
						else:
							lazy_references[h_vert] = None

						break

		lazy_vertices = list(lazy_vertices.keys())
		lazy_border = list(lazy_border.keys())
		lazy_references = list(lazy_references.keys())
		for lazy_vertex in lazy_references:
			lazy_triangle = (lazy_vertex,lazy_border[1],lazy_border[0])
			if q == 5:
				lazy_triangle = (lazy_vertex,lazy_border[0],lazy_border[1])
			triangles.append(lazy_triangle)

	new_triangles = []
	for triangle in triangles:
		triangle = Scaling(triangle,radius)
		tv = []
		for v in triangle:
			v = translation_transformation(v,center)
			tv.append(v)
		triangle = tv
		new_triangles.append(triangle)
	generators = {}
	all_sites = {}
	planes_gen = [] #planes per generator triangle
	points_gens = [] #Points per generator triangle
	for triangle in triangles:
		generator,relevant_points,triangle_type,perimeter_points,planes = Irreducible_Volume(triangle_parameters,triangle)
		generator = np.round(generator,2)
		generators[tuple(generator)] = None
		planes_gen.append(planes)
		points_gens.append(relevant_points)
		for point_r in relevant_points:
			point_r = np.round(point_r,2)
			all_sites[tuple(point_r)] = None
##################################################################################################################################################################
	fig_b = plt.figure()
	bx = fig_b.gca(projection='3d')
	u = np.linspace(0,  2*np.pi, 25)
	v = np.linspace(0, np.pi, 25)
	centerO = np.array((0,0,0))
	radiusO = 1.0
	xS = radiusO * np.outer(np.cos(u), np.sin(v)) + centerO[0]
	yS = radiusO * np.outer(np.sin(u), np.sin(v)) + centerO[1]
	zS = radiusO * np.outer(np.ones(np.size(u)), np.cos(v)) + centerO[2]
	t = np.linspace(0,1,50)
	bx.plot_surface(xS, yS, zS, cmap = 'autumn',alpha = 0.1)
	bx.scatter(centerO[0],centerO[1],centerO[2],color='r',alpha = 1,label = 'Center')
	for triangle in triangles:
		P = np.array(triangle[0])
		Q = np.array(triangle[1])
		R = np.array(triangle[2])
		xpq,ypq,zpq = Slerp(P,Q,t)
		xpr,ypr,zpr = Slerp(P,R,t)
		xqr,yqr,zqr = Slerp(Q,R,t)
		bx.scatter(P[0],P[1],P[2],color = 'g')
		bx.scatter(Q[0],Q[1],Q[2],color = 'g')
		bx.scatter(R[0],R[1],R[2],color = 'g')
		bx.plot(xpq,ypq,zpq,color = 'k')
		bx.plot(xpr,ypr,zpr,color = 'k')
		bx.plot(xqr,yqr,zqr,color = 'k')

	gx = []
	gy = []
	gz = []
	for generator in generators:
		gx.append(generator[0])
		gy.append(generator[1])
		gz.append(generator[2])

	bx.scatter(gx,gy,gz,color = 'k',alpha = 0.75,label = 'Generator vertices')
	plt.title('Sphere tiling with triangle:' + '\n' + str([p,q,r]))
	plt.axis('off')
	plt.legend()

	fig_c = plt.figure(2)
	cx = fig_c.gca(projection='3d')
	u = np.linspace(0,  2*np.pi, 25)
	v = np.linspace(0, np.pi, 25)
	centerO = np.array((0,0,0))
	cx.plot_surface(xS, yS, zS, cmap = 'spring',alpha = 0.1)
	cx.scatter(centerO[0],centerO[1],centerO[2],color='r',alpha = 1,label = 'Center')
	for triangle in triangles:
		P = np.array(triangle[0])
		Q = np.array(triangle[1])
		R = np.array(triangle[2])
		xpq,ypq,zpq = Slerp(P,Q,t)
		xpr,ypr,zpr = Slerp(P,R,t)
		xqr,yqr,zqr = Slerp(Q,R,t)
		cx.plot(xpq,ypq,zpq,color = 'k')
		cx.plot(xpr,ypr,zpr,color = 'k')
		cx.plot(xqr,yqr,zqr,color = 'k')

	sx = []
	sy = []
	sz = []
	for site in all_sites:
		sx.append(site[0])
		sy.append(site[1])
		sz.append(site[2])

	cx.scatter(sx,sy,sz,color = 'b',alpha = 0.5,label = 'Wythoff vertices')	
	cx.scatter(gx,gy,gz,color = 'k',alpha = 0.75,label = 'Generator vertices')

	for planes in planes_gen:
		for plane in planes:
			cx.add_collection3d(plane)

	plt.title('Wythoff construction of a:' + '\n' + Polyhedron)
	plt.axis('off')
	plt.legend()
	plt.show()
##################################################################################################################################################################

	generators = list(generators.keys())
	generators = Scaling(generators,radius)
	gens = []
	for generator in generators:
		generator = translation_transformation(generator,center)
		gens.append(generator)

	generators = gens

	all_sites = list(all_sites.keys())
	all_sites = Scaling(all_sites,radius)
	a_s = []
	for site in all_sites:
		site = translation_transformation(site,center)
		a_s.append(site)

	all_sites = a_s

	return(new_triangles,planes_gen,generators,all_sites,points_gens)
if __name__ == "__Wythoff__":
	Wythoff()  


#Wythoff("Tetrahedron",1,[0,0,0])

