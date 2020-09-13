from math import *
import numpy as np

def General_Rotation(point,vector,angle,obj):

	center = np.array([0,0,0])
	vector = np.append(vector,1)
	
	#We wish to rotate an object about an arbitrary axis defined by a 'point' and a direction vector 'vector'.
	#A rotation of 'angle' about that axis can be defined by a concatening the following transformations
	#First a traslation matrix is used to relocate 'point' to the origin of coordinates

	T0 = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[center[0],center[1],center[2],1]])
	#vector = np.append(vector,1)
	#v_0 = np.dot(vector,T0)
	
	#Then the vector is rotated until it matches with the yz plane using a rotation matrix, R_y 
	rotation_angle = -(atan2(radians(degrees(vector[0])),radians(degrees(vector[2]))))

	R_y = np.matrix([[cos(rotation_angle),0,-sin(rotation_angle),0],[0,1,0,0],[sin(rotation_angle),0,cos(rotation_angle),0],[0,0,0,1]])
	
	#v_yz = np.dot(v_0,R_y)

	#After that the vector is rotated with respect to X axis until it matches with Z axis by using the rotation matrix R_x
	rotation_angle = atan2(radians(degrees(vector[1])),sqrt(radians(degrees(vector[0]))**2 + radians(degrees(vector[2]))**2))
	R_x = np.matrix([[1,0,0,0],[0,cos(rotation_angle),sin(rotation_angle),0],[0,-sin(rotation_angle),cos(rotation_angle),0],[0,0,0,1]])
	#v_z = np.dot(v_yz,R_x)

	#Rotation of 'angle' radians  through Z-axis with R_z as rotation matrix

	rotation_angle = angle
	R_z = np.matrix([[cos(rotation_angle),sin(rotation_angle),0,0],[-sin(rotation_angle),cos(rotation_angle),0,0],[0,0,1,0],[0,0,0,1]])
	#v_z_a = np.dot(v_z,R_z)

	#Rotation with R_x to revert the last process

	rotation_angle = -atan2(radians(degrees(vector[1])),sqrt(radians(degrees(vector[0]))**2 + radians(degrees(vector[2]))**2))
	R_x_b = np.matrix([[1,0,0,0],[0,cos(rotation_angle),sin(rotation_angle),0],[0,-sin(rotation_angle),cos(rotation_angle),0],[0,0,0,1]])
	#v_x_b = np.dot(v_z_a,R_x_b)

	#Rotation back to the starting point using R_y

	rotation_angle = (atan2(vector[0],vector[2])) 
	R_y_b = np.matrix([[cos(rotation_angle),0,-sin(rotation_angle),0],[0,1,0,0],[sin(rotation_angle),0,cos(rotation_angle),0],[0,0,0,1]])
	#v_y_b = np.dot(v_x_b,R_y_b)


	#Finally apply a traslation back to the origin 
	T_p = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[point[0],point[1],point[2],1]])
	#v_t = np.dot(v_y_b,T_p)

	#Define the total transformation matrix as a multiplication of all the matrices used throug out the process

	T = T0*R_y*R_x*R_z*R_x_b*R_y_b*T_p

	obj = np.append(obj,1)
	obj_out = np.dot(obj,T)
	obj_rot = []
	obj_rot.append(obj_out.item(0))
	obj_rot.append(obj_out.item(1))
	obj_rot.append(obj_out.item(2))
	obj_rot = np.asarray(obj_rot)
	return(obj_rot)






