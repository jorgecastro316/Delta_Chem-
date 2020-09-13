from math import *
import numpy as np 
import deltachem.DC_library
from deltachem.DC_library import *

def Tetrahedron(radius,center):
	faces = np.array([[0,1,2],[0,3,1],[0,2,3],[3,2,1]])
	vertices = [[0.8164156395068652, 0.0, 0.5774647206268071],[-0.8164156395068652, 0.0, 0.5774647206268071],[-0.0, 0.8164156395068652, -0.5774647206268071],[0.0, -0.8164156395068652, -0.5774647206268071]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.asarray(vertices)
	#Given edge lenght 'a' calculate the volume of the octahedron
	Volume = (a*a*a)/(6*sqrt(2))
	num_faces = 4 
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Tetrahedron__":
    Tetrahedron()

def Cube(radius,center):
	faces = np.asarray([[3,7,1,5],[1,7,6,4],[4,1,5,2],[5,2,0,3],[6,0,3,7],[4,6,0,2]])
	vertices = [[-0.469,0.000,-0.664],[0.469,0.000,0.664],[-0.469,0.664,0.000],[-0.469,-0.664,0.000],[0.469,0.664,0.000],[-0.469,0.000,0.664],[0.469,0.000,-0.664],[0.469,-0.664,0.000]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.asarray(vertices)
	num_faces = 6 
	#Given edge lenght 'a' calculate the volume of the octahedron
	Volume = (a*a*a)
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Cube__":
    Cube()

def Cuboctahedron(radius,center):
	center_of_sphere = np.array(center)
	faces = [[0,7,1,3],[1,11,7],[2,3,8,10],[2,0,3],[11,9,4],[8,9,10],[4,6,5],[7,6,0],[9,10,5,4],[5,2,0,6],[11,7,6,4],[11,9,8,1],[1,3,8],[10,2,5]]
	vertices = [[0,0,0.813],[0,-0.813,0],[0.575,0.407,0.407],[0.575,-0.407,0.407],[-0.575,0.407,-0.407],[0,0.813,0],[-0.575,0.407,0.407],[-0.575,-0.407,0.407],[0.575,-0.407,-0.407],[0.000,0.000,-0.813],[0.575,0.407,-0.407],[-0.575,-0.407,-0.407]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	num_faces = 14
	#Given the edge lenght 'a' calculate the volume of the Cuboctahedron
	Volume = (5*sqrt(2)*(a*a*a))/3
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Cuboctahedron__":
    Cuboctahedron()

def Dodecahedron(radius,center):
	s = 1/sqrt(3)
	t = sqrt((3-sqrt(5))/6)
	u = sqrt((3+sqrt(5))/6)
	faces = [[0,8,9,4,16],[0,12,13,1,8],[0,16,17,2,12],[8,1,18,5,9],[12,2,10,3,13],[16,4,14,6,17],[9,5,15,14,4],[6,11,10,2,17],[3,19,18,1,13],[7,15,5,18,19],[7,11,6,14,15],[7,19,3,10,11]]
	vertices = [[s,s,s],[s,s,-s],[s,-s,s],[s,-s,-s],[-s,s,s],[-s,s,-s],[-s,-s,s],[-s,-s,-s],[t,u,0],[-t,u,0],[t,-u,0],[-t,-u,0],[u,0,t],[u,0,-t],[-u,0,t],[-u,0,-t],[0,t,u],[0,-t,u],[0,t,-u],[0,-t,-u]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.array(vertices)
	num_faces = 12
	#Given the edge lenght 'a' calculate the volume of the dodecahedron
	Volume = ((15 + 7*(sqrt(5)))/4)*(a*a*a)
	area = 3*(a*a)*sqrt(25 + 10 * sqrt(5))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Dodecahedron__":
    Dodecahedron()

def Great_Rhombicosidodecahedron(radius,center):
	poly_name = "Great_Rhombicosidodecahedron"
	faces = [[47, 23, 88, 54, 32, 97, 1, 83, 53, 21],[84, 111, 4, 3, 9, 7],[76, 6, 54, 32],[24, 9, 3, 17],[18, 62, 96, 69, 52, 74, 81, 89, 12, 55],[35, 39, 87, 75, 95, 90, 44, 80, 79, 19],[8, 16, 78, 117],[8, 42, 45, 65, 10, 117],[7, 60, 13, 84],[114, 90, 44, 14, 59, 61],[5, 0, 73, 48, 107, 86, 57, 29, 118, 63],[77, 94, 104, 50, 101, 34, 11, 72, 99, 103],[4, 22, 31, 111],[49, 37, 55, 12, 99, 72],[6, 76, 13, 84, 111, 31, 105, 51, 98, 15],[100, 33, 38, 40],[18, 55, 37, 116],[43, 70, 41, 85],[58, 42, 45, 119],[105, 31, 22, 27, 96, 69],[21, 47, 112, 93],[105, 51, 52, 69],[28, 71, 65, 45, 119, 40, 38, 68, 61, 114],[15, 86, 107, 98],[68, 61, 59, 36],[93, 92, 102, 46, 53, 21],[114, 28, 95, 90],[65, 71, 66, 10],[28, 95, 75, 30, 66, 71],[113, 62, 96, 27],[74, 81, 73, 48],[89, 103, 77, 0, 73, 81],[76, 13, 60, 25, 97, 32],[14, 91, 80, 44],[16, 8, 42, 58, 24, 9, 7, 60, 25, 64],[23, 29, 118, 109, 112, 47],[30, 106, 87, 75],[43, 85, 33, 100, 17, 3, 4, 22, 27, 113],[82, 49, 37, 116, 70, 41, 36, 59, 14, 91],[0, 5 , 94, 77],[93, 112, 109, 2, 67, 115, 108, 56, 20, 92],[74, 52, 51, 98, 107, 48],[79, 19, 101, 34],[25, 64, 1, 97],[12, 89, 103, 99],[39, 35, 108, 56],[19, 101, 50, 115, 108, 35],[64, 1, 83, 26, 78, 16],[104, 50, 115, 67],[110, 102, 46, 26, 78, 117, 10, 66, 30, 106],[113, 62, 18, 116, 70, 43],[118, 63, 2, 109],[5, 94, 104, 67, 2, 63],[100, 17, 24, 58, 119, 40],[49, 82, 11, 72],[91, 80, 79, 34, 11, 82],[83, 53, 46, 26],[102, 110, 20, 92],[106, 87, 39, 56, 20, 110],[85, 41, 36, 68, 38, 33],[88, 23, 29, 57],[6, 54, 88, 57, 86, 15]]
	vertices = [[0.519,     0.280,    -0.560],[0.214,    -0.733,     0.280],[0.280,   -0.214,    -0.733],[0.107,     0.107,     0.799],[0.280,     0.214,     0.733],[0.453,     0.107,    -0.667],[0.733,    -0.280,     0.214],[0.280,    -0.214,     0.733],[-0.280,    -0.560,     0.519],[0.107,    -0.107,     0.799],[-0.560,    -0.519,     0.280],[-0.280,     0.560,    -0.519],[0.214,     0.733,    -0.280],[0.519,    -0.280,     0.560],[-0.667,     0.453,    -0.107],[0.799,    -0.107,     0.107],[-0.107,    -0.667,     0.453],[-0.107,     0.107,     0.799],[0.107,     0.799,     0.107],[-0.453,     0.107,    -0.667],[-0.280,    -0.560,    -0.519],[0.214,    -0.733,    -0.280],[0.346,     0.387,     0.626],[0.560,    -0.519,    -0.280],[-0.107,    -0.107,     0.799],[0.280,    -0.560,     0.519],[-0.107,    -0.799,     0.107],[0.280,     0.560,     0.519],[-0.799,    -0.107,     0.107],[0.626,    -0.346,    -0.387],[-0.667,    -0.453,    -0.107],[0.519,     0.280,     0.560],[0.560,    -0.519,     0.280],[-0.346,     0.387,     0.626],[-0.346,     0.387,    -0.626],[-0.453,    -0.107,    -0.667],[-0.560,     0.519,     0.280],[-0.107,     0.799,    -0.107],[-0.519,     0.280,     0.560],[-0.519,    -0.280,    -0.560],[-0.453,     0.107,     0.667],[-0.387,     0.626,     0.346],[-0.346,    -0.387,     0.626],[-0.107,     0.667,     0.453],[-0.733,     0.280,    -0.214],[-0.519,    -0.280,     0.560],[-0.107,   -0.799,    -0.107],[0.387,    -0.626,    -0.346],[0.733,     0.280,    -0.214],[-0.214,     0.733,    -0.280],[-0.107,     0.107,    -0.799],[0.733,     0.280,     0.214],[0.667,     0.453,     0.107],[0.107,    -0.799,    -0.107],[0.667,    -0.453,     0.107],[0.107,     0.799,    -0.107],[-0.346,    -0.387,    -0.626],[0.733,    -0.280,    -0.214],[-0.280,    -0.214,     0.733],[-0.667,     0.453,     0.107],[0.346,    -0.387,     0.626],[-0.733,     0.280,     0.214],[0.214,     0.733,     0.280],[0.453,    -0.107,    -0.667],[0.107,    -0.667,     0.453],[-0.626,    -0.346,     0.387],[-0.667,    -0.453,     0.107],[0.107,    -0.107,    -0.799],[-0.626,     0.346,     0.387],[0.560,     0.519,     0.280],[-0.214,     0.733,     0.280],[-0.733,    -0.280,     0.214],[-0.107,     0.667,    -0.453],[0.626,     0.346,    -0.387],[0.667,     0.453,    -0.107],[-0.733,    -0.280,    -0.214],[0.626,    -0.346,     0.387],[0.346,     0.387,    -0.626],[-0.214,    -0.733,     0.280],[-0.519,     0.280,    -0.560],[-0.626,     0.346,    -0.387],[0.560,     0.519,    -0.280],[-0.387,     0.626,    -0.346],[0.107,    -0.799,     0.107],[0.453,    -0.107,     0.667],[-0.280,     0.560,     0.519],[0.799,    -0.107,    -0.107],[-0.626,    -0.346,    -0.387],[0.667,    -0.453,    -0.107],[0.387,     0.626,    -0.346],[-0.799,     0.107,    -0.107],[-0.560,     0.519,    -0.280],[-0.107,    -0.667,    -0.453],[0.107,    -0.667,    -0.453],[0.280,     0.214,   -0.733],[-0.799,    -0.107,    -0.107],[0.387,     0.626,     0.346],[0.387,    -0.626,     0.346],[0.799,     0.107,     0.107],[0.107,     0.667,    -0.453],[-0.280,     0.214,     0.733],[-0.280,     0.214,    -0.733],[-0.214,    -0.733,    -0.280],[0.280,     0.560,    -0.519],[0.107,     0.107,    -0.799],[0.626,     0.346,     0.387],[-0.560,    -0.519,    -0.280],[0.799,     0.107,    -0.107],[-0.280,    -0.214,    -0.733],[0.346,    -0.387,    -0.626],[-0.387,    -0.626,    -0.346],[0.453,     0.107,     0.667],[0.280,    -0.560,    -0.519],[0.107,     0.667,     0.453],[-0.799,     0.107,     0.107],[-0.107,    -0.107,    -0.799],[-0.107,     0.799,     0.107],[-0.387,    -0.626,     0.346],[0.519,    -0.280,    -0.560],[-0.453,    -0.107,     0.667]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.array(vertices)
	num_faces = 62
	#Given the edge lenght 'a' calculate the volume of the Great Rhombicosidodecahedron
	Volume = 5*(a**3)*(19 + 10*sqrt(5))
	area = 30*(a*a)*(1 + sqrt(3) + sqrt(5 + 2*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Great_Rhombicosidodecahedron__":
    Great_Rhombicosidodecahedron()

def Great_Rhombicuboctahedron(radius,center):
	faces = [[36, 13, 23, 40, 15, 31, 32, 28],[40, 15, 19, 3, 30, 7],[37, 34, 45, 26, 2, 38, 21, 47],[24, 30, 3, 17],[42, 9, 35, 19, 3, 17, 18, 44],[6, 10, 0, 46, 13, 23],[44, 38, 21, 27, 1, 42],[37, 25, 20, 34],[10, 39, 14, 0],[47, 16, 29, 33, 25, 37],[28, 5, 22, 43, 41, 36],[9, 35, 31, 32, 11, 4],[19, 35, 31, 15],[18, 17, 24, 12, 26, 2],[39, 8, 45, 34, 20, 14],[18, 44, 38, 2],[0, 46, 41, 43, 33, 25, 20, 14],[9, 42, 1, 4],[26, 12, 8, 45],[40, 7, 6, 23],[22, 29, 33, 43],[32, 28, 5, 11],[8, 12, 24, 30, 7, 6, 10, 39],[29, 16, 27, 1, 4, 11, 5, 22],[21, 47, 16, 27],[46, 13, 36, 41]]
	vertices = [[0.175,-0.774,-0.175],[-0.175,0.774,-0.175],[-0.671,0.175,0.423],[0.175,0.175,0.774],[0.175,0.774,-0.175],[0.423,0.351,-0.599],[0.423,-0.599,0.351],[0.423,-0.351,0.599],[-0.423,-0.599,0.351],[0.175,0.774,0.175],[0.175,-0.774,0.175],[0.423,0.599,-0.351],[-0.423,-0.351,0.599],[0.671,-0.423,-0.175],[-0.175,-0.774,-0.175],[0.671,0.175,0.423],[-0.423,0.351,-0.599],[-0.175,0.175,0.774],[-0.423,0.351,0.599],[0.423,0.351,0.599],[-0.423,-0.599,-0.351],[-0.671,0.423,-0.175],[0.175,0.175,-0.774],[0.671,-0.423,0.175],[-0.175,-0.175,0.774],[-0.423,-0.351,-0.599],[-0.671,-0.175,0.423],[-0.423,0.599,-0.351],[0.671,0.175,-0.423],[-0.175,0.175,-0.774],[0.175,-0.175,0.774],[0.671,0.423,0.175],[0.671,0.423,-0.175],[-0.175,-0.175,-0.774],[-0.671,-0.423,-0.175],[0.423,0.599,0.351],[0.671,-0.175,-0.423],[-0.671,-0.175,-0.423],[-0.671,0.423,0.175],[-0.175,-0.774,0.175],[0.671,-0.175,0.423],[0.423,-0.351,-0.599],[-0.175,0.774,0.175],[0.175,-0.175,-0.774],[-0.423,0.599,0.351],[-0.671,-0.423,0.175],[0.423,-0.599,-0.351],[-0.671,0.175,-0.423]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 26
	#Given the edge lenght 'a' calculate the volume of the Great Rhombicuboctahedron
	Volume = 2 * (a*a*a) * (11 + 7*sqrt(2))
	area = 12*(a*a)*(2 + sqrt(2) + sqrt(3))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Great_Rhombicuboctahedron__":
    Great_Rhombicuboctahedron()

def Icosahedron(radius,center):
	s = (1+sqrt(5))/2
	t = sqrt(1+s*s)
	s = s/t
	t = 1/t
	faces = [[0,8,4],[0,5,10],[2,4,9],[2,11,5],[1,6,8],[1,10,7],[3,9,6],[3,7,11],[0,10,8],[1,8,10],[2,9,11],[3,11,9],[4,2,0],[5,0,2],[6,1,3],[7,3,1],[8,6,4],[9,4,6],[10,5,7],[11,7,5]]
	vertices = [[s,t,0],[-s,t,0],[s,-t,0],[-s,-t,0],[t,0,s],[t,0,-s],[-t,0,s],[-t,0,-s],[0,s,t],[0,-s,t],[0,s,-t],[0,-s,-t]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces= 20
	#Given the edge lenght 'a' calculate the volume of the Icosidodecahedron
	Volume = ((a*a*a)/6)*(45 + 17*sqrt(5))
	area =  5*(a*a)*sqrt(3)
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Icosahedron__":
    Icosahedron()

def Icosidodecahedron(radius,center):
	faces = [[15,4,17],[10,2,22,8,27],[22,2,16],[28,0,24],[12,26,6,5,3],[20,0,24,12,11],[29,21,25],[7,0,20],[28,0,7,1,10],[19,27,8],[26,16,2,28,24],[14,17,9],[23,29,21,8,19],[10,2,28],[16,6,26],[12,24,26],[1,18,7],[29,25,5,15,4],[25,21,22,16,6],[11,14,20],[5,25,6],[17,4,23,13,9],[13,23,19],[3,11,12],[18,7,20,14,9],[1,27,10],[19,27,1,18,13],[13,18,9],[5,15,3],[15,3,11,14,17],[4,29,23],[21,22,8]]
	vertices = [[0.000,0.000,0.813],[-0.251,0.658,0.407],[0.658,0.407,0.251],[0.000,-0.813,0.000],[-0.407,-0.251,-0.658],[0.251,-0.658,-0.407],[0.658,-0.407,-0.251],[-0.407,0.251,0.658],[0.251,0.658,-0.407],[-0.813,0.000,0.000],[0.251,0.658,0.407],[-0.251,-0.658,0.407],[0.251,-0.658,0.407],[-0.658,0.407,-0.251],[-0.658,-0.407,0.251],[-0.251,-0.658,-0.407],[0.813,0.000,0.000],[-0.658,-0.407,-0.251],[-0.658,0.407,0.251],[-0.251,0.658,-0.407],[-0.407,-0.251,0.658],[0.407,0.251,-0.658],[0.658,0.407,-0.251],[-0.407,0.251,-0.658],[0.407,-0.251,0.658],[0.407,-0.251,-0.658],[0.658,-0.407,0.251],[0.000,0.813,0.000],[0.407,0.251,0.658],[0.000,0.000,-0.813]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces= 32
	#Given the edge lenght 'a' calculate the volume of the Icosidodecahedron
	Volume = ((a*a*a)/6)*(45 + 17*sqrt(5))
	area = (a*a)*(5*sqrt(3)+3*sqrt(25 + 10*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Icosidodecahedron__":
    Icosidodecahedron()

def Octahedron(radius,center):
	faces = np.asarray([[4,0,2],[4,2,1],[4,1,3],[4,3,0],[5,2,0],[5,1,2],[5,3,1],[5,0,3]])
	vertices = [[1, 0, 0], [-1, 0,0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.asarray(vertices)
	num_faces = 8 
	#Given edge lenght 'a' calculate the volume of the octahedron
	Volume = (a*a*a)/(3*sqrt(2))
	h = sqrt((a*a) - (a/2)*(a/2))
	area = 2*(a*a)*sqrt(3)
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Octahedron__":
    Octahedron()

def Small_Rhombicosidodecahedron(radius,center):
	faces = [[9,0,2,16],[31,22,20],[57,6,58],[24,7,3,17],[3,7,53],[35,32,1,36],[14,52,44],[59,13,34],[26,50,15,59,33],[11,52,20,30,4],[7,8,9,53],[38,28,45,43,13],[3,22,31,53],[1,27,36],[10,21,57,58],[17,54,38,28],[20,22,46,30],[54,17,3,22,46],[24,32,45,28],[35,45,32],[37,56,34,29,23],[31,44,52,20],[21,10,48],[11,19,12],[38,13,59,15],[44,0,49,14],[21,48,36,27],[45,43,39,35],[38,54,15],[8,7,24,32,1],[52,11,19,14],[50,30,46],[44,31,53,9,0],[59,33,29,34],[25,5,41,10,58],[21,57,2,16,27],[51,37,56,39],[34,13,43,56],[16,8,1,27],[19,18,47,12],[0,2,49],[17,24,28],[23,55,42,29],[4,11,12,40],[6,49,14,19,18],[37,23,5,41],[23,55,5],[43,39,56],[47,55,5,25],[15,54,46,50],[16,9,8],[6,18,25,58],[18,47,25],[48,36,35,39,51],[26,33,42,40],[33,29,42],[47,12,40,42,55],[48,51,41,10],[51,37,41],[30,4,26,50],[57,6,49,2],[26,4,40]]
	vertices = [[0.771,-0.182,0.182],[0.000,-0.659,0.477],[0.659,-0.477,0.000],[0.182,0.182,0.771],[0.182,0.771,-0.182],[-0.182,-0.182,-0.771],[0.589,-0.295,-0.477],[0.182,-0.182,0.771],[0.295,-0.477,0.589],[0.589,-0.295,0.477],[0.000,-0.659,-0.477],[0.477,0.589,-0.295],[0.295,0.477,-0.589],[-0.771,0.182,0.182],[0.771,0.182,-0.182],[-0.477,0.589,0.295],[0.477,-0.589,0.295],[-0.182,0.182,0.771],[0.477,0.000,-0.659],[0.589,0.295,-0.477],[0.477,0.589,0.295],[0.182,-0.771,-0.182],[0.295,0.477,0.589],[-0.477,0.000,-0.659],[-0.182,-0.182,0.771],[0.182,-0.182,-0.771],[-0.182,0.771,-0.182],[0.182,-0.771,0.182],[-0.477,0.000,0.659],[-0.589,0.295,-0.477],[0.182,0.771,0.182],[0.589,0.295,0.477],[-0.295,-0.477,0.589],[-0.477,0.589,-0.295],[-0.771,0.182,-0.182],[-0.477,-0.589,0.295],[-0.182,-0.771,0.182],[-0.589,-0.295,-0.477],[-0.589,0.295,0.477],[-0.659,-0.477,0.000],[0.000,0.659,-0.477],[-0.295,-0.477,-0.589],[-0.295,0.477,-0.589],[-0.771,-0.182,0.182],[0.771,0.182,0.182],[-0.589,-0.295,0.477],[0.000,0.659,0.477],[0.182,0.182,-0.771],[-0.182,-0.771,-0.182],[0.771,-0.182,-0.182],[-0.182,0.771,0.182],[-0.477,-0.589,-0.295],[0.659,0.477,0.000],[0.477,0.000,0.659],[-0.295,0.477,0.589],[-0.182,0.182,-0.771],[-0.771,-0.182,-0.182],[0.477,-0.589,-0.295],[0.295,-0.477,-0.589],[-0.659,0.477,0.000]]	
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 62
	#Given the edge lenght 'a' calculate the volume of the Small Rhombicosidodecahedron
	Volume = ((a*a*a)/3)*(60 + 29*sqrt(5))
	#area = (a*a)*(30 + 5*sqrt(3) + 3*sqrt(25 + 10*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Small_Rhombicosidodecahedron__":
    Small_Rhombicosidodecahedron()

def Small_Rhombicuboctahedron(radius,center):
	faces = [[9,3,17,20],[8,23,15,12],[6,13,11,15],[22,7,3,17],[4,7,8,19],[4,1,18,19],[13,16,11],[8,12,19],[0,17,22],[20,21,2],[23,6,15],[9,20,2,10],[3,9,1,4],[11,5,12,15],[18,14,5],[17,20,21,0],[1,9,10],[12,19,18,5],[22,23,6,0],[8,7,22,23],[6,0,21,13],[1,18,14,10],[3,7,4],[16,2,10,14],[13,16,2,21],[16,11,5,14]]
	vertices = [[-0.702,0.000,0.411],[0.702,0.411,0.000],[-0.291,0.702,-0.291],[0.291,0.291,0.702],[0.702,0.000,0.411],[0.291,-0.291,-0.702],[-0.702,-0.411,0.000],[0.291,-0.291,0.702],[0.291,-0.702,0.291],[0.291,0.702,0.291],[0.291,0.702,-0.291],[-0.291,-0.291,-0.702],[0.291,-0.702,-0.291],[-0.702,0.000,-0.411],[0.291,0.291,-0.702],[-0.291,-0.702,-0.291],[-0.291,0.291,-0.702],[-0.291,0.291,0.702],[0.702,0.000,-0.411],[0.702,-0.411,0.000],[-0.291,0.702,0.291],[-0.702,0.411,0.000],[-0.291,-0.291,0.702],[-0.291,-0.702,0.291]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 26
	#Given the edge lenght 'a' calculate the volume of the Small Rhombicuboctahedron
	Volume = 2/3 * (a*a*a) * ( 6 + 5*sqrt(2))
	#area = 2*(a*a)*(9 + sqrt(3))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Small_Rhombicuboctahedron__":
    Small_Rhombicuboctahedron()

def Snub_Cube(radius,center):
	faces = [[4,20,16,3],[5,14,11,0],[15,1,2,23],[14,3,16],[14,5,3],[3,5,1],[3,1,4],[4,1,15],[5,2,1],[5,0,2],[10,21,13,12],[19,6,22,7],[17,18,8,9],[7,22,21],[7,21,10],[8,7,10],[8,10,9],[9,10,12],[13,0,12],[12,0,11],[23,2,13],[21,23,13],[13,2,0],[11,14,17],[20,6,19],[6,20,4],[22,15,23],[6,15,22],[6,4,15],[14,16,17],[18,17,16],[18,16,20],[12,11,9],[8,19,7],[18,19,8],[20,19,18],[11,17,9],[21,22,23]]
	vertices = [[0.204,-0.691,-0.376],[0.691,-0.204,0.376],[0.691,-0.376,-0.204], [0.204,-0.376,0.691],[0.376,0.204,0.691],[0.376,-0.691,0.204],[0.204,0.691,0.376],[-0.204,0.691,-0.376],[-0.691,0.376,-0.204],[-0.691,-0.204,-0.376],[-0.376,0.204,-0.691],[-0.376,-0.691,-0.204],[-0.204,-0.376,-0.691],[0.376,-0.204,-0.691],[-0.204,-0.691,0.376],[0.691,0.376,0.204],[-0.376,-0.204,0.691],[-0.691,-0.376,0.204],[-0.691,0.204,0.376],[-0.376,0.691,0.204],[-0.204,0.376,0.691],[0.204,0.376,-0.691],[0.376,0.691,-0.204],[0.691,0.204,-0.376]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.array(vertices)
	num_faces = 38
	#Given the edge lenght 'a' calculate the volume of the Snub Cube
	#Calculate first the Tribonacci constant
	t = ( 1 +  (19 + 3*sqrt(33))**(1/3) + (19 - 3*sqrt(33))**(1/3))/3
	Volume = ((a*a*a)*(3*sqrt(t - 1) + 4*sqrt(t + 1)))/(3*sqrt(2 - t))
	#area = 2*(a*a)*(3 + 4*sqrt(3))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Snub_Cube__":
    Snub_Cube()

def Snub_Dodecahedron(radius,center):
	faces = [[16, 15, 13, 36, 10],[18, 1, 2, 12, 4],[51, 35, 48, 58, 40],[39, 11, 3, 46, 45],[47, 49, 0, 14, 50],[43, 38, 37, 34, 32],[26, 29, 20, 53, 52],[33, 56, 23, 27, 19],[9,17, 5, 6, 28],[59, 42, 57, 21, 24],[22, 25, 8, 41, 7],[31, 54, 55, 30, 44],[10, 35, 37],[36, 11, 48],[13, 0, 3],[15, 2, 14],[16, 38, 12],[51, 33, 34],[39, 41, 58],[49, 42, 46],[1,5, 50],[43, 44, 4],[40, 8, 56],[45, 59, 7],[47, 17, 57],[18, 30, 6],[32, 19, 31],[53, 28, 55],[20, 21, 9],[29, 22, 24],[26, 23, 25],[52, 54, 27],[10, 36, 35],[36, 13, 11],[13, 15,0],[15, 16, 2],[16, 10, 38],[48, 35, 36],[3, 11, 13],[14, 0, 15],[12, 2, 16],[37, 38, 10],[35, 51, 37],[11, 39, 48],[0, 49, 3],[2, 1, 14],[38, 43, 12],[51, 40, 33],[39, 45, 41],[49, 47, 42],[1, 18, 5],[43, 32, 44],[40, 58, 8],[45, 46, 59],[47, 50, 17],[18, 4, 30],[32, 34, 19],[58, 48, 39],[46, 3, 49],[50, 14, 1],[4,12, 43],[34, 37, 51],[53, 20, 28],[20, 29, 21],[29, 26, 22],[26, 52, 23],[52, 53, 54],[9, 28, 20],[24, 21, 29],[25, 22, 26],[27, 23, 52],[55, 54, 53],[28, 6, 55],[21, 57, 9],[22, 7, 24],[23, 56, 25],[54, 31, 27],[6, 5, 18],[57, 42, 47],[7, 41, 45], [56, 33, 40],[31, 44, 32],[5, 17, 50],[42, 59, 46],[41, 8, 58],[33, 19, 34],[44, 30, 4],[17, 9, 57],[59, 24, 7],[8, 25,56],[19, 27, 31],[30, 55, 6]]
	vertices = 	[[0.408,-0.371,0.597],[-0.146,-0.733,0.321],[-0.227,-0.503,0.597],[0.630,-0.051,0.512],[-0.652,-0.450,0.183],[-0.081,-0.806,-0.065],[-0.313,-0.678,-0.321],[0.652,0.365,-0.321],[0.227,0.759,-0.183],[0.243,-0.583,-0.512],[-0.243,0.197,0.750],[0.479,0.274,0.597],[-0.540,-0.329,0.512],[0.313,-0.016,0.750],[0.146,-0.615,0.512],[0.081,-0.302,0.750],[-0.262,-0.171,0.750],[0.262,-0.747,-0.183],[-0.408,-0.700,0.065],[-0.652,0.450,-0.183],[0.112,-0.292,-0.750],[0.479,-0.274,-0.597],[0.408,0.371,-0.597],[-0.227,0.503,-0.597],[0.630,0.051,-0.512],[0.146,0.615,-0.512],[0.081,0.302,-0.750],[-0.540,0.329,-0.512],[-0.112,-0.540,-0.597],[0.313,0.016,-0.750],[-0.630,-0.479,-0.183],[-0.742,0.088,-0.321],[-0.792,0.172,0.065],[-0.408,0.700,-0.065],[-0.630,0.479,0.183],[-0.112,0.540,0.597],[0.112,0.292,0.750],[-0.479,0.412,0.512],[-0.549,0.060,0.597],[0.549,0.507,0.321],[-0.081,0.806,0.065],[0.540,0.605,-0.065],[0.742,-0.326,-0.065],[-0.742,-0.088,0.321],[-0.792,-0.171,-0.065],[0.742,0.326,0.065],[0.792,-0.019,0.183],[0.540,-0.605,0.065],[0.243,0.583,0.512],[0.652,-0.365,0.321],[0.227,-0.759,0.183],[-0.313,0.678,0.321],[-0.262,0.172,-0.750],[-0.243,-0.197,-0.750],[-0.549,-0.060,-0.597],[-0.479,-0.412,-0.512],[-0.146,0.733,-0.321],[0.549,-0.507,-0.321],[0.262,0.747,0.183],[0.792,0.019,-0.183]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	#Given the edge lenght 'a' calculate the volume of the Snub Dodecahedron
	#Golden ratio
	phi = 1.618033988749895
	#Epsilon constant
	eps = ((phi/2) + (sqrt(phi - 5/27))/2)**(1/3) + ((phi/2) - (sqrt(phi - 5/27))/2)**(1/3)
	num_faces = 92 
	Volume = (a*a*a)*((12*(eps**2)*(3*phi+1))-(eps*(36*phi+7))-(53*phi+6))/(6*((sqrt(3-(eps*eps)))**3))
	#area = (a*a)*(20*sqrt(3) + 3*sqrt(25 + 10*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Snub_Dodecahedron__":
    Snub_Dodecahedron()

def Truncated_Cube(radius,center):
	faces = [[11,23,0,22,8,5,18,15],[1,14,6],[9,23,0],[4,3,21],[19,18,5],[22,10,8],[13,21,3,17,20,6,1,7],[12,17,20],[15,11,2],[15,18,19,12,20,6,14,2],[16,9,0,22,10,4,21,13],[7,16,13],[10,4,3,17,12,19,5,8],[14,2,11,23,9,16,7,1]]
	vertices = [[0.552,0.229,-0.552],[-0.229,-0.780,0.000],[-0.552,-0.229,-0.552],[0.229,0.000,0.780],[0.552,0.229,0.552],[-0.229,0.780,0.000],[-0.552,-0.552,0.229],[0.229,-0.780,0.000],[0.229,0.780,0.000],[0.552,-0.229,-0.552],[0.552,0.552,0.229],[-0.229,0.000,-0.780],[-0.552,0.229,0.552],[0.552,-0.552,0.229],[-0.552,-0.552,-0.229],[-0.552,0.229,-0.552],[0.552,-0.552,-0.229],[-0.229,0.000,0.780],[-0.552,0.552,-0.229],[-0.552,0.552,0.229],[-0.552,-0.229,0.552],[0.552,-0.229,0.552],[0.552,0.552,-0.229],[0.229,0.000,-0.780]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	#Given the edge lenght 'a' calculate the volume of the Truncated Cube
	num_faces = 14
	Volume = ((a*a*a)/3)*(21 + 14*sqrt(2))
	#area = 2*(a*a)*(6 + 6*sqrt(2) + sqrt(3))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1.0:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Truncated_Cube__":
    Truncated_Cube()

def Truncated_Dodecahedron(radius,center):
	faces = [[53, 56, 58, 13, 32, 1, 30, 59, 50, 28],[5, 13, 32],[23, 37, 26, 27, 12, 34, 19, 15, 8, 18],[42, 3, 21],[19, 34, 36],[51, 48, 2],[37, 26, 49],[24, 43, 17],[57, 27, 12],[46, 21, 42, 31, 0, 52, 14, 25, 29, 4],[10, 54, 11, 47, 45, 24, 43, 38, 30, 59],[5,55, 0, 52, 22, 20, 33, 9, 58, 13],[41, 48, 51, 7, 4, 29, 16, 35, 18, 8],[32, 5, 55, 31, 42, 3, 17, 43, 38, 1],[45, 24, 17, 3, 21, 46, 7, 51, 2, 6],[15, 41, 8],[30, 38, 1],[28, 40, 53],[39, 36, 34, 12, 57, 40, 28, 50, 10, 54],[6, 47, 45],[9, 56, 58],[35, 23, 18],[29, 16, 25],[20, 44, 33],[22, 52, 14],[7, 46, 4],[55, 0, 31],[50, 59, 10],[39, 11, 54],[16, 35, 23, 37, 49, 44, 20, 22, 14, 25],[48, 41, 15, 19, 36, 39, 11, 47, 6, 2],[44, 49, 26, 27, 57, 40, 53, 56, 9, 33]]
	vertices = [[0.580,0.443,0.358],[-0.358,0.580,0.443],[-0.137,-0.717,0.358],[0.137,0.000,0.801],[0.580,-0.443,0.358],[0.137,0.717,0.358],[-0.358,-0.580,0.443],[0.358,-0.580,0.443],[0.137,-0.717,-0.358],[0.137,0.717,-0.358],[-0.801,0.137,0.000],[-0.717,-0.358,0.137],[-0.358,-0.137,-0.717],[0.000,0.801,0.137],[0.801,0.137,0.000],[-0.137,-0.717,-0.358],[0.717,-0.358,-0.137],[-0.137,0.000,0.801],[0.358,-0.580,-0.443],[-0.358,-0.580,-0.443],[0.580,0.443,-0.358],[0.358,-0.137,0.717],[0.717,0.358,-0.137],[0.443,-0.358,-0.580],[-0.358,-0.137,0.717],[0.801,-0.137,0.000],[0.137,0.000,-0.801],[-0.137,0.000,-0.801],[-0.580,0.443,-0.358],[0.717,-0.358,0.137],[-0.580,0.443,0.358],[0.443,0.358,0.580],[-0.137,0.717,0.358],[0.358,0.580,-0.443],[-0.443,-0.358,-0.580],[0.580,-0.443,-0.358],[-0.580,-0.443,-0.358],[0.358,-0.137,-0.717],[-0.443,0.358,0.580],[-0.717,-0.358,-0.137],[-0.443,0.358,-0.580],[0.000,-0.801,-0.137],[0.358,0.137,0.717],[-0.358,0.137,0.717],[0.443,0.358,-0.580],[-0.443,-0.358,0.580],[0.443,-0.358,0.580],[-0.580,-0.443,0.358],[0.000,-0.801,0.137],[0.358,0.137,-0.717],[-0.717,0.358,-0.137],[0.137,-0.717,0.358],[0.717,0.358,0.137],[-0.358,0.580,-0.443],[-0.801,-0.137,0.000],[0.358,0.580,0.443],[-0.137,0.717,-0.358],[-0.358,0.137,-0.717],[0.000,0.801,-0.137],[-0.717,0.358,0.137]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 32
	#Given the edge lenght 'a' calculate the volume of the Truncated Dodecahedron
	Volume = (5/12) * (a**3) * (99 + 47*sqrt(5))
	#area = 5*(a*a)*(sqrt(3) + 6*sqrt(5 + 2*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Truncated_Dodecahedron__":
    Truncated_Dodecahedron()

def Truncated_Icosahedron(radius,center):
	faces = [[7, 0, 54, 23, 13, 30],[51, 25, 42, 18, 16, 44],[18, 59, 11, 27, 16],[31, 32, 3, 24, 21, 49],[48, 17, 37, 22, 56],[28, 1, 53, 12, 47, 57],[21, 24, 40, 4, 6],[45, 38, 43, 3, 24, 40],[44, 9, 33, 52, 51],[37, 17, 33, 52, 15, 20],[18, 59, 39, 38, 43, 42],[28, 1, 58, 4, 6, 36],[5,39, 38, 45, 14],[25, 51, 52, 15, 31, 32],[43, 3, 32, 25, 42],[13, 23, 35, 47, 12],[11, 26, 7, 0, 50, 27],[49, 31, 15, 20, 29],[34, 19, 30, 13, 12, 53],[10, 14, 45, 40, 4, 58],[41, 22, 56, 35, 47, 57],[19, 30, 7, 26, 55],[2, 36, 6, 21, 49, 29],[44, 9, 8, 50, 27,16],[57, 28, 36, 2, 41],[20, 37, 22, 41, 2, 29],[46, 54, 0, 50, 8],[26, 11, 59, 39, 5, 55],[5, 55, 19, 34, 10, 14],[53, 1, 58, 10, 34],[48, 46, 54, 23, 35, 56],[9, 8, 46, 48, 17, 33]]
	vertices = [[-0.265,0.328,-0.695],[-0.164,-0.796,0.000],[0.594,-0.531,0.164],[0.000,0.164,0.796],[-0.164,-0.594,0.531],[-0.796,0.000,0.164],[0.164,-0.594,0.531],[-0.531,0.164,-0.594],[0.164,0.594,-0.531],[0.328,0.695,-0.265],[-0.594,-0.531,0.164],[-0.594,0.531,-0.164],[-0.164,-0.594,-0.531],[-0.265,-0.328,-0.695],[-0.695,-0.265,0.328],[0.695,0.265,0.328],[-0.164,0.796,0.000],[0.695,0.265,-0.328],[-0.328,0.695,0.265],[-0.695,-0.265,-0.328],[0.796,0.000,0.164],[0.265,-0.328,0.695],[0.695,-0.265,-0.328],[0.000,-0.164,-0.796],[0.000,-0.164,0.796],[0.164,0.594,0.531],[-0.695,0.265,-0.328],[-0.328,0.695,-0.265],[0.164,-0.796,0.000],[0.695,-0.265,0.328],[-0.531,-0.164,-0.594],[0.531,0.164,0.594],[0.265,0.328,0.695],[0.594,0.531,-0.164],[-0.594,-0.531,-0.164],[0.265,-0.328,-0.695],[0.328,-0.695,0.265],[0.796,0.000,-0.164],[-0.531,0.164,0.594],[-0.695,0.265,0.328],[-0.265,-0.328,0.695],[0.594,-0.531,-0.164],[-0.164,0.594,0.531],[-0.265,0.328,0.695],[0.164,0.796,0.000],[-0.531,-0.164,0.594],[0.265,0.328,-0.695],[0.164,-0.594,-0.531],[0.531,0.164,-0.594],[0.531,-0.164,0.594],[-0.164,0.594,-0.531],[0.328,0.695,0.265],[0.594,0.531,0.164],[-0.328,-0.695,-0.265],[0.000,0.164,-0.796],[-0.796,0.000,-0.164],[0.531,-0.164,-0.594],[0.328,-0.695,-0.265],[-0.328,-0.695,0.265],[-0.594,0.531,0.164]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices = np.array(vertices)
	num_faces = 32
	#Given the edge lenght 'a' calculate the volume of the Truncated Icosahedron
	Volume = ((a*a*a)/4) * (125 + 43*sqrt(5))
	#area =  3*(a*a)*(10*sqrt(3) + sqrt(25 + 10*sqrt(5)))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Truncated_Icosahedron__":
    Truncated_Icosahedron()

def Truncated_Octahedron(radius,center):
	faces = [[4,14,19,13,17,23],[5,16,7,22,6,18],[16,7,14,19],[10,8,3,0,21,12],[1,2,11,20,15,9],[10,1,9,12],[11,2,17,13],[14,4,3,0,22,7],[20,5,18,15],[15,9,12,21,6,18],[16,5,20,11,13,19],[2,1,10,8,23,17],[8,3,4,23],[6,21,0,22]]
	vertices = [[0.000,-0.257,0.771],[0.727,0.257,-0.257],[0.364,0.514,-0.514],[0.000,0.257,0.771],[-0.364,0.514,0.514],[-0.364,-0.514,-0.514],[0.000,-0.771,0.257],[-0.727,-0.257,0.257],[0.364,0.514,0.514],[0.727,-0.257,-0.257],[0.727,0.257,0.257],[0.000,0.257,-0.771],[0.727,-0.257,0.257],[-0.364,0.514,-0.514],[-0.727,0.257,0.257],[0.364,-0.514,-0.514],[-0.727,-0.257,-0.257],[0.000,0.771,-0.257],[0.000,-0.771,-0.257],[-0.727,0.257,-0.257],[0.000,-0.257,-0.771],[0.364,-0.514,0.514],[-0.364,-0.514,0.514],[0.000,0.771,0.257]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 14
	#Given the edge lenght 'a' calculate the volume of the Truncated Octahedron
	Volume = 8*(a*a*a)*sqrt(2)
	#area = 6*(a*a)*(1 + 2*sqrt(3))
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Truncated_Octahedron__":
    Truncated_Octahedron()

def Truncated_Tetrahedron(radius,center):
	faces = [[9,5,7],[6,4,8],[10,11,0,6,4,2],[4,2,3,9,7,8],[1,11,0],[10,3,2],[6,0,1,5,7,8],[5,1,11,10,3,9]]
	vertices = [[ -0.347,0,-0.735],[-0.693,-0.347,-0.245],[ 0.347,0.693,0.245],[0,0.347,0.735],[ 0.693,0.347,-0.245],[-0.347,-0.693,0.245],[0.347,0.000,-0.735],[0.347,-0.693,0.245],[0.693,-0.347,-0.245],[0.000,-0.347,0.735],[-0.347,0.693,0.245],[ -0.693,0.347,-0.245]]
	vertices,a = Scaling_Translation(vertices,radius,center,faces)
	vertices=np.array(vertices)
	num_faces = 8
	#Given the edge lenght 'a' calculate the volume of the Truncated Tetrahedron
	Volume = (23*sqrt(2)*(a*a*a))/12
	area = 7*sqrt(3)*(a*a)
	faces_areas = {}
	Area = 0
	for face in faces:
		face_verts = []
		area_face = 0
		for vert_ind in face:
			face_verts.append(vertices[vert_ind])
		f_baricenter,f_radius = Baricenter(face_verts)
		for i in range((len(face))):
			v_1 = np.array(vertices[face[i-1]])
			v_2 = np.array(vertices[face[i]])
			v1_b = np.subtract(v_1,f_baricenter)
			v2_b = np.subtract(v_2,f_baricenter)
			t_area = np.linalg.norm(np.cross(v1_b,v2_b))/2
			area_face += t_area
		faces_areas[tuple(face)] = area_face
		Area += area_face

	areas = list(faces_areas.values())
	areas_faces = {}
	for area in areas:
		area_faces = []
		for face in faces_areas:
			if abs(faces_areas[face] - area) <= 1:
				area_faces.append(face)
		areas_faces[area] = area_faces
	max_area = max(areas)
	min_area = min(areas)
	return(faces,vertices,num_faces,Volume,Area,max_area,min_area,areas_faces,faces_areas)
if __name__ == "__Truncated_Tetrahedron__":
    Truncated_Octahedron()

def Scaling_Translation(vertices,radius,center,faces):
	#This function creates two 4X4 matrices for scaling and traslation 
	#of a set of coordinates to fit on a sphere of an arbitrary radius and center
	#In order to obtain the scaling factor, the maximum value of the set
	#of coordinates is obtained.
	max_coords = []
	poly_baricenter,max_distance = Baricenter(vertices)
	#print("vertices", vertices)
	for i in vertices:
		max_coords.append(max(i))
		i.append(1)
	scale_factor = radius/max_distance
	#Built the scaling matrix 4x4
	S = np.matrix([[scale_factor,0,0,0],[0,scale_factor,0,0],[0,0,scale_factor,0],[0,0,0,1]])
	#Built the translation matrix using the coordinates of the center of the sphere(baricenter
	#of the set of coordinates)
	T = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[center[0],center[1],center[2],1]])
	#take their matrix product to create one 4 × 4 matrix that gives the transformation
	#print(T)
	ST = np.dot(S,T)
	new_verts = []
	#multiply each vertor by the trasnformation matrix ST
	for j in vertices:
		j = np.asarray(j)
		new_verts.append(np.dot(j,ST))

	transformed_verts = []

	for k in range(0,len(new_verts)):
		vert = (new_verts[k].item(0),new_verts[k].item(1),new_verts[k].item(2))
		transformed_verts.append(vert)

	#Compute the lenght of the edge of the polyhedron
	base = faces[0]
	edge_lenght = np.linalg.norm(np.asarray(transformed_verts[base[1]]) - np.asarray(transformed_verts[base[0]]))

	return(transformed_verts,edge_lenght)
if __name__ == "__Scaling_Translation__":
    Scaling_Translation()
