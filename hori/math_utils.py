import math

def distance(a1, a2):
	dx = a1.x - a2.x
	dy = a1.y - a2.y
	dz = a1.z - a2.z
	return math.sqrt(dx*dx + dy*dy + dz*dz)	

def angle_between_vectors(v1, v2):
	dotp = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
	m1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
	m2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
	if m1*m2 < 1e-12:
		return 0.0
	c = dotp / (m1*m2)
	c = max(min(c, 1.0), -1.0)
	return math.degrees(math.acos(c))

def cross_product(v1, v2):
	return (v1[1]*v2[2] - v1[2]*v2[1],
			v1[2]*v2[0] - v1[0]*v2[2],
			v1[0]*v2[1] - v1[1]*v2[0])

def dist_3d(a, b):
	dx = a[0] - b[0]
	dy = a[1] - b[1]
	dz = a[2] - b[2]
	return math.sqrt(dx*dx + dy*dy + dz*dz)
