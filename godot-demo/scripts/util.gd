extends Node

# if number is negative (a 64bit integer is a negative 32 integer if 
# it is bigger than 2^31 - 1 (msb = 1)
func from_32_to_64(i : int) -> int:
	if i >= (1 << 31): return i - (1 << 32)
	return i


# open a file at given path and fill octree and 
# gravity_values arrays with data
func get_universal_octree(path : String, octree : PackedInt32Array, gravity_values : PackedFloat32Array) -> Array:
	var min : Vector3 = Vector3()
	var edge : float = 0.0
	
	var file : FileAccess = FileAccess.open(path, FileAccess.READ)
	
	var type : int = file.get_32()
	
	# read min and edge
	min.x = file.get_float()
	min.y = file.get_float()
	min.z = file.get_float()
	edge = file.get_float()
	
	# read octree size
	var octree_size = file.get_32()
	
	# fill octree
	for i in octree_size:
		octree.append(file.get_32())
	
	if(type == 0):
		# read gravity_values size
		var gravity_values_size = file.get_32()
		
		# fill gravity_values
		for j in gravity_values_size:
			var x = file.get_float()
			var y = file.get_float()
			var z = file.get_float()
			gravity_values.append(x)
			gravity_values.append(y)
			gravity_values.append(z)
			
	return [min, edge]


func read_octree(p : Vector3, octree : PackedInt32Array, gravity_values : PackedFloat32Array, min : Vector3, edge : float) -> Vector3:
	var current_min = min
	var current_edge = edge
	
	var output = Vector3(0, 0, 0)
	
	var i : int = 0
	while true:
		if octree[i] <= 0:
			if gravity_values.size() > 0:
				var values : Array[Vector3] = []
				for j in 8:
					values.append(
						Vector3(
							gravity_values[(-octree[i + j])*3],
							gravity_values[(-octree[i + j])*3 + 1],
							gravity_values[(-octree[i + j])*3 + 2]))
				output = interpolate(p, get_box(current_min, current_edge), values)
				break
			else:
				var values: Array[float] = []
				for j in range(8):
					var sl = octree.slice(i+j, i+j+1)
					var sl_to_bytes = sl.to_byte_array()
					var sl_float = sl_to_bytes.to_float32_array()
					values.append(sl_float[0])
				output = get_gradient(p, get_box(current_min, current_edge), values)
				break
				
		else:
			var k : int = 0;
			if p.x > current_min.x + current_edge/2.0: k += 1
			if p.y > current_min.y + current_edge/2.0: k += 2
			if p.z > current_min.z + current_edge/2.0: k += 4
			i = octree[i + k]
			current_edge /= 2.0
			current_min = get_box(current_min, current_edge)[k]
	return output


func interpolate(p : Vector3, box : Array[Vector3], values : Array[Vector3]) -> Vector3:
	var weights : Array[float] = trilinear_coordinates(p, box)
	var output = Vector3()
	for i in 8: 
		output += values[i]*weights[i]
	return output;


func get_box(min : Vector3, edge : float) -> Array[Vector3]:
	return [
		Vector3(min.x, min.y, min.z),
		Vector3(min.x + edge, min.y, min.z),
		Vector3(min.x, min.y + edge, min.z),
		Vector3(min.x + edge, min.y + edge, min.z),
		Vector3(min.x, min.y, min.z + edge),
		Vector3(min.x + edge, min.y, min.z + edge),
		Vector3(min.x, min.y + edge, min.z + edge),
		Vector3(min.x + edge, min.y + edge, min.z + edge)
	]


func trilinear_coordinates(p : Vector3, box : Array[Vector3]) -> Array[float]:
	var volume_per_corner : Array[float] = []
	for i in 8:
		volume_per_corner.append(abs(box[i].x - p.x)*abs(box[i].y - p.y)*abs(box[i].z - p.z))
		
	var coords : Array[float] = []
	var volume = abs(box[7].x - box[0].x)*abs(box[7].y - box[0].y)*abs(box[7].z - box[0].z)
	for i in 8:
		coords.append(volume_per_corner[7 - i] / volume)
		
	return coords
	
func get_gradient(p: Vector3, locations: Array[Vector3], values: Array[float]) -> Vector3:
	var edge: float = locations[0].distance_to(locations[1])
	var gradient := Vector3(0, 0, 0)
	
	# x axis; interpolazione fatta sui valori y e z
	# 1 - 0;  3 - 2; 5 - 4; 7 - 6
	var x_values: Array[float] = [
		values[1] - values[0], 
		values[3] - values[2], 
		values[5] - values[4], 
		values[7] - values[6]
		]
	var first_c: float = p.y - locations[0].y;
	var second_c: float = p.z - locations[0].z;
	var weights: Array[float] = [
		(edge - first_c) * (edge - second_c),
		first_c * (edge - second_c),
		(edge - first_c) * second_c,
		first_c * second_c
		]
	for i in range(4):
		gradient.x -= x_values[i] * weights[i] / pow(edge, 3)

	# y axis x e z
	# 2 - 0; 3 - 1; 6 - 4; 7 - 5
	var y_values: Array[float] = [
		values[2] - values[0], 
		values[3] - values[1], 
		values[6] - values[4], 
		values[7] - values[5]
		]
	first_c = p.x - locations[0].x
	second_c = p.z - locations[0].z
	weights = [
		(edge - first_c) * (edge - second_c),
		first_c * (edge - second_c),
		(edge - first_c) * second_c,
		first_c * second_c
		]
	for i in range(4):
		gradient.y -= y_values[i] * weights[i] / pow(edge, 3)
		
	# z axis
	# 4 - 0; 5 - 1; 6 - 2; 7 - 3
	var z_values: Array[float] = [
		values[4] - values[0], 
		values[5] - values[1], 
		values[6] - values[2], 
		values[7] - values[3]
		]
	first_c = p.x - locations[0].x;
	second_c = p.y - locations[0].y;
	weights = [
		(edge - first_c) * (edge - second_c),
		first_c * (edge - second_c),
		(edge - first_c) * second_c,
		first_c * second_c
		]
	for i in range(4):
		gradient.z -= z_values[i] * weights[i] / pow(edge, 3)
	
	return gradient
	
