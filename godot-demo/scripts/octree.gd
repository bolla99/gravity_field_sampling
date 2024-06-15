extends RefCounted
class_name Octree

var min : Vector3 = Vector3(0, 0, 0)
var edge : float = 0
var octree : PackedInt32Array = []
var gravity_values : PackedFloat32Array = []

enum {TYPE_GRAVITY, TYPE_POTENTIAL}
	
func get_type() -> int:
	if gravity_values.size() > 0:
		return TYPE_GRAVITY
	else:
		return TYPE_POTENTIAL

func loadFromFile(path: String) -> void: 
	var _output = util.get_universal_octree(path, octree, gravity_values)
	min = _output[0]
	edge = _output[1]
	
func get_gravity(p: Vector3) -> Vector3:
	if get_type() == TYPE_GRAVITY:
		return util.read_octree(p, octree, gravity_values, min, edge)
	else:
		return util.read_octree(p, octree, [], min, edge)
