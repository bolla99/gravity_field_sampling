extends CharacterBody3D

var gravity: Vector3 = Vector3(1, 0, 0)
var flat: bool = false

func _process(delta):
	# normalize gravity 
	gravity = gravity.normalized()
	
	# get collision data
	var collision = move_and_collide(gravity * 0.1, true)
	
	if collision and flat:
		var normal = collision.get_normal()
		var flat_gravity = (gravity - gravity.dot(normal)*normal).normalized()
		# neglect if flat gravity magnitude is less than 1% of gravity magnitude
		#if flat_gravity.length() < gravity.length()*0.02:
		#	flat_gravity = gravity
		rotate(basis.y.cross(flat_gravity).normalized(), basis.y.angle_to(flat_gravity))
	else:
		rotate(basis.y.cross(gravity.normalized()).normalized(), basis.y.angle_to(gravity.normalized()))
