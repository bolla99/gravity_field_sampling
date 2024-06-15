extends Node3D
class_name Trackball

# h and v are meant to be set by input

# horizontal rotation (around up axis)
@export var h : float = 0

# vertical movement (aroung right axis)
@export var v : float = 0

# rotation around front axis
@export var tilt : float = 0

# distance along z axis
@export var distance : float = 10

# camera position when distance = 0; can represents the position of the object tracked by camera
@export var tracking_pos = Vector3(0, 0, 0)

func _ready():
	position = tracking_pos + distance*transform.basis.z

# Called when the node enters the scene tree for the first time.
func _physics_process(delta):
	# reset rotation
	transform = transform.orthonormalized()
	rotation = Vector3.ZERO
	
	# rotate
	rotate(transform.basis.y, h) # around right axis
	rotate(transform.basis.x, v) # around up axis
	rotate(transform.basis.z, tilt) # aroung front axis
	
	# set position
	position = tracking_pos + distance*transform.basis.z
