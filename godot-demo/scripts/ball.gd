extends RigidBody3D


var state: PhysicsDirectBodyState3D

var _last_velocity: Vector3

@export var trail_on: bool = false
@export var trail_frequency: float = 1

@onready var trail: PackedScene = preload("res://scenes/BallTrail.tscn")

func _integrate_forces(state):
	self.state = state

func apply_gravity(gravity: Vector3) -> void:
		self.apply_force(gravity)

func cache_velocity():
	_last_velocity = linear_velocity
	self.linear_velocity = Vector3(0, 0, 0)
	
func reset_velocity(new_velocity):
	if not new_velocity:
		self.linear_velocity = _last_velocity
	else:
		self.linear_velocity = new_velocity

func _physics_process(delta):
	if trail_on:
		_draw_trail()


func _draw_trail():
	await get_tree().create_timer(trail_frequency).timeout
	var _trail = trail.instantiate()
	_trail.position = self.position
	_trail.visible = self.visible
	get_tree().get_root().add_child(_trail)
