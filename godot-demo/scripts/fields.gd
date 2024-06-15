extends Node


@onready var _fields = $FieldsList.get_children()
@onready var _timer = $FieldTimer

var mat: ShaderMaterial = preload("res://materials/mat_particle.tres")

var paused := false

var wait_time: float

var _next_field = 1
var _p_per_edge: int = 10
var _spawn_edge: float = 4.0


# Called when the node enters the scene tree for the first time.
func _ready():
	# set timer 
	_timer.wait_time = _fields[0].lifetime / _fields.size()
	_timer.timeout.connect(_timer_timeout)
	
	# set material
	for f in _fields:
		f.process_material = mat
		f.sub.process_material = mat
	
	mat.set_shader_parameter("positions", _get_positions(_p_per_edge, _spawn_edge))
	
	for f in _fields:
		f.amount = _get_positions_size(_p_per_edge)
		#f.sub.amount = _get_positions_size(_p_per_edge)
		f.sub.amount = _get_positions_size(_p_per_edge) * 10


func update_material_from_octree(ot: Octree):
	mat.set_shader_parameter("_min", ot.min)
	mat.set_shader_parameter("edge", ot.edge)
	mat.set_shader_parameter("octree_size", ot.octree.size())
	mat.set_shader_parameter("gravity_values_size", ot.gravity_values.size())
	mat.set_shader_parameter("octree", ot.octree)
	mat.set_shader_parameter("gravity_values", ot.gravity_values)


func update_material_from_mesh(mesh: Mesh):
	mat.set_shader_parameter("aabb_min", mesh.get_aabb().position)
	mat.set_shader_parameter("aabb_size", mesh.get_aabb().size)


func update_material_spawn_edge(spawn_edge: float) -> void:
	_spawn_edge = spawn_edge
	mat.set_shader_parameter("positions", _get_positions(_p_per_edge, _spawn_edge))


func update_material_p_per_edge(p_per_edge: int) -> void:
	_p_per_edge = p_per_edge
	var new_amount = _get_positions_size(_p_per_edge)
	for f in _fields:
		f.amount = new_amount
		#f.sub.amount = new_amount;
	mat.set_shader_parameter("positions", _get_positions(_p_per_edge, _spawn_edge))

func get_fields():
	return _fields

# start: restart firs field and start timer
func start():
	if _fields.size() > 0:
		_next_field = 1
		_fields[0].restart()
		_fields[0].sub.restart()
		_fields[0].sub.emitting = false
		_timer.start()
		
func _process(delta):
	_timer.wait_time = wait_time
	
func reset() -> void:
	for f in _fields:
		f.restart()
		f.sub.restart()
		f.emitting = false
		f.sub.emitting = false
		_timer.stop()


# pause and unpause: set speed scale = 0 and pause timer 
func pause() -> void:
	if !paused:
		paused = true
		for f in _fields:
			f.speed_scale = 0
		_timer.paused = true


func unpause() -> void:
	if paused:
		paused = false
		for f in _fields:
			f.speed_scale = 1
		_timer.paused = false


# timout -> start next field
func _timer_timeout():
	if _next_field < _fields.size():
		_fields[_next_field].restart()
		_fields[_next_field].sub.restart()
		_fields[_next_field].sub.emitting = false
		_next_field = (_next_field + 1) % _fields.size()


func _get_positions_size(p_per_edge: int):
	return p_per_edge * p_per_edge * 2 + ((p_per_edge - 2)*4 + 4)*(p_per_edge - 2)


func _get_positions(p_per_edge: int, spawn_edge: float) -> PackedVector3Array:
	var positions: PackedVector3Array = []
	var step: float = spawn_edge/float(p_per_edge - 1)
	for i in p_per_edge:
		for j in p_per_edge:
			for k in p_per_edge:
				if i == 0 or i == (p_per_edge - 1) or j == 0 or j == (p_per_edge - 1) or k == 0 or k == (p_per_edge - 1):
					positions.append(Vector3(
						-spawn_edge/2.0 + float(i)*step,
						-spawn_edge/2.0 + float(j)*step,
						-spawn_edge/2.0 + float(k)*step))
	return positions
