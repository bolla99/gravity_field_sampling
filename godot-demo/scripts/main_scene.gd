extends Node3D

var is_gravity_on : bool = false

var octrees: Array[Octree] = []

var raycast_location := Vector3(0, 0, 0)
var goomba_raycast_location := Vector3(0, 0, 0)

var current_gravity: Vector3

# multithreading
var octree_thread := Thread.new()
var octrees_mutex := Mutex.new()
var fields_mutex := Mutex.new()

#region CHILDREN
@onready var camera = $Camera
@onready var cam_distance_slider = $UI/VBoxContainerLeft/CamDistance
@onready var gravity_force_slider = $UI/VBoxContainerLeft/GravityForce
@onready var ball_pos_x_input = $UI/VBoxContainerLeft/BallPosX
@onready var ball_pos_y_input = $UI/VBoxContainerLeft/BallPosY
@onready var ball_pos_z_input = $UI/VBoxContainerLeft/BallPosZ
@onready var gravity_checkbox = $UI/VBoxContainerLeft/GravityCheckbox
@onready var force_label = $UI/VBoxContainerLeft/ForceLabel
@onready var fields = $Fields
@onready var choose_octree_button = $UI/VBoxContainerLeft/ChooseOctreeButton
@onready var choose_mesh_button = $UI/VBoxContainerLeft/ChooseMeshButton
@onready var octree_browser = $UI/VBoxContainerLeft/OctreeBrowser
@onready var mesh_browser = $UI/VBoxContainerLeft/MeshBrowser
@onready var ball = $Ball
@onready var octree_selection = $UI/VBoxContainerLeft/OctreeSelection
@onready var start_particles_button = $UI/VBoxContainerLeft/StartParticlesButton
@onready var pause_particles_button = $UI/VBoxContainerLeft/PauseParticlesButton
@onready var stop_particles_button = $UI/VBoxContainerLeft/StopParticlesButton
@onready var trail_duration = $UI/VBoxContainerLeft/TrailDuration
@onready var wait_time = $UI/VBoxContainerLeft/WaitTime

@onready var fps_label = $UI/VBoxContainerRight/Label
@onready var gravity_label = $UI/GravityLabel

@onready var ball_velocity_x = $UI/VBoxContainerRight/BallVelocityX
@onready var ball_velocity_y = $UI/VBoxContainerRight/BallVelocityY
@onready var ball_velocity_z = $UI/VBoxContainerRight/BallVelocityZ

@onready var spawn_edge = $UI/VBoxContainerRight/SpawnEdge
@onready var p_per_edge = $UI/VBoxContainerRight/PPerEdge

@onready var add_ball_button = $UI/VBoxContainerRight/AddBallButton
@onready var add_pos_x = $UI/VBoxContainerRight/AddPosX
@onready var add_pos_y = $UI/VBoxContainerRight/AddPosY
@onready var add_pos_z = $UI/VBoxContainerRight/AddPosZ
@onready var open_balls_list_button = $UI/VBoxContainerRight/OpenBallsListButton
@onready var balls_list_window = $UI/VBoxContainerRight/BallsListWindow
@onready var balls_list : ItemList = $UI/VBoxContainerRight/BallsListWindow/BallsList
@onready var remove_goombas_button : Button = $UI/VBoxContainerRight/RemoveGoombasButton
@onready var remove_balls_button : Button = $UI/VBoxContainerRight/RemoveBallsButton
@onready var goombas_gravity_flat = $UI/VBoxContainerRight/GoombasGravityFlat
@onready var balls_trace = $UI/VBoxContainerRight/BallsTrace

# SCENES PRELOAD
@onready var ball_scene : PackedScene = preload("res://scenes/Ball.tscn")
@onready var goomba_scene : PackedScene = preload("res://scenes/Goomba.tscn")
@onready var goomba = goomba_scene.instantiate()



#endregion


# READY
func _ready():
	# set starting value for camera distance from UI value
	camera.distance = cam_distance_slider.value
	# connect camera distance value to ui slider
	cam_distance_slider.value_changed.connect(func(value): camera.distance = value)
	
	balls_list_window.close_requested.connect(
		func():
			balls_list_window.visible = false
	)
	open_balls_list_button.pressed.connect(
		func():
			balls_list_window.visible = true
	)
	
	# set starting values for spawn_edge and p_per_edge
	fields.update_material_spawn_edge(spawn_edge.value)
	fields.update_material_p_per_edge(p_per_edge.value)

	#fields material update when spawn edge changes
	# changing it make fields to reset and restart
	spawn_edge.value_changed.connect(
			func(value): 
				fields.reset()
				fields.update_material_spawn_edge(value)
				fields.start()
				)
	
	# fields material update when p per edge changes 
	# changing it make fields to reset and restart
	p_per_edge.value_changed.connect(
			func(value): 
				fields.reset()
				fields.update_material_p_per_edge(value)
				fields.start()
				)
	
	# connect is_gravity_on to ui checkbox
	gravity_checkbox.toggled.connect(
		func(toggled_on):
			# manage main ball
			if not toggled_on:
				ball.cache_velocity() 
			else:
				ball.reset_velocity(Vector3(
						ball_velocity_x.value, 
						ball_velocity_y.value, 
						ball_velocity_z.value))
			# manage other balls
			for c in get_tree().get_root().get_children():
				if c.has_method("cache_velocity") and not toggled_on:
					c.cache_velocity()
				elif c.has_method("reset_velocity") and toggled_on:
					c.reset_velocity(null)
				
			is_gravity_on = toggled_on
	)
	
	# set octree browser (open on button pressed and selected handler)
	choose_octree_button.pressed.connect(
		func(): 
			octree_browser.visible = true
			octree_browser.set_current_dir("res://octree_data"))
	octree_browser.file_selected.connect(_on_octree_browser_file_selected)
	
	# set mesh browser (open on button pressed and selected handler)
	choose_mesh_button.pressed.connect(
		func():
			mesh_browser.visible = true
			mesh_browser.set_current_dir("res://meshes")
	)
	mesh_browser.file_selected.connect(_on_mesh_browser_file_selected)
	
	# set selected octree to update material with current octree
	octree_selection.item_selected.connect(
		func(id):
			octree_thread.wait_to_finish()
			octree_thread.start(
				func():
					fields.update_material_from_octree(octrees[id])
			)
	)
	
	# update material from mesh (mesh bounding box)
	if $Planet/Mesh.mesh:
		fields.update_material_from_mesh($Planet/Mesh.mesh)
	
	# set start particle button to restart particles
	start_particles_button.pressed.connect(
		func():
			fields.reset()
			fields.start()
	)
	
	# set particles pause button
	pause_particles_button.pressed.connect(
		func():
			if fields.paused:
				fields.unpause()
			else:
				fields.pause()
	)
	
	# set particles stop button
	stop_particles_button.pressed.connect(
		func():
			fields.reset()
	)
	
	# set starting values for trail lifetime
	for f in fields.get_fields():
		f.trail_lifetime = trail_duration.value;
		f.sub.trail_lifetime = trail_duration.value;
	
	# connect change trail lifetime 
	trail_duration.value_changed.connect(
		func(value):
			for f in fields.get_fields():
				f.trail_lifetime = value
				f.sub.trail_lifetime = value
	)
	
	# set wait_time starting value and connect 
	fields.wait_time = wait_time.value
	wait_time.value_changed.connect(
		func(value):
			fields.wait_time = value
	)
	
	add_ball_button.pressed.connect(
		func():
			add_ball(ball_scene, Vector3(add_pos_x.value, add_pos_y.value, add_pos_z.value), balls_list)
	)
	
	balls_list.item_clicked.connect(
		func(id, _0, _1):
			var ball = balls_list.get_item_metadata(id)
			balls_list.remove_item(id)
			ball.queue_free()
	)
	
	remove_goombas_button.pressed.connect(
		func():
			var goombas = get_tree().get_nodes_in_group("goomba")
			for goomba in goombas:
				goomba.queue_free()
	)
	
	remove_balls_button.pressed.connect(
		func():
			balls_list.clear()
			var balls = get_tree().get_nodes_in_group("gravitable")
			balls.erase(ball)
			for ball in balls:
				ball.queue_free()
	)
	
	goombas_gravity_flat.toggled.connect(
		func(value):
			var goombas = get_tree().get_nodes_in_group("goomba")
			for goomba in goombas:
				if "flat" in goomba:
					goomba.flat = value
	)
	balls_trace.toggled.connect(
		func(value):
			var balls = get_tree().get_nodes_in_group("gravitable")
			for ball in balls:
				ball.trail_on = value
	)
	
	


# PROCESS (labels update)
func _process(delta):
	# update gravity force slider label
	force_label.text = "gravity force: " + str(gravity_force_slider.value)
	
	# update fps label
	fps_label.text = "FPS: " + str(Engine.get_frames_per_second())
	
	# reset outline
	var balls = get_tree().get_nodes_in_group("gravitable")
	for b in balls:
		if b.has_node("Outline"):
			b.get_node("Outline").visible = false
	
	# set outline on to hovered ball
	var id = balls_list.get_item_at_position(
		get_viewport().get_mouse_position() - (balls_list_window.position as Vector2), true
		)
	if id >= 0:
		balls_list.get_item_metadata(id).get_node("Outline").visible = true

# PHYSICS PROCESS (ball processing + raycasting)
func _physics_process(delta: float):
	# get current selected octree id
	var id = octree_selection.get_selected_id()
	
	if octrees.size() > 0:
		# update goombas gravity
		var goombas = get_tree().get_nodes_in_group("goomba")
		for goomba in goombas:
			goomba.gravity = octrees[id].get_gravity(goomba.position)
	
	if is_gravity_on and octrees.size() > 0:
		# apply to gravitable objects
		var gravitables = get_tree().get_nodes_in_group("gravitable")
		for c in gravitables:
			if c.has_method("apply_gravity"):
				c.apply_gravity(octrees[id].get_gravity(c.position) * gravity_force_slider.value)
				
		# update current gravity label
		gravity_label.text = str(octrees[id].get_gravity(ball.position))
		
		ball_velocity_x.value = ball.linear_velocity.x
		ball_velocity_y.value = ball.linear_velocity.y
		ball_velocity_z.value = ball.linear_velocity.z
		
		ball_pos_x_input.value = ball.position.x
		ball_pos_y_input.value = ball.position.y
		ball_pos_z_input.value = ball.position.z
	else:
		#ball.linear_velocity.x = ball_velocity_x.value
		#ball.linear_velocity.y = ball_velocity_y.value
		#ball.linear_velocity.z = ball_velocity_z.value
		
		# this prevents ball to move even if its
		# linear velocity is greater than 0
		ball.position.x = ball_pos_x_input.value
		ball.position.y = ball_pos_y_input.value
		ball.position.z = ball_pos_z_input.value
		
		
	var to: Vector3 = camera.position + camera.project_ray_normal(get_viewport().get_mouse_position()) * 100
	var space_state: PhysicsDirectSpaceState3D = get_world_3d().direct_space_state
	var query := PhysicsRayQueryParameters3D.create(camera.position, to, 0x00000001)
	var result = space_state.intersect_ray(query)
	if result:
		var r: float = 0.01
		DebugDraw3D.draw_sphere(result.position + result.normal * r, r)
		raycast_location = result.position + result.normal * ball.get_node("Collider").shape.radius*1.1
		goomba_raycast_location = result.position + result.normal * goomba.get_node("Collider").shape.radius
	
	#if ball.state and ball.state.get_contact_count() > 0:
	#	var norm: Vector3 = ball.state.get_contact_local_normal(0).normalized()
	#	var flat_gravity = current_gravity - (current_gravity.dot(norm))*norm
	#	var _sc = DebugDraw3D.new_scoped_config().set_thickness(0.01)
	#	DebugDraw3D.draw_arrow(ball.position, ball.position + 0.2*flat_gravity, Color(0, 1, 0), 0.01)



# UNHANDLED INPUT
func _unhandled_input(event):
	# move trackball camera with mouse
	if event is InputEventMouseMotion:
		if Input.is_mouse_button_pressed(MOUSE_BUTTON_LEFT):
			camera.h += -event.relative.x / 100
			camera.v += -event.relative.y / 100
	elif event is InputEventKey:
		if event.keycode == KEY_K:
			ball_pos_x_input.value = raycast_location.x
			ball_pos_y_input.value = raycast_location.y
			ball_pos_z_input.value = raycast_location.z
			ball.linear_velocity = Vector3(0, 0, 0)
		elif event.keycode == KEY_A and event.is_pressed():
			add_ball(ball_scene, raycast_location, balls_list)
		elif event.keycode == KEY_G and event.is_pressed():
			var goomba = goomba_scene.instantiate()
			goomba.position = raycast_location
			get_tree().get_root().add_child(goomba)

# instantiate ball_scene at pos
# if itemlist is provided add item to it with node as metadata
func add_ball(ball_scene, pos, itemlist):
		var ball = ball_scene.instantiate()
		ball.position = pos
		get_tree().get_root().add_child(ball)
		if itemlist:
			var id = itemlist.add_item(ball.name)
			itemlist.set_item_metadata(id, ball)
	

# OCTREE LOAD
func _on_octree_browser_file_selected(path):
	octree_thread.wait_to_finish()
	octree_thread.start(
		func():
			# load and append new octree
			var octree = Octree.new()
			octree.loadFromFile(path)
			octrees_mutex.lock()
			octrees.append(octree)
			octrees_mutex.unlock()
			# if it is the first octree to be added, set material, then when others octrees are 
			# added, material will be updated by changing selected octree by option button callback 
			# then start the timer and start the first field; following fields will be started
			# by timer callback
			if octrees.size() == 1:
				fields.update_material_from_octree(octrees[0])

			# add new element to octree selection i.e. OptionButton
			octree_selection.add_item(path.get_file())
	)

func _on_mesh_browser_file_selected(path):
	print_debug("changing mesh to: ", path)
	var new_mesh = load(path)
	$Planet/Mesh.mesh = new_mesh
	if new_mesh.has_method('create_trimesh_shape'):
		$Planet/Collider.shape = new_mesh.create_trimesh_shape()
	fields.update_material_from_mesh($Planet/Mesh.mesh)
	var mesh_name: String = path.get_file().get_slice(".", 0)
	var particle_collider_path = "res://other_data/" + mesh_name + "_particle_collider.exr"
	$Planet/GPUParticlesCollisionSDF3D.texture = load(particle_collider_path)
