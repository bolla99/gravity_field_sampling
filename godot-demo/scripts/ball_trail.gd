extends MeshInstance3D

var _time: float = 2.5

# Called when the node enters the scene tree for the first time.
func _ready():
	get_tree().create_timer(_time).timeout.connect(queue_free)
	
func _process(delta):
	transparency += 1.0/(_time*60)
