extends Button

@onready var output_label = $output_label

var times : int = 0

@export var next_scene : PackedScene

func _on_pressed():
	if(next_scene): get_tree().change_scene_to_packed(next_scene)
	print_debug("next scene is null")

