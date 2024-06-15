shader_type spatial;
render_mode particle_trails;

varying vec4 c;
varying float d;


void vertex() {
	c = INSTANCE_CUSTOM;
	d = distance(VERTEX.xyz, c.xyz);
}

void fragment() {
	// Called for every pixel the material is visible on.
	ALBEDO = vec3(1, 0, 0);
	ALPHA = 1.0 - (d/10.0);
}

//void light() {
	// Called for every pixel for every light affecting the material.
	// Uncomment to replace the default light processing function with this one.
//}
