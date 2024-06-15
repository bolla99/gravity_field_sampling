shader_type particles;

float rand_from_seed(in uint seed) {
  int k;
  int s = int(seed);
  if (s == 0)
    s = 305420679;
  k = s / 127773;
  s = 16807 * (s - k * 127773) - 2836 * k;
  if (s < 0)
    s += 2147483647;
  seed = uint(s);
  return float(seed % uint(65536)) / 65535.0;
}

uint hash(uint x) {
  x = ((x >> uint(16)) ^ x) * uint(73244475);
  x = ((x >> uint(16)) ^ x) * uint(73244475);
  x = (x >> uint(16)) ^ x;
  return x;
}

vec3[8] get_box(vec3 __min, float _edge) {
	return vec3[8](
		vec3(__min.x, __min.y, __min.z),
		vec3(__min.x + _edge, __min.y, __min.z),
		vec3(__min.x, __min.y + _edge, __min.z),
		vec3(__min.x + _edge, __min.y + _edge, __min.z),
		vec3(__min.x, __min.y, __min.z + _edge),
		vec3(__min.x + _edge, __min.y, __min.z + _edge),
		vec3(__min.x, __min.y + _edge, __min.z + _edge),
		vec3(__min.x + _edge, __min.y + _edge, __min.z + _edge)
	);
}

float[8] trilinear_coordinates(vec3 p, vec3 box[8]) {
	float volume_per_corner[8] = float[8];
	for(int i = 0; i < 8; i++) {
		volume_per_corner[i] = abs(box[i].x - p.x) * abs(box[i].y - p.y) * abs(box[i].z - p.z);
	}
		
	float coords[8] = float[8];
	float volume = abs(box[7].x - box[0].x)*abs(box[7].y - box[0].y)*abs(box[7].z - box[0].z);
	for(int i = 0; i < 8; i++) {
		coords[i] = volume_per_corner[7 - i] / volume;
	}
	return coords;
}

vec3 interpolate(vec3 p, vec3[8] box, vec3[8] values) {
	float weights[8] = trilinear_coordinates(p, box);
	vec3 output = vec3(0.f, 0.f, 0.f);
	for(int i = 0; i < 8; i++) {
		output += values[i]*weights[i];
	}
	return output;
}

vec3 get_gradient(vec3 p, vec3[8] locations, float[8] values) {
    float gradient_edge = distance(locations[0], locations[1]);

   highp vec3 gradient = vec3(0, 0, 0);
    // x axis; interpolazione fatta sui valori y e z
    // 1 - 0;  3 - 2; 5 - 4; 7 - 6
    float x_values[4] = {values[1] - values[0], values[3] - values[2], values[5] - values[4], values[7] - values[6]};
    float first_c = p.y - locations[0].y;
    float second_c = p.z - locations[0].z;
    float[4] weights = {
            (gradient_edge - first_c) * (gradient_edge - second_c),
            first_c * (gradient_edge - second_c),
            (gradient_edge - first_c) * second_c,
            first_c * second_c
    };

    for(int i = 0; i < 4; i++) {
        gradient.x -= x_values[i] * weights[i] / pow(gradient_edge, 3);
    }

    // y axis x e z
    // 2 - 0; 3 - 1; 6 - 4; 7 - 5
    float[4] y_values = {values[2] - values[0], values[3] - values[1], values[6] - values[4], values[7] - values[5]};
    first_c = p.x - locations[0].x;
    second_c = p.z - locations[0].z;
    weights = {
            (gradient_edge - first_c) * (gradient_edge - second_c),
            first_c * (gradient_edge - second_c),
            (gradient_edge - first_c) * second_c,
            first_c * second_c
    };

    for(int i = 0; i < 4; i++) {
        gradient.y -= y_values[i] * weights[i] / pow(gradient_edge, 3);
    }
    // z axis
    // 4 - 0; 5 - 1; 6 - 2; 7 - 3
    float[4] z_values = {values[4] - values[0], values[5] - values[1], values[6] - values[2], values[7] - values[3]};
    first_c = p.x - locations[0].x;
    second_c = p.y - locations[0].y;
    weights = {
            (gradient_edge - first_c) * (gradient_edge - second_c),
            first_c * (gradient_edge - second_c),
            (gradient_edge - first_c) * second_c,
            first_c * second_c
    };

    for(int i = 0; i < 4; i++) {
        gradient.z -= z_values[i] * weights[i] / pow(gradient_edge, 3);
    }
    return gradient;
}

uniform int octree[30000000];
uniform int octree_size = 0;
uniform float gravity_values[30000000];
uniform int gravity_values_size;
uniform vec3 _min;
uniform float edge;
uniform vec3 positions[15000];
uniform vec3 aabb_min;
uniform vec3 aabb_size;

highp vec3 read_octree(vec3 p) {

	vec3 current_min = _min;
	float current_edge = edge;

	vec3 output = vec3(0.0, 0.0, 0.0);

	int counter = -1;
	int i = 0;
	while(true) {
		if(octree[i] <= 0) {
			if(gravity_values_size > 0) {
				vec3 values[8] = vec3[8];
				for(int j = 0; j < 8; j++) {
					values[j] = vec3(
							gravity_values[(-octree[i + j])*3],
							gravity_values[(-octree[i + j])*3 + 1],
							gravity_values[(-octree[i + j])*3 + 2]);
				}
				return interpolate(p, get_box(current_min, current_edge), values);
			} else {
				float values[8];
				for(int j = 0; j < 8; j++) {
					values[j] = intBitsToFloat(octree[i+j]);
				}
				output = get_gradient(p, get_box(current_min, current_edge), values);
				break;
			}
		} else {
			counter++;
			int k = 0;
			int x = 0;
			int y = 0;
			int z = 0;
			if(p.x > (current_min.x + current_edge/2.0)) k += 1;
			if(p.y > (current_min.y + current_edge/2.0)) k += 2;
			if(p.z > (current_min.z + current_edge/2.0)) k += 4;
			if(p.x > (current_min.x + current_edge/2.0)) x += 1;
			if(p.y > (current_min.y + current_edge/2.0)) y += 1;
			if(p.z > (current_min.z + current_edge/2.0)) z += 1;
			i = octree[i + k];
			current_edge = current_edge / 2.0;
			current_min = vec3(current_min.x + current_edge*float(x), current_min.y + current_edge*float(y), current_min.z + current_edge*float(z));
			//current_min = get_box(current_min, current_edge)[k];
			//if(counter == 3) return vec3(current_edge, 0, 0);
		}
	}
	return output;
}


void start() {
	/*
	int positions_size = p_per_edge * p_per_edge * 2 + ((p_per_edge - 2)*4 + 4)*(p_per_edge - 2);
	vec3 positions[2168];
	int next_i = 0;
	float step = spawn_edge/float(p_per_edge - 1);
	for(int i = 0; i < p_per_edge; i++) {
		for(int j = 0; j < p_per_edge; j++) {
			for(int k = 0; k < p_per_edge; k++) {
				bool should_add = false;
			if(i == 0 || i == (p_per_edge - 1) || j == 0 || j == (p_per_edge - 1) || k == 0 || k == (p_per_edge - 1)) {
					positions[next_i] = vec3(-spawn_edge/2.0 + float(i)*step, -spawn_edge/2.0 + float(j)*step, -spawn_edge/2.0 + float(k)*step);
					next_i++;
				}
			}
		}
	}*/
	
	
	/*
	float step = 4.0/9.0;
	vec3 positions[1000];
	for(int i = 0; i < 10; i++) {
		for(int j = 0; j < 10; j++) {
			for(int k = 0; k < 10; k++) {
				positions[i*100 + j*10 + k] = vec3(-2.0 + float(i)*step, -2.0 + float(j)*step, -2.0 + float(k)*step);
			}
		}
	}
	ACTIVE = true;
	if((positions[INDEX].x > aabb_min.x && positions[INDEX].x < aabb_min.x + aabb_size.x)
	&& (positions[INDEX].y > aabb_min.y && positions[INDEX].y < aabb_min.y + aabb_size.y)
	&& (positions[INDEX].z > aabb_min.z && positions[INDEX].z < aabb_min.z + aabb_size.z)) {
		ACTIVE = false;
	}
	*/
	if(CUSTOM.z < 0.5 && CUSTOM.x < 0.5) {
		TRANSFORM[3].xyz = positions[INDEX];
		
		float rx = rand_from_seed(hash(NUMBER + uint(1) + RANDOM_SEED));
		float ry = rand_from_seed(hash(NUMBER + uint(27) + RANDOM_SEED));
		float rz = rand_from_seed(hash(NUMBER + uint(43) + RANDOM_SEED));
		TRANSFORM[3].xyz += vec3(rx - 0.5, ry - 0.5, rz - 0.5);
	}
	

	
	// timer
	CUSTOM.y = 0.0;
	// max fade timeout
	CUSTOM.w = 0.0;
	USERDATA1.x = 700.0;
	USERDATA1.y = 0.0;
	
	// CUSTOM SEMANTIC:
	// x: has collided at least one time; to be set at the end of collision handling
	// y: timer, it starts when particle has collided and it is clamped from USERDATA1.x
	// z: 1 if subparticle
	// w: value from 0 to 1 (it starts changing after first collision); duration depends on
	// USERDATA1.x
}

void process() {

	if(octree_size > 0) {
		VELOCITY += read_octree(TRANSFORM[3].xyz) * DELTA * DELTA;
	}
	
	if(CUSTOM.x > 0.0 && CUSTOM.y < USERDATA1.x) {
		CUSTOM.y += 1.0;
		CUSTOM.w = CUSTOM.y / USERDATA1.x;
	}
	
	if(COLLIDED) {
		// execute only on the first collision, since right after CUSTOM.x is set to 1.0
		//VELOCITY = COLLISION_NORMAL * 0.1;
		VELOCITY = (VELOCITY - COLLISION_NORMAL*dot(VELOCITY, COLLISION_NORMAL))*0.2 + COLLISION_NORMAL*0.1;
		if(CUSTOM.x < 0.5 && USERDATA1.y < 0.5) {
			CUSTOM.x = 10.0;
			mat4 custom_transform = TRANSFORM;
			//custom_transform[3].xyz += COLLISION_NORMAL*0.2;
			//custom_transform[3].xyz = vec3(1.0, 1.0, 1.0);
			custom_transform[3].xyz = TRANSFORM[3].xyz + COLLISION_NORMAL * 0.05;
			bool success = emit_subparticle(
				custom_transform, 
				VELOCITY, 
				COLOR, 
				vec4(1.0, 0.0, 1.0, 0.0), 
				FLAG_EMIT_POSITION | FLAG_EMIT_CUSTOM | FLAG_EMIT_VELOCITY
			);
		}
	}
}