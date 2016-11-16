/* 
 * File:   raytrace.c
 * Author: Matthew
 *
 * Created on November 13, 2016, 1:54 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define SPHERE 0
#define PLANE 1
#define MAX_COLOR 255
#define PI acos(-1.0)
#define MAX_RECURSION_DEPTH 7

// struct for camera data
typedef struct {
    double position[3];
    double width;
    double height;
} Camera;

// struct that stores object data
typedef struct {
    int kind;
    double position[3];
    double diffuse_color[3];
    double specular_color[3];
    double reflectivity;
    double refractivity;
    double ior;
    union {
        struct {
            double radius;
        } sphere;
        struct {
            double normal[3];
        } plane;
    };
} Object;

// struct that stores light data
typedef struct {
    double color[3];
    double position[3];
    double direction[3];
    double theta;
    double radial_a2;
    double radial_a1;
    double radial_a0;
    double angular_a0;
} Light;

// struct that stores color data
typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} Pixel;

// inline function for square
static inline double sqr(double v) {
    return v*v;
}

// inline function for normalizing a vector
static inline void normalize(double* v) {
    double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
    
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
}

typedef double* v3;
// add vector "a" with vector "b" and puts result in vector "c"
static inline void v3_add(v3 a, v3 b, v3 c) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

// subtract vector "b" from vector "a" and puts result in vector "c"
static inline void v3_subtract(v3 a, v3 b, v3 c) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

// scale vector "a" by value "s" and puts result in vector "c"
static inline void v3_scale(v3 a, double s, v3 c) {
    c[0] = s * a[0];
    c[1] = s * a[1];
    c[2] = s * a[2];
}

// dot product between vector "a" and vector "b"
static inline double v3_dot(v3 a, v3 b) {
    return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

void read_scene(FILE*);
void raycast(FILE*);
double* illuminate(double*, double*, Object*);
double* shade(Object*, double*, double*, int, double);
double* direct_shade(Object*, double*, double*, double*, double*);
double* refraction(double, double, double*, double*);
int compare_objects(Object*, Object*);
double diffuse(double, double, double);
double specular(double, double, double, double);
double frad(double, double, double, double);
double fang(double*, double*, double, double);
void parse_camera(FILE*);
void parse_sphere(FILE*, Object*);
void parse_plane(FILE*, Object*);
void parse_light(FILE*, Light*);
double sphere_intersect(double*, double*, double*, double);
double plane_intersect(double*, double*, double*, double*);
void skip_ws(FILE*);
void expect_c(FILE*, int);
int next_c(FILE*);
char* next_string(FILE*);
double next_number(FILE*);
double* next_vector(FILE*);
double clamp(double, double, double);
void output_p6(FILE*);

// initialize input file line counter
int line = 1;
int w;
int h;

// create arrays for storing objects and pixels
Object** objects;
Light** lights;
Pixel* pixmap;

// default camera
Camera camera;

int main(int argc, char** argv) {
    // check for correct number of inputs
    if (argc != 5) {
        fprintf(stderr, "Error: Arguments should be in format: 'width' 'height' 'source' 'dest'.\n");
        exit(1);
    }
    
    // check that 'width' is a non-zero number
    w = atoi(argv[1]);
    if (w == 0) {
        fprintf(stderr, "Error: Argument 1, 'width' must be a non-zero integer.\n");
        exit(1);
    }
    
    // check that 'height' is a non-zero number
    h = atoi(argv[2]);
    if (h == 0) {
        fprintf(stderr, "Error: Argument 2, 'height' must be a non-zero integer.\n");
        exit(1);
    }
    
    // open and check input file
    FILE* json = fopen(argv[3], "r");
    if (json == NULL) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", argv[3]);
        exit(1);
    }
    
    // allocate space for 128 objects/lights
    objects = malloc(sizeof(Object*)*129);
    lights = malloc(sizeof(Light*)*129);
    
    read_scene(json);
    
    raycast(json);
    
    // open and check output file location
    FILE* output = fopen(argv[4], "w");
    if (output == NULL) {
        fprintf(stderr, "Error: Could not create file '%s'.\n", argv[4]);
        exit(1);
    }
    
    // write pixel data to output file then close it
    output_p6(output);
    fclose(output);
    
    for (int i=0; objects[i] != NULL; i++)
        free(objects[i]);
    free(objects);
    
    for (int i=0; lights[i] != NULL; i++)
        free(lights[i]);
    free(lights);
    
    return (EXIT_SUCCESS);
}

// reads object data from a json file
void read_scene(FILE* json) {
    // next char in file
    int c;
    
    // look for start of json file ([)
    skip_ws(json);
    expect_c(json, '[');
    skip_ws(json);
    
    // initialize object array index
    int object_index = 0;
    int light_index = 0;
    // read all objects found in file
    while (1) {
        c = next_c(json);
        // check for empty scene
        if (c == ']') {
            fprintf(stderr, "Error: This scene is empty.\n");
            fclose(json);
            objects[object_index] = NULL;
            return;
        }
        
        // checks for start of an object
        if (c == '{') {
            skip_ws(json);
            
            // check that 'type' field is first
            char* key = next_string(json);
            if (strcmp(key, "type") != 0) {
                fprintf(stderr, "Error: Expected 'type' key. (Line %d)\n", line);
                exit(1);
            }
            
            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            
            // check what object is being read and calls a function depending on what it is
            char* value = next_string(json);
            if (strcmp(value, "camera") == 0) {
                parse_camera(json);
            } else if (strcmp(value, "sphere") == 0) {
                // allocate space for an object
                objects[object_index] = malloc(sizeof(Object));
                parse_sphere(json, objects[object_index]);
                object_index++;
            } else if (strcmp(value, "plane") == 0) {
                // allocate space for an object
                objects[object_index] = malloc(sizeof(Object));
                parse_plane(json, objects[object_index]);
                object_index++;
            } else if (strcmp(value, "light") == 0) {
                lights[light_index] = malloc(sizeof(Light));
                parse_light(json, lights[light_index]);
                light_index++;
            } else {
                fprintf(stderr, "Error: Unknown type '%s'. (Line %d)\n", value, line);
                exit(1);
            }
            
            // check for more objects
            skip_ws(json);
            c = next_c(json);
            if (c == ',') {
                skip_ws(json);
            } else if (c == ']') {
                objects[object_index] = NULL;
                lights[light_index] = NULL;
                fclose(json);
                return;
            } else {
                fprintf(stderr, "Error: Expecting ',' or ']'. (Line %d)\n", line);
                exit(1);
            }
            
            // check if scene has too many objects in it
            if ((object_index + light_index) == 129) {
                objects[object_index] = NULL;
                lights[light_index] = NULL;
                fclose(json);
                fprintf(stderr, "Error: Too many objects in file.\n");
                return;
            }
        } else {
            fprintf(stderr, "Error: Expecting '{'. (Line %d)\n", line);
            exit(1);
        }
    }
}

// raycast objects found after reading json file
void raycast(FILE* json) {
    double cx = camera.position[0];
    double cy = camera.position[1];
    double cz = camera.position[2];
    
    // calculate pixel height, width
    double pixheight = camera.height/h;
    double pixwidth = camera.width/w;
    
    // allocate space for number of pixels needed
    pixmap = malloc(sizeof(Pixel)*h*w);
    // initialize pixmap index
    int pixmap_index = 0;
    
    // got through each spot pixel by pixel to see what color it should be
    for (int y=h; y>0; y--) {
        for (int x=0; x<w; x++) {
            // ray origin
            double Ro[3] = {cx, cy, cz};
            // ray destination
            double Rd[3] = {cx - (camera.width/2) + pixwidth*(x + 0.5),
                            cy - (camera.height/2) + pixheight*(y + 0.5),
                            1};
            normalize(Rd);
            
            double best_t = INFINITY;
            Object* closest_object;
            // look for intersection of an object 
            for (int i=0; objects[i] != NULL; i++) {
                double t = 0;
                switch (objects[i]->kind) {
                    case SPHERE:
                        t = sphere_intersect(Ro, Rd, objects[i]->position, objects[i]->sphere.radius);
                        break;
                    case PLANE:
                        t = plane_intersect(Ro, Rd, objects[i]->position, objects[i]->plane.normal);
                        break;
                    default:
                        fprintf(stderr, "Error: Unknown object.\n");
                        exit(1);
                }
                
                // save object if it intersects closer to the camera
                if (t > 0 && t < best_t) {
                    best_t = t;
                    closest_object = objects[i];
                }
            }
                
            if (best_t > 0 && best_t != INFINITY) {
                double Ron[3];
                v3_scale(Rd, best_t, Ron);
                v3_add(Ron, Ro, Ron);
                
                double* color = shade(closest_object, Rd, Ron, 0, 1.0);
                
                // save pixel to pixmap buffer
                pixmap[pixmap_index].r = (unsigned char) (clamp(color[0], 0, 1)*MAX_COLOR);
                pixmap[pixmap_index].g = (unsigned char) (clamp(color[1], 0, 1)*MAX_COLOR);
                pixmap[pixmap_index].b = (unsigned char) (clamp(color[2], 0, 1)*MAX_COLOR);
            } else {
                pixmap[pixmap_index].r = 0;
                pixmap[pixmap_index].g = 0;
                pixmap[pixmap_index].b = 0;
            }

            pixmap_index++;
        }
    }
}

// adds lighting to the pixel found above at best_t from the object closest to the camera
double* illuminate(double* Rd, double* Ro, Object* closest_object) {
    double* current_color = malloc(sizeof(double)*3);
    current_color[0] = 0;
    current_color[1] = 0;
    current_color[2] = 0;

    double pixel_position[3];
    pixel_position[0] = Ro[0];
    pixel_position[1] = Ro[1];
    pixel_position[2] = Ro[2];
    
    // find position of pixel in space and get vector from it to the camera
    double obj_to_cam[3];
    v3_subtract(camera.position, pixel_position, obj_to_cam);
    normalize(obj_to_cam);

    // find and save object normal
    double* N = malloc(sizeof(double)*3);
    if (closest_object->kind == SPHERE)
        v3_subtract(pixel_position, closest_object->position, N);
    else
        N = closest_object->plane.normal;
    normalize(N);

    // look through light to find the ones that influence this pixel
    Light* current_light;
    Object* current_object;
    double new_t = 0;
    for (int j=0; lights[j] != NULL; j++) {
        current_light = lights[j];

        // find vector from light to the object
        double light_to_obj[3];
        v3_subtract(pixel_position, current_light->position, light_to_obj);
        normalize(light_to_obj);

        // find vector from the object to the light
        double obj_to_light[3];
        v3_subtract(current_light->position, pixel_position, obj_to_light);
        normalize(obj_to_light);
        
        // find distance from the light to the object
        double dl = sqrt(sqr(pixel_position[0]-current_light->position[0]) + sqr(pixel_position[1]-current_light->position[1]) + sqr(pixel_position[2]-current_light->position[2]));

        // boolean that tells if object is in a shadow
        int in_shadow = 0;
        for (int k=0; objects[k] != NULL; k++) {
            current_object = objects[k];
            // skip checking for intersection with the object already being looked at
            if (compare_objects(current_object, closest_object))
                continue;

            double* offset = malloc(sizeof(double)*3);
            v3_scale(obj_to_light, 0.00001, offset);
            v3_add(offset, pixel_position, offset);
            // find new intersection between light and each object looking for one that is closer to the light 
            // and casts a shadow on the closest object to the camera at this pixel
            switch (current_object->kind) {
                case SPHERE:
                    new_t = sphere_intersect(offset, obj_to_light, current_object->position, current_object->sphere.radius);
                    break;
                case PLANE:
                    new_t = plane_intersect(offset, obj_to_light, current_object->position, current_object->plane.normal);
                    break;
                default:
                    fprintf(stderr, "Error: Unknown object.\n");
                    exit(1);
            }

            // if there is a closer object to the light, then there is a shadow
            if (new_t > 0 && new_t <= dl) {
                in_shadow = 1;
                break;
            }
        }

        // if not in a shadow then modify the color of the pixel using diffuse/specular components, and the color of the light
        if (in_shadow == 0) {
            double R[3];
            v3_scale(N, 2*v3_dot(N, light_to_obj), R);
            v3_subtract(light_to_obj, R, R);
            normalize(R);

            // calculate (N*L)
            double diffuse_component = v3_dot(N, obj_to_light);
            // calculate (R*V)
            double specular_component = v3_dot(R, obj_to_cam);

            // normalize light direction
            double* ld = current_light->direction;
            normalize(ld);

            // find diffuse color
            double diffuse_color[3];
            diffuse_color[0] = diffuse(current_light->color[0], closest_object->diffuse_color[0], diffuse_component);
            diffuse_color[1] = diffuse(current_light->color[1], closest_object->diffuse_color[1], diffuse_component);
            diffuse_color[2] = diffuse(current_light->color[2], closest_object->diffuse_color[2], diffuse_component);

            // find specular color
            double specular_color[3];
            specular_color[0] = specular(current_light->color[0], closest_object->specular_color[0], diffuse_component, specular_component);
            specular_color[1] = specular(current_light->color[1], closest_object->specular_color[1], diffuse_component, specular_component);
            specular_color[2] = specular(current_light->color[2], closest_object->specular_color[2], diffuse_component, specular_component);

            // get attenuation values
            double rad = frad(dl, current_light->radial_a0, current_light->radial_a1, current_light->radial_a2);
            double ang = fang(ld, light_to_obj, current_light->angular_a0, current_light->theta);

            // modify color if pixel to reflect changes from illumination
            current_color[0] += rad * ang * (diffuse_color[0] + specular_color[0]);
            current_color[1] += rad * ang * (diffuse_color[1] + specular_color[1]);
            current_color[2] += rad * ang * (diffuse_color[2] + specular_color[2]);
        }
    }
    
    return current_color;
}

double* shade(Object* current_object, double* Rd, double* Ro, int level, double ior) {
    double* color = malloc(sizeof(double)*3);
    color[0] = 0;
    color[1] = 0;
    color[2] = 0;
    
    if (level == MAX_RECURSION_DEPTH) {
        return color;
    }
    
    v3_add(color, illuminate(Rd, Ro, current_object), color);
    
    double* N = malloc(sizeof(double)*3);
    if (current_object->kind == SPHERE)
        v3_subtract(Ro, current_object->position, N);
    else
        N = current_object->plane.normal;
    normalize(N);

    if (current_object->reflectivity > 0) {
        double reflected_ray[3];
        v3_scale(N, 2*v3_dot(N, Rd), reflected_ray);
        v3_subtract(Rd, reflected_ray, reflected_ray);
        normalize(reflected_ray);

        double best_t = INFINITY;
        Object* closest_object;
        // look for intersection of an object 
        for (int i=0; objects[i] != NULL; i++) {
            // skip checking for intersection with the object already being looked at
            if (compare_objects(current_object, objects[i]))
                continue;

            double* offset = malloc(sizeof(double)*3);
            v3_scale(reflected_ray, 0.00001, offset);
            v3_add(offset, Ro, offset);
            
            double t = 0;
            switch (objects[i]->kind) {
                case SPHERE:
                    t = sphere_intersect(offset, reflected_ray, objects[i]->position, objects[i]->sphere.radius);
                    break;
                case PLANE:
                    t = plane_intersect(offset, reflected_ray, objects[i]->position, objects[i]->plane.normal);
                    break;
                default:
                    fprintf(stderr, "Error: Unknown object.\n");
                    exit(1);
            }

            // save object if it intersects closer to the camera
            if (t > 0 && t < best_t) {
                best_t = t;
                closest_object = objects[i];
            }
        }

        if (best_t > 0 && best_t != INFINITY) {
            double* Ron = malloc(sizeof(double)*3);
            v3_scale(reflected_ray, best_t, Ron);
            v3_add(Ro, Ron, Ron);

            double* reflected_color = malloc(sizeof(double)*3);
            v3_scale(shade(closest_object, reflected_ray, Ron, level+1, current_object->ior), current_object->reflectivity, reflected_color);
            
            v3_add(color, direct_shade(closest_object, Ron, reflected_color, reflected_ray, Ron), color);
            v3_add(color, reflected_color, color);
        }
    }
    
    if (current_object->refractivity > 0) {
        double* refracted_ray = refraction(ior, current_object->ior, Rd, N);
        
        double* offset = malloc(sizeof(double)*3);
        v3_scale(refracted_ray, 0.00001, offset);
        v3_add(offset, Ro, offset);
        
        double d;
        switch (current_object->kind) {
            case SPHERE:
                d = sphere_intersect(offset, refracted_ray, current_object->position, current_object->sphere.radius);
                break;
            case PLANE:
                d = plane_intersect(offset, refracted_ray, current_object->position, current_object->plane.normal);
                break;
            default:
                fprintf(stderr, "Error: Unknown object.\n");
                exit(1);
        }
        
        double* back = malloc(sizeof(double)*3);
        v3_scale(refracted_ray, d, back);
        v3_add(Ro, back, back);
        
        double* back_normal = malloc(sizeof(double)*3);
        if (current_object->kind == SPHERE) {
            v3_subtract(back, current_object->position, back_normal);
        } else if (current_object->kind == PLANE) {
            back_normal = current_object->plane.normal;
            v3_scale(back_normal, -1, back_normal);
        }
        normalize(back_normal);
        
        double* next_ray = refraction(current_object->ior, 1.0, refracted_ray, back_normal);
        
        double best_t = INFINITY;
        Object* closest_object;
        // look for intersection of an object 
        for (int i=0; objects[i] != NULL; i++) {
            // skip checking for intersection with the object already being looked at
            if (compare_objects(current_object, objects[i]))
                continue;

            double* offset = malloc(sizeof(double)*3);
            v3_scale(refracted_ray, 0.00001, offset);
            v3_add(offset, Ro, offset);
            
            double t = 0;
            switch (objects[i]->kind) {
                case SPHERE:
                    t = sphere_intersect(offset, next_ray, objects[i]->position, objects[i]->sphere.radius);
                    break;
                case PLANE:
                    t = plane_intersect(offset, next_ray, objects[i]->position, objects[i]->plane.normal);
                    break;
                default:
                    fprintf(stderr, "Error: Unknown object.\n");
                    exit(1);
            }

            // save object if it intersects closer to the camera
            if (t > 0 && t < best_t) {
                best_t = t;
                closest_object = objects[i];
            }
        }
        
        if (best_t > 0 && best_t != INFINITY) {
            double* Ron = malloc(sizeof(double)*3);
            v3_scale(next_ray, best_t, Ron);
            v3_add(Ro, Ron, Ron);

            double* refracted_color = malloc(sizeof(double)*3);
            v3_scale(shade(closest_object, next_ray, Ron, level+1, current_object->ior), current_object->reflectivity, refracted_color);
            
            v3_add(color, direct_shade(closest_object, Ron, refracted_color, next_ray, back), color);
            v3_add(color, refracted_color, color);
        }
    }
    
    return color;
}

double* direct_shade(Object* current_object, double* pixel_position, double* light_color, double* light_direction, double* light_position) {
    double* color = malloc(sizeof(double)*3);
    color[0] = 0;
    color[1] = 0;
    color[2] = 0;
    
    double* N = malloc(sizeof(double)*3);
    if (current_object->kind == SPHERE)
        v3_subtract(pixel_position, current_object->position, N);
    else
        N = current_object->plane.normal;
    normalize(N);
    
    // find vector from the object to the light
    double obj_to_light[3];
    v3_subtract(light_position, pixel_position, obj_to_light);
    normalize(obj_to_light);
    
    double obj_to_cam[3];
    v3_subtract(camera.position, pixel_position, obj_to_cam);
    normalize(obj_to_cam);
    
    double R[3];
    v3_scale(N, 2*v3_dot(N, light_direction), R);
    v3_subtract(light_direction, R, R);
    normalize(R);
    
    // calculate (N*L)
    double diffuse_component = v3_dot(N, obj_to_light);
    // calculate (R*V)
    double specular_component = v3_dot(R, obj_to_cam);
    
    // find diffuse color
    double diffuse_color[3];
    diffuse_color[0] = diffuse(light_color[0], current_object->diffuse_color[0], diffuse_component);
    diffuse_color[1] = diffuse(light_color[1], current_object->diffuse_color[1], diffuse_component);
    diffuse_color[2] = diffuse(light_color[2], current_object->diffuse_color[2], diffuse_component);

    // find specular color
    double specular_color[3];
    specular_color[0] = specular(light_color[0], current_object->specular_color[0], diffuse_component, specular_component);
    specular_color[1] = specular(light_color[1], current_object->specular_color[1], diffuse_component, specular_component);
    specular_color[2] = specular(light_color[2], current_object->specular_color[2], diffuse_component, specular_component);

    // modify color if pixel to reflect changes from illumination
    color[0] += (diffuse_color[0] + specular_color[0]);
    color[1] += (diffuse_color[1] + specular_color[1]);
    color[2] += (diffuse_color[2] + specular_color[2]);
    
    return color;
}

double* refraction(double outer_ior, double inner_ior, double* ray_in, double* N) {
    double* ray_out = malloc(sizeof(double)*3);
    
    double eta, c1, cs2; 
    eta = outer_ior / inner_ior;
    c1 = v3_dot(ray_in, N) * -1; // cos(theta)
    cs2 = 1 - eta * eta * (1 - c1 * c1); // cos^2(phi)

    if (cs2 < 0) { // total internal reflection
      ray_out[0] = 0.0;
      ray_out[1] = 0.0;
      ray_out[2] = 0.0;
      return ray_out;
    }

    // Linear combination:
    // transmittedRay = eta * incomingRay + (eta * c1 - sqrt(cs2)) * normal
    double* temp = malloc(3 * sizeof(double));
    v3_scale(ray_in, eta, temp);
    v3_scale(N, eta * c1 - sqrt(cs2), ray_out);
    v3_add(temp, ray_out, ray_out);

    normalize(ray_out);

    return ray_out;
}

// compares two objects based on field values because c cant compare objects directly
int compare_objects(Object* a, Object* b) {
    if (a->kind == b->kind && 
        a->diffuse_color[0] == b->diffuse_color[0] &&
        a->diffuse_color[1] == b->diffuse_color[1] &&
        a->diffuse_color[2] == b->diffuse_color[2] &&
        a->specular_color[0] == b->specular_color[0] &&
        a->specular_color[1] == b->specular_color[1] &&
        a->specular_color[2] == b->specular_color[2] &&
        a->position[0] == b->position[0] &&
        a->position[1] == b->position[1] &&
        a->position[2] == b->position[2]) {
        if (a->kind == SPHERE && 
            a->sphere.radius == b->sphere.radius) {
            return 1;
        } else if (a->kind == PLANE &&
                   a->plane.normal[0] == b->plane.normal[0] &&
                   a->plane.normal[1] == b->plane.normal[1] &&
                   a->plane.normal[2] == b->plane.normal[2]) {
            return 1;
        }
    }
    
    return 0;
}

// calculates diffuse component of illumination formula
double diffuse(double light_value, double object_value, double diffuse_component) {
    if (diffuse_component > 0)
        return light_value * object_value * diffuse_component;
    
    return 0;
}

// calculates specular component of illumination formula
double specular(double light_value, double object_value, double diffuse_component, double specular_component) {
    if (specular_component > 0 && diffuse_component > 0)
        return light_value * object_value * pow(specular_component, 20);
    
    return 0;
}

// calculates radial attenuation
double frad(double dl, double a0, double a1, double a2) {
    if (dl == INFINITY)
        return 1;
    
    return 1 / (a2*sqr(dl) + (a1*dl) + a0);
}

// calculates angular attenuation
double fang(double* light_direction, double* light_to_obj, double a0, double theta) {
    if (light_direction[0] == 0 && light_direction[1] == 0 && light_direction[2] == 0)
        return 1;
    
    if (v3_dot(light_to_obj, light_direction) < (cos(theta / 180 * PI)) * (180 / PI))
        return 0;
    
    return pow(v3_dot(light_to_obj, light_direction), a0);
}

// reads camera data from json file
void parse_camera(FILE* json) {
    int c;
    skip_ws(json);
    
    camera.position[0] = 0;
    camera.position[1] = 0;
    camera.position[2] = 0;
    
    camera.width = 1;
    camera.height = 1;
    
    // checks fields in camera
    while (1) {
        c = next_c(json);
        if (c == '}') {
            break;
        } else if (c == ',') {
            skip_ws(json);
            char* key = next_string(json);
            
            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            
            double value = next_number(json);
            // checks that field can be in camera and sets 'w' or 'h' if found
            if (strcmp(key, "width") == 0) {
                camera.width = value;
            } else if (strcmp(key, "height") == 0) {
                camera.height = value;
            } else {
                fprintf(stderr, "Error: Unknown property '%s' for 'camera'. (Line %d)\n", key, line);
                exit(1);
            }
        }
    }
}

// gets sphere information and stores it into an object
void parse_sphere(FILE* json, Object* object) {
    int c;
    
    // used to check that all fields for a sphere are present
    int hasradius = 0;
    int hascolor = 0;
    int hasposition = 0;
    
    // set object kind to sphere
    object->kind = SPHERE;
    object->reflectivity = 0;
    object->refractivity = 0;
    object->ior = 1;
    
    skip_ws(json);
    
    // check fields for this sphere
    while(1) {
        c = next_c(json);
        if (c == '}') {
            break;
        } else if (c == ',') {
            skip_ws(json);
            char* key = next_string(json);

            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            
            // set values for this sphere depending on what key was read and sets its 'boolean' to reflect the found field
            if (strcmp(key, "radius") == 0) {
                object->sphere.radius = next_number(json);
                hasradius = 1;
            } else if (strcmp(key, "reflectivity") == 0) {
                object->reflectivity = clamp(next_number(json), 0, 1);
            } else if (strcmp(key, "refractivity") == 0) {
                object->refractivity = clamp(next_number(json), 0, 1);
            } else if (strcmp(key, "ior") == 0) {
                object->ior = next_number(json);
            } else if (strcmp(key, "diffuse_color") == 0) {
                double* value = next_vector(json);
                object->diffuse_color[0] = value[0];
                object->diffuse_color[1] = value[1];
                object->diffuse_color[2] = value[2];
                hascolor = 1;
            } else if (strcmp(key, "specular_color") == 0) {
                double* value = next_vector(json);
                object->specular_color[0] = value[0];
                object->specular_color[1] = value[1];
                object->specular_color[2] = value[2];
            } else if (strcmp(key, "position") == 0) {
                double* value = next_vector(json);
                object->position[0] = value[0];
                object->position[1] = value[1];
                object->position[2] = value[2];
                hasposition = 1;
            } else {
                fprintf(stderr, "Error: Unknown property '%s' for 'sphere'. (Line %d)\n", key, line);
                exit(1);
            }
        }
    }
    
    // check for missing fields
    if (!hasradius) {
        fprintf(stderr, "Error: Sphere missing 'radius' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hascolor) {
        fprintf(stderr, "Error: Sphere missing 'color' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hasposition) {
        fprintf(stderr, "Error: Sphere missing 'position' field. (Line %d)\n", line);
        exit(1);
    }
}

// gets plane information and stores it into an object
void parse_plane(FILE* json, Object* object) {
    int c;
    
    // used to check that all fields for a plane are present
    int hasnormal = 0;
    int hascolor = 0;
    int hasposition = 0;
    
    // set object kind to plane
    object->kind = PLANE;
    object->reflectivity = 0;
    object->refractivity = 0;
    object->ior = 1;
    
    skip_ws(json);
    
    // check fields for this plane
    while(1) {
        c = next_c(json);
        if (c == '}') {
            break;
        } else if (c == ',') {
            skip_ws(json);
            char* key = next_string(json);

            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            
            // set values for this plane depending on what key was read and sets its 'boolean' to reflect the found field
            if (strcmp(key, "normal") == 0) {
                double* value = next_vector(json);
                object->plane.normal[0] = value[0];
                object->plane.normal[1] = value[1];
                object->plane.normal[2] = value[2];
                hasnormal = 1;
            } else if (strcmp(key, "diffuse_color") == 0) {
                double* value = next_vector(json);
                object->diffuse_color[0] = value[0];
                object->diffuse_color[1] = value[1];
                object->diffuse_color[2] = value[2];
                hascolor = 1;
            } else if (strcmp(key, "specular_color") == 0) {
                double* value = next_vector(json);
                object->specular_color[0] = value[0];
                object->specular_color[1] = value[1];
                object->specular_color[2] = value[2];
            } else if (strcmp(key, "position") == 0) {
                double* value = next_vector(json);
                object->position[0] = value[0];
                object->position[1] = value[1];
                object->position[2] = value[2];
                hasposition = 1;}
            else if (strcmp(key, "reflectivity") == 0) {
                object->reflectivity = clamp(next_number(json), 0, 1);
            } else if (strcmp(key, "refractivity") == 0) {
                object->refractivity = clamp(next_number(json), 0, 1);
            } else if (strcmp(key, "ior") == 0) {
                object->ior = next_number(json);
            } else {
                fprintf(stderr, "Error: Unknown property '%s' for 'sphere'. (Line %d)\n", key, line);
                exit(1);
            }
        }
    }
    
    // check for missing fields
    if (!hasnormal) {
        fprintf(stderr, "Error: Plane missing 'normal' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hascolor) {
        fprintf(stderr, "Error: Plane missing 'color' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hasposition) {
        fprintf(stderr, "Error: Plane missing 'position' field. (Line %d)\n", line);
        exit(1);
    }
}

// gets light information and stores it into a light object
void parse_light(FILE* json, Light* light) {
    int c;
    
    // used to check that all fields for a light are present
    int hascolor = 0;
    int hasposition = 0;
    int hasr2 = 0;
    int hasr1 = 0;
    int hasr0 = 0;
    
    // default theta to 0, gets changed later if one is found
    light->theta = 0;
    
    skip_ws(json);
    
    while(1) {
        c = next_c(json);
        if (c == '}') {
            break;
        } else if (c == ',') {
            skip_ws(json);
            char* key = next_string(json);

            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            
            // set values for this light depending on what key was read and sets its 'boolean' to reflect the found field
            if (strcmp(key, "color") == 0) {
                double* value = next_vector(json);
                light->color[0] = value[0];
                light->color[1] = value[1];
                light->color[2] = value[2];
                hascolor = 1;
            } else if (strcmp(key, "position") == 0) {
                double* value = next_vector(json);
                light->position[0] = value[0];
                light->position[1] = value[1];
                light->position[2] = value[2];
                hasposition = 1;
            } else if (strcmp(key, "radial-a2") == 0) {
                light->radial_a2 = next_number(json);
                hasr2 = 1;
            } else if (strcmp(key, "radial-a1") == 0) {
                light->radial_a1 = next_number(json);
                hasr1 = 1;
            } else if (strcmp(key, "radial-a0") == 0) {
                light->radial_a0 = next_number(json);
                hasr0 = 1;
            } else if (strcmp(key, "theta") == 0) {
                light->theta = next_number(json);
            } else if (strcmp(key, "direction") == 0) {
                double* value = next_vector(json);
                light->direction[0] = value[0];
                light->direction[1] = value[1];
                light->direction[2] = value[2];
            } else if (strcmp(key, "angular-a0") == 0) {
                light->angular_a0 = next_number(json);
            } else {
                fprintf(stderr, "Error: Unknown property '%s' for 'light'. (Line %d)\n", key, line);
                exit(1);
            }
        }
    }
    
    // check for missing fields
    if (!hascolor) {
        fprintf(stderr, "Error: Light missing 'color' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hasposition) {
        fprintf(stderr, "Error: Light missing 'position' field. (Line %d)\n", line);
        exit(1);
    }
    
    if (!hasr2 || !hasr1 || !hasr0) {
        fprintf(stderr, "Error: Light missing one or more radial field(s). (Line %d)\n", line);
        exit(1);
    }
}

// calculate the sphere intersect
double sphere_intersect(double* Ro, double* Rd, double* C, double r) {
    double a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    double b = 2*(Rd[0]*(Ro[0]-C[0]) + Rd[1]*(Ro[1]-C[1]) + Rd[2]*(Ro[2]-C[2]));
    double c = sqr(Ro[0]-C[0]) + sqr(Ro[1]-C[1]) + sqr(Ro[2]-C[2]) - sqr(r);
    
    // check determinant
    double det = sqr(b) - 4*a*c;
    if (det < 0)
        return -1;
    
    det = sqrt(det);
    
    // return t value if an intersect was found, otherwise return -1
    double t0 = (-b - det) / (2*a);
    if (t0 > 0)
        return t0;
    
    double t1 = (-b + det) / (2*a);
    if (t1 > 0)
        return t1;
    
    return -1;
}

// calculate the plane intersect
double plane_intersect(double* Ro, double* Rd, double* P, double* N) {
    // dot product of normal and position to find distance 
    double d = -(N[0]*P[0] + N[1]*P[1] + N[2]*P[2]); 
    double t = -(N[0]*Ro[0] + N[1]*Ro[1] + N[2]*Ro[2] + d) / (N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2]);
    
    // return t value if an intersect was found, otherwise return -1
    if (t > 0)
        return t;
    
    return -1;
}

// skips white space in file
void skip_ws(FILE* json) {
    int c = next_c(json);
    
    while(isspace(c))
        c = next_c(json);
    
    ungetc(c, json);
}

// check that a certain character is next in file
void expect_c(FILE* json, int d) {
    int c = next_c(json);
    
    if (c == d)
        return;
    
    // error if the character found was not what was expected
    fprintf(stderr, "Error: Expected '%c'. (Line %d)\n", d, line);
    exit(1);
}

// get the next character in file
int next_c(FILE* json) {
    int c = fgetc(json);
    
    // increment line if newline was found in file
    if (c == '\n')
        line++;
    
    // error if end of file found when another character was expected
    if (c == EOF) {
        fprintf(stderr, "Error: Unexpected end of file. (Line %d)\n", line);
        exit(1);
    }
    
    return c;
}

// get next string in file
char* next_string(FILE* json) {
    char buffer[129];
    int c = next_c(json);
    
    // look for start of a string
    if (c != '"') {
        fprintf(stderr, "Error: Expected string. (Line %d)\n", line);
        exit(1);
    }
    
    c = next_c(json);
    
    int i = 0;
    // get characters until the end of the string is found or string becomes bigger than 128 characters
    while (c != '"') {
        // checks length
        if (i >= 128) {
            fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported. (Line %d)\n", line);
            exit(1);
        }
        
        // checks for escape codes
        if (c == '\\') {
            fprintf(stderr, "Error: Strings with escape codes are not supported. (Line %d)\n", line);
            exit(1);
        }
        
        // checks that all characters in file are ASCII
        if (c < 32 || c > 126) {
            fprintf(stderr, "Error: Strings may only contain ASCII characters. (Line %d)\n", line);
            exit(1);
        }
        
        // saves character to buffer and increment index
        buffer[i] = c;
        i++;
        c = next_c(json);
    }
    
    // null terminate string
    buffer[i] = 0;
    return strdup(buffer);
}

// gets next number in file
double next_number(FILE* json) {
    double value;
    // look for number, error if one is not found
    if (fscanf(json, "%lf", &value) != 1) {
        fprintf(stderr, "Error: Number value not found. (Line %d)\n", line);
        exit(1);
    }
    
    return value;
}

// get next 3 number vector in file
double* next_vector(FILE* json) {
    double* v = malloc(3*sizeof(double));
    
    // check for start of vector and first number
    expect_c(json, '[');
    skip_ws(json);
    v[0] = next_number(json);
    skip_ws(json);
    
    // check for second number
    expect_c(json, ',');
    skip_ws(json);
    v[1] = next_number(json);
    skip_ws(json);
    
    // check for third number
    expect_c(json, ',');
    skip_ws(json);
    v[2] = next_number(json);
    skip_ws(json);
    
    // check for end of vector
    expect_c(json, ']');
    return v;
}

// if number is greater than max, lower it to max, if less than min, raise it to min
double clamp(double number, double min, double max) {
    // clamps number
    if (number < min)
        return min;

    if (number > max)
        return max;

    return number;
}

// outputs data in buffer to output file
void output_p6(FILE* outputfp) {
    // create header
    fprintf(outputfp, "P6\n%d %d\n%d\n", h, w, MAX_COLOR);
    // writes buffer to output Pixel by Pixel
    fwrite(pixmap, sizeof(Pixel), h*w, outputfp);
}