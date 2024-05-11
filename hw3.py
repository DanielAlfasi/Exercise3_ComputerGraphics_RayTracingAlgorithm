from helper_classes import *
import matplotlib.pyplot as plt

def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom

    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            # screen is on origin
            pixel = np.array([x, y, 0])
            origin = camera
            direction = normalize(pixel - origin)
            ray = Ray(origin, direction)

            result = find_intersection_with_ray(ray, objects)   
            
            if result is not None:
                nearest_intersection_obj, intersection_point = result
                hit = elevate_point(nearest_intersection_obj, intersection_point)

                color = get_color(ambient, lights, nearest_intersection_obj, hit, objects, max_depth, ray, 0, camera)

            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color,0,1)

    return image

def get_color(ambient, lights, nearest_object, point_of_intersection, objects, max_depth, ray, depth, camera):
    color = np.zeros(3)
    color += calc_ambient_color(ambient, nearest_object)

    for light in lights:
        if is_object_obscured(light, point_of_intersection, objects):
            continue
        color += calc_diffuse_color(nearest_object, light, point_of_intersection)
        color += calc_specular_color(nearest_object, light, point_of_intersection, camera)

    
    depth += 1
    if max_depth <= depth:
        return color
    
    reflective_ray = construct_reflective_ray(ray, nearest_object, point_of_intersection)
    result = find_intersection_with_ray(reflective_ray, objects)
    if result is not None:
        nearest_intersection_obj, intersection_point = result
        hit = elevate_point(nearest_intersection_obj, intersection_point)
        color += nearest_object.reflection * get_color(ambient, lights, nearest_intersection_obj, hit, objects, max_depth, reflective_ray, depth, camera)
    return color
    
def construct_reflective_ray(ray, object, hit):
    L = normalize(ray.direction)
    N = normalize(compute_normal_in_hit_point(object, hit))

    return Ray(hit, L - 2 * (np.dot(L, N)) * N)
    
def compute_normal_in_hit_point(object, point):
    if isinstance(object, Sphere):
        return normalize(point - object.center)
    
    return object.normal

def find_intersection_with_ray(ray, objects):
    nearest_object, distance = ray.nearest_intersected_object(objects)
    if nearest_object is None:
        return None
    
    intersection_point = ray.origin + ray.direction * distance
    return nearest_object, intersection_point

def calc_ambient_color(ambient, object):
    return ambient * object.ambient

def calc_diffuse_color(object, light_object, hit):
    object_normal = normalize(compute_normal_in_hit_point(object, hit))
    light_direction_normalized = normalize(light_object.get_light_ray(hit).direction)

    dot_prod = np.dot(object_normal, light_direction_normalized)
    return object.diffuse * light_object.get_intensity(hit) * dot_prod

def calc_specular_color(object, light_object, hit, view_point):
    object_normal = normalize(compute_normal_in_hit_point(object, hit))
    v = normalize(view_point - hit)

    light_direction_normalized = normalize((-1) * light_object.get_light_ray(hit).direction)
    light_direction_reflection_normalized = reflected(light_direction_normalized, object_normal)
    dot_prod = np.dot(v, light_direction_reflection_normalized)**object.shininess

    return object.specular * light_object.get_intensity(hit) * dot_prod

def is_object_obscured(light, hit, objects):
    ray = light.get_light_ray(hit)

    nearest_object , distance = ray.nearest_intersected_object(objects)

    if nearest_object is None:
        return False
    
    if distance < light.get_distance_from_light(hit):
        return True
 
    return False

def elevate_point(object, point):
    return (compute_normal_in_hit_point(object, point) * 0.01) + point

def your_own_scene():
    # Define the sphere
    sphere = Sphere([-0.7, 0, -1], 1)
    sphere.set_material([0.5, 0.3, 0.8], [1, 1, 1], [0.5, 0.5, 0.5], 50, 0.9)

    # Define the diamond using a custom set of vertices for a pyramid-like shape
    v_list = np.array([
        [0.5, 0, -2],  # Base vertex A
        [1.5, 0, -2],  # Base vertex B
        [1.5, 1, -2],  # Base vertex C
        [0.5, 1, -2],  # Base vertex D
        [1, 0.5, -1]  # Top vertex E (apex)
    ])
    diamond = Pyramid(v_list)
    diamond.set_material([1, 0, 0], [1, 1, 1], [0.3, 0.3, 0.3], 100, 0.9)

    # Reflective floor
    plane = Plane([0, 1, 0], [0, -1, 0])
    plane.set_material([0.2, 0.2, 0.2], [0.2, 0.2, 0.2], [1, 1, 1], 1000, 0.5)

    # Background for reflections
    background = Plane([0, 0, 1], [0, 0, -3])
    background.set_material([0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [0.1, 0.1, 0.1], 10, 0.1)

    # Lighting
    point_light = PointLight(intensity=np.array([1, 1, 1]), position=np.array([0, 2, 0]), kc=0.1, kl=0.1, kq=0.1)
    spotlight = SpotLight(intensity=np.array([1, 0, 0]), position=np.array([0, -0.5, 0]),
                          direction=np.array([0, 0, -1]), kc=0.1, kl=0.1, kq=0.1)

    # Camera setup
    camera = np.array([0, 0.5, 1])

    objects = [sphere, diamond, plane, background]
    lights = [point_light, spotlight]

    return camera, lights, objects
