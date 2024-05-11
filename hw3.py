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

                color = get_color(ambient, lights, nearest_intersection_obj, hit, objects, max_depth, ray, 0)

            
            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color,0,1)

    return image

def get_color(ambient, lights, nearest_object, point_of_intersection, objects, max_depth, ray, depth):
    color = np.zeros(3)
    color += calc_ambient_color(ambient, nearest_object)

    for light in lights:
        if is_object_obscured(light, point_of_intersection, objects):
            continue
        color += calc_diffuse_color(nearest_object, light, point_of_intersection)
        color += calc_specular_color(nearest_object, light, point_of_intersection, ray.origin)

    return color
    depth += 1
    if max_depth < depth:
        return color
    
def find_intersection_with_ray(ray, objects):
    nearest_object, distance = ray.nearest_intersected_object(objects)
    if nearest_object is None:
        return None
    
    intersection_point = ray.origin + ray.direction * distance
    return nearest_object, intersection_point

def calc_ambient_color(ambient, object):
    return ambient * object.ambient

def calc_diffuse_color(object, light_object, hit):
    object_normal = normalize(object.normal)
    light_direction_normalized = normalize(light_object.get_light_ray(hit).direction)

    dot_prod = np.dot(object_normal, light_direction_normalized)
    return object.diffuse * light_object.get_intensity(hit) * dot_prod

def calc_specular_color(object, light_object, hit, view_point):
    object_normal = normalize(object.normal)
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
    if isinstance(object, Sphere):
        return None
    
    return (object.normal * 0.01) + point
# Write your own objects and lights
# TODO
def your_own_scene():
    camera = np.array([0,0,1])
    lights = []
    objects = []
    return camera, lights, objects
