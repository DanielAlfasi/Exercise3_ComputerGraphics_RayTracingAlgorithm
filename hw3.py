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

            color = np.zeros(3)
            color = get_color(ambient, lights, objects, max_depth, ray, 0, color)

            
            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color,0,1)

    return image

def get_color(ambient, lights, objects, max_depth, ray, depth, color):

    nearest_intersection_obj, distance = ray.nearest_intersected_object(objects)                
    t, _ = nearest_intersection_obj.intersect(ray)
    if t is None:
        print("Hello")
    point_of_intersection = ray.get_point_on_ray(t)

    color += calc_ambient_color(ambient, nearest_intersection_obj)

    for light in lights:
        if is_object_obscured(light, nearest_intersection_obj, point_of_intersection, objects):
            continue
        color += calc_diffuse_color(nearest_intersection_obj, light, point_of_intersection)
        color += calc_specular_color(nearest_intersection_obj, light, point_of_intersection, ray.origin)

    return color
    depth += 1
    if max_depth < depth:
        return color
    

        
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

def is_object_obscured(light, object, hit, objects):
    ray = light.get_light_ray(hit)

    nearest_object , _ = ray.nearest_intersected_object(objects)

    if nearest_object == object:
        return False
    
    return True


# Write your own objects and lights
# TODO
def your_own_scene():
    camera = np.array([0,0,1])
    lights = []
    objects = []
    return camera, lights, objects
