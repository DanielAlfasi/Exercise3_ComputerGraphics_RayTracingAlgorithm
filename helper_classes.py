import numpy as np


# This function gets a vector and returns its normalized form.
def normalize(vector):
    return vector / np.linalg.norm(vector)


# This function gets a vector and the normal of the surface it hit
# This function returns the vector that reflects from the surface
def reflected(vector, axis):
    axis_normalized = normalize(axis)
    v = vector - 2 * np.dot(np.dot(vector, axis_normalized), axis_normalized)
    return v

## Lights


class LightSource:
    def __init__(self, intensity):
        self.intensity = intensity


class DirectionalLight(LightSource):

    def __init__(self, intensity, direction):
        super().__init__(intensity)
        self.direction = direction

    # This function returns the ray that goes from a point to the light source
    def get_light_ray(self,intersection_point):
        return Ray(intersection_point, self.direction)

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.inf

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        return self.intensity


class PointLight(LightSource):
    def __init__(self, intensity, position, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from a point to the light source
    def get_light_ray(self,intersection):
        return Ray(intersection, normalize(self.position - intersection))

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self,intersection):
        return np.linalg.norm(intersection - self.position)

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        d = self.get_distance_from_light(intersection)
        return self.intensity / (self.kc + self.kl*d + self.kq * (d**2))


class SpotLight(LightSource):
    def __init__(self, intensity, position, direction, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.kc = kc
        self.kl = kl
        self.kq = kq
        self.direction = direction

    # This function returns the ray that goes from a point to the light source
    def get_light_ray(self, intersection):
        return Ray(intersection ,self.position - intersection)

    def get_distance_from_light(self, intersection):
        return np.linalg.norm(self.position - intersection)

    def get_intensity(self, intersection):
        distance = self.get_distance_from_light(intersection)
        nominator = self.intensity
        denominator = self.kc + (self.kl * distance) + (self.kq * (distance**2))
        return nominator / denominator


class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

    # The function is getting the collection of objects in the scene and looks for the one with minimum distance.
    # The function should return the nearest object and its distance (in two different arguments)
    def nearest_intersected_object(self, objects):
        nearest_object = None
        min_distance = np.inf

        for object in objects:
            result = object.intersect(self)
            if not result:
                continue

            t, _ = result
            point_of_intersection = self.get_point_on_ray(t)
            distance = np.linalg.norm(self.origin, point_of_intersection)

            if distance < min_distance:
                min_distance = distance
                nearest_object = object

        return nearest_object, min_distance
    
    def get_point_on_ray(self, t):
        return self.origin + t * self.direction


class Object3D:
    def set_material(self, ambient, diffuse, specular, shininess, reflection):
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflection = reflection


class Plane(Object3D):
    def __init__(self, normal, point):
        self.normal = np.array(normal)
        self.point = np.array(point)

    def intersect(self, ray: Ray):
        v = self.point - ray.origin
        t = np.dot(v, self.normal) / (np.dot(self.normal, ray.direction) + 1e-6)
        if t > 0:
            return t, self
        else:
            return None


class Triangle(Object3D):
    """
        C
        /\
       /  \
    A /____\ B

    The fornt face of the triangle is A -> B -> C.
    
    """
    def __init__(self, a, b, c):
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.normal = self.compute_normal()

    # computes normal to the trainagle surface. Pay attention to its direction!
    def compute_normal(self):
        vector_ab = self.b - self.a
        vector_ac = self.c - self.a
        return np.cross(vector_ab, vector_ac)

    def intersect(self, ray: Ray):
        triangle_plane = Plane(self.normal, self.a)
        point_of_intersection_tuple = triangle_plane.intersect(ray)
        if point_of_intersection_tuple == None:
            return None
        
        t, _ = point_of_intersection_tuple
        point_of_intersection = ray.get_point_on_ray(t)
        triangle_area = np.linalg.norm(self.normal) / 2
        vector_pb = self.b - point_of_intersection
        vector_pc = self.c - point_of_intersection
        vector_pa = self.a - point_of_intersection

        alpha = np.linalg.norm(np.cross(vector_pb, vector_pc)) / (triangle_area * 2)
        beta = np.linalg.norm(np.cross(vector_pc, vector_pa)) / (triangle_area * 2)
        gamma = 1 - alpha - beta

        if (self.is_valid_barycentric_coordinates(alpha, beta, gamma)):
            return t, self
        
        return None

    def is_valid_barycentric_coordinates(self, alpha, beta, gamma, epsilon=1e-10):
        if not (0 - epsilon <= alpha <= 1 + epsilon and
                0 - epsilon <= beta <= 1 + epsilon and
                0 - epsilon <= gamma <= 1 + epsilon):
            return False
    
        if not (1 - epsilon <= alpha + beta + gamma <= 1 + epsilon):
            return False

        return True

        


class Pyramid(Object3D):
    """     
            D
            /\*\
           /==\**\
         /======\***\
       /==========\***\
     /==============\****\
   /==================\*****\
A /&&&&&&&&&&&&&&&&&&&&\ B &&&/ C
   \==================/****/
     \==============/****/
       \==========/****/
         \======/***/
           \==/**/
            \/*/
             E 
    
    Similar to Traingle, every from face of the diamond's faces are:
        A -> B -> D
        B -> C -> D
        A -> C -> B
        E -> B -> A
        E -> C -> B
        C -> E -> A
    """
    def __init__(self, v_list):
        self.v_list = v_list
        self.triangle_list = self.create_triangle_list()

    def create_triangle_list(self):
        l = []
        t_idx = [
                [0,1,3],
                [1,2,3],
                [0,3,2],
                 [4,1,0],
                 [4,2,1],
                 [2,4,0]]
        for triangle_indices in t_idx:
            a, b, c = triangle_indices
            l.append(Triangle(self.v_list[a], self.v_list[b], self.v_list[c]))
            
        return l

    def apply_materials_to_triangles(self):
        for triangle in self.triangle_list:
            triangle.set_material(self.ambient, self.diffuse, self.specular, self.shininess, self.reflection)

    def intersect(self, ray: Ray):
        object, _ = ray.nearest_intersected_object(self.triangle_list)
        return object.intersect(ray)


class Sphere(Object3D):
    def __init__(self, center, radius: float):
        self.center = center
        self.radius = radius

    def intersect(self, ray: Ray):
        #TODO
        pass

