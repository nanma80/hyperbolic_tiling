GlowScript 2.7 VPython

# hyperboloid.py
phi = (sqrt(5) + 1) / 2
hyperbolic_signature = [+1, -1, -1, -1]

def match(v1, v2):
  return abs(v1 - v2) < 10 ** (-10)

def inner(v1, v2, sig=hyperbolic_signature):
  result = 0.0
  for index in range(len(v1)):
    result += sig[index] * v1[index] * v2[index]
  return result

def norm(v):
  return sqrt(inner(v, v))

def scale(v, factor):
  return [el * factor for el in v]

def join(original, vertices):
  for v in vertices:
    matched = False
    for vo in original:
      if match(inner(v, vo), inner(v, v)):
        matched = True
        break
    if not matched:
      original.append(v)
  return original

def dedup(vertices):
  output = []
  return join(output, vertices)

def space_rotation(angle, vertex):
  matrix = [
    [cos(angle), -sin(angle)],
    [sin(angle), cos(angle)]
  ]
  new_vertex = [vertex[0]]
  new_vertex.append(vertex[1] * cos(angle) - vertex[2] * sin(angle))
  new_vertex.append(vertex[1] * sin(angle) + vertex[2] * cos(angle))
  for index in range(3, len(vertex)):
    new_vertex.append(vertex[index])
  return new_vertex


# glowscript specific
def initialize():
    scene.fov = 0.001 # simulate orthographic projection

def draw_wireframe(vertices, edges):
    initialize()
    vertex_size = 0.2
    edge_size = vertex_size / 3
    center_plot = [4, 0, 0]
    vectors = [vec(v[0] - center_plot[0], v[1] - center_plot[1], v[2] - center_plot[2]) for v in vertices] # offset the center
    for v in vectors:
        sphere(pos = v, radius = vertex_size)
    
    for edge in edges:
        cylinder(pos = vectors[edge[0]], axis = vectors[edge[1]] - vectors[edge[0]], radius = edge_size)

# H73_hyperboloid.py

def get_heptagon_next_vertex(v1, v2, v3):
  ratio = 1 + 2 * cos( 2 * pi / 7)
  return [(v3[index] - v2[index]) * ratio + v1[index] for index in range(len(v1))]

def get_heptagon_vertices(v1, v2, v3):
  # given three vertices v1, v2, v3, generate all 7 vertices such that all vertices are on the same affine regular heptagon
  vertices = [v1, v2, v3]
  for index in range(4):
    vertices.append(get_heptagon_next_vertex(vertices[-3], vertices[-2], vertices[-1]))
  return vertices

def get_edges(vertices, target_inner_product):
  edges = []
  for i in range(len(vertices)):
    for j in range(i+1, len(vertices)):
      inner_prod = inner(vertices[i], vertices[j])
      if abs(inner_prod - target_inner_product) < 0.01:
        edges.append([i, j])
  return edges

def vertex_first():
  ch2phi = ((8./3 * cos(pi / 7)**2) - 1)
  rsquare = 1/(ch2phi ** 2 - 1)
  x0square = rsquare + 1
  r = sqrt(rsquare)
  x0 = sqrt(x0square)
  costheta = ((x0 - r) ** 2 - 0.5) / ((x0 - r) ** 2 + 1.0)
  theta = acos(costheta)
  print('Angle at the vertex of the vertex first {7,3} is: ' + theta * 180 / pi)

def get_vertices_face_first():
  rotation_angle = 2 * pi / 7
  cosine = cos(2 * pi / 7)
  ch2phi = ((8./3 * cos(pi / 7)**2) - 1)
  rsquare = (1 - cosine)/(ch2phi - 1)
  x0square = rsquare + 1
  r = sqrt(rsquare)
  x0 = sqrt(x0square)
  v00 = [x0, 1, 0]
  v01 = [x0, cos(2*pi/7), sin(2*pi/7)]
  v06 = [x0, cos(2*pi/7), - sin(2*pi/7)]

  b = 1 + 2 * cosine
  a = sqrt(rsquare + b ** 2)
  v10 = [a, b, 0]

  vertices = get_heptagon_vertices(v06, v00, v01)
  vertices = dedup(vertices)
  join(vertices, [v10])
  join(vertices, get_heptagon_vertices(v10, v00, v01))

  for i in range(6):
    rotated_vertices = [space_rotation(rotation_angle, v) for v in vertices]
    join(vertices, rotated_vertices)
  return vertices

cosine = cos(2 * pi / 7)
ch2phi = ((8./3 * cos(pi / 7)**2) - 1)
rsquare = (1 - cosine)/(ch2phi - 1)
target_inner_product = ch2phi * rsquare

vertices = get_vertices_face_first()
vertices = dedup(vertices)
edges = get_edges(vertices, target_inner_product)
print('Vertex count: ' + str(len(vertices)))
print('Edge count: ' + str(len(edges)))




draw_wireframe(vertices, edges)

