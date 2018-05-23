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
    scene.fov = 1 # simulate orthographic projection

def draw_wireframe(vertices, edges, color=vec(1,1,1)):
    initialize()
    vertex_size = 0.2
    edge_size = vertex_size / 3
    center_plot = [8, 0, 0]
    # offset the center, adjust orientation
    vectors = [vec(v[1] - center_plot[1], v[2] - center_plot[2], - v[0] + center_plot[0]) for v in vertices]
    for v in vectors:
        sphere(pos = v, radius = vertex_size, color = color)
    
    for edge in edges:
        cylinder(pos = vectors[edge[0]], axis = vectors[edge[1]] - vectors[edge[0]], radius = edge_size, color = color)

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

def extend_by_rotation(vertices):
  rotation_angle = 2 * pi / 7
  for i in range(6):
    rotated_vertices = [space_rotation(rotation_angle, v) for v in vertices]
    join(vertices, rotated_vertices)


def get_edges(vertices, target_inner_product = None):
  if target_inner_product == None:
    target_inner_product = min([inner(vertices[0], vertices[i]) for i in range(1, len(vertices))])
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

  vertices = get_heptagon_vertices(v06, v00, v01) # central heptagon
  vertices = dedup(vertices)
  join(vertices, [v10])
  join(vertices, get_heptagon_vertices(v10, v00, v01)) # 7 neighbors of center

  extend_by_rotation(vertices)
  join(vertices, get_heptagon_vertices(vertices[9], vertices[8], vertices[15])) # 7 neighbors of previous
  extend_by_rotation(vertices)

  join(vertices, get_heptagon_vertices(vertices[33], vertices[32], vertices[62])) # 7 neighbors of previous
  join(vertices, get_heptagon_vertices(vertices[33], vertices[34], vertices[55])) # mirror image of above
  extend_by_rotation(vertices)

  return vertices

def extend_edges(ratio, vertices, edges):
  inner_prod_extended_edge = None
  new_vertices = []
  for edge in edges:
    v1 = vertices[edge[0]]
    v2 = vertices[edge[1]]
    v1_extended = [ - (v2[index] - v1[index]) * ratio + v1[index] for index in range(len(v1))]
    new_vertices.append(v1_extended)
    v2_extended = [ - (v1[index] - v2[index]) * ratio + v2[index] for index in range(len(v1))]
    new_vertices.append(v2_extended)
    if inner_prod_extended_edge == None:
      inner_prod_extended_edge = inner(v1_extended, v2_extended)
  return dedup(new_vertices), inner_prod_extended_edge

vertices = get_vertices_face_first()
vertices = dedup(vertices)
edges = get_edges(vertices)

ratio = 2/(1/sin(pi/14)-2)

stellated_vertices, inner_prod_extended_edge = extend_edges(ratio, vertices, edges)
stellated_edges = get_edges(stellated_vertices, inner_prod_extended_edge)

print('Vertex count: ' + str(len(vertices)))
print('Edge count: ' + str(len(edges)))

print('Stellated vertex count: ' + str(len(stellated_vertices)))
print('Stellated edge count: ' + str(len(stellated_edges)))







draw_wireframe(stellated_vertices, stellated_edges, vec(1,1,0))
draw_wireframe(vertices, edges, vec(1, 1, 2))
