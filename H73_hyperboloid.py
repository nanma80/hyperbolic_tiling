from math import *
from hyperboloid import *

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

vertices73 = get_vertices_face_first()
vertices73 = dedup(vertices73)
edges73 = get_edges(vertices73)

extend_ratio = 2/(1/sin(pi/14)-2)

vertices727, inner_prod_extended_edge727 = extend_edges(extend_ratio, vertices73, edges73)
edges727 = get_edges(vertices727, inner_prod_extended_edge727)

vertices37 = vertices727
edges37 = get_edges(vertices37)

print('{7, 3} vertex count: ' + str(len(vertices73)))
print('{7, 3} edge count: ' + str(len(edges73)))

print('{7/2, 7} vertex count: ' + str(len(vertices727)))
print('{7/2, 7} edge count: ' + str(len(edges727)))

print('{3, 7} vertex count: ' + str(len(vertices37)))
print('{3, 7} edge count: ' + str(len(edges37)))

csv_write('data_73', vertices73, edges73)
csv_write('data_727', vertices727, edges727)
csv_write('data_37', vertices37, edges37)