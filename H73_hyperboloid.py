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
