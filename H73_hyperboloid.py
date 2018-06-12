from math import *
from hyperboloid import *


def get_edges_by_distance(vertices, target_distance = None):
  if target_distance == None:
    target_distance = min([distance_square(vertices[0], vertices[i]) for i in range(1, len(vertices))])
  edges = []
  for i in range(len(vertices)):
    for j in range(i+1, len(vertices)):
      each_distance_square = distance_square(vertices[i], vertices[j])
      if abs(each_distance_square - target_distance) < 0.01:
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

def get_73_vertices_face_first():
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
  join(vertices, get_heptagon_vertices(v10, v00, v01)) # 7 neighbors of center
  extend_by_rotation(vertices, 7)

  join(vertices, get_heptagon_vertices(vertices[9], vertices[8], vertices[15])) # 7 neighbors of previous
  extend_by_rotation(vertices, 7)

  join(vertices, get_heptagon_vertices(vertices[33], vertices[32], vertices[62])) # 7 neighbors of previous
  join(vertices, get_heptagon_vertices(vertices[33], vertices[34], vertices[55])) # mirror image of above
  extend_by_rotation(vertices, 7)

  return vertices

def get_73_one_face():
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

  vertices = get_heptagon_vertices(v10, v00, v01) # central heptagon
  vertices = dedup(vertices)

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




vertices73 = get_73_vertices_face_first()
vertices73 = dedup(vertices73)
edges73 = get_edges(vertices73)

extend_ratio = 2/(1/sin(pi/14)-2)

vertices727, inner_prod_extended_edge727 = extend_edges(extend_ratio, vertices73, edges73)
edges727 = get_edges(vertices727, inner_prod_extended_edge727)

vertices37 = vertices727
edges37 = get_edges(vertices37)

rectified_vertices73 = rectify(vertices73, edges73)
rectified_edges73 = get_edges(rectified_vertices73)

# rectified_vertices37 = rectify(vertices37, edges37)
# rectified_edges37 = get_edges(rectified_vertices37)


# dual_vertices73, inner_prod_dual = dual_edges_to_points(vertices73, edges73)
# dual_edges73 = get_edges(dual_vertices73, inner_prod_dual)
# dual_vertices73, dual_edges73 = dualize(vertices73, edges73)
dual_rectified_vertices73, dual_rectified_edges73 = dualize(rectified_vertices73, rectified_edges73)



# print('{7, 3} vertex count: ' + str(len(vertices73)))
# print('{7, 3} edge count: ' + str(len(edges73)))

# print('{7/2, 7} vertex count: ' + str(len(vertices727)))
# print('{7/2, 7} edge count: ' + str(len(edges727)))

# print('{3, 7} vertex count: ' + str(len(vertices37)))
# print('{3, 7} edge count: ' + str(len(edges37)))


# dual_vertices37, dual_edges37 = dualize(vertices37, edges37)

# print('Dual {3, 7} vertex count: ' + str(len(dual_vertices37)))
# print('Dual {3, 7} edge count: ' + str(len(dual_edges37)))


# print('Dual {7, 3} vertex count: ' + str(len(dual_vertices73)))
# print('Dual {7, 3} edge count: ' + str(len(dual_edges73)))

# print('r{7, 3} vertex count: ' + str(len(rectified_vertices73)))
# print('r{7, 3} edge count: ' + str(len(rectified_edges73)))

# print('r{3, 7} vertex count: ' + str(len(rectified_vertices37)))
# print('r{3, 7} edge count: ' + str(len(rectified_edges37)))

# print('Dual r{7, 3} vertex count: ' + str(len(dual_rectified_vertices73)))
# print('Dual r{7, 3} edge count: ' + str(len(dual_rectified_edges73)))


other_end_of_zeroth_vertex = [edge[1] for edge in dual_rectified_edges73 if edge[0] == 0]
highlighted_edges = [edge for edge in dual_rectified_edges73 if edge[0] == 0]
highlighted_vertices = [dual_rectified_vertices73[i] for i in other_end_of_zeroth_vertex]

# csv_write('data_73', vertices73, edges73)
# csv_write('data_727', vertices727, edges727)
# csv_write('data_37', vertices37, edges37)

# print(dual_rectified_vertices73[0])
# print(dual_rectified_vertices73[2])
# print other_end_of_zeroth_vertex

for third_index in [5, 7]:
  highlighted_faces_vertices = get_heptagon_vertices(dual_rectified_vertices73[2], dual_rectified_vertices73[0], dual_rectified_vertices73[third_index])
  get_edges(highlighted_faces_vertices)

