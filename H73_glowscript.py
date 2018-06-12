GlowScript 2.7 VPython

# hyperboloid.py






phi = (sqrt(5) + 1) / 2
hyperbolic_signature = [+1, -1, -1, -1]
folder = "data"


def match(v1, v2):
  return abs(v1 - v2) < 10 ** (-10)


def distance_square(v1, v2):
  difference = [v1[index] - v2[index] for index in range(len(v1))]
  return inner(difference, difference)


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
  new_vertex = [vertex[0]]
  new_vertex.append(vertex[1] * cos(angle) - vertex[2] * sin(angle))
  new_vertex.append(vertex[1] * sin(angle) + vertex[2] * cos(angle))
  for index in range(3, len(vertex)):
    new_vertex.append(vertex[index])
  return new_vertex


def csv_write_vertices(filename, vertices):
  with open(filename + '.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for v in vertices:
      writer.writerow(v)


def csv_write_edges(filename, vertices, edges):
  with open(filename + '.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for e in edges:
      writer.writerow(vertices[e[0]] + vertices[e[1]])


def csv_write(filename, vertices, edges):
  csv_write_vertices(folder + '/' + filename + '_vertices', vertices)
  csv_write_edges(folder + '/' + filename + '_edges', vertices, edges)


def cross(a, b):
  c = [a[1]*b[2] - a[2]*b[1],
       a[2]*b[0] - a[0]*b[2],
       a[0]*b[1] - a[1]*b[0]]
  return c


def dual_edge_to_point(edge_coordinates):
  euclidean_cross_prod = cross(edge_coordinates[0], edge_coordinates[1])
  euclidean_cross_prod[0] *= -1
  return euclidean_cross_prod


def dual_edges_to_points(vertices, edges):
  dual_vertices = []
  zero_vertex_index = 0
  dual_vertices_edges_around_zero = []
  for edge in edges:
    edge_coordinates = [vertices[index] for index in edge]
    dual_point = dual_edge_to_point(edge_coordinates)
    dual_vertices.append(dual_point)
    dual_vertices.append(scale(dual_point, -1))
    if edge[0] == zero_vertex_index or edge[1] == zero_vertex_index:
      dual_vertices_edges_around_zero.append(dual_point)
  inner_prod = -abs(inner(dual_vertices_edges_around_zero[0], dual_vertices_edges_around_zero[1]))
  return dedup(dual_vertices), inner_prod


def dualize(vertices, edges):
  dual_vertices = []
  zero_vertex_index = 0
  dual_vertices_edges_around_zero = []
  for edge in edges:
    edge_coordinates = [vertices[index] for index in edge]
    dual_point = dual_edge_to_point(edge_coordinates)
    dual_vertices.append((dual_point, edge))
    dual_vertices.append((scale(dual_point, -1), edge))
    if edge[0] == zero_vertex_index or edge[1] == zero_vertex_index:
      dual_vertices_edges_around_zero.append(dual_point)
  target_inner_prod = - abs(inner(dual_vertices_edges_around_zero[0], dual_vertices_edges_around_zero[1]))
  dual_edges = []
  for i in range(len(dual_vertices)):
    for j in range(i+1, len(dual_vertices)):
      inner_prod = inner(dual_vertices[i][0], dual_vertices[j][0])
      if abs(inner_prod - target_inner_prod) < 0.01:
        edge_i = dual_vertices[i][1]
        edge_j = dual_vertices[j][1]
        overlap_indices = [index for index in edge_i if index in edge_j]
        if len(overlap_indices) > 0:
          dual_edges.append([i, j])
  dual_vertices_coordinates = [component[0] for component in dual_vertices]
  return dedup_polytope(dual_vertices_coordinates, dual_edges)


def dedup_polytope(vertices, edges):
  mapping_vertices = {}
  new_vertices = []
  for v_index, v in enumerate(vertices):
    matched = False
    for vo_index, vo in enumerate(new_vertices):
      if match(inner(v, vo), inner(v, v)):
        matched = True
        break
    if matched:
      mapping_vertices[v_index] = vo_index
    else:
      new_index = len(new_vertices)
      new_vertices.append(v)
      mapping_vertices[v_index] = new_index
  new_edges_set = set()
  for edge in edges:
    new_edge = [mapping_vertices[index] for index in edge]
    sorted_edge = sorted(new_edge)
    new_edges_set.add((sorted_edge[0], sorted_edge[1]))
  new_edges = [list(edge_tuple) for edge_tuple in list(new_edges_set)]
  return new_vertices, new_edges


def rectify(vertices, edges):
  output_vertices = []
  for edge in edges:
    v1 = vertices[edge[0]]
    v2 = vertices[edge[1]]
    v_mid = [(v1[index] + v2[index])/2 for index in range(len(v1))]
    output_vertices.append(v_mid)
  return output_vertices


def extend_by_rotation(vertices, degree):
  rotation_angle = 2 * pi / degree
  for i in range(degree - 1):
    rotated_vertices = [space_rotation(rotation_angle, v) for v in vertices]
    join(vertices, rotated_vertices)


def get_edges(vertices, target_inner_product = None, dual = False):
  if target_inner_product == None:
    target_inner_product = min([inner(vertices[0], vertices[i]) for i in range(1, len(vertices))])
  edges = []
  for i in range(len(vertices)):
    for j in range(i+1, len(vertices)):
      inner_prod = inner(vertices[i], vertices[j])
      if abs(inner_prod - target_inner_product) < 0.01:
        edges.append([i, j])
  return edges


def get_heptagon_next_vertex(v1, v2, v3):
  ratio = 1 + 2 * cos( 2 * pi / 7)
  return [(v3[index] - v2[index]) * ratio + v1[index] for index in range(len(v1))]


def get_heptagon_vertices(v1, v2, v3):
  # given three vertices v1, v2, v3, generate all 7 vertices such that all vertices are on the same affine regular heptagon
  vertices = [v1, v2, v3]
  for index in range(4):
    vertices.append(get_heptagon_next_vertex(vertices[-3], vertices[-2], vertices[-1]))
  return vertices












# glowscript specific
def initialize():
    scene.fov = 1 # simulate orthographic projection

def draw_wireframe(vertices, edges, color=vec(1,1,1), vertex_size = 0.2):
    initialize()
    edge_size = vertex_size / 3
    center_plot = [8, 0, 0]
    center_plot = [0, 0, 0]
    # offset the center, adjust orientation
    vectors = [vec(v[1] - center_plot[1], v[2] - center_plot[2], - v[0] + center_plot[0]) for v in vertices]
    for v in vectors:
        sphere(pos = v, radius = vertex_size, color = color)
    
    for edge in edges:
        cylinder(pos = vectors[edge[0]], axis = vectors[edge[1]] - vectors[edge[0]], radius = edge_size, color = color)

# H73_hyperboloid.py















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
  draw_wireframe(highlighted_faces_vertices, get_edges(highlighted_faces_vertices), vec(1, 1, 0), 0.21)




# csv_write('data_73', vertices73, edges73)
# csv_write('data_727', vertices727, edges727)
# csv_write('data_37', vertices37, edges37)




# draw_wireframe(dual_vertices73, dual_edges73, vec(1, 1, 1), 0.2)
# draw_wireframe(dual_vertices37, dual_edges37, vec(1, 1, 1), 0.2)

# draw_wireframe(rectified_vertices73, rectified_edges73, vec(1, 1, 1), 0.2)
# draw_wireframe(rectified_vertices37, rectified_edges37, vec(1, 1, 1), 0.2)

draw_wireframe(dual_rectified_vertices73, dual_rectified_edges73, vec(1, 1, 1), 0.2)

# draw_wireframe(vertices727, edges727, vec(1, 1, 0), 0.18)
# draw_wireframe(vertices73, edges73, vec(1, 1, 1), 0.2)
# draw_wireframe(vertices37, edges37, vec(1, 1, 1), 0.2)


# draw_wireframe(highlighted_vertices, [], vec(1, 1, 0), 0.21)
# draw_wireframe([dual_rectified_vertices73[0]], [], vec(1, 0, 0), 0.21)
