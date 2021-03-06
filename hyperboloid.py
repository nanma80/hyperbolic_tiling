from math import *
import csv
import matplotlib.pyplot as plt
import matplotlib.collections as pltcol


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


def plot_klein_model(vertices, edges):
  vertices_klein = [[v[1]/v[0], v[2]/v[0]] for v in vertices]
  vertices_klein_0 = [v[0] for v in vertices_klein]
  vertices_klein_1 = [v[1] for v in vertices_klein]
  plt.plot(vertices_klein_0, vertices_klein_1, 'o')
  lines = [[tuple(vertices_klein[j]) for j in i] for i in edges]
  lc = pltcol.LineCollection(lines)
  plt.axes().add_collection(lc)
  plt.axes().set_aspect('equal')
  plt.show()

