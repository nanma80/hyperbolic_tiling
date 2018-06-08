from math import *
import csv

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
  inner_prod = inner(dual_vertices_edges_around_zero[0], dual_vertices_edges_around_zero[1])
  # print(inner(dual_vertices_edges_around_zero[0], dual_vertices_edges_around_zero[2]))
  # print(inner(dual_vertices_edges_around_zero[1], dual_vertices_edges_around_zero[2]))
  # print(inner(dual_vertices_edges_around_zero[0], scale(dual_vertices_edges_around_zero[1], -1)))
  return dedup(dual_vertices), abs(inner_prod)


def rectify(vertices, edges):
  output_vertices = []
  for edge in edges:
    v1 = vertices[edge[0]]
    v2 = vertices[edge[1]]
    v_mid = [(v1[index] + v2[index])/2 for index in range(len(v1))]
    output_vertices.append(v_mid)
  return output_vertices
