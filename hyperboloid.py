from math import *

phi = (sqrt(5) + 1) / 2
hyperbolic_signature = [+1, -1, -1, -1]

def match(v1, v2):
  return abs(v1 - v2) < 10 ** (-10)

def inner(v1, v2, sig=hyperbolic_signature):
  result = 0.0
  for index in xrange(len(v1)):
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
