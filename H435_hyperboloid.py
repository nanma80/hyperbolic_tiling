import math
from hyperboloid import *

phi = (math.sqrt(5) + 1) / 2
radius = phi ** (-1./2)
scaling_factor = phi ** 2

def get_x0(v3d, radius):
  return math.sqrt(inner(v3d, v3d, [1, 1, 1]) + radius ** 2)

def normalize(v3d):
  v3d_scaled_down = scale(v3d, 1.0/ scaling_factor)
  x0 = get_x0(v3d_scaled_down, radius)
  return [x0] + v3d_scaled_down

center = [phi ** (-1./2),0,0,0]
dist1 = [phi**(3./2), phi,1,0]
dist2 = [phi**(5./2), 2*phi,0,0]

print "target inner prod:", inner(dist1, center)
print "target norm:", norm(center)

vzome_dist1 = normalize([1+phi, 0, -1-2*phi])
vzome_opposite = normalize([4+6*phi, 0, -6-10*phi])

vzome_face_before_scaling = [
  [1+phi, 0, -1-2*phi],
  [2+3*phi, -1-phi, -1-2*phi],
  [3+5*phi, 0, -1-2*phi],
  [2+3*phi, 1+phi, -1-2*phi]
]

vzome_face = [normalize(v) for v in vzome_face_before_scaling]
inner_point = normalize([2+4*phi, 0, 0])
outer_point = normalize([6+10*phi,0,-4-6*phi])

for x in vzome_face:
  print match(inner(inner_point, x), inner(outer_point, x))
