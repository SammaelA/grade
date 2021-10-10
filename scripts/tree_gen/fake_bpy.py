from math import atan2, pi

from chturtle import Vector
import leaf_shapes as leaf_geom
#split_stem = self.branch_curves[stem.depth].splines.new('BEZIER')
class Point(object):
    """Class to store data for each leaf in the system"""

    __slots__ = (
        'co','radius','handle_left','handle_right'
    )
    def __init__(self):
        self.co = Vector((0,0,0))
        self.handle_left = Vector((0,0,0))
        self.handle_right = Vector((0,0,0))
        self.radius = 1.0
class Spline(object):
    """Class to store data for each leaf in the system"""

    __slots__ = (
        'bezier_points','radius_interpolation','resolution_u'
    )
    def __init__(self):
        self.bezier_points = [Point()]
        self.radius_interpolation = 'error'
        self.resolution_u = 1
    def add(self):
        self.bezier_points.append(Point())
class Splines(object):
    """Class to store data for each leaf in the system"""

    __slots__ = (
        'data'
    )
    def __init__(self):
        self.data = []
    def new(self, str):
        self.data.append(Spline())
        return self.data[-1]
class Curve(object):
    """Class to store data for each leaf in the system"""

    __slots__ = (
        'resolution_u', 'splines'
    )
    def __init__(self):
        self.resolution_u = 1
        self.splines = Splines()