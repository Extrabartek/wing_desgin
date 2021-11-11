import scipy.integrate as integrate
import math

# global constant
g = 9.81


# class definition list
class Material:
    # if any parameters need to be added follow the method for the already added parameters
    def __init__(self, rho, sig_yld, sig_ult, E, G):
        """
        Initiate variable of type material
        :param rho: Density of the material [kg/m^3]
        :type rho: float
        :param sig_yld: Yield strength of the material [N/m^2]
        :type sig_yld: float
        :param sig_ult: Ultimate strength of the material [N/m^2]
        :type sig_ult: float
        :param E: The elastic modules [N/m^2]
        :type E: float
        :param G: Shear modules [N/m^2]
        :type G: float
        """
        self.rho = rho
        self.sig_yld = sig_yld
        self.sig_ult = sig_ult
        self.E = E
        self.G = G

    def mass(self, volume):
        """
        This function will return the mass (in kilograms) of an object given its volume (in meters cubed)
        :param volume: The object's volume [m^3]
        :type volume: float
        :return: The mass of the object [kg]
        :rtype: float
        """
        return self.rho * volume


class Planform:
    # if any parameters need to be added follow the method for the already added parameters
    def __init__(self, b, cr, ct, sweep_le, spar_rear, spar_front):
        """
        Initiate variable of type planform
        :param b: Wingspan of the planform [m]
        :type b: float
        :param cr: Cord at the root [m]
        :type cr: float
        :param ct: Cord at the tip [m]
        :type ct: float
        :param sweep_le: Leading edge sweep [rad]
        :param spar_rear: The percentage of the cord where the rear spar is
        :type spar_rear: float
        :param spar_front: he percentage of the cord where the front spar is
        :type spar_front: float
        """
        self.b = b
        self.cr = cr
        self.ct = ct
        self.sweep_le = sweep_le
        self.spar_dif = spar_rear - spar_front

    def area(self):
        """
        This function returns the area (in m^2) of the planform
        :return: Area of the planform [m^2]
        :rtype: float
        """
        return self.b * (self.cr + self.ct) * 0.5

    def sweep(self, c_percent):
        """
        This function returns the sweep of the planform at any given percentage of the cord
        :param c_percent: The percent of the cord (range from 0 to 1)
        :type c_percent: float
        :return: The sweep of the planform at a percentage of the cord [rad]
        :rtype: float
        """
        # need to actually write the function
        sweep = 2 * c_percent
        return sweep

    def chord(self, y):
        """
        This function returns the cord of the planform at any given distance away form the root
        :param y: Distance away from the root [m]
        :type y: float
        :return: Cord length [m]
        :rtype: float
        """
        return self.cr - 2 * y * (self.cr - self.ct) / self.b


class Stringer:
    def __init__(self, t, w, h, x_stop, material):
        """
        Initiate variable of type stringer.
        :param t: Thickness of stringer in [m]
        :type t: float
        :param w: Width of stringer in [m]
        :type w: float
        :param h: Height of stringer in [m]
        :type h: float
        :param material: The material of the stringer
        :type material: Material
        :param x_stop: The distance form the root at which the stringer ends
        :type x_stop: float
        """
        self.t = t
        self.x_stop = x_stop
        self.w = w
        self.h = h
        self.material = material

    def area(self):
        """
        This function calculates the cross-sectional area of the stringer.
        :return: Cross-sectional area of the stringer [m^2]
        :rtype: float
        """
        return self.h * self.t + self.w * self.t - self.t * self.t

    def mass(self):
        """
        :return: Mass per unit length of the stringer [kg/m]
        :rtype: float
        """
        return self.material.rho * self.area()


# same as the planform and material such for wingbox geometry
class WingBox:
    def __init__(self, t, stringers, material):
        """
        Initiate variable of type wingbox
        :param t: Thickness of the sheets used [m]
        :type t: float
        :param n: Number of stringer (min is 4)
        :type n: int
        :param stringers: The list of stringers used
        :type stringers: list of Stringer
        :param material: material of wingbox
        :type material: Material
        """
        self.a = 0
        self.thickness = t
        self.stringers = stringers
        self.planemass = (10220.5 / 1.3)
        self.material = material

    def width(self, planform, y):
        """
        This function returns the wingbox width from the wing geometry at any given distance from the root
        :param planform: The planform used
        :type planform: Planform
        :param y: Distance away from root [m]
        :type y: float
        :return: Wingbox width [m]
        :rtype: float
        """
        return planform.chord(y) * planform.spar_dif

    def height(self, planform, y):
        """
        This function returns the wingbox height from the wing geometry at any given distance from the root
        :param planform: The planform used
        :type planform: Planform
        :param y: The distance away from the root [m]
        :type y: float
        :return: Wingbox height [m]
        :rtype: float
        """
        return 0.1296 * planform.chord(y)

    def moment_of_inertia(self, y):
        """
        This function returns the moment of inertia for the wingbox at given distance from root with given stringers
        :param y: The distance away from the root [m]
        :type y: float
        :return: Wingbox moment of inertia [m^4]
        :rtype: float
        """
        stringer_moment = 0
        for stringer in self.stringers:
            stringer_moment += (y < stringer.x_stop) * (stringer.area() * self.height(y) ** 2)
        return 2 * self.width(y) * self.thickness * self.height(y) ** 2 + self.thickness * 8 * self.height(
            y) ** 3 / 6 + stringer_moment

    def torsional_constant(self, y):
        """
        This function returns the torsional constant J for the wingbox at given distance from root
        :param y: The distance away from the root [m]
        :type y: float
        :return: The torsional constant [m^4]
        :rtype: float
        """
        return 4 * (2 * self.width(y) * self.height(y)) ** 2 / (
                2 * self.width(y) / self.thickness + 4 * self.height(y) / self.thickness)

    def cross_section(self, planform, x):
        """
        This function calculates the cross-section of the wingbox
        :param planform: The planform used
        :type planform: Planform
        :param x: distance from the root [m]
        :type x: float
        :return: cross-sectional area of the wingbox in [m^2]
        :rtype: float
        """
        stringer_area = 0
        for stringer in self.stringers:
            stringer_area += (x < stringer.x_stop) * stringer.area()
        return 2 * self.height(planform, x) * self.thickness + 2 * self.width(planform, x) * self.thickness \
               + stringer_area + (2 * self.height(planform, x) + 2 * self.width(planform, x)) \
               * 1.816 * (8 / 10000)

    def mass_distribution(self, planform, x):
        """
        This function returns the mass per unit length (kilograms per meter) a given distance (in meters) away from the root
        This also includes the weight of the rest of the aircraft at the root of the wing
        :param planform: wing planform dimensions
        :type planform: Planform
        :param x: Distance away from the root [m]
        :type x: float
        :return: Mass per unit length [kg/m]
        :rtype: float
        """
        return self.cross_section(planform, x) * self.material.rho + (x == 0) * self.planemass + (x == 2) * 1000 #+ (4 < x < 5) * 10000


# function definition list

def lift_distribution(x):
    """
    This function returns the lift per unit length (in newton per meter) a given distance (in meters) away from the root
    :param x: Distance away from the root [m]
    :type x: float
    :return: A lift force per unit length [N/m]
    :rtype: float
    """

    # This function is to be writen by the team responsible for the data collection
    lift = 1000*math.sqrt(1-(2*x/34)**2)
    # lift = 7319.5
    return lift


def shear_force(x, wingbox, planform):
    shear = \
        integrate.quad(lambda b: lift_distribution(b) - g * wingbox.mass_distribution(planform, b), x, planform.b / 2)[
            0]
    return shear


def bending_moment(x, wingbox, planform):
    """
    This function returns the bending moment (in newton meter) a given distance (in meters) away from the root
    :param x: Distance away from the root [m]
    :type x: float
    :return: Bending moment [Nm]
    :rtype: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    """

    # According to the Mechanics of Materials the bending moment should be a double integral of the distributed load
    # intensity (lift_distribution - g*(mass_distribution)). Should discuss that in the session

    moment = integrate.quad(lambda a: shear_force(a, wingbox, planform), x, planform.b / 2)[0]
    return -moment  # minus sign is included for coordinates


# deflection profiles
def lateral_deflection(y, moment, material, wingbox):  # dv/dy , E modulus is for one material (can be improved later)
    """
    This function returns the lateral deflection (v) at y distance away from the root chord,
    :param y: Distance away from the root [m]
    :type y: float
    :param moment: The bending moment function
    :type moment: function
    :param material: Type of material used
    :type: Material
    :param wingbox: The wingbox used
    :type: WingBox
    :return: Lateral deflection [m]
    :type: float
    """
    return -1 * integrate.dblquad(moment / (material.E * wingbox.moment_of_inertia), 0, y, 0, y)


def twist_angle(y, torsion, wingbox, material):  # lower limit must be set for the fuselage
    """
    This function returns the twist angle at y distance away from the root chord
    :param y: Distance away from the root [m]
    :type y: float
    :param torsion: The torsion function
    :type torsion: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param material: The material used
    :type material: Material
    :return: Twist angle [rad]
    :rtype: float
    """
    return integrate.quad(torsion / (wingbox.torsional_constant * material.G), 0, y)

# finish this
# def second_moment_of_inertia(x):
#    moment_of_inertia = 0
#    return moment_of_inertia

# def bending_moment(L, x, a):
#   M = L * (x - a)
#  return M
