import scipy.integrate as integrate

# global constant
g = 9.81


# class definition list
class Material:
    # if any parameters need to be added follow the method for the already added parameters
    def __init__(self, rho, sig_yld, sig_ult):
        """
        Initiate variable of type material.
        :param rho: Density of the material [kg/m^3]
        :type rho: float
        :param sig_yld: Yield strength of the material [N/m^2]
        :type sig_yld: float
        :param sig_ult: Ultimate strength of the material [N/m^2]
        :type sig_ult: float
        """
        self.rho = rho
        self.sig_yld = sig_yld
        self.sig_ult = sig_ult

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
    def __init__(self, t, w, h, material):
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
        """
        self.t = t
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
    def __init__(self, t, n, stringer):
        """
        Initiate variable of type wingbox
        :param t: Thickness of the sheets used [m]
        :type t: float
        :param n: Number of stringer (min is 4)
        :type n: int
        :param stringer: The stringer type used
        :type stringer: Stringer
        """
        self.a = 0
        self.thickness = t
        self.n_stringers = n  # Minimum possible number is 4
        self.stringer = stringer
        self.Wp = (10220.5 / 1.3) * g

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

    def moment_of_inertia(self, stringer, y):
        """
        This function returns the moment of inertia for the wingbox at given distance from root with given stringers
        :param stringer: The stringer used
        :type stringer: Stringer
        :param y: The distance away from the root [m]
        :type y: float
        :return: Wingbox moment of inertia [m^4]
        :rtype: float
        """
        return 2 * self.width(y) * self.thickness * self.height(y) ** 2 + self.thickness * 8 * self.height(
            y) ** 3 / 6 + self.n_stringers * stringer.area() * self.height(y) ** 2

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
        return 2 * self.height(planform, x) * self.thickness + 2 * self.width(planform,
                                                                              x) * self.thickness + self.stringer.area() * self.n_stringers

    def mass_distribution(self, x):
        """
        This function returns the mass per unit length (kilograms per meter) a given distance (in meters) away from the root
        This also includes the weight of the rest of the aircraft at the root of the wing
        :param x: Distance away from the root [m]
        :type x: float
        :return: Mass per unit length [kg/m]
        :rtype: float
        """

        # This function has to be decided on, but the "sliced" approached should make pretty easy to combine with the
        # lift distribution

        mass = 0
        return mass


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

    lift = 0
    return lift


def bending_moment(x):
    """
    This function returns the bending moment (in newton meter) a given distance (in meters) away from the root
    :param x: Distance away from the root [m]
    :type x: float
    :return: Bending moment [Nm]
    :rtype: float
    """

    # According to the Mechanics of Materials the bending moment should be a double integral of the distributed load
    # intensity (lift_distribution - g*(mass_distribution)). Should discuss that in the session

    moment = 0
    return moment


# deflection profiles
def lateral_deflection(moment):  # dv/dy
    integrate.dblquad()

# finish this
# def second_moment_of_inertia(x):
#    moment_of_inertia = 0
#    return moment_of_inertia

# def bending_moment(L, x, a):
#   M = L * (x - a)
#  return M
