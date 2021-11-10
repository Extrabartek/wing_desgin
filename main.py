import scipy.integrate as integrate

# global constant
g = 9.81


# class definition list
class Material:
    """
    This is a class that will make the base for all the materials that will be used in the project
    :param self.rho: Density of the material [kg/m^3]
    :param self.sig_yld: Yield strength of the material [N/m^2]
    :param self.sig_ult: Ultimate strength of the material [N/m^2]
    """

    # if any parameters need to be added follow the method for the already added parameters
    def __init__(self, rho, sig_yld, sig_ult):
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
    """
    This is a class that will make the base for all the planforms that will be used in the project
    :param self.b: The wingspan [m^2]
    :param self.cr: The cord at the root [m]
    :param self.ct: The cord at the tip [m]
    :param self.sweep_le: Sweep at the leading edge [rad]
    """

    # if any parameters need to be added follow the method for the already added parameters
    def __init__(self, b, cr, ct, sweep_le):
        self.b = b
        self.cr = cr
        self.ct = ct
        self.sweep_le = sweep_le

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
        return self.cr - 2 * y * (self.cr - self.ct) / self.b


# same as the planform and material such for wingbox geometry
class WingBox:
    def __init__(self,t,n,spar_dif):
        self.a = 0
        self.thickness = t
        self.n_stringers = n

    def box_width(self,):
        Planform.chord(y) * spar_dif

    box_height = 0.1296 * Planform.chord(y)

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


def mass_distribution(x):
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


# finish this
def second_moment_of_inertia(x):
    moment_of_inertia = 0
    return moment_of_inertia

# def bending_moment(L, x, a):
#   M = L * (x - a)
#  return M
