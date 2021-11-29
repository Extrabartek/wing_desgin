import numpy as np
import scipy.integrate as integrate
import math
import WP_41 as WP41

# global constant
g = 9.80665
AoA = 10  # Degrees
load_factor = 2.5
rho_fuel = 800  # kg/m^3
fuel = 0
weight_final = 23528


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
        :param spar_front: The percentage of the cord where the front spar is
        :type spar_front: float
        """
        self.b = b
        self.cr = cr
        self.ct = ct
        self.sweep_le = sweep_le
        self.spar_rear = spar_rear
        self.spar_front = spar_front
        self.spar_dif = spar_rear - spar_front

    def area(self):
        """
        This function returns the area (in m^2) of the planform

        :return: Area of the planform [m^2]
        :rtype: float
        """
        return self.b * (self.cr + self.ct) * 0.5

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
    def __init__(self, t, base, flange, x_stop, material):
        """
        Initiate variable of type stringer.

        :param t: Thickness of stringer in [m]
        :type t: float
        :param base: Width of stringer in [m]
        :type base: float
        :param flange: Height of stringer in [m]
        :type flange: float
        :param material: The material of the stringer
        :type material: Material
        :param x_stop: The distance form the root at which the stringer ends
        :type x_stop: float
        """
        self.t = t
        self.x_stop = x_stop
        self.base = base
        self.flange = flange
        self.material = material

    def area(self):
        """
        This function calculates the cross-sectional area of the stringer.

        :return: Cross-sectional area of the stringer [m^2]
        :rtype: float
        """
        return self.flange * self.t + self.base * self.t - self.t * self.t

    def mass(self):
        """
        This function calculates the mass per unit (in kilograms per meter) length of a stringer

        :return: Mass per unit length of the stringer [kg/m]
        :rtype: float
        """
        return self.material.rho * self.area()

    def centroid(self):
        """
        This function calculates the vertical centroid position of an L-stringer

        :return: Distance of the centroid from the attached flange
        :rtype: float
        """
        return ((self.flange / 2) * self.flange * self.t + (self.t / 2) * self.base * self.t) / (
                self.t * (self.flange + self.base))

    def moment_inertia(self):
        """
        This function returns the moment of inertia of the stringers
        :return: moment of inertia of the stringer in m^4
        :rtype: float
        """

        return (self.t * self.flange ** 3) / 12 + self.t * self.flange * (
                self.flange / 2 - self.centroid()) ** 2 + self.base * self.t * self.centroid() ** 2


# same as the planform and material such for wingbox geometry
class WingBox:
    def __init__(self, t_list, stringers, material):
        """
        Initiate variable of type wingbox

        :param t_list: List of the thicknesses of the sheets used [m]
        # first two are the top and bottom plate
        # third and fourth are front and rear spar
        :type t_list: list
        :param stringers: The list of stringers used
        :type stringers: list of Stringer
        :param material: material of wingbox
        :type material: Material
        """
        self.a = 0
        self.thickness = 0.005
        self.thickness_list = t_list
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
        This function returns the wingbox half-height from the wing geometry at any given distance from the root

        :param planform: The planform used
        :type planform: Planform
        :param y: The distance away from the root [m]
        :type y: float
        :return: Wingbox half-height [m]
        :rtype: float
        """
        return 0.1296 * self.width(planform, y)

    def moment_of_inertia(self, planform, y):
        """
        This function returns the moment of inertia for the wingbox at given distance from root with given stringers

        :param y: The distance away from the root [m]
        :type y: float
        :param planform: Planform used
        :type planform: Planform
        :return: Wingbox moment of inertia [m^4]
        :rtype: float
        """
        stringer_moment = 0
        wingbox_moment = 0
        for stringer in self.stringers:
            stringer_moment += (y < stringer.x_stop) * ((stringer.area() * (
                    self.height(planform, y) - stringer.centroid()) ** 2) + stringer.moment_inertia())
        # top plate
        wingbox_moment += self.width(planform, y) * self.thickness_list[0] * self.height(planform, y) ** 2
        # bottom plate
        wingbox_moment += self.width(planform, y) * self.thickness_list[1] * self.height(planform, y) ** 2
        # front spar
        wingbox_moment += self.width(planform, y) * self.thickness_list[2] * self.height(planform, y) ** 2
        # rear spar
        wingbox_moment += self.width(planform, y) * self.thickness_list[3] * self.height(planform, y) ** 2

        return stringer_moment + wingbox_moment

    def torsional_constant(self, y, planform):
        """
        This function returns the torsional constant J for the wingbox at given distance from root

        :param planform: The planform used
        :type planform: Planform
        :param y: The distance away from the root [m]
        :type y: float
        :return: The torsional constant [m^4]
        :rtype: float
        """

        return (4 * (self.width(planform, y) * 2 * self.height(planform, y)) ** 2) / (
                (2 * self.width(planform, y) / self.thickness) + (4 * self.height(planform, y) / self.thickness))

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
        return 4 * self.height(planform, x) * self.thickness + 2 * self.width(planform,
                                                                              x) * self.thickness + stringer_area + (
                           4 * self.height(planform, x) + 2 * self.width(planform, x)) * 1.816 * (1 / 1000)

    def mass_distribution(self, planform, x):
        """
        This function returns the mass per unit length (kilograms per meter) a given distance (in meters) away from the
        root

        :param planform: wing planform dimensions
        :type planform: Planform
        :param x: Distance away from the root [m]
        :type x: float
        :return: Mass per unit length [kg/m]
        :rtype: float
        """
        mass = self.cross_section(planform, x) * self.material.rho + (fuel == 1) * 2 * self.height(planform,
                                                                                                   x) * self.width(
            planform, x) * rho_fuel
        return mass

    def Q(self, planform, x):
        """
        This function calculates the first moment of area

        :param planform: The planform used
        :type planform: Planform
        :param x: The distance from the root [m]
        :type x: float
        :return: This function returns the first moment of area in [m^3]
        :rtype: float
        """
        q = self.height(planform, x) ** 2 * self.thickness + self.height(planform, x) * self.width(planform,
                                                                                                   x) * self.thickness
        for stringer in self.stringers:
            q += stringer.area() * (self.height(planform, x) - stringer.centroid()) * (stringer.x_stop > x)
        return q

    def total_weight(self, planform):
        return integrate.quad(lambda a: self.mass_distribution(planform, a), 0, planform.b / 2, epsrel=0.3, epsabs=0.3)[
            0]


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
    # lift = 1000 * math.sqrt(1 - (2 * x / 34) ** 2)
    lift = 7319.5
    return lift


def shear_force(x, wingbox, planform):
    """
    This function returns the values of the shear per unit length (in newton per meter) at a given distance (in meters)
    from the root

    :param x: The distance away from the root [m]
    :type x: float
    :param wingbox: The wingbox used for the calculation
    :type wingbox: WingBox
    :param planform: The planform used for the calculation
    :type planform: Planform
    :return: The shear force per unit length [N/m]
    :rtype: float
    """
    shear = integrate.quad(lambda a: (WP41.Ndis0(a, AoA)), x, planform.b / 2, epsabs=0.3, epsrel=0.3)[0] - \
            g * np.cos(np.radians(AoA)) * load_factor * \
            integrate.quad(lambda a: wingbox.mass_distribution(planform, a), x, planform.b / 2, epsabs=0.3, epsrel=0.3)[
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

    moment = integrate.quad(lambda a: shear_force(a, wingbox, planform), x, planform.b / 2, epsrel=0.5, epsabs=0.5)[0]
    return -moment  # minus sign is included for coordinates


def torsion(x, planform):
    """
    This function returns the torsion per unit area at any distance from the root chord

    :param x: Distance away from the root [m]
    :type x: float
    :return: Torsion per unit area
    :rtype: float
    :param planform: The planform used
    """
    return (1 / 8) * planform.chord(x) * WP41.Ndis0(x, AoA) + WP41.Mdis(x, AoA)


# deflection profiles
def slope_deflection(y, material, wingbox, planform):  # dv/dy , E modulus is for one material (can be improved later)
    """
    This function returns the lateral deflection (v) at y distance away from the root chord,

    :param planform: Planform used
    :type planform: Planform
    :param y: Distance away from the root [m]
    :type y: float
    :param material: Type of material used
    :type: Material
    :param wingbox: The wingbox used
    :type: WingBox
    :return: Lateral deflection [m]
    :type: float
    """
    return -1 * integrate.quad(
        lambda b: bending_moment(b, wingbox, planform) / wingbox.moment_of_inertia(planform, b), 0, y,
        epsabs=0.5, epsrel=0.5)[0]


def vertical_deflection(y, material, wingbox, planform):
    """
    This function calculates the vertical deflection in m at a certain distance from the root.

    :param y: distance from root [m]
    :type y: float
    :param material: the material the thing is made of
    :type material: Material
    :param wingbox: the wingbox used
    :type wingbox: WingBox
    :param planform: the planform used
    :type planform: Planform
    :return: the vertical deflection in [m]
    :rtype: float
    """

    return -1 / material.E * integrate.quad(
        lambda b:
        integrate.quad(lambda a: bending_moment(a, wingbox, planform) / wingbox.moment_of_inertia(planform, a), 0, b,
                       epsabs=0.5, epsrel=0.5)[0], 0, y,
        epsabs=0.5, epsrel=0.5)[0]


def twist_angle(x, wingbox, material, planform):  # lower limit must be set for the fuselage
    """
    This function returns the twist angle at y distance away from the root chord

    :param planform: The planform used for the calculations
    :type planform: Planform
    :param x: Distance away from the root [m]
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param material: The material used
    :type material: Material
    :return: Twist angle [rad]
    :rtype: float
    """
    return \
        integrate.quad(lambda a: (torsion(a, planform) / (wingbox.torsional_constant(a, planform) * material.G)), 0, x)[
            0]


def normal_stress(x, wingbox, planform):
    """
    This function returns the maximum normal stress (in Pascal) at a given distance away from the root (in meters)

    :param x: Distance from the root [m]
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The maximum normal stress at a given distance [Pa]
    """
    return bending_moment(x, wingbox, planform) * wingbox.height(planform, x) / (wingbox.moment_of_inertia(planform, x))


def shear_stress(x, wingbox, planform):
    """
    This function calculates the shear stress due to torsion

    :param x:distance from root
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: the planform used
    :type planform: Planform
    :return: The shear stress in four places[N/m^2]
    :rtype: list
    """
    point_edge_1 = (shear_force(x, wingbox, planform) * wingbox.width(planform, x) * wingbox.height(planform, x)) / \
                   wingbox.moment_of_inertia(planform, x) - torsion(x, planform) / (
                           2 * wingbox.thickness * (2 * wingbox.height(planform, x) * wingbox.width(planform, x)))
    point_edge_2 = (shear_force(x, wingbox, planform) * wingbox.width(planform, x) * wingbox.height(planform, x)) / \
                   wingbox.moment_of_inertia(planform, x) + torsion(x, planform) / (
                           2 * wingbox.thickness * (2 * wingbox.height(planform, x) * wingbox.width(planform, x)))
    point_mid_1 = -1 * torsion(x, planform) / (
            2 * wingbox.thickness * (2 * wingbox.height(planform, x) * wingbox.width(planform, x))) \
                  + (shear_force(x, wingbox, planform) * wingbox.Q(planform, x)) / (
                          2 * wingbox.moment_of_inertia(planform, x) * wingbox.thickness)
    point_mid_2 = torsion(x, planform) / (
            2 * wingbox.thickness * (2 * wingbox.height(planform, x) * wingbox.width(planform, x))) \
                  + (shear_force(x, wingbox, planform) * wingbox.Q(planform, x)) / (
                          2 * wingbox.moment_of_inertia(planform, x) * wingbox.thickness)
    return point_mid_1, point_mid_2, point_edge_1, point_edge_2


def tau_max(x, wingbox, planform):
    """
    This function will return the maximum shear stress (according to the Mohr's circle) at a given distance away from
    the root (in meters)

    :param x: Distance away from the root [m]
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The maximum shear stress [N/m^2]
    :rtype: float
    """
    # zero line
    point1 = math.sqrt(shear_stress(x, wingbox, planform)[2] ** 2)
    point2 = math.sqrt(shear_stress(x, wingbox, planform)[3] ** 2)

    # edge
    point3 = math.sqrt(shear_stress(x, wingbox, planform)[1] ** 2 + (normal_stress(x, wingbox, planform) / 2) ** 2)
    point4 = math.sqrt(shear_stress(x, wingbox, planform)[1] ** 2 + (normal_stress(x, wingbox, planform) / 2) ** 2)

    return max(point1, point2, point3, point4)


def optimize_stringers(x, wingbox1, planform1):
    """
    This function will optimize the number of stringers (find the minimal safe number) at a given distance away from the
    root (in meters)

    :param x: Distance away from the root [m]
    :type x: float
    :param wingbox1: The wingbox used
    :type wingbox1: WingBox
    :param planform1: The planform used
    :type planform1: Planform
    :return: Returns the minimal safe number of stringers
    :rtype: list
    """

    stringer_used = wingbox1.stringers[0]
    stringer_list = []
    wingbox_cal = WingBox(wingbox1.thickness_list, wingbox1.stringers, wingbox1.material)
    planform_cal = Planform(planform1.b, planform1.cr, planform1.ct, planform1.sweep_le, planform1.spar_rear,
                            planform1.spar_front)
    # The following part is to make sure that only 4 full length stringer exist
    for a in range(4):
        stringer_list.append(stringer_used)

    wingbox_cal.stringers = stringer_list
    done = True

    stress = tau_max(x, wingbox_cal, planform_cal)
    while stress > (wingbox_cal.material.sig_yld / 2) * 0.8:
        if len(stringer_list) * stringer_used.base > wingbox_cal.width(planform_cal, 0) or len(stringer_list) > 60:
            done = False
            break
        for b in range(2):
            stringer_list.append(stringer_used)
        wingbox_cal.stringers = stringer_list
        stress = tau_max(x, wingbox_cal, planform_cal)

    return len(wingbox_cal.stringers), done


def stringer_length_conversion(number_of_stringer_list, stringer, step_size, rangy):
    """
    This function will create a list of stringers that can be used to crate a wingbox

    :param rangy: The rage used for the calculations
    :type rangy: list
    :param number_of_stringer_list: List of the number of stringers at a given distance
    :type number_of_stringer_list: list of float
    :param stringer: The stringers used
    :type stringer: Stringer
    :param step_size: A step size for the calculations
    :type step_size: float
    :return: A list of stringers
    :rtype: list of Stringer
    """

    list_stringers = []

    for a in range(int(max(number_of_stringer_list))):
        list_stringers.append(Stringer(stringer.t, stringer.base, stringer.flange, 0, stringer.material))

    for x in rangy:
        for a in range(int(number_of_stringer_list[int(x / step_size)])):
            list_stringers[a].x_stop = x + step_size

    return list_stringers


def dynamic_pressure(wingbox, planform):
    """
    This function calculates the dynamic pressure needed for a certain load factor
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The dynamic pressure in [kg/m*s^2]
    """
    q = load_factor * g * ((weight_final / 2)
                           + integrate.quad(lambda a: wingbox.mass_distribution(planform, a), 0, planform.b / 2,
                                            epsrel=0.1, epsabs=0.1)[0]) / \
        (integrate.quad(lambda a: WP41.c(a) * (
                (AoA == 0) * WP41.cly1(a) + (AoA == 10) * WP41.cly2(a) + (AoA == -10) * WP41.cly3(a)), 0,
                        planform.b / 2)[0])

    return q


def column_buckling(x, stringer):
    """
    This function returns the critical stress of the stringer column at a distance X

    :param x: Distance away from the root [m]
    :type x: float
    :param stringer: The type of stringer used
    :type stringer: Stringer
    :return: Critical stress for column buckling [Pa]
    :rtype: float
    """
    return ((math.pi ** 2) * stringer.material.E * stringer.moment_inertia()) / (x ** 2 * stringer.area())
