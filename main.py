import numpy as np
import scipy.integrate as integrate
import math
import WP_41 as WP41

# Global constant
g = 9.80665
AoA = 10  # Degrees
load_factor = 2.5
rho_fuel = 800  # kg/m^3
fuel = 1
weight_final = 23528


class Material:
    def __init__(self, rho, sig_yld, sig_ult, E, G, nu):
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
        :param nu: The Poisson's ratio of the material
        :type nu: float
        """
        self.rho = rho
        self.sig_yld = sig_yld

        self.sig_ult = sig_ult
        self.E = E
        self.G = G
        self.nu = nu

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
    def __init__(self, t, a, x_stop, material):
        """
        Initiate variable of type stringer.

        :param t: Thickness of stringer in [m]
        :type t: float
        :param a: The perimeter of the stringer [m]
        :type a: float
        :param material: The material of the stringer
        :type material: Material
        :param x_stop: The distance form the root at which the stringer ends
        :type x_stop: float
        """
        self.t = t
        self.x_stop = x_stop
        self.a = a
        self.material = material

    def area(self):
        """
        This function calculates the cross-sectional area of the stringer.

        :return: Cross-sectional area of the stringer [m^2]
        :rtype: float
        """
        return self.a * self.t

    '''
    def variable_area(self, h):
        """
        This function gives the area (in meters squared) if the stringer was cut at given distance (in meter)

        :param h: The distance from the mounted part of the stringer [m]
        :type h: float
        :return: The area of the stringer is it was cut at a given distance h [m].
        :rtype: float
        """
        area = 0
        if h < 0:
            return 0

        if h < self.a / 4:

            if h < self.t:
                area = h * self.a / 4
            elif self.t < h < (self.a / 4 - self.t):
                area = self.t * self.a / 4 + 2 * self.t * h
            elif (self.a / 4 - self.t) < h < self.a / 4:
                area = self.t * self.a / 4 + 2 * self.t * (self.a / 4) + (self.a / 4 - h) * self.a / 4
        else:
            area = self.area

        return area
    '''

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

        :return: Distance of the centroid from the attached base
        :rtype: float
        """
        return self.a / 8

    # def variable_centroid(self, h):

    def moment_inertia(self):
        """
        This function returns the moment of inertia of the stringers

        :return: moment of inertia of the stringer in m^4
        :rtype: float
        """

        # return (self.a / 4 * (self.a / 4) ** 3 - (self.a / 4 - 2 * self.t) * (self.a / 4 - 2 * self.t) ** 3) / 12

        return (self.a ** 3 * self.t) / 96


class WingBox:
    def __init__(self, t_list, stringers_top, stringers_bottom, rib_list, material):
        """
        Initiate variable of type wingbox

        :param t_list: List of the thicknesses of the sheets used [m]
        # first two are the top and bottom plate
        # third is the thickness of the spars
        :type t_list: list
        :param stringers_top: The list of stringers used at the top of the wingbox
        :type stringers_top: list of Stringer
        :param stringers_bottom: The list of stringers used at the bottom of the wingbox
        :type stringers_bottom: list of Stringer
        :param rib_list: The list with the locations of the ribs (distance from the root in meter)
        :type rib_list: list
        :param material: material of wingbox
        :type material: Material
        """

        self.t_top = t_list[0]
        self.t_bottom = t_list[1]
        self.t_spar = t_list[2]
        self.stringers_top = stringers_top
        self.stringers_bottom = stringers_bottom
        self.rib_list = rib_list
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
        This function returns the wingbox half-height (in meters) from the wing geometry at any given
        distance from the root (in meters)

        :param planform: The planform used
        :type planform: Planform
        :param y: The distance away from the root [m]
        :type y: float
        :return: Wingbox half-height [m]
        :rtype: float
        """
        return 0.1296 * self.width(planform, y)

    def centroid(self, planform, y):
        """
        This function returns the position of the centroid (distance from the bottom sheet in meters)
        at a given distance from the root (in meters)

        :param y: Distance from the root [m]
        :type y: float
        :param planform: Planform used
        :type planform: Planform
        :return: Distance of the centroid from the bottom sheet [m]
        :rtype: float
        """

        total_area = 0
        total_product = 0

        for stringer in self.stringers_bottom:
            total_area += stringer.area() * (y < stringer.x_stop)
        for stringer in self.stringers_top:
            total_area += stringer.area() * (y < stringer.x_stop)

        total_area += (
                self.width(planform, y) * (self.t_top + self.t_bottom) + 4 * self.height(planform, y) * self.t_spar)

        for stringer in self.stringers_bottom:
            total_product += stringer.area() * (stringer.centroid() + self.t_bottom) * (y < stringer.x_stop)
        for stringer in self.stringers_top:
            total_product += stringer.area() * (2 * self.height(planform, y) - stringer.centroid() - self.t_top) * (
                    y < stringer.x_stop)

        total_product += (2 * self.height(planform, y) * (
                2 * self.height(planform, y) * self.t_spar + self.width(planform, y) * self.t_top))

        return total_product / total_area

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
        for stringer in self.stringers_top:
            stringer_moment += (y < stringer.x_stop) * (stringer.moment_inertia() + stringer.area() * (
                    2 * self.height(planform, y) - self.centroid(planform,
                                                                 y) - stringer.centroid() - self.t_top) ** 2)
        for stringer in self.stringers_bottom:
            stringer_moment += (y < stringer.x_stop) * (stringer.moment_inertia() + stringer.area() * (
                    self.centroid(planform, y) - stringer.centroid() - self.t_bottom) ** 2)
        # top plate
        wingbox_moment += self.width(planform, y) * self.t_top * (
                2 * self.height(planform, y) - self.centroid(planform, y)) ** 2
        # bottom plate
        wingbox_moment += self.width(planform, y) * self.t_bottom * self.centroid(planform, y) ** 2
        # front spars
        wingbox_moment += (4 / 3) * self.height(planform, y) ** 3 * self.t_spar + 4 * self.height(planform, y) * \
                          self.t_spar * (self.centroid(planform, y) - self.height(planform, y)) ** 2

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
                (self.width(planform, y) / self.t_top) + (self.width(planform, y) / self.t_bottom) + (
                4 * self.height(planform, y) / self.t_spar))

    def cross_section(self, planform, y):
        """
        This function calculates the cross-section of the wingbox

        :param planform: The planform used
        :type planform: Planform
        :param y: distance from the root [m]
        :type y: float
        :return: cross-sectional area of the wingbox in [m^2]
        :rtype: float
        """
        total_area = 0

        for stringer in self.stringers_bottom:
            total_area += stringer.area() * (y < stringer.x_stop)
        for stringer in self.stringers_top:
            total_area += stringer.area() * (y < stringer.x_stop)

        total_area += (
                self.width(planform, y) * (self.t_top + self.t_bottom) + 4 * self.height(planform, y) * self.t_spar)

        return total_area + 2 * (self.width(planform, y) + 2 * self.height(planform, y)) * 1.816 * (1 / 1000)

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

    def Q(self, planform, h, y):
        """
        This function calculates the first moment of area

        :param planform: The planform used
        :type planform: Planform
        :param y: The distance from the root [m]
        :type y: float
        :param h: The chosen distance from the bottom plate [m]
        :type h: float
        :return: This function returns the first moment of area in [m^3]
        :rtype: float
        """
        yprime1 = abs(self.centroid(planform, y) - (self.height(planform, y) - (h / 2)))
        yprime2 = ((2 * self.height(planform, y) - self.centroid(planform,
                                                                 y)) - self.t_top / 2)  # centroid of the top plate
        yprime3 = (2 * self.height(planform, y) - self.stringers_top[0].centroid() - self.t_top - self.centroid(
            planform, y))  # centroid of one of the top stringers
        product_area = len(self.stringers_top) * yprime3 * self.stringers_top[
            0].area() + yprime2 * self.t_top * self.width(planform, y)
        # spar contribution
        q_total = 2 * yprime1 * (2 * self.height(planform, y) - h) * self.t_spar + product_area
        return q_total

    def total_weight(self, planform, material):
        total_rib = 0
        for rib in self.rib_list:
            total_rib += self.rib_mass(rib, planform, material)

        return integrate.quad(lambda a: self.mass_distribution(planform, a),
                              0, planform.b / 2, epsrel=0.3, epsabs=0.3)[0] + total_rib

    def rib_mass(self, y, planform, material):
        return 2 * self.height(planform, y) * self.width(planform, y) * 0.002 * material.rho


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
    shear = integrate.quad(lambda a: (WP41.Ndis0(a, AoA)), x, planform.b / 2,
                           epsabs=0.3, epsrel=0.3)[0] - \
            g * np.cos(np.radians(AoA)) * load_factor * \
            integrate.quad(lambda a: wingbox.mass_distribution(planform, a),
                           x, planform.b / 2, epsabs=0.3, epsrel=0.3)[
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

    moment = integrate.quad(lambda a: shear_force(a, wingbox,
                                                  planform), x, planform.b / 2, epsrel=0.5, epsabs=0.5)[0]
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
    # return (1 / 8) * planform.chord(x) * WP41.Ndis0(x, AoA) + WP41.Mdis(x, AoA)
    return \
        integrate.quad(lambda a: (1 / 8) * planform.chord(a) * WP41.Ndis0(a, AoA) + WP41.Mdis(a, AoA), x,
                       planform.b / 2,
                       epsabs=0.5, epsrel=0.5)[0]


def slope_deflection(y, wingbox, planform):
    """
    This function returns the lateral deflection (v) at y distance away from the root chord,

    :param planform: Planform used
    :type planform: Planform
    :param y: Distance away from the root [m]
    :type y: float
    :param wingbox: The wingbox used
    :type: WingBox
    :return: Lateral deflection [m]
    :type: float
    """
    return -1 * integrate.quad(
        lambda b: bending_moment(b, wingbox, planform) \
                  / wingbox.moment_of_inertia(planform, b), 0, y,
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
        integrate.quad(lambda a: bending_moment(a, wingbox,
                                                planform) / wingbox.moment_of_inertia(planform, a), 0, b,
                       epsabs=0.5, epsrel=0.5)[0], 0, y,
        epsabs=0.5, epsrel=0.5)[0]


def twist_angle(x, wingbox, material, planform):
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


def normal_stress(x, wingbox, planform, h):
    """
    This function returns the maximum normal stress (in Pascal) at a given distance away from the root (in meters)

    :param x: Distance from the root [m]
    :type x: float
    :param h: height from the bottom plate
    :type h: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The maximum normal stress at a given distance [Pa]
    """
    return bending_moment(x, wingbox, planform) * (h - wingbox.centroid(planform, x)) / (
        wingbox.moment_of_inertia(planform, x))


def shear_stress(x, h, wingbox, planform):
    """
    This function calculates the shear stress due to torsion

    :param x:distance from root
    :type x: float
    :param h: height from the bottom plate
    :type h: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: the planform used
    :type planform: Planform
    :return: The shear stress in four places[N/m^2]
    :rtype: float
    """
    shear_stress_result = abs((shear_force(x, wingbox, planform) * wingbox.Q(planform, h, x)) / (
            2 * wingbox.moment_of_inertia(planform, x) * wingbox.t_spar)) + abs(
        torsion(x, planform) / (4 * wingbox.height(planform, x) * wingbox.width(planform, x) * wingbox.t_spar))
    return shear_stress_result


def tau_max(x, wingbox, planform):
    """
    This function will return the maximum shear stress (according to
    the Mohr's circle) at a given distance away from the root (in meters)

    :param x: Distance away from the root [m]
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The maximum shear stress [N/m^2]
    :rtype: list
    """
    max_list = []
    step_size = wingbox.height(planform, x) / 10

    for h in np.arange(0, (2 * wingbox.height(planform, x) + step_size * 0.95), step_size):
        max_list.append(math.sqrt(
            shear_stress(x, h, wingbox, planform) ** 2 + (1 / 2 * normal_stress(x, wingbox, planform, h)) ** 2))

    return max(max_list), np.argmax(max_list)


def optimize_stringers(x, wingbox1, planform1):
    """
    This function will optimize the number of stringers (find the minimal safe number)
    at a given distance away from the root (in meters)

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
    planform_cal = Planform(planform1.b, planform1.cr, planform1.ct,
                            planform1.sweep_le, planform1.spar_rear, planform1.spar_front)
    # The following part is to make sure that only 4 full length stringer exist
    for a in range(4):
        stringer_list.append(stringer_used)

    wingbox_cal.stringers = stringer_list
    done = True

    stress = tau_max(x, wingbox_cal, planform_cal)
    while stress > (wingbox_cal.material.sig_yld / 2) * 0.8:
        if len(stringer_list) * stringer_used.base > wingbox_cal.width(planform_cal, 0) \
                or len(stringer_list) > 60:
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

    for a in range(math.ceil(max(number_of_stringer_list))):
        list_stringers.append(Stringer(stringer.t, stringer.a, 0, stringer.material))

    for x in rangy:
        for a in range(math.ceil(number_of_stringer_list[int(x / step_size)])):
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
                           + integrate.quad(lambda a: \
                                                wingbox.mass_distribution(planform, a), 0,
                                            planform.b / 2, epsrel=0.1, epsabs=0.1)[0]) / \
        (integrate.quad(lambda a: WP41.c(a) * (
                (AoA == 0) * WP41.cly1(a) + (AoA == 10) * WP41.cly2(a) +
                (AoA == -10) * WP41.cly3(a)), 0, planform.b / 2)[0])

    return q


def rib_spacing_column(x, wingbox, planform):
    """
    This function returns the rib spacing to meet a certain critical stress


    :param x: Distance away from the root [m]
    :type x: float
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :param planform: The planform used
    :type planform: Planform
    :return: The spacing of the next rib
    :rtype: float
    """

    K = 1
    E = wingbox.stringers_top[0].material.E
    # a = stringer.a  # Stringer a [m]
    # t = stringer.t # Stringer thickness [m]
    interval = 0.5
    stepsize = 0.1

    height = (load_factor == 2.5) * wingbox.height(planform, x) * 2

    stress_list = []

    for x in np.arange(x, x + interval, stepsize):
        stress_list.append(abs(normal_stress(x, wingbox, planform, height)))
    stress = max(stress_list)

    sigma_cr = 1.1 * stress  # Critical stress defined [Pa]
    A = wingbox.stringers_top[0].area()  # Area of the stringer [m^2]
    I = wingbox.stringers_top[0].moment_inertia()
    # I = 2 * a / 8 * t * a ** 2 / 64 + a / 4 * t * a ** 2 / 64 + 2 / 12 * a ** 3 / 64 * t  # Moment of inertia of the stringer [m^4]

    M_max = 1 * 10 ** 6
    # sigma_M = M * y / I
    L = np.sqrt(K * np.pi ** 2 * E * I / (sigma_cr * A))  # Spacing of the ribs [m]
    print('Rib spacing is', L)
    return L


def vertstringer_spacing_web(material, wingbox, planform, x):
    '''
    This function gives the vertical stringer spacing required to account for web buckling
    :param material: material used
    :param wingbox: wingbox used
    :param planform: planform used
    :param x: distance along the wing
    :return:
    '''
    ks = 4.5
    interval = 0.5
    stepsize = 0.1

    stress_list = []

    for x in np.arange(x, x + interval, stepsize):
        stress_list.append(shear_stress(x, wingbox.height(planform, x), wingbox, planform))
    stress = max(stress_list)
    b = np.pi * wingbox.t_spar * np.sqrt((ks * material.E) / (12 * (1 - material.nu ** 2) * 1.1 * stress))
    a = wingbox.height(planform, 0) * 2
    output = b
    if a / b < 1:
        b = wingbox.height(planform, x) * 2
        a = np.pi * wingbox.t_spar * np.sqrt((ks * material.E) / (12 * (1 - material.nu ** 2) * 1.1 * stress))
        output = a
    return output


def column_buckling(x, stringer, riblist, planform):
    """
    This function returns the critical stress of the stringer column at a distance X

    :param x: Distance away from the root [m]
    :type x: float
    :param stringer: The type of stringer used
    :type stringer: Stringer
    :return: Critical stress for column buckling [Pa]
    :rtype: float
    """

    K = 1.  # Factor for end conditions of stringers
    L = 0

    if x < riblist[0]:
        L = riblist[0]
    elif x > riblist[-1]:
        L = planform.b / 2 - riblist[-1]
    else:
        for i in range(len(riblist)):
            if x <= riblist[i] and x >= riblist[i - 1]:
                L = riblist[i] - riblist[i - 1]

    stress = (K * (np.pi ** 2) * stringer.material.E * stringer.moment_inertia()) / (L ** 2 * stringer.area())
    return stress


def skin_buckling_stringer_count(y, wingbox, planform):
    """
    This function returns the critical stress for skin buckling to occur

    :param y: Distance from the root [m]h
    :type y: float
    :param planform: The planform used
    :type planform: Planform
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :return: The critical stress for skin buckling to occur [Pa]
    :rtype: float and bool
    """

    k_c = 4  # assumed to be constant

    number_of_stringers = len(wingbox.stringers_top) * (load_factor == 2.5) + len(wingbox.stringers_bottom) * (load_factor == -1)
    b_basic = (wingbox.width(planform, y) - number_of_stringers * (wingbox.stringers_top[0].a / 2)) / (
            number_of_stringers - 1)
    b_basic_seq = [b_basic]
    b_real = 0

    for a in range(15):
        b_basic_seq.append(b_basic_seq[a] * 2 + wingbox.stringers_top[0].a / 2)

    if load_factor >= 0:
        b_min = wingbox.t_top * math.sqrt((math.pi * k_c * wingbox.material.E) / (
                abs(normal_stress(y, wingbox, planform, 2 * wingbox.height(planform, y))) * 12 * (
                1 - wingbox.material.nu ** 2)))

    else:
        b_min = wingbox.t_bottom * math.sqrt((math.pi * k_c * wingbox.material.E) / (
                abs(normal_stress(y, wingbox, planform, 0)) * 12 * (1 - wingbox.material.nu ** 2)))

    if b_basic > b_min:
        return -1, False

    if b_min > wingbox.width(planform, y):
        return wingbox.width(planform, y), True

    i = 0
    while b_basic_seq[i] < b_min:
        b_real = b_basic_seq[i]
        i += 1

    # return ((math.pi ** 2) * k_c * wingbox.material.E) / (12 * (1 - wingbox.material.nu ** 2)) * (
    # wingbox.t_top / b) ** 2

    # return b_min, b_real return (2*(b_real - wingbox.width(planform, y))) / (wingbox.stringers_top[0].a + 2 *
    # b_real - 2 * wingbox.width(planform, y)), True
    return round(2 * (b_real + wingbox.width(planform, y) / (wingbox.stringers_top[0].a + 2 * b_real))), True


def web_buckling(y, wingbox, planform, material, stringers_list):
    """
    This function returns the critical stress for web buckling to occur

    :param y: Distance from the root [m]
    :type y: float
    :param planform: The planform used
    :type planform: Planform
    :param wingbox: The wingbox used
    :type wingbox: WingBox
    :return: The critical stress for web buckling to occur [Pa]
    :rtype: float
    """
    '''
    rangrang = np.linspace(0, 15.50, 500)
    tau = []
    for i in range(rangrang):
        tau.append(tau_max(i, wingbox, planform))
    taumax = max(tau)
    print('Tau max is', taumax)
    '''
    k_s = 4.5
    a = wingbox.height(planform, y) * 2  # Not relevant (only short side used in calculations for tau_cr)

    """
    dummy = 1
    for i in range(len(stringers_list) - 1):
        if y < stringers_list[i]:
            if i == 0:
                y = 0
            else:
                y = stringers_list[i - 1]
                dummy = i
                break

    b = stringers_list[dummy] - stringers_list[dummy - 1]
    output = b
    """
    if y < stringers_list[0]:
        b = stringers_list[0]
    elif y > stringers_list[-1]:
        b = planform.b / 2 - stringers_list[-1]
    else:
        for i in range(len(stringers_list)):
            if y <= stringers_list[i] and y >= stringers_list[i - 1]:
                b = stringers_list[i] - stringers_list[i - 1]
    output = b
    '''
    b = np.pi * wingbox.t_spar * np.sqrt(
        (k_s * material.E) / (12 * (1 - material.nu ** 2) * 1.1 * shear_stress(y, wingbox.height(planform, y), wingbox, planform)))
    output = b
    '''
    if a / b < 1:
        a = b
        b = wingbox.height(planform, y) * 2
        output = b

    return ((np.pi ** 2) * k_s * material.E) / (12 * (1 - material.nu ** 2)) * (
            wingbox.t_spar / output) ** 2
