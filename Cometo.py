# Thatchapon Unprasert
# 58882220
# Section 1
# Project name : Cometo

from visual import *
from visual.graph import *
from Tkinter import *
import threading
import re
import random
import wx

# Sidenote
# Tkinter is not considered external library, so it's not included here as package

# -constants-
# general constants
G = 6.67e-11
one_au = 1.5e11
one_year = 365.25 * 24 * 3600
dE = 1.496e11
vE = 30e3
CANNOT_DETERMINE = 'fatal error (trajectory behavior beyond approximation)'
NO_OBJECT_FOUND = 'error finding corresponding object'

# scales
planet_true_scale = false
comet_true_scale = false
sun_true_scale = false
default_scale_sun = 10
scale_sun = 1
default_scale_planet = 160
if not sun_true_scale:
    scale_sun = default_scale_sun
scale_planet = 1
if not planet_true_scale:
    scale_planet = default_scale_planet
default_scale_comet = 3e6
scale_comet = 1
if not comet_true_scale:
    scale_comet = default_scale_comet

# objects
earth_mass = 5.973e24
earth_radii = 6371000
sun_mass = 333000 * earth_mass
sun_radii = 6.597e8

# simulation properties
fps = 500
bounds = 6e12
dt = 500
running = true
comet_trail = true
planet_trail = false
create_stat = false
planet_label = false
comet_label = false
auto_clean_comet = false
any_direction_comet = false
max_comet_speed = 1e5
precision_rate = 1.2
rate_rate = 1.2

# scene configuration
scene.width = 1920
scene.height = 660
scene.range = (9e11, 10e11, 10e11)
scene.center = (0, 0, 0)
scene.autoscale = 0
scene.ambient = color.black
scene.title = 'Cometo'

# construct solar system
# assume zero-plane as ecliptic
print 'init sun (' + str(scale_sun) + 'x)'
sun = sphere(pos=(0, 0, 0), color=color.orange, radius=sun_radii * scale_sun, material=materials.emissive)
sun.defaultRadii = sun_radii
sun.comet = false
sun.sun = true
sun.reject = false
sun.mass = sun_mass
sun.name = 'sun'
# let sun fixed object
sun.v = vector(0, 0, 0)
# mimic sun light
scene.lights = [local_light(pos=(0, 0, 0), color=color.white)]
# default center is sun
center_object = sun

# total number of comets
n_comet = 0

# string representation for fps and simulation speed
def fpsString():
    return 'Desire FPS: {:6}'.format(int(fps))

def simulationSpeedString():
    return 'Simulation Speed: {:6.1f}x'.format(float(dt))

# fps label
fps_xoffset = 180
fps_yoffset = 50
fps_label = label(pos=sun.pos, text=fpsString(), xoffset=scene.width * 0.5 - fps_xoffset, yoffset=-scene.height * 0.5 + fps_yoffset, border=3, line=0)
# simulation speed label
ss_xoffset = fps_xoffset + 240
ss_yoffset = fps_yoffset
ss_label = label(pos=sun.pos, text=simulationSpeedString(), xoffset=scene.width * 0.5 - ss_xoffset, yoffset=-scene.height * 0.5 + ss_yoffset,
                 border=3, line=0)

# convert au to meter
def au(meter):
    return meter * one_au

# class for storing graph skeletons
class PreDots():
    title, ytitle, display = "", "", None

    def __init__(self, title, ytitle):
        self.title = title
        self.ytitle = ytitle

    def exitGraph(self, event):
        # destroy its gdisplay
        self.display.display._destroy()
        # remove itself from predot list
        predots.remove(self)

# list for all planets, comets, and comets' stats
planets = []
comets = []
stats = []

# list for graph elements
predots = []

# add a planet to simulation
def addPlanet(name, pos, v, requireColor, color, radii, mass):
    print 'init ' + name + ' (' + str(scale_planet) + 'x)'
    planet = sphere(pos=pos, radius=radii * scale_planet, material=materials.rough)
    planet.mass = mass
    planet.v = v
    planet.comet = false
    planet.sun = false
    planet.defaultRadii = radii
    planet.name = name
    planet.reject = false
    # set planet color if needed
    if requireColor:
        planet.color = color
    # planet orbit
    planet.trail = curve(color=planet.color, visible=planet_trail)
    # planet label
    _label = label(pos=planet.pos, color=planet.color, linecolor=planet.color, text=planet.name, xoffset=0, yoffset=30, line=1, opacity=0.4)
    planet.label = _label
    planet.label.visible = planet_label
    planets.append(planet)
    return planet

# add planets to simulation
mercury = addPlanet('mercury', (0.387 * dE, 0, 0), vector(0, 0, -1.607 * vE), true, color.gray(0.7), 0.383 * earth_radii, 0.0553 * earth_mass)
venus = addPlanet('venus', (0.723 * dE, 0, 0), vector(0, 0, -1.174 * vE), true, (1, 1, 0.5), 0.97 * earth_radii, 0.815 * earth_mass)
earth = addPlanet('earth', (dE, 0, 0), vector(0, 0, -vE), false, None, earth_radii, earth_mass)
earth.material = materials.earth
mars = addPlanet('mars', (1.52 * dE, 0, 0), vector(0, 0, -0.802 * vE), true, color.red, 0.532 * earth_radii, 0.107 * earth_mass)
jupiter = addPlanet('jupiter', (-5.2 * dE, 0, 0), vector(0, 0, 0.434 * vE), false, color.orange, earth_radii * 11.209, 317.83 * earth_mass)
saturn = addPlanet('saturn', (0, 0, -9.58 * dE), vector(-0.323 * vE, 0, 0), true, color.yellow, earth_radii * 9.499, 95.16 * earth_mass)
uranus = addPlanet('uranus', (19.20 * dE, 0, 0), vector(0, 0, -0.228 * vE), true, color.cyan, earth_radii * 4.007, 14.54 * earth_mass)
neptune = addPlanet('neptune', (30.05 * dE, 0, 0), vector(0, 0, -0.182 * vE), true, color.blue, earth_radii * 3.883, 17.15 * earth_mass)

# update perihelion of comet, update as current distance < current perihelion
def updatePerihelion(obj):
    d = distanceFromSun(obj)
    if d < obj.perihelion:
        obj.perihelion = d

# string representation of comet's perihelion
def perihelionString(comet):
    return 'perihelion: ' + '{:.3f}'.format(comet.perihelion / one_au) + ' AU'

# string representation of comet's distance from sun
def distanceString(comet):
    return 'distance: ' + '{:.3f}'.format(distanceFromSun(comet) / one_au) + ' AU'

# string representation of comet's speed
def speedString(comet):
    return 'speed: ' + '{:.3f}'.format(mag(comet.v) / 1000.) + ' km/s'

# string representation of comet
def cometString(comet):
    return comet.name + ' ' + speedString(comet) + ', ' + distanceString(comet) + ', ' + perihelionString(comet)

# get random color
def randomColor():
    rr = random.uniform(0.7, 0.9) + random.uniform(-0.5, 0.1)
    rg = random.uniform(0.7, 0.9) + random.uniform(-0.5, 0.1)
    rb = random.uniform(0.7, 0.9) + random.uniform(-0.5, 0.1)
    return rr, rg, rb

# calculate distance vector between two objects
def distanceVector(obj1, obj2):
    return vector(obj1.pos - obj2.pos)

# calculate distance between two objects
def distance(obj1, obj2):
    return mag(distanceVector(obj1, obj2))

# calculate distance between object and sun
def distanceFromSun(obj):
    return distance(obj, sun)

# calculate kinetic energy of the current object
def kinetic_energy(obj):
    return 0.5 * obj.mass * obj.v.mag2

# calculate potential energy of the current object
def potential_energy(obj):
    energy = 0
    # have potential energy with respect to any other object
    for _planet in planets:
        if obj != _planet:
            energy -= G * obj.mass * _planet.mass / distance(obj, _planet)
    for _comet in comets:
        if obj != _comet:
            energy -= G * obj.mass * _comet.mass / distance(obj, _comet)
    return energy

# calculate total energy of the current object
def totalEnergy(obj):
    return kinetic_energy(obj) + potential_energy(obj)

# add a comet to simulation
def addComet(pos, v, color, radii, mass):
    global n_comet
    _comet = sphere(pos=pos, color=color, radius=radii * scale_comet)
    # comet generic properties
    _comet.name = 'comet' + str(n_comet + 1)
    _comet.reject = false
    _comet.toBeReject = false
    _comet.comet = true
    _comet.sun = false
    _comet.defaultRadii = radii
    _comet.mass = mass
    _comet.v = v
    # comet orbit
    _comet.trail = curve(color=_comet.color, visible=comet_trail, reject=false)
    _comet.stat = None
    _comet.perihelion = distanceFromSun(_comet)
    print 'init ' + _comet.name + ' with mass ' + str(mass) + ', radii ' + str(radii) + ', velocity ' + str(v) + ', origin ' + str(pos)
    _comet.index = size(comets)
    comets.append(_comet)
    # construct stat for each comet, and place right in appropriate scene locations
    x_offset = scene.width * 0.5 - 610
    _n_comet = n_comet
    if n_comet > 18:
        x_offset = -x_offset
        _n_comet -= 19
    # comet stat
    _stat = label(pos=sun.pos, color=_comet.color, text=cometString(_comet), xoffset=x_offset, yoffset=260 - 6 * (1.2 + 4.4 * _n_comet), line=0,
                  opacity=0.4)
    _comet.stat = _stat
    _comet.stat.visible = create_stat
    _label = label(pos=_comet.pos, color=_comet.color, linecolor=_comet.color, text=_comet.name, xoffset=0, yoffset=30, line=1, opacity=0.4)
    _comet.label = _label
    _comet.label.visible = comet_label
    stats.append(_stat)
    n_comet += 1

# check if two objects hit each other (distance <= r1 + r2)
def tooClose(obj1, obj2):
    return mag(obj1.pos - obj2.pos) <= obj1.defaultRadii + obj2.defaultRadii

# calculate gravitational force exerted by other object, then return the acceleration term
def interaction(obj1, obj2):
    vec = obj1.pos - obj2.pos
    accel = -G * obj2.mass * vec.norm() / vec.mag2
    return accel

# update velocity according to external force
def updateInteraction(obj1, obj2):
    da = interaction(obj1, obj2)
    # velocity changes
    obj1.v += da * dt
    # comet hit by something, will be rejected
    if obj1.comet:
        if tooClose(obj1, obj2):
            obj1.toBeReject = true
            obj1.hitSun = obj2.sun
            if obj2.comet:
                obj2.toBeReject = true
    return da

# add predefined comets to simulation
print 'comet scale: ' + str(scale_comet) + 'x'
addComet((-au(6), 1.2e11, 10e8), vector(10e2, -12e2, -22e2), color.orange, 300., 5000.)
addComet((au(2), 10e9, -11e8), vector(5e2, 5.3e3, 6e3), color.green, 200., 150.)
addComet((-1e10, au(4), 3e7), vector(4e3, -10e4, 4e3), color.cyan, 125., 200.)
addComet((au(1.1), 300, 300), vector(-3e3, 5e3, -9e3), color.blue, 350., 600.)
addComet((5.2e8, -au(1.9), -1e7), vector(-12e3, 27e3, -9e3), color.red, 200., 100.)
addComet((-au(5.18), 1e4, 6e5), vector(-10.5e2, -9e3, 3e3), color.magenta, 100., 187.)

# get comet object by index
def comet_object(comet_idx):
    if comet_idx >= len(comets) or comet_idx < 0:
        print 'comet index out of bounds'
        return None
    return comets[comet_idx]

# get planet object by name
def planet_object(target):
    if target == 'sun':
        return sun
    elif target == 'mercury':
        return mercury
    elif target == 'venus':
        return venus
    elif target == 'earth':
        return earth
    elif target == 'mars':
        return mars
    elif target == 'jupiter':
        return jupiter
    elif target == 'saturn':
        return saturn
    elif target == 'uranus':
        return uranus
    elif target == 'neptune':
        return neptune
    return None

m_3tuple_str = '\((\-*\d+\.*\d*),(\-*\d+\.*\d*),(\-*\d+\.*\d*)\)'  # tuple of three doubles
m_3tuple = re.compile(m_3tuple_str)
m_comet = re.compile("comet(\d+)")

# get any object (comet or planet or sun) by name
def what_object(target):
    o = planet_object(target)
    if not o:
        m_comet_check = m_comet.match(target)
        if m_comet_check:
            o = comet_object(int(m_comet_check.group(1)) - 1)
    return o

# deal with all console commands
def process(command):
    if len(command) == 0:
        return
    global running, comets
    size_comets = len(comets)
    if (command == 'pause' and running) or (command == 'resume' and not running):
        running = not running
        return
    # add comet feature
    # example: add mass=100 o=(1,-1,1) v=(2,-3,6) r=156 oau
    m_addComet = re.match('add (.+)', command)
    if m_addComet:
        r, mass = 100., 100.
        v, o = vector(1e3, 1e3, 1e3), (1e11, 1e11, 1e11)
        oau, specify_o = false, false
        add_args = m_addComet.group(1)
        for add_arg in add_args.split(' '):
            sub_add_arg = add_arg.split('=')
            opt = sub_add_arg[0]
            if len(sub_add_arg) == 2:
                _val = sub_add_arg[1]
                if opt == 'mass':
                    mass = float(_val)
                elif opt == 'r':
                    r = float(_val)
                elif opt == 'v':
                    m_v = m_3tuple.match(_val)
                    if m_v:
                        vx, vy, vz = float(m_v.group(1)), float(m_v.group(2)), float(m_v.group(3))
                        v = vector(vx, vy, vz)
                elif opt == 'o':
                    m_o = m_3tuple.match(_val)
                    if m_o:
                        x, y, z = float(m_o.group(1)), float(m_o.group(2)), float(m_o.group(3))
                        o = (x, y, z)
                        specify_o = true
            else:
                if opt == 'oau':
                    oau = true

        if oau and specify_o:
            o = tuple(x * one_au for x in o)
        addComet(o, v, randomColor(), r, mass)
        return

    # delete all outdated comets
    if command == 'clean':
        for _comet in comets:
            if _comet.reject and not _comet.trail.reject:
                print 'delete ' + _comet.name
                _comet.trail.visible = false
                _comet.trail.reject = true
        return

    # adjust precision and fps rate
    m_set_rate = re.match('rate (\w+) (\d+\.*\d*)', command)
    if m_set_rate:
        global precision_rate, rate_rate
        rate_type = m_set_rate.group(1)
        rate_value = float(m_set_rate.group(2))
        if rate_value == 0:
            print 'ignore zero rate'
            return
        if rate_type == 'precision':
            precision_rate = rate_value
        elif rate_type == 'fps' or rate_type == 'rate':
            rate_rate = rate_value
        return

    # follow object
    # example: follow mercury
    m_follow = re.match('(follow|center) (\w+)', command)
    if m_follow:
        global center_object
        target = m_follow.group(2)
        center_object = planet_object(target)
        if not center_object:
            target_comet = what_object(target)
            if target_comet:
                # we won't follow rejected comet
                if target_comet.reject:
                    print 'cannot follow rejected comet'
                    return
                center_object = target_comet
            else:
                print 'unknown target'
        if center_object:
            scene.center = center_object.pos
            # the proper range is still mysterious
            scene.range = (8e10, 7e10, 1e11)
        return

    # velocity multiplication for comets
    # example: speedmult comet2 2.5
    # example: speedmult comet6 (-1,2,3.2)
    m_speed = re.match('speedmult (\w+) (\d+\.*\d*)', command)
    if not m_speed:
        m_speed = re.match('speedmult (\w+) ' + m_3tuple_str, command)
    if m_speed:
        target_comet = what_object(m_speed.group(1))
        if target_comet:
            # we won't modify rejected comet
            if target_comet.reject:
                print "cannot modify rejected comet's speed"
                return
            m_speed_group_size = len(m_speed.groups())
            if m_speed_group_size == 4:
                # every components are being modified individually
                vx, vy, vz = float(m_speed.group(2)), float(m_speed.group(3)), float(m_speed.group(4))
                target_comet.v = vector(target_comet.v.x * vx, target_comet.v.y * vy, target_comet.v.z * vz)
            elif m_speed_group_size == 2:
                # every components are being modified equally
                target_comet.v *= float(m_speed.group(2))
        return

    # set scale for objects
    # example: scale sun 1.2
    m_changeScale = re.match('scale (\w+) (\d+\.*\d*)', command)
    if m_changeScale:
        scale = float(m_changeScale.group(2))
        # zero-scaling is meaningless
        if scale == 0:
            print "scale = 0, ignoring"
            return
        target = what_object(m_changeScale.group(1))
        if target:
            target.radius = target.defaultRadii * scale
        return

    # graphing
    # example: graph d(istance)-t(ime) comet1 comet2
    # example: graph d(istance)-t(ime) comet1 sun
    # example: graph s(peed)-t(ime) comet1
    # example: graph a(ccel)-t(ime) comet1
    # example: graph e(nergy)-t(ime) comet1
    # example: graph analyze comet1
    m_graph_dt_pattern = re.match('graph (d|distance)\-(t|time) (\w+) (\w+)', command)
    if m_graph_dt_pattern:
        object1 = what_object(m_graph_dt_pattern.group(3))
        object2 = what_object(m_graph_dt_pattern.group(4))
        if object1 and object2:
            if object1 == object2:
                print "won't plot graph for the same object"
                return
            if object1.reject or object2.reject:
                print 'error, ' + object1.name + ' or ' + object2.name + ' is rejected'
                return
            # graph skeletons for d-t are being added
            title = 'Distance between ' + object1.name + ' and ' + object2.name + ' (AU) vs Time (year)'
            _predot = PreDots(title, 'Distance (AU)')
            _predot.color, _predot.size, _predot.type = color.black, 4, 'dt'
            _predot.object1, _predot.object2 = object1, object2
            predots.append(_predot)
        else:
            print 'error finding corresponding objects'
        return

    m_graph_energy_pattern = re.match('graph (e|energy)\-(t|time) (\w+)', command)
    if m_graph_energy_pattern:
        _object = what_object(m_graph_energy_pattern.group(3))
        if _object:
            if _object.reject:
                print _object.name + ' is rejected, ignoring'
                return
            if _object.sun:
                print 'exception for sun'
                return
            # graph skeletons for e-t are being added
            title = 'Kinetic Energy (Blue) and Potential Energy (Red) for ' + _object.name + ' (J) vs Time (year)'
            _predot = PreDots(title, 'Energy (J)')
            _predot.kcolor, _predot.pcolor, _predot.size, _predot.type = color.blue, color.red, 4, 'et'
            _predot.object = _object
            predots.append(_predot)
        else:
            print NO_OBJECT_FOUND
        return

    m_graph_analyze_pattern = re.match('graph analyze (\w+)', command)
    if m_graph_analyze_pattern:
        _object = what_object(m_graph_analyze_pattern.group(1))
        if _object:
            if _object.reject:
                print _object.name + ' is rejected, ignoring'
                return
            if _object.sun:
                print 'exception for sun'
                return
            process('graph d-t ' + _object.name + ' sun')
            process('graph v-t ' + _object.name)
            process('graph a-t ' + _object.name)
            process('graph e-t ' + _object.name)
        return

    m_graph_vt_at_pattern = re.match('graph (v|speed|a|accel)\-(t|time) (\w+)', command)
    if m_graph_vt_at_pattern:
        _object = what_object(m_graph_vt_at_pattern.group(3))
        if _object:
            if _object.reject:
                print _object.name + ' is rejected, ignoring'
                return
            if _object.sun:
                print 'exception for sun'
                return
            ytitle, title, mtype = '', '', ''
            graph_type = m_graph_vt_at_pattern.group(1)
            if graph_type == 'v' or graph_type == 'speed':
                ytitle = 'Speed (km/s)'
                title = 'Speed of ' + _object.name + ' (km/s) vs Time (year)'
                mtype = 'vt'
            else:
                ytitle = 'Acceleration (m/s^2)'
                title = 'Acceleration of ' + _object.name + ' (m/s^2) vs Time (year)'
                mtype = 'at'
            _predot = PreDots(title, ytitle)
            _predot.type, _predot.size, _predot.object = mtype, 4, _object
            if _object.name == 'earth':
                _predot.color = color.blue
            else:
                _predot.color = _object.color
            predots.append(_predot)
        else:
            print NO_OBJECT_FOUND
        return

    # approximate comet orbit type
    # example: approx comet3
    m_orbit_approximation = re.match('(approx|approximate) (\w+)', command)
    if m_orbit_approximation:
        _object = what_object(m_orbit_approximation.group(2))
        if _object:
            if not _object.comet:
                print 'support only comet'
                return
            if _object.reject:
                print _object.name + ' is rejected, ignoring'
                return
            Ek = 0.5 * _object.mass * (_object.v.mag ** 2)
            d_sun = distanceFromSun(_object)
            Ep = - G * sun.mass * _object.mass / d_sun
            Ep_sun = -Ep
            # here, we find the massive body with highest gravitational force
            # since it would be the origin of this comet
            Ep_max = -Ep / d_sun
            Ep_max_host = sun
            for _planet in planets:
                # gather potential energy by every planets
                r = distance(_object, _planet)
                d = G * _planet.mass * _object.mass / r
                if d / r > Ep_max:
                    Ep_max = d / r
                    Ep_max_host = _planet
                Ep -= d
            for _comet in comets:
                if _comet.reject or _comet == _object:
                    continue
                # gather potential energy by other comets
                r = distance(_object, _comet)
                d = G * _comet.mass * _object.mass / r
                if d / r > Ep_max:
                    Ep_max = d / r
                    Ep_max_host = _comet
                Ep -= d
            # total energy = kinetic energy + potential energy
            Et = Ek + Ep
            print 'approx. orbit type:',
            if Et < 0:
                print 'elliptic'
            elif Et < 1e-25:
                print 'parabolic'
            elif Et > 0:
                print 'hyperbolic'
            print 'approx. total energy: ' + str(Et) + ' J'
            # angular momentum
            M = Ep_max_host.mass
            L = _object.mass * mag(distanceVector(_object, Ep_max_host).cross(_object.v))
            if Ep_max_host != sun:
                Et += Ep_sun
            print 'approx. center: ' + Ep_max_host.name
            numerator_e = 2 * Et * ((M + _object.mass) ** 3) * L * L
            denominator_e = G * G * (M ** 5) * (_object.mass ** 3)
            # eccentricity
            _val_ = 1. + numerator_e / denominator_e
            if _val_ < 0:
                print CANNOT_DETERMINE
                return
            e = sqrt(_val_)
            print 'approx. orbital eccentricity: ' + str(e)
            # semi-major axis
            a = - G * M * _object.mass / (2. * Et)
            print 'approx. semi-major axis: ' + str(a / one_au) + ' AU'
            # perihelion
            rp = a * (1 - e)
            print 'approx. perihelion: ' + str(rp / one_au) + ' AU'
            if Et < 0:
                # aphelion for elliptic
                ra = a * (1 + e)
                print 'approx. aphelion: ' + str(ra / one_au) + ' AU'
            # max speed
            vp = sqrt(G * M * (1 + e) / a / (1 - e))
            print 'approx. maximum speed: ' + str(vp / 1000) + ' km/s'
            if Et < 0:
                # min speed for elliptic
                va = sqrt(G * M * (1 - e) / a / (1 + e))
                print 'approx. minimum speed: ' + str(va / 1000) + ' km/s'
                # orbital period for elliptic
                T = 2 * pi * sqrt(pow(a, 3) / G / M)
                print 'approx. orbital period: ' + str(T / one_year) + ' year'
            elif Et > 0:
                # terminal speed for hyperbolic
                vinf = sqrt(G * M / -a)
                print 'approx. terminal speed: ' + str(vinf / 1000) + ' km/s'
        return

    # print all adjustable values
    if command == 'allvars':
        print 'default_scaled_sun: ' + str(default_scale_sun)
        print 'scale_sun: ' + str(scale_sun)
        print 'default_scaled_planet: ' + str(default_scale_planet)
        print 'scale_planet: ' + str(scale_planet)
        print 'planet_orbit: ' + str(planet_trail)
        print 'planet_label: ' + str(planet_label)
        print 'default_scale_comet: ' + str(default_scale_comet)
        print 'scale_comet: ' + str(scale_comet)
        print 'comet_orbit: ' + str(comet_trail)
        print 'comet_label: ' + str(comet_label)
        print 'auto_clean_comet: ' + str(auto_clean_comet)
        print 'any_direction_comet: ' + str(any_direction_comet)
        print 'create_comet_stat: ' + str(create_stat)
        print 'max_comet_speed: ' + str(max_comet_speed)
        print 'desire_fps: ' + str(fps)
        print 'fps_rate: ' + str(rate_rate)
        print 'simulation_speed ' + str(dt)
        print 'precision_rate: ' + str(precision_rate)
        print 'simulation_bounds: ' + str(bounds)
        print 'center_object:',
        if center_object:
            print center_object.name
        else:
            print 'null'
        return

    if command != 'pause' and command != 'resume':
        print 'unrecognized command: ' + command

# control window
class ControlWindow(threading.Thread):
    def accept(self, event):
        user_input = self.input.get()
        user_input_len = len(user_input)
        if user_input_len > 0:
            # we separate input by |
            for sub_input in user_input.split('|'):
                sub_input = sub_input.strip()
                print "-> " + sub_input
                process(sub_input)
        self.input.delete(0, user_input_len)

    def __init__(self):
        threading.Thread.__init__(self)
        self.root = Tk()
        self.root.geometry('800x30+0+' + str(int(scene.height)))
        self.root.title('Console')
        self.label = Label(self.root, text='Enter command')
        self.label.pack(side=LEFT)
        self.input = Entry(self.root, width=80)
        self.input.bind('<Return>', self.accept)
        self.input.pack(side=LEFT)
        self.acceptBtn = Button(self.root, text='Run', command=self.accept)
        self.acceptBtn.pack(side=RIGHT)

    def run(self):
        self.root.mainloop()

def adjustRate(type):
    global fps, rate_rate
    if type == 1:
        fps *= rate_rate
    elif type == -1:
        if fps <= rate_rate:
            return
        fps /= rate_rate
    if fps < 1:
        fps = 1
    fps_label.text = fpsString()

def adjustPrecision(type):
    global dt, precision_rate
    if type == 1:
        if dt <= precision_rate:
            return
        dt /= precision_rate
    elif type == -1:
        dt *= precision_rate
    if dt < 1:
        dt = 1
    ss_label.text = simulationSpeedString()

def simulationState(isrun):
    if isrun:
        return 'pause'
    return 'resume'

def showOrbitState(show, object):
    if show:
        return 'hide ' + object + ' orbit'
    return 'show ' + object + ' orbit'

class SimulationWindow(threading.Thread):
    def increaseRate(self):
        adjustRate(1)

    def decreaseRate(self):
        adjustRate(-1)

    def increasePrecision(self):
        adjustPrecision(1)

    def decreasePrecision(self):
        adjustPrecision(-1)

    def toggleCometScale(self):
        global scale_comet, default_scale_comet
        self.comet_scale.set(not self.comet_scale.get())
        if self.comet_scale.get():
            scale_comet = 1
        else:
            scale_comet = default_scale_comet
        for _comet in comets:
            # skip rejected comet
            if _comet.reject:
                continue
            _comet.radius = _comet.defaultRadii * scale_comet

    def toggleCometOrbit(self):
        global comet_trail
        self.comet_trail.set(not self.comet_trail.get())
        comet_trail = self.comet_trail.get()
        for _comet in comets:
            # skip comet with rejected orbit
            if _comet.trail.reject:
                continue
            _comet.trail.visible = comet_trail

    def toggleCometLabel(self):
        global comet_label
        self.comet_label.set(not self.comet_label.get())
        comet_label = self.comet_label.get()
        for _comet in comets:
            # skip rejected comet in both ways
            if _comet.reject or _comet.trail.reject:
                continue
            _comet.label.visible = comet_label

    def togglePlanetOrbit(self):
        global planet_trail
        self.planet_trail.set(not self.planet_trail.get())
        planet_trail = self.planet_trail.get()
        for _planet in planets:
            _planet.trail.visible = planet_trail

    def togglePlanetScale(self):
        global scale_planet, default_scale_planet
        self.planet_scale.set(not self.planet_scale.get())
        if self.planet_scale.get():
            scale_planet = 1
        else:
            scale_planet = default_scale_planet
        for _planet in planets:
            _planet.radius = _planet.defaultRadii * scale_planet

    def togglePlanetLabel(self):
        global planet_label
        self.planet_label.set(not self.planet_label.get())
        planet_label = self.planet_label.get()
        for _planet in planets:
            _planet.label.visible = planet_label

    def toggleStat(self):
        global create_stat
        self.create_stat.set(not self.create_stat.get())
        create_stat = self.create_stat.get()
        for _stat in stats:
            _stat.visible = create_stat
            _stat.pos = center_object.pos

    def toggleRedCyan(self):
        self.redcyan.set(not self.redcyan.get())
        if self.redcyan.get():
            scene.stereo = 'redcyan'
        else:
            scene.stereo = 'nostereo'

    def toggleSunScale(self):
        global scale_sun, default_scale_sun, sun
        self.sun_scale.set(not self.sun_scale.get())
        if self.sun_scale.get():
            scale_sun = 1
        else:
            scale_sun = default_scale_sun
        sun.radius = sun.defaultRadii * scale_sun

    def toggleSimulation(self):
        global running
        running = not running
        self.toggleSimulationBtn.config(text=simulationState(running))

    def cleanUp(self):
        process('clean')

    def autoCleanup(self):
        global auto_clean_comet
        self.auto_clean.set(not self.auto_clean.get())
        auto_clean_comet = self.auto_clean.get()
        process('clean')

    def anyDirection(self):
        global any_direction_comet
        self.any_direction_comet.set(not self.any_direction_comet.get())
        any_direction_comet = self.any_direction_comet.get()

    def randomComet(self):
        # random mass and radii
        mass, r = random.uniform(10, 100000), random.uniform(200, 900)
        # random origin
        ox, oy, oz = random.uniform(-6, 6), random.uniform(-6, 6), random.uniform(-6, 6)
        # random velocity
        vx, vy, vz, ov = random.uniform(-5, 5), random.uniform(-5, 5), random.uniform(-5, 5), random.uniform(2, 4.5)
        if not any_direction_comet:
            # constrain velocity so object moves towards sun
            pre_o = (ox, oy, oz)
            pre_v_mag = vector(vx * (10 ** ov), vy * (10 ** ov), vz * (10 ** ov)).mag
            to_sun = sun.pos - pre_o
            proper_v = vector(to_sun.x + random.uniform(-1, 1), to_sun.y + random.uniform(-1, 1), to_sun.z + random.uniform(-1, 1))
            proper_v *= pre_v_mag
            if proper_v.mag > max_comet_speed:
                proper_v *= max_comet_speed * random.uniform(0.96, 1) / proper_v.mag
            vx, vy, vz = proper_v.x, proper_v.y, proper_v.z
        else:
            vx, vy, vz = vx * (10 ** ov), vy * (10 ** ov), vz * (10 ** ov)
            # workaround: too low velocity magnitude
            if vx < 2000:
                vx = 2000
            if vy < 2000:
                vy = 2000
            if vz < 2000:
                vz = 2000
        process('add mass={:f} o=({:f},{:f},{:f}) v=({:f},{:f},{:f}) r={:f} oau'.format(mass, ox, oy, oz, vx, vy, vz, r))

    def __init__(self):
        threading.Thread.__init__(self)
        global planet_trail, comet_trail, create_stat, scale_planet, scale_comet, auto_clean_comet
        self.planet_trail = BooleanVar()
        self.comet_trail = BooleanVar()
        self.create_stat = BooleanVar()
        self.planet_scale = BooleanVar()
        self.comet_scale = BooleanVar()
        self.sun_scale = BooleanVar()
        self.planet_label = BooleanVar()
        self.comet_label = BooleanVar()
        self.redcyan = BooleanVar()
        self.auto_clean = BooleanVar()
        self.any_direction_comet = BooleanVar()
        self.root = Tk()
        self.root.geometry('460x300+' + str(int(scene.width) - 460 - 50) + '+' + str(int(scene.height)))
        self.root.title('Simulation control')
        simulation_frame = LabelFrame(self.root, text='Simulation rate')
        simulation_frame.pack()
        self.increaseRateBtn = Button(simulation_frame, text='rate +', command=self.increaseRate)
        self.increaseRateBtn.pack(side=LEFT)
        self.decreaseRateBtn = Button(simulation_frame, text='rate -', command=self.decreaseRate)
        self.decreaseRateBtn.pack(side=LEFT)
        self.increasePrecisionBtn = Button(simulation_frame, text='precision +', command=self.increasePrecision)
        self.increasePrecisionBtn.pack(side=LEFT)
        self.decreasePrecisionBtn = Button(simulation_frame, text='precision -', command=self.decreasePrecision)
        self.decreasePrecisionBtn.pack(side=LEFT)
        self.toggleSimulationBtn = Button(simulation_frame, text=simulationState(running), command=self.toggleSimulation)
        self.toggleSimulationBtn.pack(side=LEFT)

        comet_frame = LabelFrame(self.root, text='Comet properties')
        comet_frame.pack()
        self.toggleCometOrbitBtn = Checkbutton(comet_frame, text='orbit', variable=self.comet_trail, command=self.toggleCometOrbit)
        if comet_trail:
            self.comet_trail.set(true)
            self.toggleCometOrbitBtn.select()
        self.toggleCometOrbitBtn.pack(side=LEFT)
        self.toggleStatBtn = Checkbutton(comet_frame, text='stat', variable=self.create_stat, command=self.toggleStat)
        if create_stat:
            self.create_stat.set(true)
            self.toggleStatBtn.select()
        self.toggleStatBtn.pack(side=LEFT)
        self.toggleCometScaleBtn = Checkbutton(comet_frame, text='true scale', variable=self.comet_scale, command=self.toggleCometScale)
        if scale_comet == 1:
            self.comet_scale.set(true)
            self.toggleCometScaleBtn.select()
        self.toggleCometScaleBtn.pack(side=LEFT)
        self.toggleCometLabelBtn = Checkbutton(comet_frame, text='label', variable=self.comet_label, command=self.toggleCometLabel)
        if comet_label:
            self.comet_label.set(true)
            self.toggleCometLabelBtn.select()
        self.toggleCometLabelBtn.pack(side=LEFT)

        comet_other_frame = LabelFrame(self.root, text='Comet control')
        comet_other_frame.pack()
        self.randomCometBtn = Button(comet_other_frame, text='randomize', command=self.randomComet)
        self.randomCometBtn.pack(side=LEFT)
        self.anyDirectionBtn = Checkbutton(comet_other_frame, text='any direction', variable=self.any_direction_comet, command=self.anyDirection)
        if any_direction_comet:
            self.any_direction_comet.set(true)
            self.anyDirectionBtn.select()
        self.anyDirectionBtn.pack(side=LEFT)
        self.cleanUpBtn = Button(comet_other_frame, text='clean up', command=self.cleanUp)
        self.cleanUpBtn.pack(side=LEFT)
        self.autoCleanupBtn = Checkbutton(comet_other_frame, text='auto clean', variable=self.auto_clean, command=self.autoCleanup)
        if auto_clean_comet:
            self.auto_clean.set(true)
            self.autoCleanup.select()
        self.autoCleanupBtn.pack(side=LEFT)

        planet_frame = LabelFrame(self.root, text='Planet properties')
        planet_frame.pack()
        self.togglePlanetOrbitBtn = Checkbutton(planet_frame, text='orbit', variable=self.planet_trail, command=self.togglePlanetOrbit)
        if planet_trail:
            self.planet_trail.set(true)
            self.togglePlanetOrbitBtn.select()
        self.togglePlanetOrbitBtn.pack(side=LEFT)
        self.togglePlanetScaleBtn = Checkbutton(planet_frame, text='true scale', variable=self.planet_scale, command=self.togglePlanetScale)
        if scale_planet == 1:
            self.planet_scale.set(true)
            self.togglePlanetScaleBtn.select()
        self.togglePlanetScaleBtn.pack(side=LEFT)
        self.togglePlanetLabelBtn = Checkbutton(planet_frame, text='label', variable=self.planet_label, command=self.togglePlanetLabel)
        if planet_label:
            self.planet_label.set(true)
            self.togglePlanetLabelBtn.select()
        self.togglePlanetLabelBtn.pack(side=LEFT)

        scene_frame = LabelFrame(self.root, text='Miscellaneous')
        scene_frame.pack()
        self.redCyanBtn = Checkbutton(scene_frame, text='red-cyan 3D', variable=self.redcyan, command=self.toggleRedCyan)
        if scene.stereo == 'redcyan':
            self.redcyan.set(true)
            self.redCyanBtn.select()
        self.redCyanBtn.pack(side=LEFT)
        self.toggleSunScaleBtn = Checkbutton(scene_frame, text='true sun scale', variable=self.sun_scale, command=self.toggleSunScale)
        if scale_sun == 1:
            self.sun_scale.set(true)
            self.toggleSunScaleBtn.select()
        self.toggleSunScaleBtn.pack(side=LEFT)

    def run(self):
        self.root.mainloop()

# might be the messy way for multi-threading, at least it works
# open "console" window
def run_control_window():
    cwindow = ControlWindow()
    # cwindow.start()

# open speed window
def run_simulation_window():
    swindow = SimulationWindow()
    swindow.start()

run_control_window()
run_simulation_window()

# simulation
t = 0
while true:
    rate(fps)
    if running:
        # update interaction between planets
        for planet1 in planets:
            # acceleration due to sun
            planet1.a = updateInteraction(planet1, sun)
            for planet2 in planets:
                if planet1 != planet2:
                    # acceleration due to other planets
                    planet1.a += updateInteraction(planet1, planet2)
            # update positions of planets
            # append trail of planets
            planet1.pos += planet1.v * dt
            planet1.label.pos = planet1.pos
            planet1.trail.append(pos=planet1.pos)
        # update interaction between comets and their surroundings
        # update interaction between planets and sun
        for comet1 in comets:
            if comet1.reject:
                continue
            # acceleration due to sun
            comet1.a = updateInteraction(comet1, sun)
            for planet in planets:
                # acceleration due to planets
                comet1.a += updateInteraction(comet1, planet)
            for comet2 in comets:
                if comet2.reject:
                    continue
                if comet1 != comet2:  # exert itself is meaningless
                    # acceleration due to other comets
                    comet1.a += updateInteraction(comet1, comet2)
                    if create_stat:
                        updatePerihelion(comet1)
                        updatePerihelion(comet2)
                        comet1.stat.text = cometString(comet1)
                        comet2.stat.text = cometString(comet2)
            # update positions of comets
            # append trails of comets
            comet1.pos += comet1.v * dt
            if comet1.trail.reject:
                continue
            comet1.label.pos = comet1.pos
            comet1.trail.append(pos=comet1.pos)
            # remove any comet if any of its following conditions met
            # out of bounds situation:
            # exclude comet from calculation to prevent disappear bug of VPython
            # this is also useful because we would not care any distant thing
            out_of_bounds = distanceFromSun(comet1) > bounds
            # collision situation (toBeReject=true):
            # comet either hit sun, planet, or another comet; this comet will be permanently rejected
            if out_of_bounds or comet1.toBeReject:
                comet1.stat.text += ' (X)'
                print 'disable ' + comet1.name + ' simulation',
                if out_of_bounds:
                    print '(bounds limit)'
                    if center_object == comet1:
                        # reset center object to sun
                        center_object = sun
                        print 'reset center object to sun'
                if comet1.toBeReject:
                    if comet1.hitSun:
                        print '(within sun)'
                    else:
                        print '(collision)'
                # hide comet label
                comet1.label.visible = false
                # hide comet
                comet1.visible = false
                # comet is rejected
                comet1.reject = true
                # comet is already rejected
                comet1.toBeReject = false
                if auto_clean_comet:
                    print 'delete ' + comet1.name
                    # hide comet trajectory
                    comet1.trail.visible = false
                    # trajectory is rejected
                    comet1.trail.reject = true
        # plot data of every graph defined in graphs list
        for predot in predots:
            # construct corresponding graph display if necessary
            if not predot.display:
                # first create window for gdisplay
                _window = window(width=900, height=300, x=len(predots) * 20, y=720, title=predot.title)
                # then assign gdisplay for each predot
                predot.display = gdisplay(window=_window, width=900, height=300, background=color.white, foreground=color.black,
                                          xtitle='t(yr)', ytitle=predot.ytitle)
                try:
                    # override calling method when "X" button is clicked
                    predot.display.display.window.win.Bind(wx.EVT_CLOSE, predot.exitGraph)
                except:
                    pass
                if predot.type == 'dt':
                    gdot = gdots(display=predot.display, size=predot.size)
                    gdot.color = predot.color
                    gdot.object1, gdot.object2 = predot.object1, predot.object2
                    predot.gdot = gdot
                elif predot.type == 'et':
                    gdot1 = gdots(display=predot.display, size=predot.size)
                    gdot2 = gdots(display=predot.display, size=predot.size)
                    gdot1.color, gdot2.color = predot.kcolor, predot.pcolor
                    predot.gdot1 = gdot1
                    predot.gdot2 = gdot2
                else:
                    gdot = gdots(display=predot.display, size=predot.size)
                    gdot.color = predot.color
                    gdot.object = predot.object
                    predot.gdot = gdot
            # check graph type, then plot
            if predot.type == 'dt':
                predot.gdot.plot(pos=(t / one_year, distance(predot.object1, predot.object2) / one_au))
            elif predot.type == 'vt':
                predot.gdot.plot(pos=(t / one_year, predot.object.v.mag / 1000.))
            elif predot.type == 'at':
                predot.gdot.plot(pos=(t / one_year, predot.object.a.mag))
            elif predot.type == 'et':
                predot.gdot1.plot(pos=(t / one_year, kinetic_energy(predot.object) / 1e6))
                predot.gdot2.plot(pos=(t / one_year, potential_energy(predot.object)))
        t += dt

    # have some key shortcuts for simulation
    if scene.kb.keys:
        k = scene.kb.getkey()
        if k == 'p':  # run or not run simulation
            running = not running
        elif k == 's':  # toggle comets' stat
            create_stat = not create_stat
            for stat in stats:
                stat.visible = create_stat
        elif k == 'e':  # exit entirely
            exit()

    # scene always has center at specific target object
    if center_object:
        scene.center = center_object.pos
        # also update fps, simulation speed label, and comets stat positions
        fps_label.pos = ss_label.pos = scene.center
        for stat in stats:
            stat.pos = scene.center
