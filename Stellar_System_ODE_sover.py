import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib.animation as animation
plt.style.use('seaborn-pastel')
import random

# Define parameters
G = 1
delta_t = 0.025
t_end = 1000
time = np.arange(0,t_end,delta_t)

# Planet/Star class
class Planet:

    def __init__(self, m, x0, y0, vx0, vy0, name, isStar):
        self.m = m
        self.x0 = x0
        self.y0 = y0
        self.vx0 = vx0
        self.vy0 = vy0
        self.params = [x0, y0, vx0, vy0]
        self.isStar = isStar
        self.name = name

# Equations
def equations(z0,t, planets):
    f = []
    planid = 0
    thisx = 0
    thisy = 0
    for i in range(0,4*len(planets)):
        if i in np.arange(0,len(planets) * 4,4):                                    # If i is x
            f.append(z0[i+2])
            thisx = z0[i]
        elif i in np.arange(1,len(planets) * 4, 4):                                 # If i is y
            f.append(z0[i+2])
            thisy = z0[i]
        elif i in np.arange(2,len(planets) * 4,4):                                  # If i is vx
            dvxdt = 0
            for k in range(len(planets)):
                if planid != k:                                         # Calculating force for each other planet
                    otherx = z0[k*4]
                    othery = z0[k*4+1]
                    rvec = [thisx - otherx, thisy - othery]
                    r = np.linalg.norm(rvec)
                    rhat = rvec / r
                    acc = -G * planets[k].m / pow(r,2) * rhat[0]
                    dvxdt += acc
            f.append(dvxdt)

        elif i in np.arange(3,len(planets) * 4,4):                                  # If i is vy
            dvydt = 0
            for k in range(len(planets)):                               # Calculating for on each other planet
                if planid != k:
                    otherx = z0[k * 4]
                    othery = z0[k * 4 + 1]
                    rvec = [thisx - otherx, thisy - othery]
                    r = np.linalg.norm(rvec)
                    rhat = rvec / r
                    acc = -G * planets[k].m / pow(r,2) * rhat[1]
                    dvydt += acc
            f.append(dvydt)
            planid += 1
    return f

# Creating planets
star = Planet(1000, 0, 0, 0, 0, '4DMadsen', True)
planet1 = Planet(3, 75, 0, 0, -2.5, 'Holleupiter', False)
planet2 = Planet(2, 100, 0, 0, 2.3, 'Gormgantua', False)
planet3 = Planet(10, 150, 40, -2, 1.3, 'FUMIO', False)
planet4 = Planet(5, 0, 80, -2.5, 0, 'Looprevil96', False)

# List of planets
planets = [planet1, planet2, planet3,planet4, star]

# Getting all start-values
p0 = []
for planet in planets:
    p0 += planet.params

# Solving equations
res = odeint(equations, p0, time, rtol = 10 ** (-12), args = (planets,))


# Animation Function
def animate_stellar_system(res,planets):

    # Creating figure and setting up limits
    fig = plt.figure()
    ax = plt.axes(xlim=(-100, 150), ylim=(-150, 100))


    # Possible colors
    colors = ['g','b','r','c','orange', 'purple', 'brown', 'lightsalmon', 'navy']
    def random_color():
        return random.choice(colors)


    # Lists / Sizes to expand for each planet
    m_tot = 0
    circles = []
    lines = []
    names = []

    # Setting up all objects assiciated with each planet
    for i in range(0,len(planets)):
        if planets[i].isStar:
            col = 'yellow'
            size = 15
        else:
            col = random_color()
            size = planets[i].m * 0.4
            colors.remove(col)
        circle = plt.Circle((i,0), size, fc = col)
        circles.append(circle)
        ax.add_artist(circle)
        line, = ax.plot([],[], color = col, lw = 0.3)
        lines.append(line)
        name = ax.annotate(planets[i].name, xy = [0,0])
        name.set_animated('True')
        names.append(name)
        m_tot += planets[i].m

    # Creating barh-plot
    y = [138,123,108]
    width = [1, 1, 1]
    height = 10
    rects = plt.barh(y,width, height, left= -50)

    # Creating texts needed
    pot_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    kin_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    COM_text = ax.text(0.02, 0.80, '', transform = ax.transAxes)
    energies = [pot_text, kin_text, energy_text, COM_text]

    # Lists for lines
    rows = len(res)
    cols = len(planets)
    foosx = [[0 for z in range(rows)] for j in range(cols)]
    foosy = [[0 for z in range(rows)] for j in range(cols)]
    plt.grid(True)
    plt.axis('equal')

    # Init-setup
    def init():
        for rect, yi in zip(rects,energy(0)):
            rect.set_width(0)
        for i in range(0,len(circles)):
            circles[i].center = (0,0)
            lines[i].set_data([],[])

        axes = plt.gca()
        axes.set_ylim([-200, 200])
        axes.set_xlim([-200, 200])
        COM_text.set_text('')
        pot_text.set_text('')
        kin_text.set_text('')
        energy_text.set_text('')
        return list(rects) + lines + circles + names + energies

    # Each step in the animation
    def animate(i):
        for rect, yi in zip(rects,energy(i)):
            rect.set_width(yi * 0.1)

        for k in range(0,len(planets)):
            x = res[i][k * 4]
            y = res[i][k * 4 + 1]
            circles[k].center = (x,y)
            foosx[k][i] = x
            foosy[k][i] = y
            lines[k].set_data(foosx[k][:i], foosy[k][:i])
            names[k].set_position((x + 3,y + 3 ))
        pot_text.set_text('U = %.5f' % energy(i)[0])
        kin_text.set_text('K = %.5f' % energy(i)[1])
        energy_text.set_text('E = %.5f' % energy(i)[2])
        (comx, comy) = center_of_mass(i)
        COM_text.set_text("COM = ( %.3f" % comx + ", %.3f" % comy + ")")
        return list(rects) + lines + circles + names  + energies

    # Calculating the potential, kinetic and total mechanic energy
    def energy(i):
        K = 0
        U = 0
        for k in range(0,len(planets)):
            K += 0.5 * planets[k].m * (res[i][k * 4 + 2] **2 + res[i] [k * 4 + 3] **2 )
            thisx = res[i][k * 4]
            thisy = res[i][k * 4 + 1]
            for z in range(0, len(planets)):
                if k != z:
                    otherx = res[i][z * 4]
                    othery = res[i][z * 4 + 1]
                    rvec = [thisx - otherx, thisy - othery]
                    r = np.linalg.norm(rvec)
                    Uz = -G * planets[k].m * planets[z].m / r
                    U += Uz
        U  = U * 0.5     # So that each planet isn't counted twice
        return U, K, U+K

    # Finding the center of mass for the entire system
    def center_of_mass(i):
        COM_X = 0
        COM_Y = 0
        for k in range(0,len(planets)):
            COM_X += res[i][k * 4] / m_tot
            COM_Y += res[i][k * 4 + 1] / m_tot
        return (COM_X, COM_Y)

    # Animating the function
    anim = animation.FuncAnimation(fig, animate, init_func = init, frames = 12000 , interval = 5, blit = True)
    #plt.show()
    anim.save('stellar_sytem.mp4', writer= 'ffmpeg'); print('mp4 created')

# Call of function
animate_stellar_system(res,planets)
