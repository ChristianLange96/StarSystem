import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib.animation
import time
plt.style.use('seaborn-pastel')

# Define the equations
def equations (y0, t):
    theta, theta_dot = y0
    f = [theta_dot, -b*theta_dot - g/L * np.sin(theta)]
    return f

def plot_result(thetas,time):
    y = thetas[:,0]
    fig, ax = plt.subplots()
    line, = ax.plot(time, thetas[:,0], color='k')
    plt.show()
    for n in range(len(time)):
        line.set_data(time[:n], y[:n])
        ax.axis([0, 20, -2, 2])
        fig.canvas.draw()




def animate_result(thetas,time):
    ys = thetas[:, 0]
    fig = plt.figure()
    ax = plt.axes(xlim=(0,20), ylim=(-2,2))
    line, = ax.plot([], [])


    def animate(i):
        thisx = time[0:i]
        thisy = ys[0:i]
        line.set_data(thisx, thisy)
        return line,

    def init():
        line.set_data([], [])
        return line,


    anim = animation.FuncAnimation(fig,animate, init_func=init, frames = 800, interval = 20, blit = True)
    plt.show()

def animate_pendulum(thetas):
    ys = -L*np.cos(thetas[:,0])
    vs = L*(thetas[:,1])
    xs = L*np.sin(thetas[:,0])
    fig = plt.figure()
    ax = plt.axes(xlim = (-2,2), ylim = (-2,2))
    line, = ax.plot([],[],'o-')
    plt.grid(True)
    plt.axis('equal')
    pot_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    kin_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    texts = [pot_text, kin_text, energy_text]

    def animate(i):
        line.set_data([0,xs[i]],[0, ys[i]])
        pot_text.set_text('U = %.5f' % energy(i)[0])
        kin_text.set_text('K = %.5f' % energy(i)[1])
        energy_text.set_text('E = %.5f' % energy(i)[2])
        return line, pot_text, kin_text, energy_text

    def init():
        line.set_data([],[])
        pot_text.set_text('')
        kin_text.set_text('')
        energy_text.set_text('')
        return line, pot_text, kin_text, energy_text

    def energy(i):
        U = m * g * L * ys[i]
        K = 0.5 * m * vs[i]**2
        E = U + K
        return U, K , E


    anim = animation.FuncAnimation(fig,animate, init_func=init, frames = 800, interval = 20, blit = True)
    plt.show()



# Parameters
m = 1
g = 9.82
L = 1
b = 0.2
delta_t = 0.025
t_end = 20
time = np.arange(0,t_end,delta_t)

# Initial Conditions
THETA_0 = np.pi/2
THETA_DOT_0 = 0

# Solver
thetas = odeint(equations, [THETA_0, THETA_DOT_0], time)

# Plot the results
animate_pendulum(thetas)


