import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib.animation
import time
plt.style.use('seaborn-pastel')


# Define parameters
g = 9.82
delta_t = 0.025
t_end = 200
time = np.arange(0,t_end, delta_t)


# Double Pendulum class

class DoublePendulum:

    def __init__(self,m1,l1,th1_0,w1_0,m2,l2,th2_0,w2_0):
        self.m1 = m1
        self.l1 = l1
        self.th1 = th1_0
        self.w1 = w1_0
        self.m2 = m2
        self.l2 = l2
        self.th2 = th2_0
        self.w2 = w2_0
        self.origin = (0,0)
        self.params = np.array([th1_0,w1_0,th2_0,w2_0])

# Initial Conditions

pendulum  = DoublePendulum(m1 = 1000, l1 = 2, th1_0 = 0 + 0.0001, w1_0 =0, m2 = 1, l2 = 2, th2_0 = np.pi, w2_0 = 0)
pendulum2 = DoublePendulum(m1 = 1, l1 = 1, th1_0 = np.pi/3 - 0.0001, w1_0 =0, m2 = 1, l2 = 2, th2_0 = np.pi, w2_0 = 0)
pendulum3 = DoublePendulum(m1 = 1, l1 = 1, th1_0 = np.pi/3, w1_0 =0, m2 = 1, l2 = 2, th2_0 = np.pi, w2_0 = 0)
pendulum4 = DoublePendulum(m1 = 1, l1 = 1, th1_0 = np.pi/3, w1_0 =0, m2 = 1, l2 = 2, th2_0 = np.pi + 0.001, w2_0 = 0)
pendulum5 = DoublePendulum(m1 = 1, l1 = 1, th1_0 = np.pi/3, w1_0 =0, m2 = 1, l2 = 2, th2_0 = np.pi - 0.001, w2_0 = 0)

pendulums = [pendulum, pendulum2, pendulum3, pendulum4, pendulum5]


# Plot the results

# Equations

def equations(y0,t, pendulum):
    th1, w1, th2, w2 = y0
    m1 = pendulum.m1
    m2 = pendulum.m2
    l1 = pendulum.l1
    l2 = pendulum.l2

    num1 = -m2 * l1 * w1**2 * np.sin(th1-th2) * np.cos(th1-th2) + g * m2 * np.sin(th2)* np.cos(th1-th2) - m2 * l2 * np.sin(th1-th2) * w2**2 - (m1+m2) * g * np.sin(th1)
    den1 = l1 * (m1 + m2) - m2 * l1 * np.cos(th1-th2)**2

    num2 = m2 * l2 * w2**2 * np.sin(th1-th2)* np.cos(th1-th2) + g * np.sin(th1) * np.cos(th1-th2)*(m1+m2) + l1 * w1**2 * np.sin(th1-th2) * (m1+m2) - g * np.sin(th2) * (m1+m2)
    den2 = l2 * (m1+m2) - m2 * l2 * np.cos(th1-th2)**2
    f = [w1, num1/den1, w2, num2/den2]
    return f

# Solver

res1 = odeint(equations,pendulum.params, time, args = (pendulum,) , rtol= 10** (-12))
res2 = odeint(equations,pendulum2.params, time, rtol= 10** (-12), args = (pendulum2,))
res3 = odeint(equations,pendulum3.params, time, rtol= 10** (-12), args = (pendulum3,))
res4 = odeint(equations,pendulum4.params, time, rtol= 10** (-12), args = (pendulum4,))
res5 = odeint(equations,pendulum5.params, time, rtol= 10** (-12), args = (pendulum5,))

res = [res1, res2, res3, res4, res5]

# Functions

def animate_pendulum(res):
    th1 = res[:, 0]
    w1 = res[:, 1]
    th2 = res[:, 2]
    w2 = res[:, 3]

    m1 = pendulum.m1
    m2 = pendulum.m2
    l1 = pendulum.l1
    l2 = pendulum.l2

    x1 = l1 * np.sin(th1)
    x2 = x1 + l2 * np.sin(th2)
    y1 = -l1 * np.cos(th1)
    y2 = y1 - l2 * np.cos(th2)



    fig = plt.figure()
    ax = plt.axes(xlim=(-5, 5), ylim=(-5, 5))
    line1, = ax.plot([],[],'o-')
    line2, = ax.plot([], [], 'o-')
    plt.grid(True)
    plt.axis('equal')
    pot_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    kin_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    texts = [pot_text, kin_text, energy_text]


    def init():
        line1.set_data([],[])
        line2.set_data([], [])
        pot_text.set_text('')
        kin_text.set_text('')
        energy_text.set_text('')
        return line1, line2, pot_text, kin_text, energy_text

    def energy(i):

        U = m1 * g * y1[i] + m2 * g * y2[i]
        K = 0.5 * m1 * (l1*w1[i]) ** 2 + 0.5 * m2 *( (l1*w1[i])**2 + (l2*w2[i])**2 + 2*l1*l2*w1[i]*w2[i]*np.cos(th1[i]-th2[i]))
        E = U + K
        return U, K, E

    def animate(i):
        line1.set_data([0, x1[i]], [0, y1[i]])
        line2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        pot_text.set_text('U = %.4f' % energy(i)[0])
        kin_text.set_text('T = %.4f' % energy(i)[1])
        energy_text.set_text('E = %.4f' % energy(i)[2])
        return line1, line2, pot_text, kin_text, energy_text

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=10000, interval=20, blit=True)
    plt.show()

def animiate_pendulums(res,pendulums):
    x1s = []
    y1s = []
    w1s = []
    x2s = []
    y2s = []
    w2s = []
    m1s = []
    m2s = []
    l1s = []
    l2s = []
    line1s = []
    line2s = []
    fig = plt.figure()
    ax = plt.axes(xlim = (-5,5), ylim = (-5,5))

    for i in range(0,len(pendulums)):
        pend = pendulums[i]
        result = res[i]
        x1 = pend.l1 * np.sin(result[:,0])
        x2 = x1 + pend.l2 * np.sin(result[:,2])
        y1 = -pend.l1 * np.cos(result[:,0])
        y2 = y1 - pend.l2 * np.cos(result[:,2])

        line1, = ax.plot([],[],'o-')
        line2, = ax.plot([],[],'o-')

        x1s.append(x1)
        y1s.append(y1)
        w1s.append(result[:,1])
        x2s.append(x2)
        y2s.append(y2)
        w2s.append(result[:,3])
        m1s.append(pend.m1)
        m2s.append(pend.m2)
        l1s.append(pend.l1)
        l2s.append(pend.l2)
        line1s.append(line1)
        line2s.append(line2)
    lines = line1s + line2s
    plt.grid(True)
    plt.axis('equal')

    def init():
        for line in line1s:
            line.set_data([],[])
        for line in line2s:
            line.set_data([],[])
        return lines
    def animate(i):

        for k in range(0,len(line1s)):
            line1s[k].set_data([0, x1s[k][i]], [0, y1s[k][i]])
            line2s[k].set_data([x1s[k][i], x2s[k][i]], [y1s[k][i], y2s[k][i]])
        return lines

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=2000, interval=20, blit=True)
    plt.show()

animiate_pendulums(res, pendulums)



