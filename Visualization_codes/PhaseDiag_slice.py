import matplotlib.pyplot as plt
import numpy as np
"""
Blue_x = [1e-14, 1e-13, 1e-12, 1e-11]
Blue_y = [pow(10, -15.75)]*4
Orange_x = Blue_x
Orange_y = [pow(10,-3)]*4
Purple_x = np.linspace(-10,-5,6)
Purple_y = [pow(10, -17.7), pow(10, -15.75), pow(10, -14.5), pow(10, -13.2), pow(10, -11.75), pow(10, -10.3)]
Green_x = np.linspace(-7,-5,3)
Green_y = [pow(10, -16.9), pow(10, -15.0), pow(10, -12.9)]

fig, TI = plt.subplots()

TI.scatter(Blue_x, Blue_y, s=10, c='blue', label='Blue')
TI.set_xscale('log')
TI.set_yscale('log')
TI.grid()
TI.xaxis.grid(which='minor')
plt.legend(loc='upper left')
plt.show()
"""


dt = 0.1
x = np.arange(-500.0, 500.0, dt)
y1 = np.arange(0, 1000.0, dt)
y3 = np.sin(x / 3.0)

fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4)

ax0.plot(x, y1)
ax0.set_xscale('symlog')
ax0.set_ylabel('symlogx')
ax0.grid()
ax0.xaxis.grid(which='minor')

ax1.plot(y1, x)
ax1.set_yscale('symlog')
ax1.set_ylabel('symlogy')

ax2.plot(x, y3)
ax2.set_xscale('symlog')
ax2.set_yscale('symlog', linthresh=0.015)
ax2.grid()
ax2.set_ylabel('symlog both')

ax3.plot(x, y3)
ax3.set_xscale('symlog')
ax3.set_yscale('linear')
ax3.grid()

plt.show()
