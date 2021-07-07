import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import sys
import glob


param_filename = 'param_1.5.3_2D.txt'
filename_prefix = 'step_'
interfaces_filename = 'interfaces_creep.txt'

surface_filenames = sorted(glob.glob(f'sp_surface_global_*'))

nproc = len(glob.glob('step_0-rank_new*'))

s = np.array([])
t = np.array([])

size = len(surface_filenames)
progress = 0

for f in surface_filenames:
    print(f'{100*progress/(size-1):3.0f}%', end='\r')
    progress += 1


    step = int(f.split('.')[0].split('_')[-1])

    _s = np.loadtxt(f, skiprows=2, comments='P')
    s = np.append(s, np.abs(np.max(_s)))

    if step == 0 and int(s) == 0:
        s = 1.5e5

    # read time
    with open(f'time_{step}.txt') as time_file:
        time_info = time_file.readline()
        _t = float(time_info.split(':')[-1])
        t = np.append(t, _t)

fig, ax = plt.subplots(figsize=(12, 6))

s = 1.5e5 - s

ax.plot(t/1.0e6, s, linestyle='-', color='C0', label='MANDYOC')

x_underworld, y_underworld = np.loadtxt('UNDERWORLD.txt', unpack=True, delimiter = ',')
ax.plot(x_underworld, y_underworld, color='C1', label='UNDERWORLD')

x_stagyy, y_stagyy = np.loadtxt('STAGYY.txt', unpack=True, delimiter = ',')
ax.plot(x_stagyy, y_stagyy, color='C2', label='STAGYY')

x_i2ves, y_i2ves = np.loadtxt('I2VES.txt', unpack=True, delimiter = ',')
ax.plot(x_i2ves, y_i2ves, color='C4', label='I2VIS')

plt.legend(loc='upper left')
ax.grid('on')

ax.set_xlim((0, 20))
ax.set_ylim((0, 900))

ax.set_xlabel('Time (Myr)')
ax.set_ylabel('Maximum topography (m)')

# inset axes....
axins = ax.inset_axes([0.5, 0.1, 0.45, 0.45])
axins.plot(t/1.0e6, s, linestyle='-', color='C0')
axins.plot(x_underworld, y_underworld, color='C1', label='UNDERWORLD')
axins.plot(x_stagyy, y_stagyy, color='C2', label='STAGYY')
axins.plot(x_i2ves, y_i2ves, color='C4', label='I2VIS')

# sub region of the original image
x1, x2, y1, y2 = 0, 0.2, 0, 400
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.grid('on')
ax.indicate_inset_zoom(axins)

figname = f'fig-max-surface.png'
fig.savefig(figname)
print(f'\nSaved {figname}')
plt.close('all')

# save data
with open('data-surface-vec.txt', 'w') as f:
    np.savetxt(f, np.stack((t, s), axis=1), fmt='%.6e')
