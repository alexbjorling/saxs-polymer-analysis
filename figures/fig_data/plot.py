"""
BMC biology: width 85 or 175 mm, so 3.35 or 6.89 inches
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

font = {'size': 12}
matplotlib.rc('font', **font)
plt.ion()


def read(*nums, conc=1):
    pattern = '../../data/RBNX_%d_sub.dat'
    dilutions = [1, 2, 4, 10]
    rtn = []
    for i in range(len(nums)):
        data = np.loadtxt(pattern % nums[i], skiprows=3)[:, :-2]
        data[:, 1] /= conc / dilutions[i]
        rtn.append(data)
    return rtn


# load data
mono_A = read(210, 212, conc=5.01e-3)
mono_AB = read(215, 218, conc=7.31e-3)
fimb_A = read(221, 223, 296, conc=7.14e-3)
fimb_AB = read(226, 228, 299, conc=12.2e-3)

# plot data
fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(6.89, 3))

for data in mono_A:
    ax[0].plot(data[:, 0], data[:, 1], color='royalblue')

for data in mono_AB:
    ax[0].plot(data[:, 0], data[:, 1], color='darkorange')

for data in fimb_A:
    ax[1].plot(data[:, 0], data[:, 1], color='royalblue')

for data in fimb_AB:
    ax[1].plot(data[:, 0], data[:, 1], color='darkorange')

# annotate
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_ylim([.02, 2e3])
ax[0].set_xlabel('$q = 4\\pi/\\lambda$ $\\sin\\theta$  [1/Å]').set_x(1)
ax[0].set_ylabel('$I / c$  [$\mathrm{cm}^2 / \mathrm{g}$]')
xl = ax[0].get_xlim()
# ax[0].text(xl[0] * 1.3, 6e2, 'A', weight='bold', ha='left')
# ax[1].text(xl[1] / 1.3, 6e2, 'B', weight='bold', ha='right')
plt.subplots_adjust(wspace=.05, bottom=.2, top=.99, left=.12, right=.99)

# make ideal lines
x = np.array([.01, .08])
ax[1].plot(x, x**(-1) / 20, 'k', lw=.5)
ax[1].text(.030, 60 / 20, 'rod (I $\propto q^{-1})$',
           fontsize=10, rotation=-21, ha='center', va='center')
x = np.array([.015, .18])
ax[1].plot(x, x**(-2) / 6, 'k', lw=.5)
ax[1].text(.055, 120, 'ideal coil (I $\propto q^{-2})$',
           fontsize=10, rotation=-37, ha='center', va='center')

# make a legend somehow
def legend(ax, name, units):
    h = .04
    l = .02
    w = 1.5
    ax.text(.004, h, name, va='center')
    ax.plot([l, l * w], [h, h], color='royalblue')
    ax.text(l * w * 1.1, h, 'A', va='center')
    ax.plot([l * 2 * w, l * 2 * w**2], [h, h], color='darkorange')
    ax.text(l * 2 * w**2 * 1.1, h, units, va='center')
legend(ax[0], 'Subunits', 'AB')
legend(ax[1], 'Fimbriae', 'A+B')

# save
plt.savefig('fig_data.pdf')
