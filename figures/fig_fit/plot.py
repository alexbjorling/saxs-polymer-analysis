"""
BMC biology: width 85 or 175 mm, so 3.35 or 6.89 inches
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os


font = {'size': 12}
matplotlib.rc('font', **font)
plt.ion()

# Form factor from Crysol
crysol = np.loadtxt('../../pdb/5lvy_dimer/5lvy_dimer.int', skiprows=1)[:, :2]
F2_form = crysol[:, 1]
F2_form[:] = F2_form / F2_form[0]

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,
                       figsize=(6.89, 3))
LEFT, WSPACE, RIGHT = .12, .05, .98
plt.subplots_adjust(
    left=LEFT, right=RIGHT, top=.99, bottom=.2, wspace=WSPACE
    )

samples = [
    {'name': 'fimbriae_A', 'conc': 7.14e-3, 'color': 'dodgerblue',
     'scale': 1.15, 'num': 221, 'sub': 6e-2, 'start': 20000,
     'd': 5.0},
    {'name': 'fimbriae_AB', 'conc': 12.2e-3, 'color': 'orange',
     'scale': .60, 'num': 226, 'sub': 1.5e-2, 'start': 20000,
     'd': 5.5}
]

for i, sample in enumerate(samples):
    # Read simulations
    files = os.listdir(sample['name'])
    runs = [n.split('.')[0] for n in files if n.endswith('.npz')]
    print('%s: loading' % sample['name'], runs)
    steps, bonds, angles, F2_debye = [], [], [], []
    for run in runs:
        steps_, bonds_, angles_ = np.loadtxt(
            '%s/%s.traj' % (sample['name'], run), skiprows=1
        ).T
        steps.append(steps_)
        bonds.append(bonds_)
        angles.append(angles_)
        dct = np.load('%s/%s.npz' % (sample['name'], run))
        q = dct['q']
        assert np.allclose(crysol[:, 0], q)
        F2_debye.append(
            np.mean(
                dct['I'][np.where(steps_ >= sample['start'])[0][0]:], axis=0
            )
        )

    # Experimental data in absolute units
    pattern = '../../data/RBNX_%d_sub.dat'
    data = np.loadtxt(pattern % sample['num'], skiprows=3)[:, :-2]
    data[:, 1] /= sample['conc']
    data[:, 1] -= sample['sub']
    ax[i].plot(data[:, 0], data[:, 1], '.', ms=7, color=sample['color'])

    # Model in absolute units
    drho = 2e10  # cm/g
    N_A = 6.022e23  # 1/mol
    M = 16.5e3 * 2 * 40  # g / mol
    absfactor = sample['scale'] * M * (drho**2) / N_A
    ax[i].plot(q, np.mean(F2_debye, axis=0) * F2_form * absfactor, 'k')
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].set_xlim(data[0, 0] / 2, data[-1, 0] * 1.5)
    ax[i].set_ylim(3e-2, 2e3)

    # trajectory inset
    ax_width = (1 - LEFT - WSPACE - (1 - RIGHT)) / 2
    pos = .01 + LEFT + (ax_width + WSPACE - .01) * i
    a2 = fig.add_axes(rect=[pos, .23, .2, .3])
    a2.set_xlim((0, 50000))
    for j in range(len(runs)):
        a2.plot(steps[j], angles[j] / np.pi * 180, color=[.35, .35, .35], lw=.1)

    a2.xaxis.tick_top()
    a2.xaxis.set_label_position('top')
    a2.xaxis.set_tick_params(labelsize=8, pad=0)
    a2.axvline(sample['start'], linestyle='--', color='k', lw=.5)
    a2.set_xlabel('Monte Carlo step', fontsize=8, labelpad=2, ha='left').set_x(0)

    a2.yaxis.tick_right()
    a2.yaxis.set_label_position('right')
    a2.yaxis.set_tick_params(labelsize=8, pad=0)
    a2.set_ylim((0, 70))
    a2.set_ylabel('Av. angle ($^\circ$)', fontsize=8, labelpad=1, ha='left').set_y(0)

    # cartoon inset
    ax_width = (1 - LEFT - WSPACE - (1 - RIGHT)) / 2
    pos = .22 + .01 + LEFT + (ax_width + WSPACE - .01) * i
    a2 = fig.add_axes(rect=[pos, .4, .25, .55])

    # read image and convert white to transparent
    im = plt.imread(sample['name'] + '/render.png')[300:640, 300:500]
    im = np.pad(im, ((0, 0), (0, 0), (0, 1)), constant_values=1)
    white = np.where(im.sum(axis=-1) > 2.999)
    im[white[0], white[1], -1] = 0
    a2.imshow(im)
    a2.axis('off')

    # scalebar
    # the whole chain is 597 pixels - manually measured on the rendering!
    chain_length_pix = 597
    # the whole chain is this long
    chain_length_nm = 39 * sample['d']
    pix_per_nm = chain_length_pix / chain_length_nm
    scalebar_nm = 20
    scalebar_pix = scalebar_nm * pix_per_nm
    off = i * 7
    a2.plot([80 - off, 80 - off + scalebar_pix], [20, 20], 'k', lw=4)
    a2.text(80 - off + scalebar_pix / 2, 10, '%u nm' % scalebar_nm,
            fontsize=8, ha='center', va='bottom')

ax[0].set_xlabel('$q = 4\\pi/\\lambda$ $\\sin\\theta$  [1/Å]').set_x(1.)
ax[0].set_ylabel('$I / c$  [$\mathrm{cm}^2 / \mathrm{g}$]')
ax[1].set_yticklabels([])

plt.savefig('fig_fit.pdf')
