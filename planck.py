#! /usr/bin/env python3

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def planck(wavelength, temperature):
    intensity = np.zeros(len(wavelength))
    for i, wl in enumerate(wavelength):
        if wl < 1.0e-20:
            intensity[i] = 0.0
        else:
            frequency = sc.c / wl
            energy = sc.Planck * frequency
            dist = np.exp(energy / (sc.k * temperature)) - 1
            intensity[i] = 2.0 * energy * (frequency / sc.c)**2 / dist
    return intensity

minT = 2500.0
maxT = 50000.0
temperatures = np.linspace(np.log10(minT), np.log10(maxT), 50)
for itemp, temp in enumerate(temperatures):
    temperatures[itemp] = 10**temp

wl1 = np.linspace(0.0, 3.0e-6, 1000)

# Assume R, G and B in RGB lie in middle of respective wavelength ranges
wl2 = np.array([682.6e-9, 532.5e-9, 472.5e-9])
n = 1
for temp in temperatures:
    fig = plt.figure()
    ax1 = plt.axes()
    rad = planck(wl1, temp)
    # Get the intensities of each RGB index to use give approximately the
    # curve the right colour for the given temperature
    rgb = planck(wl2, temp)
    # Normalise the colour
    rgb = rgb / max(rgb)
    plt.plot(wl1*1e6, rad/max(rad), color=rgb, label="{:6.0f}K".format(temp))
    ax1.set_facecolor("black")
    plt.imshow([[1.0, 0.0], [1.0, 0.0]],
               extent=[380e-3, 750e-3, 0.0, 1.0],
               cmap=plt.cm.gist_rainbow,
               interpolation="bicubic",
               aspect="auto")
    ax1.set_xlim(0.0, 3.0)
    ax1.set_xlabel("Wavelength ($\mu m$)")
    ax1.legend(loc="upper right")

    ax2 = ax1.twiny()
    dimfrac = 0.3
    x0, y0 = ax2.transAxes.transform((0, 0))
    x1, y1 = ax2.transAxes.transform((1, 1))
    dx = x1 - x0
    dy = y1 - y0
    maxdim = max(dx, dy)
    circ_wid = dimfrac * maxdim / dx
    circ_hgt = dimfrac * maxdim / dy
    circle = Ellipse([0.67, 0.5], circ_wid, circ_hgt, color=rgb)
    ax2.add_artist(circle)

    plt.savefig('plots/planck_{:04d}.png'.format(n), format='png')
    plt.close(fig)
    n += 1

temperatures_rev = np.flip(temperatures)
for temp in temperatures_rev:
    fig = plt.figure()
    ax1 = plt.axes()
    rad = planck(wl1, temp)
    # Get the intensities of each RGB index to use give approximately the
    # curve the right colour for the given temperature
    rgb = planck(wl2, temp)
    # Normalise the colour
    rgb = rgb / max(rgb)
    plt.plot(wl1*1e6, rad/max(rad), color=rgb, label="{:6.0f}K".format(temp))
    ax1.set_facecolor("black")
    plt.imshow([[1.0, 0.0], [1.0, 0.0]],
               extent=[380e-3, 750e-3, 0.0, 1.0],
               cmap=plt.cm.gist_rainbow,
               interpolation="bicubic",
               aspect="auto")
    ax1.set_xlim(0.0, 3.0)
    ax1.set_xlabel("Wavelength ($\mu m$)")
    ax1.legend(loc="upper right")

    ax2 = ax1.twiny()
    dimfrac = 0.3
    x0, y0 = ax2.transAxes.transform((0, 0))
    x1, y1 = ax2.transAxes.transform((1, 1))
    dx = x1 - x0
    dy = y1 - y0
    maxdim = max(dx, dy)
    circ_wid = dimfrac * maxdim / dx
    circ_hgt = dimfrac * maxdim / dy
    circle = Ellipse([0.67, 0.5], circ_wid, circ_hgt, color=rgb)
    ax2.add_artist(circle)

    plt.savefig('plots/planck_{:04d}.png'.format(n), format='png')
    plt.close(fig)
    n += 1
