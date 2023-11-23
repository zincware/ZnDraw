import numpy as np


def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


def _get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)
