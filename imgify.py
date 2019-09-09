#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-x',
                    '--line-width',
                    type=int)
parser.add_argument('-c',
                    '--color-map',
                    default='gray')
parser.add_argument('-e',
                    '--exp',
                    default=False,
                    action='store_true')
parser.add_argument('-l',
                    '--log',
                    default=False,
                    action='store_true')
parser.add_argument('-a',
                    '--abs',
                    default=False,
                    action='store_true')


def parse_img(line_width, apply_exp=False, apply_log=False, apply_abs=False):
    data = sys.stdin.read().split()
    img = []
    for i in range(0, len(data), line_width):
        line = list(map(float, data[i:i + line_width]))
        img.append(line)
    img = np.array(img)

    if apply_exp: img = np.exp(img)
    if apply_log: img = np.log(img)
    if apply_abs: img = np.abs(img)

    return img


def draw(img, cmap='gray'):
    plt.imshow(img, cmap=cmap)
    plt.gca().invert_yaxis()
    plt.show()


if __name__ == '__main__':
    args = parser.parse_args()
    img = parse_img(args.line_width);
    draw(img, args.color_map)
