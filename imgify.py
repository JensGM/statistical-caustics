#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-x',
                    '--line-width',
                    type=int)
parser.add_argument('-y',
                    '--line-count',
                    type=int)
parser.add_argument('-z',
                    '--frames',
                    type=int,
                    default=1)
parser.add_argument('-f',
                    '--framerate',
                    type=int,
                    default=30)
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


def parse_imgs(line_width,
               line_count,
               frames,
               apply_exp=False,
               apply_log=False,
               apply_abs=False):
    data = sys.stdin.read().split()
    imgs = []
    for frame in range(frames):
        img = []
        for i in range(0, line_count * line_width, line_width):
            idx = frame * line_count * line_width + i
            line = list(map(float, data[idx:idx + line_width]))
            img.append(line)
        img = np.array(img)
        if apply_exp: img = np.exp(img)
        if apply_log: img = np.log(img)
        if apply_abs: img = np.abs(img)
        imgs.append(img)
    return imgs


def draw(imgs, frames, cmap='gray'):
    if frames == 1:
        plt.imshow(imgs[0], cmap=cmap)
        plt.gca().invert_yaxis()
        plt.show()
        return

    fig = plt.figure()
    pltimgs = []
    for i in range(frames):
        img = plt.imshow(imgs[i], cmap=cmap, animated=True)
        pltimgs.append([img])
    ani = animation.ArtistAnimation(fig, pltimgs, interval=50, blit=True)
    plt.show()


if __name__ == '__main__':
    args = parser.parse_args()
    imgs = parse_imgs(args.line_width, args.line_count, args.frames);
    draw(imgs, args.frames, args.color_map)
