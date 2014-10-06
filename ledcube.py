#!/usr/bin/python
# -*- coding: utf8 -*-

'''
Python script for importing, exporting, generating and visualizing LED
sequences for the Velleman 5x5x5 LED cube
The encoding is reverse-engineered, therefore there is no guarantee for
correctness or completeness and no warranty whatsoever!

The generated .cb5 files can be uploaded to the Velleman LED-cube via their
CubeAnimator program.


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from math import *
import numpy as np
import random

s = 5 # cube size 5x5x5

# internal data structure:
# * sequence = [frame, ..., frame]
# * frame = {'array': np.array((s, s, s), dtype=bool),
#            'brightness':int, 'scene_start':bool, 'scene_end':bool}
#
# cb5 data structure:
# * raw_sequence = [raw_frame, ..., raw_frame]
# * raw_frame = [block0, block1, block2, block3, block4]
# * block = [byte0, byte1, byte2, byte3]
    
    
def encode_frame(frame):
    ar = frame['array']
    raw_frame = []
    for z in range(s):
        bytesc = [0] * int(ceil(s**2 / 8.))
        for i in range(s**2):
            if ar[s-1-z][i / s][i % s]:
                bytesc[i / 8] |= 1 << (7 - i%8)
        raw_frame.append(bytesc)
    
    # brightness 0=off, 5=max
    raw_frame[1][3] += int(2**(4 - frame['brightness']))
    
    # scenes
    if frame['scene_start']:
        raw_frame[0][3] += 1
    if frame['scene_end']:
        raw_frame[4][3] += 2
    
    return raw_frame


def encode_sequence(sequence):
    return [encode_frame(frame) for frame in sequence]


def empty_array():
    return np.zeros((s, s, s), dtype=bool)


def new_frame(array=None, bright=None, start=None, end=None):
    frame = {'brightness':5, 'scene_start':False, 'scene_end':False}
    if array != None:
        frame['array'] = np.array(array, dtype=bool)
    else:
        frame['array'] = empty_array()
    if bright != None:
        frame['brightness'] = bright
    if start != None:
        frame['scene_start'] = start
    if end != None:
        frame['scene_end'] = end
    return frame


def decode_frame(raw_frame):
    array = empty_array()
    for z in range(s):
        for i in range(s**2):
            array[s-1-z][i / s][i % s] = raw_frame[z][i / 8] >> (7 - i%8) & 1
        
    frame = {'array': array}
    frame['brightness'] = 4 - int(np.log2(0.5 + raw_frame[1][3] % 32))
    frame['scene_start'] = 0 != raw_frame[0][3] & 1
    frame['scene_end'] = 0 != raw_frame[0][3] & 1 << 1
    
    return frame


def decode_sequence(raw_sequence):
    return [decode_frame(raw_frame) for raw_frame in raw_sequence]


def add_start_end(scene):
    if len(scene) > 0:
        scene[0]['scene_start'] = True
        scene[-1]['scene_end'] = True
    return scene


def import_cb5(fname):
    raw_sequence = []
    with open(fname, 'rb') as f:
        while True:
            raw_frame = f.read(4*s)
            if len(raw_frame) != 4*s:
                break
            raw_sequence.append([[ord(raw_frame[j+4*i])
                for j in range(4)] for i in range(s)])
    return decode_sequence(raw_sequence)


def export_cb5(sequence, fname):
    with open(fname, 'wb') as f:
        print 'saving', fname
        for raw_frame in encode_sequence(sequence):
            for bytes in raw_frame:
                f.write(bytearray(bytes))


def visualize(sequence, dt=0.1):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.animation as animation
    import time
    
    fig = plt.figure()
    ax = Axes3D(fig, aspect='equal')
    fig.add_axes(ax)
    
    x = np.array([i % s for i in range(s**3)], dtype='int')
    y = np.array([(i/s) % s for i in range(s**3)], dtype='int')
    z = np.array([i / s**2 for i in range(s**3)], dtype='int')
    
    def update_leds(num):
        frame = sequence[num]
        array = frame['array']
        b = 1.0 - frame['brightness']/5.
        col = [{False:(1, 1, 1), True:(b, b, 1)}[array[i/s**2][(i/s)%s][i%s]]
            for i in range(s**3)]
        # set_facecolor() has a bug, therefore we replot completely
        plt.gca().cla()
        pts = ax.scatter(x, y, z, s=60, c=col)
        ax.set_xlim([-0.5, s-0.5])
        ax.set_ylim([-0.5, s-0.5])
        ax.set_zlim3d([-0.5, s-0.5])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('LED cube')
        return pts,
    
    ani = animation.FuncAnimation(fig, update_leds, frames=len(sequence),
        interval=1000.*dt, blit=False)
    plt.show()
        


################################# scenes ######################################


def scene_fade():
    sequence = []
    for b in [0, 1, 1, 2, 3, 4, 4, 5, 5, 5, 4, 4, 3, 2, 1, 1, 0, 0]:
        fr = empty_array()
        fr[:,:,:] = 1
        sequence.append(new_frame(fr, bright=b))
    return sequence


def scene_iteration():
    sequence = []
    for i in range(s**3):
        fr = empty_array()
        fr[i / s**2, (i/s) % s, i % s] = True
        sequence.append(new_frame(fr, i%5+1,
            (i==25 or i==0), (i==24 or i==s**3-1)))
    return sequence


def scene_randomwalk():
    sequence = []
    c = [2, 2, 2]
    for i in range(50):
        ar = empty_array()
        ar[c[0], c[1], c[2]] = True
        sequence.append(new_frame(ar, 5))
        
        r = random.randint(0, 2)
        d = -1 + 2 * random.randint(0, 1)
        c[r] = s-1 - abs(s-1 - abs(c[r] + d))
    return sequence


def scene_frame():
    sequence = []
    for j in [2, 1, 1, 0, 0, 0, 1, 1, 2, 2]:
        ar = empty_array()
        for k in range(s**3):
            x, y, z = k%s, (k/s)%s, k/s**2
            d = [abs(l-s/2) for l in x, y, z]
            ar[z, y, x] = d.count(j) >= 2 and max(d) == j
        sequence.append(new_frame(ar))
    return sequence


def scene_planeupdown():
    sequence = []
    for i in range(2 * (s-1)):
        ar = empty_array()
        ar[s-1 - abs(s-1 - i),:,:] = 1
        sequence.append(new_frame(ar))
    return sequence


def scene_cornerzoom(inv_x=False, inv_y=False, inv_z=False):
    sequence = []
    for i in range(2 * s):
        ar = empty_array()
        for k in range(s**3):
            x, y, z = k%s, (k/s)%s, k/s**2
            if inv_x: x2 = s-1 - x
            else: x2 = x
            if inv_y: y2 = s-1 - y
            else: y2 = y
            if inv_z: z2 = s-1 - z
            else: z2 = z
            ar[z, y, x] = all([j > i-s and j <= i for j in x2, y2, z2])
        sequence.append(new_frame(ar))
    return sequence


def scene_snake(n=70):
    sequence = []
    snake = [[s/2, s/2, s/2]]
    for i in range(n):
        ar = empty_array()
        for p in snake:
            ar[p[0], p[1], p[2]] = 1
        sequence.append(new_frame(ar))
        neighbours = [list(snake[0] + np.array(k)) for k in
            [[0,0,1], [0,0,-1], [0,1,0], [0,-1,0], [1,0,0], [-1,0,0]]]
        neighbours = [p for p in neighbours if
            (p not in snake) and (-1 not in p) and (s not in p)]
        snake = ([random.choice(neighbours)] + snake)[:5]
    return sequence


def scene_rain(n=80, density=0.06):
    sequence = []
    ar = empty_array()
    for i in range(n):
        ar[s-1,:,:] = np.random.random((s,s)) < density
        sequence.append(new_frame(ar))
        ar = np.roll(ar, -1, 0)
    return sequence


def scene_skyrocket(n=1, m=12, v0=0.5):
    sequence = []
    for i in range(n):
        x0 = random.randint(1, s-2)
        y0 = random.randint(1, s-2)
        x1 = random.randint(1, s-2)
        y1 = random.randint(1, s-2)
        
        # rocket launch
        for z in range(s):
            if z >= s-2: z -= 1
            ar = empty_array()
            x = int(0.5 + x0 + (x1 - x0) * z / float(s-1))
            y = int(0.5 + y0 + (y1 - y0) * z / float(s-1))
            ar[z, y, x] = True
            sequence.append(new_frame(ar))
        
        # explosion
        v = []
        for j in range(m):
            phi = 2*pi*random.random()
            vz = -1 + 2 * random.random()
            vx = sqrt(1-vz**2) * sin(phi)
            vy = sqrt(1-vz**2) * cos(phi)
            v.append([vx, vy, vz])
        ne = 2*s + 3
        for j in range(1, ne):
            ar = empty_array()
            for vi in v:
                x = int(0.5 + x1 + vi[0] * v0 * j)
                y = int(0.5 + y1 + vi[1] * v0 * j)
                z = int(0.5 + s-2 + vi[2] * v0 * j - 0.02 * j**2)
                if all([l >= 0 and l < s for l in x, y, z]):
                    ar[z, y, x] = True
            sequence.append(new_frame(ar, min(s, ne-1-j)))
    return sequence


############################## custom code ####################################


scene = (add_start_end(scene_fade())
    + add_start_end(scene_planeupdown() + scene_planeupdown())
    + add_start_end(sum([scene_cornerzoom(*inv) for inv in
      [[0,0,0],[0,0,1],[0,1,1],[1,1,1],[1,0,1],[1,0,0],[1,1,0],[0,1,0]]], []))
    + add_start_end(scene_frame()+scene_frame()+scene_frame()+scene_frame())
    + add_start_end(scene_rain())
    + add_start_end(scene_skyrocket(10))
    + add_start_end(scene_snake()))

export_cb5(scene, 'scene.cb5')

visualize(scene, dt=0.1)
