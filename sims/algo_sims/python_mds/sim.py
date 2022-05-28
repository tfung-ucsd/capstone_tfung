"""
Author: Trevor Fung
Date: 4/28/22
Filename: sim.py
Purpose: play a video of mds working with simulated data
"""

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

import math
import time
import code

METER_SCALING = 3
AVERAGING_FACTOR = 10
np.random.seed(0)
counter = 0

def generate_similarities():
    """creates simmed data to feed the visualization
    """
    num_nodes = 3
    # X_true = np.array([[0.0, 0.0],[1.0, 0.0],[0.5, 0.75**0.5]]) #create equilateral triangle
    X_true = np.array([[0.0, 0.0],[3.0/METER_SCALING, 0.0],[7.0/METER_SCALING, 0.0]]) #create straight line
    X_true *= METER_SCALING

    # Center the data
    # X_true -= X_true.mean()

    # #rotate input pts
    # theta = 180*math.pi/180
    # for i in range(len(X_true)):
    #     pt = X_true[i]
    #     temp = np.array([0.0, 0.0])
    #     temp[0] = math.cos(theta)*pt[0]-math.sin(theta)*pt[1]
    #     temp[1] = math.sin(theta)*pt[0]+math.cos(theta)*pt[1]
    #     X_true[i] = temp

    similarities = euclidean_distances(X_true)

    # Add noise to the similarities
    avg_noise = np.zeros((num_nodes, num_nodes))
    for i in range(AVERAGING_FACTOR):
        noise = np.random.normal(loc=0, scale=1.0, size=(num_nodes, num_nodes)) #NxN
        noise = noise + noise.T #make symmetric (alg limitation)
        noise[np.arange(noise.shape[0]), np.arange(noise.shape[0])] = 0 #set diagonals to zero
        avg_noise += noise
    avg_noise /= AVERAGING_FACTOR
    print(avg_noise)
    # code.interact(local=locals())

    similarities += avg_noise
    print(similarities)

    return X_true, similarities

def gen_audio_data_1():
    global counter
    values1 = [33, 20, 28, 27, 29, 28, 29, 25, 30, 24, 37, 40, 20, 34, 88, 86, 78, 86, 80, 74, 80, 79, 71, 81, 88, 88, 89, 72, 72, 74]
    out = values1[counter]
    return out

def gen_audio_data_2():
    global counter
    values2 = [25, 20, 16, 26, 24, 6, 24, 33, 6, 27, 36, 2, 11, 16, 67, 78, 83, 60, 68, 79, 67, 71, 74, 70, 75, 79, 68, 65, 58, 66]
    out = values2[counter]
    counter += 1
    return out

def single_mds_iter(X_true, similarities):
    """based on https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#sphx-glr-auto-examples-manifold-plot-mds-py
        # Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
        # License: BSD
    """

    mds = manifold.MDS(
        n_components=2,
        max_iter=3000,
        eps=1e-9,
        random_state=3,
        dissimilarity="precomputed",
        n_jobs=1,
    )
    pos = mds.fit(similarities).embedding_

    # Rescale the data
    # pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())

    # Rotate and center the output
    # clf = PCA(n_components=2)
    # pos = clf.fit_transform(pos)
    code.interact(local=locals())
    pos -= pos[0] #crappy attempt to anchor one end
    # pos -= pos.mean()

    # fig = plt.figure(1)
    ax = plt.axes([0.0, 0.0, 1.0, 1.0])

    s = 100
    plt.scatter(X_true[:, 0], X_true[:, 1], color="navy", s=s, lw=0, label="True Position")
    plt.scatter(pos[:, 0], pos[:, 1], color="turquoise", s=s, lw=0, label="MDS")
    plt.legend(scatterpoints=1, loc="best", shadow=False)

    np.fill_diagonal(similarities, 0)
    # Plot the edges
    start_idx, end_idx = np.where(pos)
    # a sequence of (*line0*, *line1*, *line2*), where::
    #            linen = (x0, y0), (x1, y1), ... (xm, ym)
    segments = [
        [X_true[i, :], X_true[j, :]] for i in range(len(pos)) for j in range(len(pos))
    ]
    values = np.abs(similarities)
    lc = LineCollection(
        segments, zorder=0, cmap=plt.cm.Blues, norm=plt.Normalize(0, values.max())
    )
    lc.set_array(similarities.flatten())
    lc.set_linewidths(np.full(len(segments), 0.5))
    ax.add_collection(lc)

    # label original pts
    for pt in X_true:
        x = pt[0]
        y = pt[1]
        plt.text(x, y+0.25, '({:.1f}, {:.1f})'.format(x, y))

    # label calculated pts
    for pt in pos:
        x = pt[0]
        y = pt[1]
        plt.text(x, y-0.25, '({:.1f}, {:.1f})'.format(x, y))

    #label distances between orig nodes
    for x in segments:
        length = euclidean_distances(x).sum()/2 #because it measures both ways idk
        if length != 0:
            text_loc = np.average(x,axis=0)
            x = text_loc[0]
            y = text_loc[1]
            plt.text(x, y, '{:.1f}'.format(length))

    #apply sensor data
    # for pt in pos:
    #     x = pt[0]
    #     y = pt[1]
    #     plt.text(x, y-0.25, '({:.1f}, {:.1f})'.format(x, y))

    bounds = 4 * METER_SCALING
    plt.xlim([-bounds, bounds])
    plt.ylim([-bounds, bounds])

    # plt.show()
    plt.draw()
    plt.pause(0.001)
    # code.interact(local=locals())

if __name__ == "__main__":
    plt.ion()
    while True:
        X_true, similarities = generate_similarities()
        single_mds_iter(X_true, similarities)
        time.sleep(0.5)