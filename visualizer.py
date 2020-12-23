import matplotlib.pyplot as plt, numpy as np
from mpl_toolkits.mplot3d import proj3d
import pandas as pd
from scipy.spatial import ConvexHull
import random
from itertools import combinations

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import animation

class Faces():
    def __init__(self, tri):
        self.tri = np.array(tri)
        norms = [self.norm(s) for s in self.tri]
        _, self.inv = np.unique(norms,return_inverse=True, axis=0)

    def norm(self,sq):
        cr = np.cross(sq[2]-sq[0],sq[1]-sq[0])
        return np.around(np.abs(cr/np.linalg.norm(cr)),5)

    def isneighbour(self, tr1,tr2):
        a = np.concatenate((tr1,tr2), axis=0)
        return len(a) == len(np.unique(a, axis=0))+2

    def triangle_intersect(self, tr1, tr2):
        pts = []
        for pt in tr1:
            if (pt == tr2).all(1).any():
                pts.append(pt)
        return(pts)

    def polytope(self):
        edge_list = []
        for t in self.tri:
            for nbr in self.tri:
                if self.isneighbour(t,nbr):
                    if not np.array_equal(self.norm(t),self.norm(nbr)):
                        edge_list.append(self.triangle_intersect(t,nbr))
                    #else:
                        #face_list.append(np.unique(np.concatenate((t,nbr), axis=0)))
        return(edge_list)

def norm(sq):
    cr = np.cross(sq[2]-sq[0],sq[1]-sq[0])
    return np.around(np.abs(cr/np.linalg.norm(cr)),5)


def plot_convex_hull(GERMS):

    CR = GERMS[GERMS['CR-lam']>.5]
    x = CR["twist_a"].tolist()
    y = CR["twist_b"].tolist()
    z = CR["twist_c"].tolist()
    verts = np.array([x,y,z]).transpose()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    hull = ConvexHull(verts)
    simplices = hull.simplices
    
    tri_vecs = [verts[s] for s in simplices]
    ext_tri_vecs = [verts[s]+0.02*verts[s] for s in simplices]
    f = Faces(ext_tri_vecs)

    norm_vec = [norm(tri) for tri in tri_vecs]
    norms = np.unique(norm_vec, axis=0)
    colours = np.random.rand(len(norms),3)
    edges = f.polytope()

    # draw the polygons of the convex hull
    for s in hull.simplices:
        tri = Poly3DCollection(verts[s])
        col = colours[(norms == norm(verts[s])).all(1)][0]
        tri.set_facecolor(col)
        tri.set_alpha(0.6)
        ax.add_collection3d(tri)

    for e in edges:
        edg = Poly3DCollection(e)
        edg.set_color('k')
        edg.set_alpha(1)
        ax.add_collection3d(edg)

    # draw the vertices
    #
    def stretch(vec, percent):
        return([(1+percent)*val for val in vec])

    # draw all completions
    x = GERMS["twist_a"].tolist()
    y = GERMS["twist_b"].tolist()
    z = GERMS["twist_c"].tolist()
    #ax.scatter(x,y,z, marker='.', color='yellow', alpha=0.2)
    
    # draw the matching twists
    TWISTMATCH = GERMS[GERMS['CR-lam']>.0]
    x = TWISTMATCH["twist_a"].tolist()
    y = TWISTMATCH["twist_b"].tolist()
    z = TWISTMATCH["twist_c"].tolist()
    ax.scatter(x,y,z, s=5, marker='o', color='blue', alpha=0.7)

    # draw the max. chain recurrent
    x = stretch(CR["twist_a"].tolist(),0.02)
    y = stretch(CR["twist_b"].tolist(),0.02)
    z = stretch(CR["twist_c"].tolist(),0.02)
    ax.scatter(x,y,z, s=15, marker='o', color='red', alpha=0.7)

    def rotate(angle):
        ax.view_init(azim=angle)

    angle = 3
    ani = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 360, angle), interval=50)
    ani.save('chamf_cube.gif', writer=animation.PillowWriter(fps=20))

    plt.show()


def visualize3DData (envelope):
    """Visualize data in 3d plot with popover next to mouse position.

    Args:
        envelope (pd.DataFrame) - dataframe: "twist_a", "twist_b" and "twist_c" are the x,y,z
            and "lamination" is the triangulation data in a string
    Returns:
        None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    x = list(envelope["twist_a"])
    y = list(envelope["twist_b"])
    z = list(envelope["twist_c"])
    X = np.array((x,y,z), dtype = float).transpose()
    print(X)

    #hull = ConvexHull(X, qhull_options = "QJ")
    hull = ConvexHull(X)

    # draw the polygons of the convex hull
    for s in hull.simplices:
        tri = Poly3DCollection(X[s])
        tri.set_color('r')
        tri.set_color([random.random(),random.random(), random.random()])
        tri.set_alpha(0.7)
        ax.add_collection3d(tri)

    ax.scatter(x, y, z, depthshade = False, picker = True)


    def distance(point, event):
        """Return distance between mouse position and given data point

        Args:
            point (np.array): np.array of shape (3,), with x,y,z in data coords
            event (MouseEvent): mouse event (which contains mouse position in .x and .xdata)
        Returns:
            distance (np.float64): distance (in screen coords) between mouse pos and data point
        """
        assert point.shape == (3,), "distance: point.shape is wrong: %s, must be (3,)" % point.shape

        # Project 3d data space to 2d data space
        x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], plt.gca().get_proj())
        # Convert 2d data space to 2d screen space
        x3, y3 = ax.transData.transform((x2, y2))

        return np.sqrt ((x3 - event.x)**2 + (y3 - event.y)**2)


    def calcClosestDatapoint(X, event):
        """"Calculate which data point is closest to the mouse position.

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            event (MouseEvent) - mouse event (containing mouse position)
        Returns:
            smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
        """
        distances = [distance (X[i, 0:3], event) for i in range(X.shape[0])]
        return np.argmin(distances)


    def annotatePlot(X, index):
        """Create popover label in 3d chart

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            index (int) - index (into points array X) of item which should be printed
        Returns:
            None
        """
        # If we have previously displayed another label, remove it first
        if hasattr(annotatePlot, 'label'):
            annotatePlot.label.remove()
        # Get data point from array of points X, at position index
        x2, y2, _ = proj3d.proj_transform(X[index, 0], X[index, 1], X[index, 2], ax.get_proj())
        annotatePlot.label = plt.annotate( "Lamination: " + list(envelope["lamination"])[index],
            xy = (x2, y2), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        fig.canvas.draw()


    def onMouseMotion(event):
        """Event that is triggered when mouse is moved. Shows text annotation over data point closest to mouse."""
        closestIndex = calcClosestDatapoint(X, event)
        annotatePlot (X, closestIndex)

    fig.canvas.mpl_connect('motion_notify_event', onMouseMotion)  # on mouse motion
    plt.show()


if __name__ == '__main__':
    envelope = pd.DataFrame({"lamination" : ["a", "b", 
        "c", "d"], "twist_a" : [1.45,2.1,2.1,2.1], "twist_b" : [-1.2,0,0,0], "twist_c" : [1.4,0,-1.1,-1]})
    print(envelope)
    plot_convex_hull(envelope)
