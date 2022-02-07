from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import time
from multiprocessing import Pool
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation as rot
from collections import defaultdict

def tan(x):
    return np.tan(x)

def ctan(x):
    return 1 / np.tan(x)

def plot_3D(my_set):
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(my_set))
    #scale = my_set.flatten()
    axes.auto_scale_xyz([1.06, 1.19], [0.26, 0.31], [-0.01, 0.01])
    #axes.set_xlim(1.06, 1.19)
    #axes.set_ylim(0.26, 0.31)
    #axes.set_zlim(-0.01, 0.01)
    axes.set_xlabel("X")
    axes.set_ylabel("Y")
    axes.set_zlabel("Z")
    #axes.set_aspect("equal")
    # Show the plot to the screen
    pyplot.show()

# Load the STL files and add the vectors to the plot
my_mesh = mesh.Mesh.from_file('./150_finale.stl')

# Reduce the mesh to a subset
#my_mesh.data = my_mesh.data[0 : int(len(my_mesh.data) / 2)]

### https://stackoverflow.com/questions/26303878/alpha-shapes-in-3d
def alpha_shape_3D(points, alpha=2.0):
    # Delaunay-triangulated surface
    tetra = Delaunay(points)
    tetrapos = np.take(points, tetra.vertices, axis=0)
    normsq = np.sum(tetrapos**2, axis=2)[:, :, None]
    ones = np.ones((tetrapos.shape[0], tetrapos.shape[1], 1))
    a = np.linalg.det(np.concatenate((tetrapos, ones), axis=2))
    Dx = np.linalg.det(np.concatenate((normsq, tetrapos[:, :, [1, 2]], ones), axis=2))
    Dy = -np.linalg.det(np.concatenate((normsq, tetrapos[:,:,[0,2]], ones), axis=2))
    Dz = np.linalg.det(np.concatenate((normsq, tetrapos[:,:,[0,1]], ones), axis=2))
    c = np.linalg.det(np.concatenate((normsq, tetrapos), axis=2))
    r = np.sqrt(Dx**2 + Dy**2 + Dz**2 - 4*a*c) / (2 * np.abs(a))
    tetras = tetra.vertices[r < alpha, :]
    TriComb = np.array([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
    Triangles = tetras[:, TriComb].reshape(-1, 3)
    Triangles = np.sort(Triangles, axis=1)
    TrianglesDict = defaultdict(int)
    for tri in Triangles:TrianglesDict[tuple(tri)] += 1
    Triangles = np.array([tri for tri in TrianglesDict if TrianglesDict[tri] ==1])
    EdgeComb = np.array([(0, 1), (0, 2), (1, 2)])
    Edges = Triangles[:,EdgeComb].reshape(-1,2)
    Edges = np.sort(Edges,axis=1)
    Edges = np.unique(Edges,axis=0)
    Vertices = np.unique(Edges)
    #return Vertices,Edges,Triangles
    return Vertices


def makeBox(origin, delta):
    points = np.array([[origin[0], origin[1], origin[2]],
                       [origin[0] + delta[0], origin[1], origin[2]],
                       [origin[0] + delta[0], origin[1] + delta[1], origin[2]],
                       [origin[0], origin[1] + delta[1], origin[2]],
                       [origin[0], origin[1], origin[2]],
                       [origin[0] + delta[0], origin[1], origin[2] + delta[2]],
                       [origin[0] + delta[0], origin[1] + delta[1], origin[2] + delta[2]],
                       [origin[0], origin[1] + delta[1], origin[2] + delta[2]]])
    return points

#points = makeBox([-100, -400, -400], [350, 800, 800])

def make_box_manually():
    L = 0.5
    r = 0.26
    R = 0.31
    delta_theta = 2.5
    delta_theta = np.radians(delta_theta)

    L = L * 1000
    r = r * 1000
    R = R * 1000

    points = []
    points.append([-L/2,  r*np.cos(delta_theta/2), -r*np.sin(delta_theta)])
    points.append([-L/2,  r*np.cos(delta_theta/2),  r*np.sin(delta_theta/2)])
    points.append([-L/2,  R*np.cos(delta_theta/2), -R*np.sin(delta_theta)])
    points.append([-L/2,  R*np.cos(delta_theta/2),  R*np.sin(delta_theta/2)])
    points.append([ L/2,  r*np.cos(delta_theta/2), -r*np.sin(delta_theta)])
    points.append([ L/2,  r*np.cos(delta_theta/2),  r*np.sin(delta_theta/2)])
    points.append([ L/2,  R*np.cos(delta_theta/2), -R*np.sin(delta_theta)])
    points.append([ L/2,  R*np.cos(delta_theta/2),  R*np.sin(delta_theta/2)])
    points = np.array(points)

    return points


Dx = -1084
Dy = 0
Dz = 0


def make_delimiting_surface(source, Dx, Dy, Dz):

    def make_circular(my_set, n_theta_steps):
        for i in range(len(my_set)):
            #v_length = np.linalg.norm(my_set[i])
            for j in range(1, n_theta_steps):
                theta = 2 * np.pi / n_theta_steps * j
                r = rot.from_rotvec([theta, 0, 0])
                #new_point = np.array([[my_set[i, 0], v_length * np.sin(theta), v_length * np.cos(theta)]])
                new_point = np.array([r.apply(my_set[i])])
                my_set = np.append(my_set, new_point, axis = 0)
        return my_set

    def apply_deltas(coordinate_set, Dx, Dy, Dz):
        if Dx != 0:
            coordinate_set[:, 0] += Dx
        if Dy != 0:
            coordinate_set[:, 1] += Dy
        if Dz != 0:
            coordinate_set[:, 2] += Dz

    if not isinstance(source, list):
        points = np.genfromtxt(source, delimiter = ",")
        apply_deltas(points, Dx, Dy, Dz)
        points = make_circular(points, 72)
    else:
        for i in range(len(source)):
            if i == 0:
                points = np.genfromtxt(source[0], delimiter = ",")
                apply_deltas(points, Dx, Dy, Dz)
                points = make_circular(points, 72)
            else:
                new_points = np.genfromtxt(source[i], delimiter = ",")
                apply_deltas(points, Dx, Dy, Dz)
                points = np.append(points, make_circular(new_points, 72), axis = 0)


    return points

#points = make_box_manually()
#print(points)
deltas_mm = [-1084, 0, 0]
deltas_m = [d  / 1000 for d in  deltas_mm]
points_inner = make_delimiting_surface("inner_wall.dat", *deltas_mm)
points_outer = make_delimiting_surface("outer_wall.dat", *deltas_mm)
#points = alpha_shape_3D(points)
#with open('test.npy', 'wb') as f:
#        np.save(f, points)


#fig = pyplot.figure()
#ax = fig.add_subplot(111, projection='3d')
#x = [hull.points[i, 0] for i in range(len(points))]
#y = [hull.points[i, 1] for i in range(len(points))]
#z = [hull.points[i, 2] for i in range(len(points))]
#ax.scatter(x, y, z)
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#pyplot.show()
#exit()

hull_inner = ConvexHull(points_inner)
hull_outer = ConvexHull(points_outer)
# Delaunay triangulation for the hull
dhull_inner = Delaunay(hull_inner.points)
dhull_outer = Delaunay(hull_outer.points)

parallel = True

t0 = time.time()

if not parallel:
    deletion_list = []

    for p in range(len(my_mesh.data)):
        if ((dhull_outer.find_simplex(my_mesh.points[p][0:3]) == -1 or
            dhull_outer.find_simplex(my_mesh.points[p][3:6]) == -1 or
            dhull_outer.find_simplex(my_mesh.points[p][6:10]) == -1) and
            (dhull_inner.find_simplex(my_mesh.points[p][0:3]) != -1 or
            dhull_inner.find_simplex(my_mesh.points[p][3:6]) != -1 or
            dhull_inner.find_simplex(my_mesh.points[p][6:10]) != -1)):

            deletion_list.append(p)
else:
    def find_bad_triangles(mesh_data):
        #global my_mesh
        global dhull_outer
        global dhull_inner
        if ((dhull_outer.find_simplex(mesh_data[1][0]) == -1 or
            dhull_outer.find_simplex(mesh_data[1][1]) == -1 or
            dhull_outer.find_simplex(mesh_data[1][2]) == -1) or
            (dhull_inner.find_simplex(mesh_data[1][0]) != -1 or
            dhull_inner.find_simplex(mesh_data[1][1]) != -1 or
            dhull_inner.find_simplex(mesh_data[1][2]) != -1)):
            return True

    pool = Pool()

    # Map the function to the STL's triangles and execute in parallel
    deletion_list = pool.map(find_bad_triangles, my_mesh.data)
    deletion_list = [i for i in range(len(deletion_list)) if deletion_list[i] == True]

print("Elapsed time: ", time.time() - t0)
print("Length of deletion list:", len(deletion_list))
print("N triangles of the original mesh ", len(my_mesh.data))
#exit()

new_mesh = mesh.Mesh(np.delete(my_mesh.data, deletion_list))

print("N triangles of the new mesh", len(new_mesh.data))

for t in new_mesh.data:
    t[1][0][0] /= 1000
    t[1][0][0] += deltas_m[0]
    t[1][0][1] /= 1000
    t[1][0][1] += deltas_m[1]
    t[1][0][2] /= 1000
    t[1][0][2] += deltas_m[2]
    t[1][1][0] /= 1000
    t[1][1][0] += deltas_m[0]
    t[1][1][1] /= 1000
    t[1][1][1] += deltas_m[1]
    t[1][1][2] /= 1000
    t[1][1][2] += deltas_m[2]
    t[1][2][0] /= 1000
    t[1][2][0] += deltas_m[0]
    t[1][2][1] /= 1000
    t[1][2][1] += deltas_m[1]
    t[1][2][2] /= 1000
    t[1][2][2] += deltas_m[2]


#exit()

# These are the triplets of coordinates
#print(my_mesh.vectors)

#for i in range(len(my_mesh.vectors)):
#    np.savetxt("data.csv", my_mesh.vectors[i], delimiter=",")


#print(np.shape(new_mesh.points))
#exit()
#np.savetxt("data.csv", new_mesh.points.reshape((np.shape(new_mesh.points)[0]*3, 3)), delimiter=",")

#with open('test.txt', 'w') as outfile:
#    # I'm writing a header here just for the sake of readability
#    # Any line starting with "#" will be ignored by numpy.loadtxt
#    outfile.write('# Array shape: {0}\n'.format(data.shape))
#    # Iterating through a ndimensional array produces slices along
#    # the last axis. This is equivalent to data[i,:,:] in this case
#    for data_slice in data:
#        # The formatting string indicates that I'm writing out
#        # the values in left-justified columns 7 characters in width
#        # with 2 decimal places.  
#        np.savetxt(outfile, data_slice, fmt='%-7.2f')
#        # Writing out a break to indicate different slices...
#        outfile.write('# New slice\n')

# These are the coordinates of the triplets of points (triangles)
#print(my_mesh.points)

#exit()

new_mesh.save("pre-diffuser.stl")

#plot_3D(new_mesh.vectors)

