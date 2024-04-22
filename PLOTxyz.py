import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_point_cloud(file_path):
    points = []
    with open(file_path, "r") as file:
        num_points = int(file.readline())
        for _ in range(num_points):
            x, y, z = map(float, file.readline().split())
            points.append((x, y, z))

    x_coords, y_coords, z_coords = zip(*points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_coords, y_coords, z_coords)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
    else:
        file_path = sys.argv[1]
        plot_point_cloud(file_path)
