from pyrite.datatypes import Node, ElementLight

from scipy.spatial import Delaunay
import numpy as np
from itertools import combinations

def try_float(string: str):

    if string.strip() == "null":
        return None 
    else:
        return float(string)

nodes = []


with open("nodes.csv", 'r') as f:

    header_skipped = False

    for line in f.readlines():
        if not header_skipped:
            header_skipped = True 
            continue

        items = line.split(",")

        try:
            id = int(items[0])
            x = float(items[1])
            y = float(items[2])
            ux = try_float(items[3])
            uy = try_float(items[4])
            fx = try_float(items[5])
            fy = try_float(items[6])
        except Exception:
            print("Error loading nodes from nodes.csv!")
            raise

        nodes.append(
            Node(
                x,
                y,
                ux,
                uy,
                fx,
                fy,
            )
        )


points = np.empty((len(nodes), 2))

for i, node in enumerate(nodes):
    points[i] = [node.x, node.y]

print(points)

tri = Delaunay(points)


elements = []

print("simplices:", tri.simplices)

    
for simp in tri.simplices:
    pairs = list(combinations(simp,2)) 
    
    for pair in pairs:

        n1 = nodes[pair[0]]
        n2 = nodes[pair[1]]
        element = ElementLight(n1, n2)

        if element not in elements:
            elements.append(element)
        else:
            print("duplicate")
    
print(elements)


import matplotlib.pyplot as plt
plt.triplot(points[:,0], points[:,1], tri.simplices)
plt.plot(points[:,0], points[:,1], 'o')
plt.show()