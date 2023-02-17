import numpy as np
from dhart.geometry import LoadOBJ, CommonRotations
from dhart.raytracer import EmbreeBVH
from dhart.graphgenerator import GenerateGraph
#from dhart.graph import GetNodeID, getNodes, ConvertToList
import dhart.spatialstructures.graph
import dhart
import time
import math
import itertools
from queue import PriorityQueue
from quad_mesh_simplify import simplify_mesh
from graph import Graph
import ctypes
import random
import csv



'''
###
Run in cmd-prompt from venv with: python  ..\..\OneDrive\Desktop\GAIDG\graphTraversa.py
###
When importing into blender, ensure that the import properties are: axis forward -Y; Up Z
###
'''
#REDUCE THE NUMBER OF NODES IN THE MESH BY THIS PERCENTAGE (0,1)
#IF THERE ARE 100 NODES IN THE MESH A PERCENTAGE VALUE OF 0.9 HAS A TARGET REDUCTION OF 10 NODES FOR 90 REMAINING
percentage = 0.9


#obj_path = dhart.get_sample_model("Weston_Analysis_z-up.obj")
obj_path = 'C:/Users/tjric/OneDrive/Desktop/GAIDG/moon.obj'

obj = LoadOBJ(obj_path, rotation=CommonRotations.Yup_to_Zup)
#obj = LoadOBJ(obj_path)
bvh = EmbreeBVH(obj)

#start_point = (-0, -20, 20) #energyblob
start_point = (-0.0226852,-11.3705,0.126856) #moon
#start_point = (-1, -6, 1623.976928)#plane
#start_point = (1.25229, 1.11708, 0.117656)
#start_point = (-6.660584, -8.460370, -0.010891)#columns
#start_point = (886.644531, 218.166840, 287.826385) #wesitin
#start_point = (964.497070,171.959366,425.407349) #westing half stair case

spacing = (10,10,10) #coverage of weston stairs
max_nodes = 500

up_step, down_step = 90,90 
up_slope, down_slope = 90,90

#up_step, down_step = 20,20
#up_slope, down_slope = 20,20
max_step_connections = 1
min_connections = 7


graph = GenerateGraph(bvh, start_point, spacing, max_nodes,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)

#graph = GenerateGraph(bvh, start_point, spacing, max_nodes, cores=-1)
print(graph)

csr_graph = graph.CompressToCSR()
#print(csr_graph)
nodes = graph.getNodes()



'''
Transfroms the csr matrix from a class object into a list object which will be used subsequnetly. 
The function will parse the csrMartix class and transform into a useable list using scipy nonzero() to 
determine the locations of non zero entries via the rows and columns. A list of tuples with the format 
(int, int, float)] where the ints are the node ID's which are connected and the float is their distance
apart is returned. 
'''
def createGraphStruct(matrix):
    rows, cols = matrix.nonzero()
    temp = [(int(i), int(j), float(matrix[i,j])) for i, j in zip(*matrix.nonzero())]
    return temp

def distanceList(csrList):
    distances = {}
    for item in csrList:
        #print(item[1])
        distances[(item[0], item[1])] = item[2]
    return distances

'''
Create a dictionary which contains the edge connections for each node. Since
the graph is a digraph, only outgoing edges are considered. 
'''
def createGraphDict(gStruct):
    graphDict = {}
    for i in gStruct:
        if (i[0] in graphDict):
            graphDict[i[0]].append(i[1])
        else:
            graphDict[i[0]] = [i[1]]
    for j in gStruct:
        if(j[1] not in graphDict):
            graphDict[j[1]] = []
        else:
            continue
    return graphDict


def structFromCSR(csrList):
    struct = {}
    for i in csrList:
        if (i[0] in struct):
            struct[i[0]].append(i[1])
        else:
            struct[i[0]] = [i[1]]
    for j in csrList:
        if(j[1] not in struct):
            struct[j[1]] = [j[0]]
        else:
            if(j[0] in struct[j[1]]):
                continue
            struct[j[1]].append(j[0])
    return struct

#-----------------------------------------------QUADMESH IMPORT FUNCTIONS-------------------------
#-------------------------------------------------------------------------------------------------

'''
Create the np.array of faces
'''
def quadMeshCreator(faces, coords):
    temp = np.array([[faces[0][0], faces[0][1], faces[0][2]]], dtype=np.uint32)
    for i in range(1,len(faces)):
        temp2 = np.array([[faces[i][0], faces[i][1], faces[i][2]]], dtype=np.uint32)
        temp = np.append(temp,temp2, axis=0)
    return temp

'''
create an array of poistions with form [[x,y,z], ...] that corresponds to the order in which the vertices appear 
in the faces
'''
def quadMeshPositions(faces, coords):
    pos = np.array([[float(coords[0][0]), float(coords[0][1]), float(coords[0][2])]])
    for k,v in coords.items():
        if(k == 0):
            continue
        else:
            pos = np.append(pos, [[float(v[0]), float(v[1]), float(v[2])]], axis = 0)
    return pos


#translate the array of new faces into a tuple with the same faces for use in the remainder of the script
def convertArrToList(inptArr):
    temp = []
    for item in inptArr:
        temp2 = (item[0], item[1], item[2])
        temp.append(temp2)
    return temp

#translate the new coordinates to a coordinate dictionary for use later 
def convertCoordDict(inptArr):
    temp = {}
    count = 0
    for item in inptArr:
        #print('item: ', item)
        temp[count] = (item[0], item[1], item[2])
        count +=1
    return temp

'''
After decimation the new faces of the mesh must be used to create a new adjacecny structure so that A* is 
not run on an older structure. 
'''
def createNewAdjDict(newFaces):
    resultDict = {}
    for item in newFaces:
        for i in item:
            temp = [j for j in item if(j!=i)]
            if(i in resultDict):
                resultDict[i].append(temp[0])
                resultDict[i].append(temp[1])
            else:
                resultDict[i] = temp
    for k,v in resultDict.items():
        temp = set(v)
        resultDict[k] = list(temp)
    return resultDict
                

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
'''
Create a dictionary of where the nodeID is the key and the value is a tuple of (x,,y,z) coordinates.
This function will be used for normal vector computation and other follow on functions
'''
def createCoordDict(nodes):
    dict1 = {}
    for item in nodes:
        dict1[item[4]] = (item[0], item[1], item[2])
    return dict1

'''
Get the angle between two vectors. This function has the intention of checking if two vectors are 
parallel. This will aide in triangle decimation. Returns Theta (the angle between the vectors). 
Theta can also be used as a threshold value so if they are within that threshold, the triangle can be 
decimated. 
'''
def computeVectorAngle(vec1, vec2):
    dot = np.dot(vec1, vec2)
    mag1 = np.sqrt(vec1.dot(vec1))
    mag2 = np.sqrt(vec2.dot(vec2))
    theta = math.acos(dot/(mag1*mag2))
    return theta

'''
compute the normal vectors of each triangle as a group or individually 
compute the nomral vector of a singular triangle face
'''
def computeNormal(faces, nodeCoord):
    norms = {}
    for item in faces:
        vec1 = np.array(list(nodeCoord[item[0]]))
        vec2 = np.array(list(nodeCoord[item[1]]))
        vec3 = np.array(list(nodeCoord[item[2]]))
        vecA = vec2 - vec1
        vecB = vec1 - vec3
        normX = (vecA[1] * vecB[2]) - (vecA[2] * vecB[1])
        normY = (vecA[2] * vecB[0]) - (vecA[0] * vecB[2])
        normZ = (vecA[0] * vecB[1]) - (vecA[1] * vecB[0])
        normal = np.array(list((normX, normY, normZ)))
        norms[item] = normal
    return norms

def computeSingleNormal(face, nodeCoord):
    vec1 = np.array(list(nodeCoord[face[0]]))
    vec2 = np.array(list(nodeCoord[face[1]]))
    vec3 = np.array(list(nodeCoord[face[2]]))
    vecA = vec2 - vec1
    vecB = vec1 - vec3
    normX = (vecA[1] * vecB[2]) - (vecA[2] * vecB[1])
    normY = (vecA[2] * vecB[0]) - (vecA[0] * vecB[2])
    normZ = (vecA[0] * vecB[1]) - (vecA[1] * vecB[0])
    normal = np.array(list((normX, normY, normZ)))
    return normal


'''
Helper function for getSharedFaces()
Determines if two tuples share two like components. Indicates a shared edge between two
triangle faces in the mesh.
'''
def checkLikeComponents(tup1, tup2):
    a = set(tup1)
    b = set(tup2)
    if(len(a.intersection(b)) == 2):
        return True
    else:
        return False

'''
Determine which faces are connected to one another. Will aide in decimation 
as triangles which are connected are subject to decimation provided they meet certain criteria. 
returns a dictionary of shared edge faces
'''
def getSharedFaces(faces):
    sh_faces = {}
    for i in range(len(faces)):
        temp = []
        for j in range( len(faces)):
            if(checkLikeComponents(faces[i], faces[j])):
                temp.append(faces[j])
            else:
                continue 
        sh_faces[faces[i]] = temp
    #[item for item in range(len(faces)) for otherItem in range(len(faces)) if checkLikeComponents(faces[item], faces[otherItem])]
    return sh_faces

'''
Helper function for fastRemoveOverlap(). Determines which two nodes have the largest separation of the four inputs. This 
distance is the hypotenuse and is used to determine which triangles should be removed to maintain a continuous mesh. Additional 
description in triangleComplement().
'''
def getLongestDist(verts,dists):
    combs = list(itertools.combinations(verts, 2))
    longest = 0
    longestKey =(0,0)
    for item in combs:
        try:
            if(dists[item] > longest):
                longest = dists[item]
                longestKey = item
            else: 
                continue
        except KeyError:
            continue
    return longestKey


'''
Helper functions for FastRemoveOverlap(). 
Given 2 triangle tuples, determine if they share a hypotenuse.
       2  ____  0
         |\  /|
         | \/ |
         | /\ |
       1 ------ 4

for traingle removal, since set operations on traingles which share four vertices, do not yeild definitve results,so
it must be determined for given a triangle tuple, find the other tuple among the remaing three which do not contain the 
unique element found in triangle tuple 1.
Example:

(4, 1, 2)
(4, 0, 1)
(4, 0, 2)
(1, 0, 2)

All four of these triangles represent the four triangles that can be made out of four vertices (0,1,2,4)

if we are using triangle (4, 1, 2) as the comparison triangle, check the other three for if they contain
the value 4. The one that does not, (1, 0, 2) is the complimentary triangle, so (4,0,1) and (4,0,2) should 
be selected for removal. 

Find the largrst distance between two vertices, this will be the hypotenuse. Then compare the raminaing two vertices 
for equality. In the example, the loargest distance is either (2,4) or (1,0). The only sets of vertices that are equal 
are: (1,0) of traingle (2,1,0) and triangle (4,1,0) or (2,4) of triangle (2,0,4) and triangle (1,2,4).
'''
def triangleComplement(lot, distances):
    verts = list(set(set(lot[0]).union(set(lot[1]))))
    marked = []
    dist = getLongestDist(verts, distances)
    for item in lot:
        if(dist[0] in item and dist[1] in item):
            marked.append(item)
    return marked

'''
FastRemoval helper. Determines if the combination of input nodes are connected.
'''
def connected(n1, n2, n3, structure):
    if((n1 in structure[n2]) and (n2 in structure[n1]) and (n1 in structure[n3]) and (n3 in structure[n1]) and (n2 in structure[n3]) and (n3 in structure[n2])):
        return True
    return False

'''
Fast removal helper
'''
def trace(key, vals, structure):
    fullConnections = []
    combs = list(itertools.combinations(vals,3))
    #print(combs)
    for item in combs:
        if(connected(item[0], item[1], item[2], structure)):
            fullConnections.append((key, item[0], item[1], item[2]))
        else:
            continue
    return fullConnections


'''
Faster Removal of overlapping triangles which share four vertices. The function creates a copy of the triangles present to 
prevent mutating the input list of faces. Using the graph structure (which vertices are adjacent to each other in a dictionary
form), check the connections of each k,v pair. If there are less than 3 connections, then a triangle is not formed, so there
is no need to check further for overlaps.This reduces the number of later operations performed. 

The connections that are greater than 2 are then 'traced' using trace() and 'connected() to get the 'squares' of 4 vertices that 
are fully connected (i.e. degree of v == 3). Becuase it is possible that the ordering of these squares is represented more than
one time, the vertexSquares are purged of duplicated before further processing. 

For each square, the 4 traingles of that square are found, and removed using the methods described in triangleComplement().

Finally, for each triangle 'marked' for removal in triangleComplement(), those triangles are removed from the copied list. The 
copied list is then returned.
'''
def fastOverlapRemoval(gstructure, faces, distances, incidence):
    print(len(gstructure))
    vertSquares = []
    noOverlaps = faces.copy()
    for k,v in gstructure.items():
        if(len(v) > 2):
            temp = trace(k, v, gstructure)
            if(len(temp) == 0):
                continue
            elif(len(temp) == 1):
                vertSquares.append(temp[0])
            else:
                for i in temp:
                    vertSquares.append(i)
                print('Type1 overlap Section 1:{:.3f}%'.format((1/(len(gstructure) * len(temp)))*100), end='\r')
        else:
            continue
    noDup = list(set(map(tuple, map(sorted, vertSquares))))
    count = 0
    temp2 = []
    for val in noDup:
        i=val[0]
        j=val[1]
        k=val[2]
        l=val[3]
        iFaces = incidence[i]
        jFaces = incidence[j]
        kFaces = incidence[k]
        lFaces = incidence[l]
        sharedFaces = iFaces + jFaces + kFaces + lFaces
        sortedTup = [tuple(sorted(val)) for val in sharedFaces]
        temp2 = list(set(sortedTup))
        temp = [item for item in temp2 if(i in item and j in item and k in item) or (i in item and k in item and l in item) or (i in item and j in item and l in item) or (j in item and k in item and l in item)]
        rmv = triangleComplement(temp, distances)
        #for each triangle selected to removal, remove that face from the copied list. 
        for tup in rmv:
            try:
                noOverlaps.remove(tup)
            except ValueError:
                continue
        count+=1
        print('Type1 overlap section 2: {:.3f}%'.format((count/len(noDup))*100), end='\r')
    return noOverlaps


#FUCTIONS FOR EDGEOVERLAP()
'''
Get the partial connected triangles. vertices that only have three components are considered tri connected.
pruning happens later.
'''
def getPartial(key, vals, struct):
    connections = []
    if(len(vals)< 3):
        vals.append(key)
        combs = list(itertools.combinations(vals, 3))
        return combs
    else:
        vals.append(key)
        combs = list(itertools.combinations(vals,3))
        return combs

'''
get the coordinates of the new vertex at the intersection of two lines n1n2, and n3n4. 
returns coordninates. 
'''
def lineIntersection(n1, n2, n3, n4, coords):
    if(len(set([n1,n2,n3,n4])) < 4):
        return 0
    else:
        v1 = np.array([coords[n1][0], coords[n1][1], coords[n1][2]])
        v2 = np.array([coords[n2][0], coords[n2][1], coords[n2][2]])
        v3 = np.array([coords[n3][0], coords[n3][1], coords[n3][2]])
        v4 = np.array([coords[n4][0], coords[n4][1], coords[n4][2]])
        a = np.array([[(v2[0] - v1[0]), (v4[0] - v3[0])],
                      [(v2[1] - v1[1]), (v4[1] - v3[1])],
                      [(v2[2] - v1[2]), (v4[2] - v3[2])]])
        b = np.array([[v3[0] - v1[0]],
                      [v3[1] - v1[1]],
                      [v3[2] - v1[2]]])
        x = np.linalg.lstsq(a,b,rcond=-1)
        t = x[0][0,0]
        intersetion_x = (v2[0]-v1[0])*t + v1[0]
        intersetion_y = (v2[1]-v1[1])*t + v1[1]
        intersetion_z = (v2[2]-v1[2])*t + v1[2]
    return (intersetion_x, intersetion_y, intersetion_z)

'''
Get the four vertices that the partial overlaps compose.
'''
def getSquares(tris, struct):
    squares = []
    for i in range(len(tris)):
        temp = []
        for j in range(1, len(tris)):
            if(tris[i] == tris[j]):
                continue
            elif(len(set(tris[i]).intersection(set(tris[j]))) == 2):
                squares.append(tuple(set(tris[i]).union(set(tris[j]))))
    noDup = [item for item in (set(tuple(k) for k in squares))]
    return noDup

'''
Determines if the tuple given is in fact a triangle. 
'''
def triangle(tup, struct):
    if((tup[0] in struct[tup[1]]) and (tup[0] in struct[tup[2]]) and (tup[1] in struct[tup[0]]) and (tup[1] in struct[tup[2]]) and (tup[2] in struct[tup[0]]) and (tup[2] in struct[tup[1]])):
        return True
    else:
        return False


'''
Function aims to get the two longest edges 'the hypotenuses' of the four vertices which have been deemed to be 
the four containing the two partially overlappign triangles. Returns a tuple containing both vertices connected 
by the two longest edges of the square. 
'''
def getHypotenuses(nodes, distances, coords):
    combs = list(itertools.combinations(nodes, 2))
    hyps = []
    for item in combs:
        try:
            hyps.append((item, distances[item]))
        except KeyError:
            try:
                temp = (item[1], item[0])
                hyps.append((temp, distances[temp]))
            except KeyError:
                continue
    longest = 0
    longestIndx = 0
    longest2 = 0
    longest2Indx = 0
    for item in hyps:
        if(item[1] > longest):
            longest = item[1]
            longestIndx = hyps.index(item)
        elif(longest>= item[1] and item[1]> longest2):
            longest2 = item[1]
            longest2Indx = hyps.index(item)
    return ((hyps[longestIndx][0][0], hyps[longestIndx][0][1]),(hyps[longest2Indx][0][0], hyps[longest2Indx][0][1]))

#update the coordinate dictionary to include the new found coordinates
def updateCoords(newCoords, coordDict):
    coordDict[len(coordDict)] = newCoords
    return coordDict

'''
Update two input faces to be the three new faces which will be returned and added to the triangle face list for 
.obj writing
'''
def updateFaces(nodes, newNode):
    temp = list(set(nodes[2]).intersection(set(nodes[3])))
    new_face1 = (temp[0], newNode, temp[1])
    temp4 = [nodes[0][0], nodes[1][1]]
    temp5 = [nodes[0][1], nodes[1][0]]
    new_face2 = (temp4[0], newNode, temp4[1])
    new_face3 = (temp5[0], newNode, temp5[1])
    return [new_face1, new_face2, new_face3]

'''
Find out if the connecting edge of two triangles is the hypotenuse or not. Explained further 
in the function edgeOverlap().
'''
def findNonPartial(knownEdge, tri1, tri2, distances):
    e1 = itertools.combinations(tri1, 2)
    e2 = itertools.combinations(tri2, 2)
    edges1 = []
    edges2 = []
    for i in e1:
        try:
            edges1.append((i,distances[i]))
        except KeyError:
            edges1.append(((i[1], i[0]), distances[(i[1], i[0])]))
            continue
    for j in e2:
        try:
            edges2.append((j, distances[j]))
        except KeyError:
            edges2.append(((j[1], j[0]), distances[(j[1], j[0])]))
            continue
    allEdges = edges1 + edges2
    longestE1 = tuple()
    longestE2 = tuple()
    longest = 0
    longest2 = 0
    for item in allEdges:
        if(item[1] > longest):
            longest = item[1]
            longestE1 = item[0]
        elif(item[1] <= longest and item[1] > longest2):
            longest2 = item[1]
            longestE2 = item[0]
        else:
            continue
    if(len(set(knownEdge).intersection(set(longestE1).intersection(set(longestE2)))) == 2):
        return True
    else:
        return False

'''
Remove overlap edges.

1_____2
 \   /|
   \/ |
   / \|
3 -----4

remove the overlap created by the above example, by creating a vertex at the intersection of 1,4 and 3,2. 
The resulting vertex creates three triangles (1,2,newV), (2,newV,4), and (3,newV,4). These are added to the
faces of the mesh and (1,2,4), (3,2,4) are removed. The new vertex is also added to the vertex list for 
.obj writing.
returns the new coordinate dictionary and the updated list of trinagle faces
'''
def edgeOverlap(gStruct, distances, faces, coords):
    startT = time.time()
    #Gets all partial overlaps 
    tris = []
    noOverlaps = faces.copy()
    print(len(coords))
    for k,v in gStruct.items():
        if(len(v) > 1 and len(v) < 4):
            #print(k,v)
            temp = getPartial(k,v, gStruct)
            #print(temp)
            for i in temp:
                if(i in faces and  triangle(i, gStruct)):
                    tris.append(i)
                    print('Type 2 Section 1: {:.3f}'.format((1/(len(gStruct)*len(temp)))*100), end='\r')
                else:
                    continue
        else:
            continue

    #Due to how the mesh can be set up depending on the number of nodes, connections, etc. it is poosible
    #for non-partially overlapping trinagles to make it into the list of partial overlaps. This section removes
    #those non-partials, by checking if the only connecting edge between two triangles is the hypotenuse, or the 
    #longest edge of the two triangles. By the nature of the partially overlappign triangles, there will be 2 
    #hypotenuses, which are not made up of the same vertices. In the case that there is only one, those two triangles
    #are removed from the list of partial triangles to be subidvied further.
    doubles = [(item, item2) for item in tris for item2 in tris if(len(set(item).intersection(set(item2))) == 2)]
    for item in doubles:
        if(findNonPartial(tuple(set(item[0]).intersection(set(item[1]))), item[0], item[1], distances)):
            try:
                tris.remove(item[0])
                tris.remove(item[1])
                print('Type 2 Section 2: {:.3f}'.format((1/len(doubles))*100), end='\r')
            except ValueError:
                continue
        else:
            continue
    squares = getSquares(tris, gStruct) 
    hypotenuses = []
    #get the hypotenuses which are part of the partial overlap, if there are 3 vertices in this instead of 
    #4, do not add that set of hypotenuses for further analysis. (i.e. if the hypotenuse makes this < shape 
    #rather than X this shape, do not add)
    for item in squares:
        hypots = getHypotenuses(item, distances, coords)
        if(len(set([hypots[0][0], hypots[0][1], hypots[1][0], hypots[1][1]])) < 4):
            continue
        else:
            hypotenuses.append(hypots)
            print('Type 2 Section 3: {:.3f}'.format((1/len(squares))*100), end='\r')
    #This section is to get the faces which are associated with the hypotenuses found. The aim is to find the 
    #faces where both ends of the hypotenuse (item[0][0], item[0][1]), are 2 of the 3 face components. The next
    #for loop find the face where both ends of the second hypotenuse (item[1][0], item[1][1]) are 2 of the 3 face 
    #components.
    faceHyp1 = []
    for item in hypotenuses:
        for item2 in faces:
            if(((item[0][0] in item2) and (item[0][1] in item2) and (len(set(item2).intersection(set(item[1]))) == 1))):
                faceHyp1.append((item, item2))
                print('Type 2 Section 4: {:.3f}'.format((1/(len(hypotenuses) * len(faces)))*100), end='\r')
    faceHyp2 =[]
    for item in faceHyp1:
        for item2 in faces:
            if(((item[0][1][0] in item2) and (item[0][1][1] in item2) and (item2 not in item) and (len(set(item2).intersection(set(item[0][0]))) == 1))): 
                faceHyp2.append((item[0][0], item[0][1], item[1], item2))
                print('Type 2 Section 5: {:.3f}'.format((1/(len(faceHyp1) * len(faces)))*100), end='\r')
    #-------------------
    #This section of the function now deals with the pruned triangles and squares. For each of the found hypotenuses, their corresponding
    #squares and partially overlappign triangles, this section gets the intersection of the two hypotenuses, updates the coordinate 
    #dictionary with the new coodirnate ID and location, and also updates the faces of the square to include the new vertex. The old faces
    #are deleted from the original triangle face list. The new coordinate dictionary and the updated faces are returned for further use in 
    #creating the object file. 
    new_cords = coords.copy()
    fRemovals =[]
    fAdditions =[]
    for item in faceHyp2:
        temp = lineIntersection(item[0][0], item[0][1], item[1][0], item[1][1], coords)
        if(temp == 0):
            continue
        else:
            new_cords = updateCoords(temp,new_cords)
            newFaces = updateFaces(item, len(new_cords)-1)
            fRemovals.append(item[2])
            fRemovals.append(item[3])
            for i in newFaces:
                fAdditions.append(i)
                print('Type 2 Section 6: {:.3f}'.format((1/(len(faceHyp2) * len(newFaces)))*100), end='\r')
    for item in fRemovals:
        if(item not in noOverlaps):
            continue
            print('Type 2 Section 7: {:.3f}'.format((1/len(fRemovals)*100)), end='\r')
        noOverlaps.remove(item)
    for item in fAdditions:
        noOverlaps.append(item)
        print('Type 2 Section 8: {:.3f}'.format((1/len(fAdditions))*100), end='\r')
    return (new_cords, noOverlaps)

'''
------------------------------------------------------DECIMATION AND DECIMATION FUNCTIONS--------------------------------------------------------------
------------------------------------------------------NOT CURRENTLY IN USE-----------------------------------------------------------------------------
decimation is influenced by the decimation heuristics and methods outlined in
Hussain M, Okada Y, Niijima K. Efficient and feature-preserving triangular mesh decimation. Journal of WSCG, 12(1): 2004, pp.167–174.
and adapted for use on the mesh created using Dhart. 
'''
'''
Computes the area of a given trangle face. Retruns the area. 
'''
def computeArea(face, coords):
    a = np.array([coords[face[0]][0], coords[face[0]][1], coords[face[0]][2]])
    b = np.array([coords[face[1]][0], coords[face[1]][1], coords[face[1]][2]])
    c = np.array([coords[face[2]][0], coords[face[2]][1], coords[face[2]][2]])
    line1 = np.subtract(b,a)
    line2 = np.subtract(c,a)
    cross = np.cross(line1, line2)
    mag = np.linalg.norm(cross)
    mag = mag * 0.5
    return mag 

'''
Creates a dictionary for use in calculating decimation heuristics. Given the coordinates and faces of of the input mesh
the the normal vector of each triangle face and the area of each triangle are determined. Returns dictionary with the 
format: {face:[normal, area]}
'''
def createDecimationDict(coords, faces):
    dict1 = {}
    for item in faces:
        normal = computeSingleNormal(item, coords)
        area = computeArea(item, coords)
        dict1[item] = [normal, area]
    return dict1

'''
For the computation of the vertex importance heuristic, the function must have all triangles incident on each vertex. 
This is stored as a dictionary with the fromat {ID: [face1, face2,...]}
'''
def getTriangleIncidence(coords, faces):
    dict1 = {}
    for k,v in coords.items():
        temp = [item for item in faces if(k in item)]
        dict1[k] = temp
    return dict1

'''
Visual Importance is the heuristic by which vertices are judged for decimation on a scale of 0<=w<=1 where w is the weight
each vertex recives based on the surrounding topography of the mesh. The hueristic utelizes the incident triangles to a vertex 
and the normals, area of each of those triangles. w is calculated as 1 - ((SUM ^n) / (SUM ^)) where ^ = area of incident triangles 
and n = nromal of incident triangle. If w == 0 then the surface is flat and will be a primary candidate for decimation. Values of w
which have a value closer to 1 add more to the topography of the mesh and should be preserved. Returns a dictionary of the vertex
visual importance: {ID: w}
'''
def visualImportance(normArea, incidence):
    #print(normArea)
    viDict = {}
    for k,v in incidence.items():
        areaSum = 0
        normAreaSum = 0
        for item in v:
            normAreaSum += normArea[item][0] * normArea[item][1]
            #print(normArea[item][0], normArea[item][1])
            areaSum += normArea[item][1]
        try:
            temp = normAreaSum / areaSum
            mag = np.linalg.norm(temp)
            w = 1-mag
            if(w < 0):
                viDict[k] = 0.0
                continue
            else:
                viDict[k] = w
                continue
        except ZeroDivisionError:
            continue
    return viDict

'''
To ensure that the mesh retains its original shape, the borders must be preserved. This function returns all of the borders of the 
input mesh, such that they may be used to determine if an edge selcted for half edge collapse is in fact a border. The borders of 
the mesh will only be found in a single triangle face of the mesh, so the function finds those edges, since they are not known prior. 
'''
def getBorderEdges(faces):
    edges = []
    for item in faces:
        possibleEdges = itertools.combinations(item, 2)
        #temp = [i for i in possibleEdges if((i[0], i[1]) not in temp and (i[1], i[0]) not in temp)]
        temp = []
        for k in possibleEdges:
            if((k[0], k[1]) not in temp and (k[1], k[0]) not in temp):
                temp.append(k)
        for j in temp:
            edges.append(j)
    borders = []
    for item in edges:
        count = 0
        for item2 in faces:
            if(len(set(item).intersection(set(item2))) == 2):
                count +=1 
        if(count == 1):
            borders.append(item)
    return borders

'''
To ensure that all half edges are chedked fro their cost associated with deciimation, all half edges of the mesh must be identified.
A hald edge is more or less a directed edge between vertex u,v where for example the if the origin is u then the head is v. Teh reverse
is also true where the origin v would have head u and the two edges make the undirected edge u,v. One important component of the decimation
is the preservation of borders, so all half edges that have head and origin on a border edge are ignored. Additionally, edges which have an 
origin as a border veterx will not be used for edge collapse, but border edges that have an origin as a the half edge vertex will be condsidered
since the half edge collapse will not alter the shape in this scenario (i.e. the head is collapsed into the origin). Retruns a dictionary with the half edges of the mesh with the form:
{ origin: [(origin, head), ...]}
'''
def getHalfEdges(gStruct, borders):
    halfE = []
    for k,v in gStruct.items():
        if(len(v) >= 2):
            for item in v:
                if((k,item) in borders or (item, k) in borders):
                    continue
                else:
                    origin = k
                    head = item
                    halfE.append((origin, head))
        else:
            continue
    order = []
    for item in halfE:
        if(item[0] not in order):
            order.append(item[0])
    hEdges = {}
    for item in order:
        temp = [item2 for item2 in halfE if(item2[0] == item)]
        hEdges[item] = temp
    return hEdges

'''
Helper for calcQ(). Swaps the head and origin vertices in the triangle for gerometric error calculation
'''
def swap(origin, head, tri):
    temp = list(tri)
    temp.remove(head)
    temp.append(origin)
    return temp

'''
calculates the Q value described in the gerometricError() description.
'''
def calcQ(origin, head, tri, coords):
    #the 'new' triangle
    tPrime = swap(origin, head, tri)
    tArea = computeArea(tri, coords)
    tPrimeArea = computeArea(tPrime, coords)
    l = 0.5 * (tArea + tPrimeArea)
    tnormal = computeSingleNormal(tri, coords)
    tPrimeNormal = computeSingleNormal(tPrime, coords)
    theta = 1 - (np.dot(tnormal, tPrimeNormal))
    q = l * theta
    return q


'''
The geometric error calculation. Determines the error asscociated with any edge collapse of the mesh. If head v is being collapsed
into origin u, then each vertex which is incident to v will be 'connected' to vertex u. This will provide a set of two triangles, the
original triangles when vertex v exist and the new triangles that would exist if the collapse were to happen. The geometric error 
(cost of collapse) is then calculated using Q = l*Theta where l = 0.5(area of original triangle + area of new triangle) and 
Theta = 1 - (normal of original triangle DOT normal of new triangle). The sum of Q for eaach triangle in the set differnce of triangles 
for the vertices (v - u) returns a dictionary of the form: {halfEdge: cost}
'''
def geometricError(halfEdges, incidence, borders, coords):
    #print('geometric error origin is index 0 head is index 1, collapse head to origin')
    geomError = {}
    #key is the origin, values are the possible heads
    for k,v in halfEdges.items():
        #print(v)
        for item in v:
            inc = incidence[item[1]]
            if(len(inc) > 0):
                origin = item[0]
                head = item[1]
                sumQ = 0
                for val in inc:
                    if(origin in val):
                        continue
                    else:
                        q = calcQ(origin, head, val, coords)
                        sumQ = sumQ + q
                        continue
                geomError[item] = sumQ
            else:
                continue
    return geomError

'''
collapses the selected edge head into the selected head origin if the head is not on a boundary to preserve the mesh boundaries.
returns a list of new faces and a list of the faces which are to be removed from the original face list. 
'''
def collapse(item, coords, faces):
    #print('-------------IN COLLAPSE FUNCTION--------------')
    inc = getTriangleIncidence(coords, faces)
    headInc = inc[item[1][1]]
    #copy to return for removal from the face list
    orig = headInc.copy()
    #triangles which are incident to both origin and head are deleted for free, but must be removed prior to the actual collapse
    temp = [i for i in headInc if(item[1][0] in i)]
    for j in temp:
        headInc.remove(j)
    #print(headInc)
    #swap all head values with origin values in the faces
    new = []
    for k in headInc:
        temp2 = swap(item[1][0], item[1][1], k)
        new.append(tuple(temp2))
    #print(new)
    return (orig, new)

'''
Checks if the head of a selected edge is on the boundary. Returns true of False
'''
def onBoundary(border, head):
    for val in border:
        if(head in val):
            return True
        else:
            continue 
    return False

'''
Gets the optimal half edge collapse cost associated with each vertex adn scales it by the visual importance. retruns a priority 
queue where the lowest cost half edge collaps is the root. Reduction is the percentage of the vertices that should be added to the 
priority queue for mesh recution, defaults to 25 percent. 
'''
def optimizeCollapse(geometricErr, visualImp, coords, faces, borders):
    #for geometricErr (origin, head):val, head is collapsed into the origin, lower val = higher priority on heap queue
    #first find the lowest edge cost collapse for each vertex's neighborhood. c_i = min(c_ij v_j in N_v_i). So the minimun cost 
    # collapse an edge into the current vertex
    #initialize the heapQueue
    hq = PriorityQueue()
    #get the ID's of the vertices which are in the VI dictionary
    idList = list(visualImp.keys())
    for item in idList:
        lowest = math.inf
        lowestItem = None
        for item2 in geometricErr:
            if(item2[1] == item ):
                temp = geometricErr[item2]
                if(temp < lowest):
                    lowest = temp
                    lowestItem = item2
        if(lowest is not None and lowestItem is not None):
            scale = visualImp[lowestItem[1]]
            cost = scale * lowest
            #push the loswest cost edge collapse for the vertex onto the PQ. Format (cost, (origin, head))
            hq.put((cost, (lowestItem[1], lowestItem[0])))
        else:
            continue
    return hq

def decimationFaceUpdate(oldF, newF, faceCopy):
    temp = faceCopy.copy()
    for i in oldF:
        temp.remove(i)
    for j in newF:
        temp.append(j)
    return temp

'''
Compute the new structure of the mesh in graph form fro use in the decmation iteration
'''
def resetStruct(triInc):
    #print(triInc)
    updatedStruct = {}
    for k,v in triInc.items():
        temp = [item2 for item in v for item2 in item if(item2 != k)]
        #print(temp)
        updatedStruct[k] = temp
    return updatedStruct

'''
Complete the reset of the VI, GE and PQ. Return the list of reset obejcts for update in the decimation function
'''
def reset(faces, coords):
    resetDict = {}
    resetDict['decDict'] = createDecimationDict(coords, faces)
    resetDict['triInc'] = getTriangleIncidence(coords, faces)
    resetDict['vi'] = visualImportance(resetDict['decDict'], resetDict['triInc'])
    resetDict['border'] = getBorderEdges(faces)
    resetDict['struct'] = resetStruct(resetDict['triInc'])
    resetDict['half'] = getHalfEdges(resetDict['struct'], resetDict['border'])
    resetDict['gERR'] = geometricError(resetDict['half'], resetDict['triInc'], resetDict['border'], coords)
    resetDict['pq'] = optimizeCollapse(resetDict['gERR'], resetDict['vi'], coords, faces, resetDict['border'])
    return resetDict['pq']


'''
Decimate the triangles. Given an input threshold, the function decimates the percentage fo the PQ vertices. Updates the 
coords, faces, VI, GE, and PQ between iterations. 
'''
def decimate(gError, visualImp, coords, faces, borders, pq, reduction=25):
    faceCopy = faces.copy()
    decLvl = pq.qsize() * (0.01*reduction)
    PriorityQ = pq
    count = 0
    while(count < math.ceil(decLvl)):
        selected = pq.get()
        #check if the head is on boundary. If so, move on to the next.
        if(not onBoundary(borders, selected[1][1])):
            #collapse the edge
            collapsedFaces = collapse(selected, coords, faceCopy)
            origFaces = collapsedFaces[0]
            newFaces = collapsedFaces[1]
            #update the coords, and faces
            #print('old faces')
            #print(faceCopy)
            #print('----------')
            faceUpdate = decimationFaceUpdate(origFaces, newFaces, faceCopy)
            faceCopy = faceUpdate
            #print('new faces')
            #print(faceCopy)
            #print('-----------')
            #resetVars = reset(faceCopy, coords)
            #print('Printing pq objects on pass {}'.format(count))
            #print(resetVars.qsize())
            #for k in range(resetVars.qsize()):
            #    print(resetVars.get())
            count +=1
        else:
            continue
    return faceCopy

#------------------------------------OLD DECIMATION FUNCTIONS END--------------------------------------------

'''
a helper function for checkTriangle(). Intersects two adjacency lists to find the common 
nodes which the 'parent' nodes connect to
'''
def intersect(list1, list2):
    list3 = [item for item in list1 if item in list2]
    return list3


'''
Check if a node in the adjacent list is connected to another node that is connected to the 
'parent' node. i.e. while in the traverse loop, the node nxt is the parent node, check if
the node designated as adjacent is connected to another node which is in the adjacency list for nxt
(use intersection of both node adjacenecy lists).Returns a triangle(s) of the from 
(node1, node2, node3) which composed of the threee adjacent nodes.
'''
def checkTriangle(gDict, pNode, adjNode):
    pNodAdj = gDict[pNode]
    aNodeAdj = gDict[adjNode]
    adjNodes = intersect(pNodAdj, aNodeAdj)
    triangles = []
    for item in adjNodes:
        #triangles.append((int(pNode),int(adjNode), int(item)))
        triangles.append((int(pNode), item, int(adjNode)))
    return triangles

'''
OLD FUNCTION SEE: dupRemoval()
Helper function for removeDuplicate(). Determines if a given set of vertices
has already been added to a list of known triangle vertices. Adds the initial 
set automatically.
'''
def duplicate(verts, dupList):
    temp = set(verts)
    if(not dupList):
        return False
    else:
        for item in dupList:
            if(set(item) == temp):
                return True
            else:
                continue
        return False


'''
OLD FUNCTION SEE: dupRemoval()
This function is meant to remove duplicate entries of the same triangle. Naturally,
this script's faceTriangluation() will create duplicate triangle faces. Example:
entry of vertices (1,2,3) and (3,2,1) is the same triangle.
'''
def removeDuplicate(triangles):
    nonD = []
    for item in triangles:
        if(not duplicate(item, nonD)):
            nonD.append(item)
        else:
            continue
    return nonD

'''
Winding order: Determine the way in which the vertices should be ordered to have a uniform 'facing' of the faces. 
Possible methods are to find the windsing order of the first triangle found then do something with the normals? to 
compare agains subsequent triangles found. 
https://gamedev.stackexchange.com/questions/159379/getting-the-winding-order-of-a-mesh

If the z component of the cross product of two of the edges of a triangle is negative, then the winding order is 
incorrect, assuming the script is to work on real world objects (z -is up). It means that the winding order is counter 
clockwise when veiwing from the positve z axis. The function should compute the cross product of each face it is given and 
then if the z componet is positive, add it back to the list of faces. If not reorder to and check. 
'''
def windingOrder(faces, nodeCoord):
    returnFaces = []
    for item in faces:
        normal = computeSingleNormal(item, nodeCoord)
        if(normal[2] < 0):
            #print('Normal less than 0 {}'.format(normal))
            returnFaces.append(item)
        else:
            temp = (item[2], item[1], item[0])
            normal2 = computeSingleNormal(temp, nodeCoord)
            #print('Normal rearranged {}'.format(normal2))
            returnFaces.append(temp)
    return returnFaces

'''
Traverse the graph using BFS traversal. Setup of visited and queue will ensure that no vertex is 
reached that has already been visted. Starts from the intital node, InitID which is a node id, 
corresponding to the graph generation.
'''
def traverseGraph(gDict, initID):
    incidence = {}
    visited = []
    queue = []
    faces = []
    visited.append(initID)
    queue.append(initID)
    while(queue):
        nxt = queue.pop(0)
        #print('Node to be dict key: ',nxt)
        incident = []
        for adjacent in gDict[nxt]:
            temp = checkTriangle(gDict, nxt, adjacent)
            incident = incident + temp
            for item in temp:
                faces.append(item)
                #print('Graph traversal completed: {:.3f}%'.format((1/(len(gDict)*len(temp)))*100), end='\r')
            if adjacent not in visited:
                visited.append(adjacent)
                queue.append(adjacent)
            else:
                continue
        sortedTup = [tuple(sorted(val)) for val in incident]
        incident2 = list(set(sortedTup))
        #print(incident2)
        incidence[nxt] = incident2
        print('graph traversal {:.3f}%'.format((nxt/len(gDict))*100), end='\r')
    #print('Incidcene in traversal: ',incidence)
    return (visited,faces, incidence)

'''
get the vertex order for writing .obj file
should be in the order of reference 
example:

shapeFaces = [(1,3,2), (2,3,0)]

The vertex here corresponds to a specific set of x,y,z values and since 
.obj files refernce the vertices in the faces to the order of appearance 
in the .obj file, the obj file would have:

(vertexID 1) v 1.0 1.0 0.0
(vertexID 3) v 1.0 2.0 0.0
(vertexID 2) v 2.0 1.0 0.0
(vertexID 0) v 0.0 0.0 0.0

(face 0)     f 1 3 2
(face 1)     f 2 3 0
'''
def getVertexOrder(faces):
    vOrder = []
    for i in faces:
        for j in i:
            if(j in vOrder):
                continue
            else:
                vOrder.append(j)
    #print(vOrder)
    return vOrder



'''
Get the coordinates of the vertices which are found to be a part of the triangulation 
of the object. The order is the order in which the vertices will be written into the 
obj file and the vertices are gathered from graph.getNodes() function from the Dhart 
library
'''

def getVertices(order, vertices):
    verts = []
    '''
    for item in order:
        for j in vertices:
            if(j[4] == item):
                verts.append(j)
                break
            else:
                continue
    '''
    for item in order:
        temp = vertices[item]
        temp2 = (temp[0], temp[1], temp[2], item)
        verts.append(temp2)
    return verts



'''
Because the ordering of nodes into the .obj file is unknown before processing, there must 
be a reorientation of the nodes which are present in the faces of each triangle to align 
with the order of the nodes when they are written as vertices.
Example:
v 0.0000 -20.0000 3.385700
v -1.0000 -21.0000 3.417600
v -1.0000 -20.0000 3.157900
v 0.0000 -21.0000 3.654700
v -1.0000 -19.0000 3.033600
v 0.0000 -19.0000 3.249600
v 1.0000 -20.0000 3.609100
v -2.0000 -21.0000 3.154200
v -2.0000 -20.0000 2.929700
v -1.0000 -22.0000 3.712600
v -2.0000 -19.0000 2.808600
f 1 2 3
f 1 2 5
f 1 3 4
f 1 3 5
f 1 3 6
f 1 4 6
f 1 5 7
f 2 3 5
f 2 3 10
f 2 3 11
f 2 5 12
f 3 4 6
f 3 4 11
f 3 4 13
There are 11 nodes used as vertices of the triangles, but the faces call on nodes 12, 13 
which 'do not exist' in this sense as their old nodeID's were used when writing to the file.
'''
def realignFaces(faces, order, dict1):
    #dict1 = {}
    newFaces = []
    #for item in verts:
    #    #the 'old' nodeID will be the value, the new nodeID for .obj writing is the key
    #    #prevents the issue of face ID misalignements while writing.
    #    dict1[verts.index(item)+1] = item[4]
    #print(dict1)
    count = 0
    for item in faces:
        temp = []
        for j in item:
            #get the key (newNodeID) which is added to a temp face and then returned to a face
            #list to have proper triangle alignement.
            temp1 = list(dict1.keys())[list(dict1.values()).index(j)]
            temp.append(temp1)
        print('Face realignment: {:.3f}%'.format((count/len(faces))*100), end='\r')
        count +=1
        newFaces.append(tuple(temp))
    return newFaces

'''
create a node face index where the the key is the realignment ID and the value is the original 
Node ID. For use when constructing the obj file and for the computation of the normal vectors.
'''
def createAlignmentDict(verts):
    dict1 = {}
    newFaces = []
    for item in verts:
        #the 'old' nodeID will be the value, the new nodeID for .obj writing is the key
        #prevents the issue of face ID misalignements while writing.
        dict1[verts.index(item)+1] = item[3]
    return dict1

'''
faster duplicate removal for duplicate trinagle faces
'''
def dupRemoval(faces):
    temp = [tuple(sorted(sub)) for sub in faces]
    new_faces = list(set(temp))
    return new_faces

'''
Create a simple object file from the vertices and the discovered faces. 
'''
#def createOBJ(vOrder, faces, vertices):
def createOBJ(verts, newFaces, name):
    OBJfile = open(name, 'w')
    for item in verts:
        OBJfile.write("v {:.4f} {:.4f} {:4f}\n".format(item[0], item[1], item[2]))
    for item in newFaces:
        OBJfile.write("f {} {} {}\n".format(item[0], item[1], item[2]))
    OBJfile.close()



'''
---------------------------------------------------------A* search on the navmesh------------------------------------------------------
'''
'''
create an adjacency list for the A* search algorithm with format {ID: [(neigborID, distance)]}
'''

def euclideanStruct(undirStruct, coords):
    edist = {}
    for k,v in undirStruct.items():
        temp = []
        for item in v:
            if(k!=item):
                distV = np.array([coords[item][0], coords[item][1], coords[item][2]]) - np.array([coords[k][0], coords[k][1], coords[k][2]])
                dist = np.linalg.norm(distV)
                temp.append((item, dist))
            else:
                continue
        edist[k] = temp
    return edist

'''
create hueristic dictionary where the distance from every node to the goal node is determined
'''
def euclideanHeuristic(coords, goal):
    hDict = {}
    for i in range(len(coords)):
        #distV = np.array([coords[goal][0], coords[goal][1], coords[goal][2]]) - np.array([coords[i][0], coords[i][1], coords[i][2]])
        distV = goal - np.array([coords[i][0], coords[i][1], coords[i][2]])
        dist = np.linalg.norm(distV)
        hDict[i] = dist
    hDict['end'] = 0.0
    return hDict

#determine if a point is in a given triangle 
def pointInTri(point, tri, coords):
    v0 = np.array([coords[tri[2]][0], coords[tri[2]][1], coords[tri[2]][2]]) - np.array([coords[tri[0]][0], coords[tri[0]][1], coords[tri[0]][2]])
    v1 = np.array([coords[tri[1]][0], coords[tri[1]][1], coords[tri[1]][2]]) - np.array([coords[tri[0]][0], coords[tri[0]][1], coords[tri[0]][2]])
    #v2 = np.array([point[0], point[1], point[2]]) - np.array([coords[tri[0]][0], coords[tri[0]][1], coords[tri[0]][2]])
    v2 = point - np.array([coords[tri[0]][0], coords[tri[0]][1], coords[tri[0]][2]])
    dot00 = np.dot(v0,v0)
    dot01 = np.dot(v0,v1)
    dot02 = np.dot(v0,v2)
    dot11 = np.dot(v1,v1)
    dot12 = np.dot(v1,v2)
    invDenom = 1/((dot00*dot11)-(dot01*dot01))
    u = ((dot11*dot02)-(dot01*dot12))*invDenom
    v = ((dot00*dot12)-(dot01*dot02))*invDenom
    if(u>=0 and v>=0 and (u+v)<1):
        return True 
    else:
        return False

#find the starting and ending triangles
def detrmineStartEnd(faces, coords,start, end):
    startFin = []
    #find starting point triangle
    for item in faces:
        if(pointInTri(start, item, coords)):
            startFin.append(item)
            break
        else:
            continue
    for item in faces:
        if(pointInTri(end, item, coords)):
            startFin.append(item)
            break
        else:
            continue
    if(len(startFin) != 2):
        #print('Start or end point is not on the mesh!')
        return 1
    elif(startFin[0] == startFin[1]):
        #print('Start and end point on the same Triangle!')
        dist = np.linalg.norm(end - start)
        #print('They have a distance of: {}'.format(dist))
        return 2
    else:
        return startFin


#find the cost of the start/end point to the points of the triangle which contains it. 
def getNextPoint(point, tri, coords):
    dist1 = np.linalg.norm(point - np.array([coords[tri[0]][0], coords[tri[0]][1], coords[tri[0]][2]]))
    dist2 = np.linalg.norm(point - np.array([coords[tri[1]][0], coords[tri[1]][1], coords[tri[1]][2]]))
    dist3 = np.linalg.norm(point - np.array([coords[tri[2]][0], coords[tri[2]][1], coords[tri[2]][2]]))
    distList = [(tri[0], dist1), (tri[1], dist2), (tri[2], dist3)]
    return distList

#generate a random start and end point on two of the randomly selected faces of the mesh for a-star testing 
def generateStartEnd(faces, coords):
    #------------------
    #random start point
    #------------------
    startTri = random.choice(faces)
    line1 = np.array([coords[startTri[1]][0], coords[startTri[1]][1],coords[startTri[1]][2]]) - np.array([coords[startTri[0]][0], coords[startTri[0]][1],coords[startTri[0]][2]])
    line2 = np.array([coords[startTri[2]][0], coords[startTri[2]][1],coords[startTri[2]][2]]) - np.array([coords[startTri[0]][0], coords[startTri[0]][1],coords[startTri[0]][2]])
    tOnLine = random.uniform(0.0,1.0)
    #paramterized lines
    #point 1 on line1
    pl1 = np.array([coords[startTri[0]][0], coords[startTri[0]][1],coords[startTri[0]][2]]) + tOnLine*line1
    #point 2 on line2
    pl2 = np.array([coords[startTri[0]][0], coords[startTri[0]][1],coords[startTri[0]][2]]) + tOnLine*line2
    #find the line between the two points of the two previous lines, this paramaterized line will be a random point on the triangle
    line3 = pl2-pl1
    t2Online = random.uniform(0.0,1.0)
    startPoint = pl1 +t2Online*line3
    #---------------
    #random endpoint
    #---------------
    endTri = random.choice(faces)
    line4 = np.array([coords[endTri[1]][0], coords[endTri[1]][1],coords[endTri[1]][2]]) - np.array([coords[endTri[0]][0], coords[endTri[0]][1],coords[endTri[0]][2]])
    line5 = np.array([coords[endTri[2]][0], coords[endTri[2]][1],coords[endTri[2]][2]]) - np.array([coords[endTri[0]][0], coords[endTri[0]][1],coords[endTri[0]][2]])
    t3OnLine = random.uniform(0.0,1.0)
    #paramterized lines
    #point 1 on line1
    pl3 = np.array([coords[endTri[0]][0], coords[endTri[0]][1],coords[endTri[0]][2]]) + t3OnLine*line4
    #point 2 on line2
    pl4 = np.array([coords[endTri[0]][0], coords[endTri[0]][1],coords[endTri[0]][2]]) + t3OnLine*line5
    #find the line between the two points of the two previous lines, this paramaterized line will be a random point on the triangle
    line6 = pl4-pl3
    t4Online = random.uniform(0.0,1.0)
    endPoint = pl3 +t2Online*line4
    return (startPoint, endPoint)

#unused function 
def checkIfGoalTriangle(goalTri, inPoint):
    if(inPoint in goalTri):
        return True 
    else:
        return False 

#add the start and end nodes of a-star test to the goodinate dictionary for testing
def addSEtoadj(adj, start, startTri, end, endTri, coords):
    strtList = getNextPoint(start, startTri, coords)
    endList = getNextPoint(end, endTri, coords)
    for item in strtList:
        try:
            adj[item[0]].append(('start', item[1]))
        except KeyError:
            adj[item[0]] = [('start', item[1])]
    adj['start'] = strtList
    for item in endList:
        try:
            adj[item[0]].append(('end', item[1]))
        except KeyError:
            adj[item[0]] = [('end', item[1])]
    adj['end'] = endList
    return adj


#conduct the a-star test
def aStarSearch(adj, heuristic, start, end):
    aStarG = Graph(adj, heuristic)
    path = aStarG.a_star_algorithm(start, end)
    return path

#write the path found in the a-star path to csv and txt for use in importing the path into blender
def writePath(path, start, end,coords):
    pathFile = open('a_star_path.txt', 'w')
    pathFile.write('{}, {}, {}\n'.format(start[0], start[1], start[2]))
    for item in path[1:len(path)-1]:
        pathFile.write('{}, {}, {}\n'.format(coords[item][0], coords[item][1], coords[item][2]))
    pathFile.write('{}, {}, {}\n'.format(end[0], end[1], end[2]))
    pathFile.close()
    #write the csv for impprt to blender
    fields = ['x', 'y', 'z']
    rows = []
    strt = [start[0], start[1], start[2]]
    finish = [end[0], end[1], end[2]]
    rows.append(strt)
    for item in path[1:len(path)-1]:
        rows.append([coords[item][0], coords[item][1], coords[item][2]])
    rows.append(finish)
    with open('a_star_path.csv', 'w') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
#------------------------------------------END A* FUNCTIONS---------------------------------------------------------------------

'''
dictionary tester for the new addition of the incidcence dictionary being created while traversing the graph checkign for triangles as well
'''
def incidencetest(control, test):
    print(type(test))
    print('Length of control {} length of test {}'.format(len(control), len(test)))
    if(len(control) != len(test)):
        print('lengths are not the same')
        print(control.keys())
        print(test.keys())
    else:
        key1 = list(control.keys())
        key2 = list(test.keys())
        for i in range(len(control)):
            if(key1[i] != key2[i]):
                return ''
            else:
                print('the keys are the same')
                for j in range(len(control)):
                    vals1 = control[j]
                    vals2 = test[j]
                    if(len(vals1) != len(vals2)):
                        print('values at index {} not the same length'.format(i))
                        return''
                    else:
                        size = len(vals1)
                        #print('------------')
                        #print(vals1)
                        #print(vals2)
                        intersect = set(vals1).intersection(set(vals2))
                        if(len(intersect) == len(vals1) and len(intersect) == len(vals2)):
                            print('vals are the same')
                        else:
                            print("VALS ARE NOT THE SAME")
    print('they are the same')
    return''




'''
make a unified function call that checks if the user wishes to create an obj along with additional 
computations mkOBJ defaults to False 

Args:
nodeList: a list of generated nodes and identifying information from Dhart
CSRmat: an adjacency matrix generated from Dhart graphToCSR function
mkOBJ: T/F Defaults to false. Make an obj File of the generated mesh
OBJpath: Defaults to  '' The path in which to write the obj file
'''
'''
def makeMesh(nodeList, CSRmat,mkOBJ = False, OBJpath ='' ):
    startTime = time.time()
    coords = createCoordDict(nodeList)
    graphStruct = createGraphStruct(CSRmat)
    gDict = createGraphDict(graphStruct)
    et = startTime - time.time()
    print('Base components created. Elapsed time: --- {:2f} seconds ---'.format(et))
    faces = traverseGraph(gDict, 0)
    noDup =  dupRemoval(faces)
    et2 = et - time.time()
    print('Graph traversed duplicate triangles removed. Total triangles in graph: {}. Elasped time: --- {:2f} seconds ---'.format(len(noDup), et2))
    if(mkOBJ):
        print('Starting OBJ file write...')
        verts = getVertices(order, nodeList)
        alignDict = createAlignmentDict(verts)
        alignedFaces = realignFaces(noDup, order, alignDict)
        createOBJ(verts, alignedFaces, OBJpath)
        et3 = et2 - time.time()
        print('OBJ file written to: {}, Elapsed time: --- {:2f} minutes---'.format(OBJpath, et3))
        normals = computeNormal(alignedFaces, alignDict, coords)
        et4 = et3 - time.time() 
        print('Normal vectors computed for triangles. Elapsed time: --- {:2f} minutes ---'.format(et4))
    else:
        print('OBJ file not written. Computing normal vectors...')
'''
'''

Function calls
'''
print('Starting graphTraversal.py...')
startTime = time.time()
buildTimeS = time.time()

print('---------------------CREATING NECESSARY STRUCTURES--------------')
#create the coordinate dictionary k=node id v=coodinates
coords = createCoordDict(nodes)
#print(coords)
print('Node coordinate dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

#create the graph structure list 
graphStruct = createGraphStruct(csr_graph)
distDict = distanceList(graphStruct)
#print(graphStruct)
#print(distDict)
print('CSR matrix converted to list')
print("--- %s seconds ---"% (time.time() - startTime))

#create undirected graph structure for use with overlap removal
waystruct = structFromCSR(graphStruct)
#print(waystruct)
print('Unidrected graph structure dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

#convert the graph structure list into a graph structure dictionary k= node id v= connected nodes
gDict = createGraphDict(graphStruct)
#print(gDict)
print('Graph structure dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))
print('-----------------------------------------------------------------------------------')

print('----------------------CREATING THE MESH--------------------------------------------')
#travrse graph using BFS, get all all triangles having form (nodeID, nodeID, nodeID)
faces = traverseGraph(waystruct, 0)
#print(faces[2])
print('Graph traversed, all trinagles found.')
print("--- %s seconds ---"% (time.time() - startTime))

#remove dulicate triangle faces 
new_faces = dupRemoval(faces[1])
#new_faces = removeDuplicate(faces[1])
#print(new_faces)
print('Duplicate triangles removed.')
print("--- %s seconds ---"% (time.time() - startTime))


'''
triangleIncidence = getTriangleIncidence(coords, new_faces)
print("Triangle incidence dictionary created.")
#print(triangleIncidence)
print("--- %s seconds ---"% (time.time() - startTime))
'''

#incidencetest(triangleIncidence, faces[2])


#remove overlapping faces (TYPE 1)
#rmvOverlap = removeOverlap(wOrderFaces, coords)
#print(faces[2])
rmvOverlap = fastOverlapRemoval(waystruct, new_faces, distDict, faces[2])
print('Type 1 overlaps Removed.')
print("--- %s seconds ---"% ((time.time() - startTime)))

#remove TYPE 2 overlap
rmvEdgeOvlp = edgeOverlap(waystruct, distDict, rmvOverlap, coords)
coords = rmvEdgeOvlp[0]
#print(rmvOverlap)
print('Type 2 Overlaps removed.')
print("--- %s seconds ---"% ((time.time() - startTime)))


#complete the winding order of the faces
wOrderFaces = windingOrder(rmvEdgeOvlp[1], coords)
#print(len(wOrderFaces))
#wOrderFaces = dupRemoval(wOrderFaces)
#print(len(wOrderFaces))
print('Winding order complete.')
print("--- %s seconds ---"% (time.time() - startTime))
buildtimeE = time.time() - buildTimeS
print('----------------INITIAL MESH CREATED----------------------------')

#---------------DECIMATION FUNCTIONS-------------------------
print('------------------------DECIMATION FUNCTIONS--------------------')
arrFaces = quadMeshCreator(wOrderFaces, coords)
#print(arrFaces)
#print('old faces: ', wOrderFaces)
print('Quad Mesh array of faces created.')
print("--- %s seconds ---"% (time.time() - startTime))

arrPos = quadMeshPositions(wOrderFaces, coords)
#print(arrPos)
print('Quad Mesh array of Positions created.')
#print('Old coords: ', coords)
print("--- %s seconds ---"% (time.time() - startTime))

#-------------------------------------DECIMATION ----------------------------------------------
numNodes = math.ceil(len(set([j for item in wOrderFaces for j in item])) * percentage)
if(numNodes < 3):
    print('NumNodes must be greater than 3 for mesh decimation.')
    exit(1)
decStart = time.time()
newPositions, newMeshFace = simplify_mesh(arrPos, arrFaces, numNodes, threshold=0.5)
decFin = time.time() - decStart
#print('New Positions: ', newPositions)
#print('New Faces: ', newMeshFace)
print('Completed quad mesh decimation.')
print("--- %s seconds ---"% (time.time() - startTime))

decFaces = convertArrToList(newMeshFace)
#print('decFaces: ',decFaces)
#print(len(decFaces))
print('Array of decimated faces converted to list of tuples.')
print("--- %s seconds ---"% (time.time() - startTime))

decCoords = convertCoordDict(newPositions)
#print('decCoords: ',decCoords)
print('Array of corrdinates translated to new coordinate dictionary.')
print("--- %s seconds ---"% (time.time() - startTime))
#------------------------------------------------------------------------------------------------

decimationDict = createDecimationDict(coords, wOrderFaces)
print('Decimation dictionary created.')
#print(decimationDict) 
print("--- %s seconds ---"% (time.time() - startTime))

triangleIncidence = getTriangleIncidence(coords, wOrderFaces)
print("Triangle incidence dictionary created.")
#print(triangleIncidence)
print("--- %s seconds ---"% (time.time() - startTime))

visualImportanceDict = visualImportance(decimationDict, triangleIncidence)
print('Visual Importance dictionary created.')
#print(visualImportanceDict)
print("--- %s seconds ---"% (time.time() - startTime))

borders = getBorderEdges(wOrderFaces)
print('Found borders of the mesh.')
#print(borders)
print("--- %s seconds ---"% (time.time() - startTime))

halfEdges = getHalfEdges(gDict, borders)
print('Gathered all half edges minus borders.')
#print(halfEdges)
print("--- %s seconds ---"% (time.time() - startTime))

gError = geometricError(halfEdges, triangleIncidence, borders, coords)
print('Calculated geometric error.')
#print(gError)
print("--- %s seconds ---"% (time.time() - startTime))

decimationPQ = optimizeCollapse(gError, visualImportanceDict, coords, wOrderFaces, borders)
print("Created PQ for halfedge collapses")
#print(decimationPQ)
print("--- %s seconds ---"% (time.time() - startTime))

decimateResult = decimate(geometricError, visualImportanceDict, coords, wOrderFaces, borders, decimationPQ)
print("decitmation complete.")
#print(decimateResult)
print("--- %s seconds ---"% (time.time() - startTime))

#complete the winding order of the faces
secondWOrder = windingOrder(decimateResult, coords)
#print(wOrderFaces)
print('2nd pass Winding order complete.')
print("--- %s seconds ---"% (time.time() - startTime))

print('--------------------------DECIMATION COMPLETE---------------------------------------')
#-------------------------------------------------------------

print('---------------------------WRITIGN OBJ FILES-----------------------------------------')
'''
write the decimated mesh to an obj file 
------------------------------------------------------------------------------------------
'''

print('Writign decimated mesh to obj file')
#get order of appearance of vertices in the faces for obj writing 
order = getVertexOrder(decFaces)
#order = getVertexOrder(wOrderFaces) #wOrderFaces
#print(order)
print('Vertex order created.')
print("--- %s seconds ---"% (time.time() - startTime))

#get cooridinates of the vertices for obj writing
#can probably update to use coords for faster look up using dictionary
#verts = getVertices(order, rmvEdgeOvlp[0]) #nodes
verts = getVertices(order, decCoords)
#print(verts)
print('Vertex coordinates found.')
print("--- %s seconds ---"% (time.time() - startTime))

#create a node face index where the the key is the realignment ID and the value is the original 
#node ID. Used for fast conversion between aligned nodes and non alinged nodes
alignDict = createAlignmentDict(verts)
#print(alignDict)
print('Alignment dictionary for .obj writing created.')
print("--- %s seconds ---"% (time.time() - startTime))

#ensure that nodes are assinged the correct id in the faces for obj writing
#alignedFaces = realignFaces(rmvEdgeOvlp[1], order, alignDict)
#alignedFaces = realignFaces(wOrderFaces, order, alignDict)
alignedFaces = realignFaces(decFaces, order, alignDict)
#print(alignedFaces)
print('Faces aligned for .obj writing.')
print("--- %s seconds ---"% (time.time() - startTime))

#write the obj 
createOBJ(verts, alignedFaces, 'decimatedMesh.obj')
print('.obj file written')
print("--- %s seconds ---"% ((time.time() - startTime)))


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#Write the original mesh before decimation to an obj file
print('writing original mesh to obj file')
#get order of appearance of vertices in the faces for obj writing 
orderOrig = getVertexOrder(wOrderFaces) #wOrderFaces
#print(order)
print('Vertex order created.')
print("--- %s seconds ---"% (time.time() - startTime))

#get cooridinates of the vertices for obj writing
#can probably update to use coords for faster look up using dictionary
vertsOrig = getVertices(orderOrig, rmvEdgeOvlp[0]) #nodes
#print(verts)
print('Vertex coordinates found.')
print("--- %s seconds ---"% (time.time() - startTime))

#create a node face index where the the key is the realignment ID and the value is the original 
#node ID. Used for fast conversion between aligned nodes and non alinged nodes
alignDictOrig = createAlignmentDict(vertsOrig)
#print(alignDict)
print('Alignment dictionary for .obj writing created.')
print("--- %s seconds ---"% (time.time() - startTime))

#ensure that nodes are assinged the correct id in the faces for obj writing
#alignedFaces = realignFaces(rmvEdgeOvlp[1], order, alignDict)
alignedFacesOrig = realignFaces(wOrderFaces, orderOrig, alignDictOrig)
#print(alignedFacesOrig)
print('Faces aligned for .obj writing.')
print("--- %s seconds ---"% (time.time() - startTime))

#write the obj 
createOBJ(vertsOrig, alignedFacesOrig, 'originalMesh.obj')
print('.obj file written')
print("--- %s seconds ---"% ((time.time() - startTime)))
print('----------------------------------------------------------------------------------------')


#compute the normals for all faces
#normals = computeNormal(rmvEdgeOvlp[1], coords)
#normals = computeNormal(wOrderFaces, coords)
#print(normals)
#print('Normal vectors of all faces found.')
#print("--- %s seconds ---"% ((time.time() - startTime)))

#get the faces which share an edge
#sharedFaces = getSharedFaces(rmvEdgeOvlp[1])
#sharedFaces = getSharedFaces(wOrderFaces)
#print(sharedFaces)
#print('Faces with shared edges determined.')
#print("--- %s seconds ---"% ((time.time() - startTime)))



####SINGLE FUNCTION 
print('#######################################################################################')
print('----------------------------A* FUNCTIONS original MESH--------------------------------') 
#create a post decimation adjacecny structure for use on A* *** NON DECIMATED MESHES MAY USE PREVIOUS ADJACECNY STRCUTRE: WAYSTRUCT ***
#postDecGStruct = createNewAdjDict(wOrderFaces)
#print(postDecGStruct)
#print('Post decimation adjacecny dict created.')
#print("--- %s seconds ---"% (time.time() - startTime))

#create a euclidean distance dictionary from the undirected graph structure
euclidS = euclideanStruct(waystruct, coords)
#print(euclidS)
print('Euclidean distance dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

beginEnd = generateStartEnd(wOrderFaces, coords)
print('Created random start and end points on the mesh. ')
#print(beginEnd)
print("--- %s seconds ---"% (time.time() - startTime))

startFinList = detrmineStartEnd(wOrderFaces, coords, beginEnd[0], beginEnd[1])
print('Found starting and ending triangles.')
#print(startFinList)
print("--- %s seconds ---"% (time.time() - startTime))

updatedAdj = addSEtoadj(euclidS, beginEnd[0], startFinList[0], beginEnd[1], startFinList[1], coords)
print('Updated adjacency list to include start and end points.')
#print(updatedAdj)
print("--- %s seconds ---"% (time.time() - startTime))

#create a euclidean distance heuristic dictionary for the distance from each node to goal node
euclidH = euclideanHeuristic(coords, beginEnd[1])
#print(euclidH)
print('Hueristic dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

#find the optimal path usign a* search 
optimalPath = aStarSearch(updatedAdj, euclidH, 'start', 'end')
print('Shortest path between start and end nodes on mesh found.')
print('Optimal path, Non-Decimated: ', optimalPath)
print("--- %s seconds ---"% (time.time() - startTime))

writePath(optimalPath, beginEnd[0], beginEnd[1], decCoords)
print('text file with path vertices written')

print('------------------------------------------------------------------------------------------------')
print('################################################################################################')
print('----------------------------A* FUNCTIONS SCALE TEST-Original MESH-------------------------------')
numTest = 200
#create a euclidean distance dictionary from the undirected graph structure
euclidS = euclideanStruct(waystruct, coords)
#print(euclidS)
print('Euclidean distance dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))
aStarDict = {}
numNonPaths = 0
for i in range(numTest):
    print('{:.3f}%'.format((i/numTest)*100), end='\r')
    beginEnd = generateStartEnd(wOrderFaces, coords)
    startFinList = detrmineStartEnd(wOrderFaces, coords, beginEnd[0], beginEnd[1])
    if(startFinList == 1):
        #one of the points is not on the mesh
        numNonPaths+=1
        continue
    elif(startFinList == 2):
        # points are on the same triangle
        aStarDict[i] = AstarEnd
        continue
    else:
        updatedAdj = addSEtoadj(euclidS, beginEnd[0], startFinList[0], beginEnd[1], startFinList[1], coords)
        #create a euclidean distance heuristic dictionary for the distance from each node to goal node
        euclidH = euclideanHeuristic(coords, beginEnd[1])
        #find the optimal path usign a* search 
        AstartTime = time.time()
        optimalPath = aStarSearch(updatedAdj, euclidH, 'start', 'end')
        AstarEnd = time.time() - AstartTime
        aStarDict[i] = AstarEnd
        continue
print('Paths found: {}'.format(len(aStarDict)))
print('Non-paths encountered: {}'.format(numNonPaths))
timeTotal = 0
for k,v in aStarDict.items():
    timeTotal += v
print('Average time for path determination: {}'.format(timeTotal/len(aStarDict)))
filewrite = open('astarTEST_ORIGINAL.txt', 'w')
filewrite.write('A* batch test on {} nodes, moon.obj\n'.format(max_nodes))
filewrite.write('Tests conducted: {}\n'.format(numTest))
filewrite.write('Paths Found: {}\n'.format(len(aStarDict)))
filewrite.write('Non-paths encountered: {}\n'.format(numNonPaths))
filewrite.write('Total test time: {}\n'.format(timeTotal))
filewrite.write('Average time for path determination: {}\n'.format(timeTotal/len(aStarDict)))

print('------------------COMPLETE BATCH TEST A* original mesh---------------------------------')
print('#######################################################################################')
print('#######################################################################################')


print('#######################################################################################')
print('#######################################################################################')
print('----------------------------A* FUNCTIONS DECIMATED MESH--------------------------------')
#create a post decimation adjacecny structure for use on A* *** NON DECIMATED MESHES MAY USE PREVIOUS ADJACECNY STRCUTRE: WAYSTRUCT ***
postDecGStruct = createNewAdjDict(decFaces)
#print(postDecGStruct)
print('Post decimation adjacecny dict created.')
print("--- %s seconds ---"% (time.time() - startTime))

#create a euclidean distance dictionary from the undirected graph structure
euclidS = euclideanStruct(postDecGStruct, decCoords)
#print(euclidS)
print('Euclidean distance dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

beginEnd = generateStartEnd(decFaces, decCoords)
print('Created random start and end points on the mesh. ')
#print(beginEnd)
print("--- %s seconds ---"% (time.time() - startTime))

startFinList = detrmineStartEnd(decFaces, decCoords, beginEnd[0], beginEnd[1])
print('Found starting and ending triangles.')
#print(startFinList)
print("--- %s seconds ---"% (time.time() - startTime))

updatedAdj = addSEtoadj(euclidS, beginEnd[0], startFinList[0], beginEnd[1], startFinList[1], decCoords)
print('Updated adjacency list to include start and end points.')
#print(updatedAdj)
print("--- %s seconds ---"% (time.time() - startTime))

#create a euclidean distance heuristic dictionary for the distance from each node to goal node
euclidH = euclideanHeuristic(decCoords, beginEnd[1])
#print(euclidH)
print('Hueristic dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))

#find the optimal path usign a* search 
optimalPath = aStarSearch(updatedAdj, euclidH, 'start', 'end')
print('Shortest path between start and end nodes on mesh found.')
print('Optimal path, Decimated Mesh: ',optimalPath)
print("--- %s seconds ---"% (time.time() - startTime))

writePath(optimalPath, beginEnd[0], beginEnd[1], decCoords)
print('text file with path vertices written')
print('---------------------------------------------------------------------------------------')



print('----------------------------A* FUNCTIONS SCALE TEST-DECIMATED MESH-------------------------------')
numTest = 10
#create a post decimation adjacecny structure for use on A* *** NON DECIMATED MESHES MAY USE PREVIOUS ADJACECNY STRCUTRE: WAYSTRUCT ***
postDecGStruct = createNewAdjDict(decFaces)
#print(postDecGStruct)
print('Post decimation adjacecny dict created.')
print("--- %s seconds ---"% (time.time() - startTime))

#create a euclidean distance dictionary from the undirected graph structure
euclidS = euclideanStruct(postDecGStruct, decCoords)
#print(euclidS)
print('Euclidean distance dictionary created.')
print("--- %s seconds ---"% (time.time() - startTime))
aStarDict = {}
pathLengths = []
numNonPaths = 0
for i in range(numTest):
    print('{:.3f}%'.format((i/numTest)*100), end='\r')

    beginEnd = generateStartEnd(decFaces, decCoords)
    startFinList = detrmineStartEnd(decFaces, decCoords, beginEnd[0], beginEnd[1])
    if(startFinList == 1):
        #one of the points is not on the mesh
        pathLengths.append(0)
        numNonPaths+=1
        continue
    elif(startFinList == 2):
        # points are on the same triangle
        #plen = np.linalg.norm(beginEnd[1] - beginEnd[0])
        #print('PLEN:' ,plen)
        aStarDict[i] = AstarEnd
        continue
    else:
        updatedAdj = addSEtoadj(euclidS, beginEnd[0], startFinList[0], beginEnd[1], startFinList[1], decCoords)
        #create a euclidean distance heuristic dictionary for the distance from each node to goal node
        euclidH = euclideanHeuristic(decCoords, beginEnd[1])
        #find the optimal path usign a* search 
        AstartTime = time.time()
        optimalPath = aStarSearch(updatedAdj, euclidH, 'start', 'end')
        AstarEnd = time.time() - AstartTime
        aStarDict[i] = AstarEnd
        continue
print('Paths found: {}'.format(len(aStarDict)))
print('Non-paths encountered: {}'.format(numNonPaths))
timeTotal = 0
for k,v in aStarDict.items():
    timeTotal += v
print('Average time for path determination decimated mesh: {}'.format(timeTotal/len(aStarDict)))
filewrite = open('astarTEST_DECIMATED.txt', 'w')
filewrite.write('A* batch test on {} nodes, moon.obj\n'.format(max_nodes))
filewrite.write('Tests conducted: {}\n'.format(numTest))
filewrite.write('Paths Found: {}\n'.format(len(aStarDict)))
filewrite.write('Non-paths encountered: {}\n'.format(numNonPaths))
filewrite.write('Total test time: {}\n'.format(timeTotal))
filewrite.write('Average time for path determination: {}\n'.format(timeTotal/len(aStarDict)))
print('------------------COMPLETE BATCH TEST A*-------------------------------------')
print('##############################################################################')
'''
print('----------------------BATCH TEST DECIMATE------------------------------------')
maxNodeList = [50, 500, 5000, 50000,100000,500000]
#REDUCE THE NUMBER OF NODES IN THE MESH BY THIS PERCENTAGE (0,1)
#IF THERE ARE 100 NODES IN THE MESH A PERCENTAGE VALUE OF 0.9 HAS A TARGET REDUCTION OF 10 NODES FOR 90 REMAINING
percentages = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
resultDict6 = decimateTest(maxNodeList, percentages, bvh)
print(resultDict6)
'''

'''
TESTING script functions here, not needed for the running of the pipeline

'''
print('###############################################################################')
print('###############################################################################')
testSumm = open('TESTSUMMARY.txt', 'w')
testSumm.write('Mesh Build time: {}\n'.format(buildtimeE))
testSumm.write('Decimation percentage: {}\n'.format(percentage))
testSumm.write('Decimation completed in: {}\n'.format(decFin))
testSumm.write('Faces prior: {}\n'.format(len(alignedFacesOrig)))
testSumm.write('Vertices Prior: {}\n'.format(len(coords)))
testSumm.write('Faces after: {}\n'.format(len(alignedFaces)))
testSumm.write('Vertices After: {}\n'.format(len(decCoords)))
testSumm.write('Spacing: {}\n'.format(spacing))
testSumm.write('Upslope/Downslope: {}/{}\n'.format(up_slope, down_slope))
testSumm.write('Upstep/Downstep: {}/{}\n'.format(up_step, down_step))
print('###############################################################################')
print('###############################################################################')
print('Mesh Build time: {}'.format(buildtimeE))
print('Decimation percentage: {}'.format(percentage))
print('Decimation completed in: {}'.format(decFin))
print('Faces prior: {}'.format(len(alignedFacesOrig)))
print('Vertices Prior: {}'.format(len(coords)))
print('Faces after: {}'.format(len(alignedFaces)))
print('Vertices After: {}'.format(len(decCoords)))
print('Spacing: {}'.format(spacing))
print('Upslope/Downslope: {}/{}'.format(up_slope, down_slope))
print('Upstep/Downstep: {}/{}'.format(up_step, down_step))
