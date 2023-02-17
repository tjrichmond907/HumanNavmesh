import numpy as np
import matplotlib.pyplot as plt
from dhart.geometry import LoadOBJ, CommonRotations
from dhart.raytracer import EmbreeBVH
from dhart.graphgenerator import GenerateGraph
#from dhart.graph import GetNodeID, getNodes, ConvertToList
import dhart.spatialstructures.graph
from graphTraversal import *
import dhart
import time
import math
from graph import Graph
from quad_mesh_simplify import simplify_mesh
import random
'''
###
Run in cmd-prompt with: python ..\..\OneDrive\Desktop\GAIDG\createTriangleTest.py
###
'''

#obj_path = dhart.get_sample_model("plane.obj")
obj_path = 'C:/Users/tjric/OneDrive/Desktop/GAIDG/moon.obj'

obj = LoadOBJ(obj_path, rotation=CommonRotations.Yup_to_Zup)
#obj = LoadOBJ(obj_path)
bvh = EmbreeBVH(obj, True)
'''
start_point = (-0, -20, 20)
#start_point = (-1, -6, 1623.976928)

spacing = (1, 1, 10)
max_nodes = 5
#up_step, down_step = 0.5,0.5 
#up_slope, down_slope = 20, 20

up_step, down_step = 90,90 
up_slope, down_slope = 90, 90
max_step_connections = 1


graph = GenerateGraph(bvh, start_point, spacing, max_nodes,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, cores=-1)

#graph = GenerateGraph(bvh, start_point, spacing, max_nodes, cores=-1)

csr_graph = graph.CompressToCSR()
nodes = graph.getNodes()
'''

def test(maxNodes, BVH):
    start_point = (-0.0226852,-11.3705,0.126856) #moon
    spacing = (10, 10, 10)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        #print(i)
        print("Begining test of BFS traversal")
        print("Creating graph from Dhart generate Graph...")
        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        print('Creating graph structure and dictionary...')
        graphStruct = createGraphStruct(csr_graph)
        gDict = createGraphDict(graphStruct)
        print('Traversing graph and creating faces of triangles...')
        startTime = time.time()
        faces = traverseGraph(gDict, 0)
        endTime = time.time() - startTime
        num_tris = len(faces[1])
        num_verts = len(faces[0])
        print('Done')
        print('Completed test of BFS taversal. Found {} triangles with {} vertices visited in {} seconds.'.format(num_tris,num_verts,endTime))
        results[i] = endTime
    return results

def plotGeneralResults(results):
    keys = list(results.keys())
    values = list(results.values())
    plt.title('Dhart Graph Generation through Traversal Test')
    plt.xlabel('Max Nodes Used')
    plt.ylabel('Time (seconds)')
    plt.bar(range(len(results)), values, tick_label=keys)
    plt.show()

def plotDuplicateResults(results):
    keys = list(results.keys())
    values = list(results.values())
    txt = 'Results for duplicate test using createTriangles.py.\n Test Cases which had no Duplicates are included with\n runtimes and nodes generated.'
    plt.title('Duplicate Traingles Test')
    plt.xlabel('Max Nodes Used\n '+txt)
    plt.ylabel('Time (seconds)')
    plt.bar(range(len(results)), values, tick_label=keys)
    plt.show()

def plotWindingOrderResults(results):
    keys = list(results.keys())
    values = list(results.values())
    plt.title('Dhart Graph Winding Order Test')
    plt.xlabel('Max Nodes Used')
    plt.ylabel('Uniform normal (z) value')
    plt.bar(range(len(results)), values, tick_label=keys)
    plt.show()

def plotOverlapResults(results):
    keys = list(results.keys())
    values = list(results.values())
    txt = 'Results for Overlap removal test using createTriangles.py.'
    plt.title('Remove Overlap Test')
    plt.xlabel('Max Nodes Used\n '+txt)
    plt.ylabel('Time (seconds)')
    plt.bar(range(len(results)), values, tick_label=keys)
    plt.show()

def plotAstarResults(results):
    keys = list(results.keys())
    values = list(results.values())
    txt = 'Results for A* search test using createTriangles.py.'
    plt.title('A* search Test')
    plt.xlabel('Max Nodes Used\n '+txt)
    plt.ylabel('Time (seconds)')
    plt.bar(range(len(results)), values, tick_label=keys)
    plt.show()

'''
faster duplicate removal for duplicate trinagle faces
'''
def dupRemove(faces):
    #new_faces = list(set(tuple(set(item)) for item in faces))
    temp = [tuple(sorted(sub)) for sub in faces]
    new_faces = list(set(temp))
    #print(new_faces)
    return new_faces

def duplicateTest(maxNodes, BVH):
    start_point = (-0.0226852,-11.3705,0.126856) #moon
    spacing = (10, 10, 10)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        #print(i)
        print("Begining test of Duplicate Triangles")
        print("Creating graph from Dhart generate Graph...")

        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        print('Creating graph structure and dictionary...')
        graphStruct = createGraphStruct(csr_graph)
        gDict = createGraphDict(graphStruct)
        print('Traversing graph and creating faces of triangles...')
        faces = traverseGraph(gDict, 0)
        trisWithDup = len(faces[0])
        num_verts = len(faces[1])
        print('Removing duplicate triangles...')
        #new_faces = removeDuplicate(faces[1])
        startTime = time.time()
        new_faces = dupRemoval(faces[1])
        endTime = time.time() - startTime
        print('Starting comparison of faces which have been deemed clear of duplicates')
        count = 0
        '''
        for j in range(len(new_faces)):
            for k in range(j+1,len(new_faces)):
                if(set(new_faces[j]) == set(new_faces[k])):
                    print("Duplicate triangles at indexes {} and {} in new_faces".format(j,k))
                    print('j = {} with {}, k = {} with {}'.format(j,new_faces[j], k,new_faces[k]))
                    count += 1
                    continue
                else:
                    print('{:.3f}%'.format((1/(len(new_faces)))*100),end='\r')
                    continue
        '''
        if(count == 0):
            print("No duplicates detetected")
            print('Completed duplicate test. Duplcates counted: {}. Completed graph creation, traversal, and duplicate removal in {} seconds'.format(count, endTime))
            results[i] = endTime
        else:
            print('Completed duplicate test. Duplcates counted: {}. Completed graph creation, traversal, and duplicate removal in {} seconds'.format(count, endTime))
    return results


def windingOrderTest(maxNodes, BVH):
    start_point = (-0.0226852,-11.3705,0.126856) #moon
    spacing = (10, 10, 10)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        #print(i)
        print("Begining test of Duplicate Triangles")
        print("Creating graph from Dhart generate Graph...")
        startTime = time.time()
        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        print('Creating node coordinate dictionary...')
        coords = createCoordDict(nodes)
        print('Creating graph structure and dictionary...')
        graphStruct = createGraphStruct(csr_graph)
        gDict = createGraphDict(graphStruct)
        print('Traversing graph and creating faces of triangles...')
        faces = traverseGraph(gDict, 0)
        trisWithDup = len(faces[0])
        num_verts = len(faces[1])
        print('Removing duplicate triangles...')
        #new_faces = removeDuplicate(faces[1])
        new_faces = dupRemoval(faces[1])
        print('Creating uniform winding order...')
        startTime = time.time()
        wOrderFaces = windingOrder(new_faces, coords)
        endTime = time.time() - startTime
        print('Unifrom winding order complete in {} seconds'.format(endTime))
        print('Verifying uniform Winding order...')
        normals = computeNormal(wOrderFaces, coords)
        count = 0
        for item in normals:
            if(item[2] < 0):
                count += 1
            else:
                continue
        if(count == 0):
            results[i] = True
        else:
            results[i] = False
        print('Found {} normals facing (-z)'.format(count))
    return results

def removeOverlapTest(maxNodes, BVH):
    start_point = (-0.0226852,-11.3705,0.126856) #moon
    spacing = (10, 10, 10)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        print("Begining test of Overlap removal")
        print("Creating graph from Dhart generate Graph...")

        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        print('Creating node coordinate dictionary...')
        coords = createCoordDict(nodes)
        print('Creating graph structure and distance dictionary, directed and undirected...')
        graphStruct = createGraphStruct(csr_graph)
        distDict = distanceList(graphStruct)
        waystruct = structFromCSR(graphStruct)
        gDict = createGraphDict(graphStruct)
        print('Traversing graph and creating faces of triangles...')
        faces = traverseGraph(gDict, 0)
        trisWithDup = len(faces[0])
        num_verts = len(faces[1])
        print('Removing duplicate triangles...')
        #new_faces = removeDuplicate(faces[1])
        new_faces = dupRemoval(faces[1])
        print('Removing overlapping triangles...')
        startTime = time.time()
        rmvOverlap = fastOverlapRemoval(waystruct, new_faces, distDict,faces[2])
        rmvOverlap2 = edgeOverlap(waystruct, distDict, rmvOverlap, coords)
        endTime = time.time() - startTime
        print('triangle removal order complete in {} seconds'.format(endTime))
        print('Creating uniform winding order...')
        wOrderFaces = windingOrder(rmvOverlap2[1], rmvOverlap2[0])
        results[i] = endTime
    return results

def aStarTest(maxNodes, BVH):
    start_point = (-0, -20, 20)
    spacing = (1, 1, 10)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        #print("Begining test of Overlap removal")
        #print("Creating graph from Dhart generate Graph...")
        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        #print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        #print('Creating node coordinate dictionary...')
        coords = createCoordDict(nodes)
        #print('Creating graph structure and distance dictionary, directed and undirected...')
        graphStruct = createGraphStruct(csr_graph)
        distDict = distanceList(graphStruct)
        waystruct = structFromCSR(graphStruct)
        gDict = createGraphDict(graphStruct)
        #print('Traversing graph and creating faces of triangles...')
        faces = traverseGraph(gDict, 0)
        trisWithDup = len(faces[0])
        num_verts = len(faces[1])
       # print('Removing duplicate triangles...')
        #new_faces = removeDuplicate(faces[1])
        new_faces = dupRemoval(faces[1])
        #print('Removing overlapping triangles...')
        rmvOverlap = fastOverlapRemoval(waystruct, new_faces, distDict)
        rmvOverlap2 = edgeOverlap(waystruct, distDict, rmvOverlap, coords)
        #print('triangle removal order complete in {} seconds'.format(endTime))
        #print('Creating uniform winding order...')
        wOrderFaces = windingOrder(rmvOverlap2[1], rmvOverlap2[0])
        print('Beginning A* search test.')
        startTime = time.time()
        euclidS = euclideanStruct(waystruct, rmvOverlap2[0])
        beginEnd = generateStartEnd(wOrderFaces, rmvOverlap2[0])
        startFinList = detrmineStartEnd(wOrderFaces, rmvOverlap2[0], beginEnd[0], beginEnd[1])
        if(startFinList == 1 or startFinList == 2):
            print('Generated start endpoints outside of graph or in the same triangle.')
            endTime = time.time() - startTime
            results[i] = endTime
            continue
        updatedAdj = addSEtoadj(euclidS, beginEnd[0], startFinList[0], beginEnd[1], startFinList[1], rmvOverlap2[0])
        euclidH = euclideanHeuristic(rmvOverlap2[0], beginEnd[1])
        optimalPath = aStarSearch(updatedAdj, euclidH, 'start', 'end')
        print(optimalPath)
        endTime = time.time() - startTime
        results[i] = endTime
    return results 

def decimateTest(maxNodes, decimateLevel, BVH):
    start_point = (-0.0226852,-11.3705,0.126856) #moon
    spacing = (20, 20, 20)
    up_step, down_step = 90,90 
    up_slope, down_slope = 90, 90
    max_step_connections = 1
    min_connections = 7
    results = {}
    for i in maxNodes:
        graph = GenerateGraph(BVH, start_point, spacing, i,
                    up_step, up_slope, down_step, down_slope,
                    max_step_connections, min_connections, cores=-1)
        #print('Compressing to CSR...')
        csr_graph = graph.CompressToCSR()
        nodes = graph.getNodes()
        #print('Creating node coordinate dictionary...')
        coords = createCoordDict(nodes)
        #print('Creating graph structure and distance dictionary, directed and undirected...')
        graphStruct = createGraphStruct(csr_graph)
        distDict = distanceList(graphStruct)
        waystruct = structFromCSR(graphStruct)
        gDict = createGraphDict(graphStruct)
        #print('Traversing graph and creating faces of triangles...')
        faces = traverseGraph(gDict, 0)
        new_faces = dupRemoval(faces[1])
        rmvOverlap = fastOverlapRemoval(waystruct, new_faces, distDict, faces[2])
        rmvOverlap2 = edgeOverlap(waystruct, distDict, rmvOverlap, coords)
        wOrderFaces = windingOrder(rmvOverlap2[1], rmvOverlap2[0])
        for j in decimateLevel:
            arrFaces = quadMeshCreator(wOrderFaces, coords)
            arrPos = quadMeshPositions(wOrderFaces, coords)
            percentage = j
            numNodes = math.ceil(len(set([k for item in wOrderFaces for k in item])) * percentage)
            if(numNodes < 3):
                print('NumNodes must be greater than 3 for mesh decimation.')
                result[i] = ['No Data, mesh too small for decimation']
                continue
            startTime = time.time()
            newPositions, newMeshFace = simplify_mesh(arrPos, arrFaces, numNodes, threshold=0.5)
            endTime = time.time() - startTime
            decFaces = convertArrToList(newMeshFace)
            decCoords = convertCoordDict(newPositions)
            print('Writign decimated mesh to obj file')
            order = getVertexOrder(decFaces)
            verts = getVertices(order, decCoords)
            alignDict = createAlignmentDict(verts)
            alignedFaces = realignFaces(decFaces, order, alignDict)
            createOBJ(verts, alignedFaces, 'decimatedMesh_{0}_{0}.obj'.format(i,j))
            results[i] = [endTime, j, len(arrFaces), len(decFaces)]
    return results



#decimate test results {maxNodes: [time, decimationPercent, numfacesPrior, numfacesAfter]}

#endTime = time.time() - startTime

#maxNodeList = [5, 250, 500, 1000, 2000, 4000, 8000, 10000]
#maxNodeList = [10, 20, 30, 40]
#maxNodeList = [50, 500, 5000, 50000,100000,500000]
maxNodeList = [10000,20000,30000,40000,50000]
#REDUCE THE NUMBER OF NODES IN THE MESH BY THIS PERCENTAGE (0,1)
#IF THERE ARE 100 NODES IN THE MESH A PERCENTAGE VALUE OF 0.9 HAS A TARGET REDUCTION OF 10 NODES FOR 90 REMAINING
#percentages = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]


resultDict = test(maxNodeList, bvh)
#print(resultDict)
plotGeneralResults(resultDict)

resultDict2 = duplicateTest(maxNodeList, bvh)
plotDuplicateResults(resultDict2)

resultDict3 = windingOrderTest(maxNodeList, bvh)
plotWindingOrderResults(resultDict3)

resultDict4 = removeOverlapTest(maxNodeList, bvh)
plotOverlapResults(resultDict4)


resultDict5 = aStarTest(maxNodeList, bvh)
plotAstarResults(resultDict5)

resultDict6 = decimateTest(maxNodeList, percentages, bvh)
print(resultDict6)

