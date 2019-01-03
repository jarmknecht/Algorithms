#!/usr/bin/python3
import math

from CS312Graph import *
import time

# Class that keeps track of Node dist and prev and index in heap for fast look up
# All functions are an amortized time complexity of O(1)
# However space is O(|V|) since there are three different dictionaries of length |V| being made
class Dictionaries(object):
    def __init__(self):
        self.dist = {}
        self.prev = {}
        self.indexDict = {}

    def getDistance(self):  #Return whole distance dictionary
        return self.dist

    def setDistance(self, node, distance):  # set node(key) new value(distance)
        self.dist[node] = distance

    def getPrev(self):  #Return whole previous dictionary
        return self.prev

    def setPrev(self, node, prev):  # Set node(key) to new value(previous)
        self.prev[node] = prev

    def getNodeDist(self, node): # get the distance for key node
        return self.dist[node]

    def getNodePrev(self, node): # get the previous for key node
        return self.prev[node]

    def getindexDict(self):  # Return whole index dictionary
        return self.indexDict

    def getNodeIndex(self, node):  # get the index for key node
        return self.indexDict[node]

    def setindexDict(self, node, index):  # set node(key) to new value(index)
        self.indexDict[node] = index


class ArrayQueue:
    def __init__(self, network, dicts):
        self.network = network
        self.dicts = dicts
        self.H = []
        self.makeQueue()

    def decreaseKey(self, index):
        # O(1) for time and space since it is not doing or storing anything
        # Don't need since deletemin finds smallest key
        pass

    def makeQueue(self):
        # O(|V|) for time and space cause it is called |V| times and makes an array of size |V|
        nodes = list(self.network.nodes)
        for x in nodes:
            self.insert(x)

    def insert(self, node):
        # O(1) for time and O(|V|) space; it eventually adds |V| things to the array and does it does it in constant time
        self.H.append(node)

    def deleteMin(self):
        # O(|V|) it looks at the distance dictionary and grabs the smallest value.
        # it can do this at most the size of H which is |V|
        # Space complexity here is O(1), since we aren't storing any new information
        d = self.dicts.getNodeDist(self.H[0]) # O(1)
        index = 0
        for i in range(len(self.H)):
            x = self.dicts.getNodeDist(self.H[i]) # O(1)
            if x < d:
                d = x
                index = i
        return self.H.pop(index) # O(1)

    def getQueue(self):
        # Time and space both O(1) since it is just returning the queue and not adding to it
        return self.H


class HeapQueue:
    def __init__(self, network, dicts):
        self.network = network
        self.dicts = dicts
        self.H = []
        self.makeQueue()

    def decreaseKey(self, node):
        # O(log|V|) is time complexity cause the binary tree's height is at most log|V| which means that bubble up can take
            # this long.
        # Space complexity is O(1) since nothing is being added to the array
        index = self.dicts.getNodeIndex(node) # O(1) time and space
        return self.bubbleUp(index)  # O(log|V|) for time

    def makeQueue(self):
        # Time complexity is O(|V|log|V|) since it runs V times and inserting into a tree is log|V|
        # The space complexity is O(|V|) since we are making a new array of size |V|
        nodes = list(self.network.nodes)  # O(|V|)
        for x in nodes:  # O(|V|)
            self.insert(x)  # O(log|V|)

    def insert(self, node):
        # O(log|V|) is time complexity cause inserting a node into a binary tree can take at most the height of the tree
            # log|V|. this means the newest addition goes from the bottom all the way to the top
        # Space complexity is O(|V|) since that is the size of the heap array made
        childIndex = len(self.H)  # O(1) grab new node index
        self.H.append(node)  # O(1) to append
        self.dicts.setindexDict(node, childIndex)  # O(1) to set new key-value pair
        self.bubbleUp(childIndex)  # O(log|V|)

    def deleteMin(self):
        # Time complexity is O(log|V|) since the furthest a node can bubble down is from the top to the bottom
            # and the height of a binary tree is log|V|
        # Space complexity here is O(1), since we aren't storing any new information
        u = self.H[0]  # Take from the front cause it is highest priority O(1)
        z = self.H[len(self.H) - 1]  #Take the right most node and put it at the top O(1)
        self.H[0] = z  # O(1)
        self.dicts.setindexDict(self.H[0], 0)  #O(1)
        self.H.pop(len(self.H) - 1)  # Pop off the rightmost node since it is a duplicate now O(1)
        self.bubbleDown(0) # O(log|V|)
        return u

    def getQueue(self):
        # Time and space both O(1) since it is just returning the queue and not adding to it
        return self.H

    def bubbleUp(self, childIndex):
        # Time complexity here is O(log|V|) since the longest time bubble up a node can be from the bottom
            # to the top of the heap which has a height of log|V|
        # Space complexity is O(1) since we aren't adding or taking away from the array
        if childIndex == 0:  # Base case we are the root node now
            return
        parentIndex = ((childIndex - 1) // 2)  # O(1)

        x = self.dicts.getNodeDist(self.H[parentIndex])  # O(1)
        y = self.dicts.getNodeDist(self.H[childIndex]) # O(1)
        if x > y:
            self.switch(parentIndex, childIndex) # O(1)
            return self.bubbleUp(parentIndex)  # O(log|V|)
        else:
            return

    def switch(self, parentIndex, childIndex):
        # Time and space complexities are O(1) since switches happen in constant time and no new memory taken up
        temp = self.H[parentIndex]  # O(1)
        self.H[parentIndex] = self.H[childIndex]  # O(1)
        self.H[childIndex] = temp  # O(1)
        self.dicts.setindexDict(self.H[parentIndex], parentIndex)  # O(1)
        self.dicts.setindexDict(self.H[childIndex], childIndex)  # O(1)
        return

    def bubbleDown(self, parentIndex):
        # Time complexity is O(log|V|) since the longest time to bubble down a node can be from the top
               # to the bottom of the heap which has a height of log|V|
        # Space complexity is O(1) since we aren't adding or taking away from the array
        if parentIndex == (len(self.H) - 1):  # Base case shifted the node all the way to the end  O(1)
            return
        leftChildIndex = (2 * parentIndex) + 1  # O(1)
        rightChildIndex = (2 * parentIndex) + 2 # O(1)

        if leftChildIndex < len(self.H): # This means it has a left child # O(1)
            if self.dicts.getNodeDist(self.H[parentIndex]) > self.dicts.getNodeDist(self.H[leftChildIndex]): # O(1)
                if rightChildIndex < len(self.H): # Means it has a right child and need to switch with smallest distance # O(1)
                    if self.dicts.getNodeDist(self.H[leftChildIndex]) <= self.dicts.getNodeDist(self.H[rightChildIndex]):
                        # This means the left child is smaller or equal to right so just switch with left kid
                        self.switch(parentIndex, leftChildIndex) # O(1)
                        return self.bubbleDown(leftChildIndex) # O(log|V|)
                    else: # Right kid is smaller
                        self.switch(parentIndex, rightChildIndex) # O(1)
                        return self.bubbleDown(rightChildIndex) # O(log|V|)

        # no left child check if there is a right one
        if rightChildIndex < len(self.H): # O(1)
            if self.dicts.getNodeDist(self.H[parentIndex]) > self.dicts.getNodeDist(self.H[rightChildIndex]): # O(1)
                self.switch(parentIndex, rightChildIndex) # O(1)
                return self.bubbleDown(rightChildIndex) # O(log|V|)

class NetworkRoutingSolver:
    def __init__(self):
        self.dicts = Dictionaries()
        pass

    def initializeNetwork( self, network ):
        assert( type(network) == CS312Graph )
        self.network = network

    # Time complexity here is O(|V|) since the for loop isn't always gone
        #  through in the while loop cause of the sparsity of the graph so we can add |E| and |V| together as seen below
    # Space complexity is O(|V|) since we are making an array of path_edges and worst case would be to have all the nodes
        # a part of the path
    def getShortestPath( self, destIndex ):
        self.dest = destIndex
        path_edges = []
        nodes = self.network.getNodes() # O(|V|)
        endNode = nodes[self.dest]  # O(1)
        total_length = self.dicts.getNodeDist(endNode)  # O(1)
        currentNode = endNode  # O(1)
        # Overall time here is O(|V| + |E|) = O(|V|)
        while self.dicts.getPrev().get(currentNode) is not None:  # O(|V|)
            nextNode = self.dicts.getPrev().get(currentNode)  # O(1)
            for edge in nextNode.neighbors:  # O(|E|) = O(|V|)
                if (edge.dest == currentNode):  # Only add edge if it is one we want O(1)
                    path_edges.append( (nextNode.loc, currentNode.loc, '{:.0f}'.format(edge.length)) )  #O(1) and space O(|V|)
            currentNode = nextNode # O(1)
        return {'cost':total_length, 'path':path_edges}

    # Compute shortest paths has a time complexity of
        # O(|V^2|) if priority queue is an array and
        # O(|V|log|V|) if the priority queue is a heap and the graph is sparse
            # otherwise the O(|V^2|log|V|) since |E| becomes equal to |V^2| not |V| anymore
    # Space complexity is O(|V| + 3|V|) = O(|V|) since this is the size of the queue(|V|) and dictionaries(3|V|) made
    def computeShortestPaths(self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        # Run Dijkstra's
            # goes through the all the nodes and sets dist for a node to be inf
            # and prev to be nil O(|V|) for time O(|V^2|) for space cause two dictionaries are being stored
        for x in self.network.nodes:  # O(|V|)
            self.dicts.setDistance(x, math.inf)  # O(1)
            self.dicts.setPrev(x, None)  # O(1)
        nodes = self.network.getNodes()  # O(1)
        self.dicts.setDistance(nodes[self.source], 0)  # O(1)

        # This makeQueue runs in O(|V|) if it is an array and O(|V|log|V|) if it is a heap
        if use_heap:
            H = HeapQueue(self.network, self.dicts)  # O(|V|log|V|)
        else:
            H = ArrayQueue(self.network, self.dicts)  # O(|V|)
        # O(|E|) = O(|V|) cause of the predefined maxiumum out-degree then:
        # If queue is an array this runs O(|V^2| + |V|) = O(|V^2|)
            # If queue is a heap runs in O(|V|log|V| + |V|log|V|) which is just O(|V|log|V|)
            # no new memory taken up so space is O(1)
        while len(H.getQueue()) > 0:  # O(|V|) for both
            u = H.deleteMin()  # O(|V|) for array O(log|V|)
            for n in u.neighbors: # O(|V|) for array O(|V|log|V|) for heap since |E| = |V|
                v = n.dest  # u and v are neighbors, they have an edge together
                if self.dicts.getNodeDist(v) > (self.dicts.getNodeDist(u) + n.length):  # O(1)
                    self.dicts.setDistance(v, (self.dicts.getNodeDist(u) + n.length))  # O(1)
                    self.dicts.setPrev(v, u)  # O(1)
                    H.decreaseKey(v)  # O(1) for array O(log|V|) for heap
        t2 = time.time()
        return (t2-t1)

