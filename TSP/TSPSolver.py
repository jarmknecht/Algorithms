#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import copy
import itertools


class TSPSolver:
	def __init__(self, gui_view):
		self._scenario = None
		self.nodeQueue = []
		self.numSols = 0
		self.totalStates = 0
		self.bssf = None
		self.currState = None
		self.prunedStates = 0
		self.maxQSize = 0
		self.numStatesAddedToQ = 0
		self.cities = 0  # Stores the array of x y coords for each city
		self.cost = 0
		self.numBSSFUpdates = 0

	def setupWithScenario(self, scenario):
		self._scenario = scenario

	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	# O(n!) time since it runs till it finds a correct tour. This is because the worst case would be the number of all city path permutations
	# which is O(n!). However the average case is more likely closer to O(n) because it most likely finds a random tour without having to go
	# through all of the permutations.
	# The space complexity here is O(n) since we are storing a route of length n + 1 to
	# get back to the start city but the one is dropped.
	def defaultRandomTour(self, time_allowance=60.0):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time() - start_time < time_allowance:  # O(n!) this may have to
			# create a random permutation
			perm = np.random.permutation(ncities)
			route = []
			# Now build the route using the random permutation
			for i in range(ncities):  # O(n)
				route.append(cities[perm[i]])
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy(self, time_allowance=60.0):
		pass

	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
	# The O(n^3 *(n-1)!) is the worst case time complexity if no pruning happened. This is because the queue can have up to
	# all the states on it except the first one made. However due to pruning this is not what the big O averages out to be
	# closer to O(n^3) because of my priority key making the search look at states that have the a lower cost and are found deeper
	# within the tree. Thus the O(n^3) comes from O(n^2 * b^n) which since I travel down one state and try to get to the bottom of the tree
	# right away becomes O(n^2 * n) since n the most children we can have on the first level and it goes down each level. So the O(n^3) comes
	# from having to make all of the child nodes at the beginning cause there is less and less work going down the tree.
	# The space complexity is O(n^2) from making the 2d array
	def branchAndBound(self, time_allowance=60.0):
		results = {}  # O(1)
		self.cities = self._scenario.getCities()[:]  # Copies over x,y 2d array for each city so O(n^2)
		self.numSols = 0 # O(1)
		self.totalStates = 0 # O(1)
		self.prunedStates = 0 # O(1)
		self.nodeQueue = [] # O(1)
		self.maxQSize = 0 # O(1)
		self.numBSSFUpdates = 0 # O(1)
		self.numStatesAddedToQ = 0 # O(1)

		defaultResults = self.defaultRandomTour(time.clock())  # Grab the default route to use as initial BSSF O(n)
		self.bssf = defaultResults['soln']  # Take what is stored in the dictionary of defaultResults at soln key and use as bssf
		print("Initial BSSF: ", self.bssf._costOfRoute()) # Checks that we got a BSSF
		start_time = time.time() # Record the start time so it can be compared with end time
		self.cost = self.bssf._costOfRoute()  # Stores initial cost of BSSF in cost
		city1 = copy.deepcopy(self.cities[0])  # Copy info of first city into city 1 O(n) deepcopy so it isn't a pointer

		city1.setCostMatrix(
			self.reduceMatrix(CostMatrix(self.cities[:])))  # O(n^2) since it copies n cities into a 2d array of size n

		# O(log n) is time complexity cause inserting a node into a binary tree can take at most the height of the tree
		# log n. Where n is the number of states pushed onto the heap.
		# Space complexity of the heap is O((n - 1)!) since that is the worst case size of the heap
		# this means that all of the states except the first city chosen are on the heap
		heapq.heappush(self.nodeQueue, city1)

		self.numStatesAddedToQ += 1  # State city1 was added to the queue so increment

		# This could at the worst case go for O(n^3 * (n - 1)!) where n is the number of states on the queue which could be all of them
		#  but the first city state. However due to pruning on average it will be much better. Closer to O(n^3)
		while (len(self.nodeQueue) > 0) and ((time.time() - start_time) < time_allowance):
			self.maxQSize = max(len(self.nodeQueue), self.maxQSize) # O(1)
			self.lookAtState(heapq.heappop(self.nodeQueue))  # Takes about O(n^3) time
			self.maxQSize = max(len(self.nodeQueue), self.maxQSize)  # O(1)

		end_time = time.time() # O(1)
		results['cost'] = self.bssf._costOfRoute()  # O(1)
		results['time'] = end_time - start_time  # O(1)
		results['count'] = self.numSols  # O(1)
		results['soln'] = self.bssf  # O(1)
		results['max'] = self.maxQSize  # O(1)
		results['total'] = self.totalStates  # O(1)
		results['pruned'] = self.prunedStates  # O(1)
		return results

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''

	def fancy(self, time_allowance=60.0):
		pass

	# As time goes on the time complexity of this averages out to be O(n^2) the space complexity here is O(1) since we
	# are only updating values
	def reduceMatrix(self, costMatrix):
		cost = costMatrix.getCost()  # O(1)
		numRows = len(self.cities)   # O(1)
		numCols = numRows  # O(1)

		# Reduces rows first worst case O(n^3) if have to go through and reduce each entry but some will
		# be inf or 0 as time goes on so it comes down to averaging to be about O(n^2)
		for row in range(0, numRows):  # O(n)
			minValue = math.inf  # Store what the min value is for each row  O(1)
			for col in range(0, numCols):  # O(n)
				matrixValue = (costMatrix.getMatrix())[col][row]  # O(1)
				minValue = min(minValue, matrixValue)  # Take the min of the two  O(1)

			if minValue > 0 and minValue < math.inf:  # Will be entered in less and less down the tree so it's for loop isn't counted in overall time complexity
				cost += minValue  # Add minvalue to the total bound cost O(1)
				# Changes each entry in a column by the minValue found
				for col in range(0, numCols):  # O(n)
					matrixValue = (costMatrix.getMatrix())[col][row]  # O(1)
					costMatrix.setMatrixValue(col, row, matrixValue - minValue)  # O(1)

		# Reduces cols second wosrtcase O(n^3) if have to go through and reduce each entry but some will
		# be inf or 0 as time goes on so it comes down to averaging to be about O(n^2)
		for col in range(0, numCols):  # O(n)
			minValue = math.inf  # stores min value for each column O(1)
			for row in range(0, numRows):  # O(n)
				matrixValue = (costMatrix.getMatrix())[col][row]  # O(1)
				minValue = min(minValue, matrixValue)  # Take min of the two  # O(1)

			if minValue > 0 and minValue < math.inf:  # Will be entered in less and less down the tree so it's for loop isn't counted in overall time complexity
				cost += minValue  # Add minvalue to the total bound cost O(1)
				for row in range(0, numRows):  # O(n)
					matrixValue = (costMatrix.getMatrix())[col][row]  # O(1)
					costMatrix.setMatrixValue(col, row, matrixValue - minValue)  # O(1)
		costMatrix.setCost(cost)  # O(1)
		return costMatrix

	# Overall the time complexity here is O(n^3) with space complexity of O(n^2) since it makes new child state matrix
	def lookAtState(self, currState):
		self.currState = currState  # O(1)

		# With the for loop being O(n) and creating a child being O(n^2) overall code here is O(n^3)
		if (self.contExploring()):  # Returns true if should still explore state's kids as time goes on will be true less and less from the bound averaging to O(1)
			for childCity in range(0, len(self.cities)):  # O(n)
				if ((self.currState.getCostMatrix()).getMatrix())[childCity][self.currState.getIndex()] < math.inf:  # O(1) this just checks that the matrix entry is less than inf
					childState = self.makeChild(childCity)  # making a child is O(n^2) cause of 2d array same with space
					self.totalStates += 1  # made a new state
					if childState != None: # will be None if the cost of the child is greater than the BSSF so it is not always entered cause of pruning
						# O(log n) is time complexity cause inserting a node into a binary tree can take at most the height of the tree
						# log n. Where n is the number of states pushed onto the heap.
						# Space complexity of the heap is O((n - 1)!) since that is the worst case size of the heap
						# this means that all of the states except the first city chosen are on the heap
						heapq.heappush(self.nodeQueue, childState)

						self.numStatesAddedToQ += 1  # inc num states added  # O(1)
						self.maxQSize = max(len(self.nodeQueue), self.maxQSize) # O(1)

	# Overall time complexity averages to be O(1) since as time goes on more and more states will be pruned as BSSF is updated
	# This makes it so the first if statement is only looked at. Space complexity is O(n^2) when we have to copy over the 2d array
	# however this also averages to O(1) for the same reason as time complexity
	def contExploring(self): # Returns true if the cost is less than BSSF otherwise it is pruned
		matrix = self.currState.getCostMatrix()  # O(1)
		currMatrix = self.reduceMatrix(matrix)  # O(n^2)

		self.currState.setCostMatrix(currMatrix)  # O(1)

		if currMatrix.getCost() > self.bssf._costOfRoute():  # O(1)
			self.prunedStates += 1  # O(1)
			return False  # Don't need to look at this states kids

		if len(self.currState.getPath()) == len(self.cities) + 1:  # we found a cycle and can update BSSF possibly O(1)
			path = self.currState.getPath()[:]  # O(n^2) for copying 2d array
			if (path[0].getName() == path[-1].getName()):  # the first and last city are the same and it is a cycle O(1)
				self.numSols += 1  # found a solution O(1)
				if (self.currState.getCostMatrix()).getCost() <= self.cost:  # O(1)
					self.numBSSFUpdates += 1  # O(1)
					self.bssf = TSPSolution(self.currState.getPath()[:-1])  # Get the path O(1)
					self.cost = (self.currState.getCostMatrix()).getCost()  # set the cost of the path O(1)
			return False  # found a leaf node
		return True

	# Overall time and space complexity are O(n^2) this comes from copying over the 2d array from the parent to the child made
	def makeChild(self, childCity):
		childName = nameForInt(childCity + 1)  # O(1)
		currLevel = len(childName)  # O(1) this checks where in the tree we are based on the length of the name

		cost = (self.currState.getCostMatrix()).getCost() + ((self.currState.getCostMatrix()).getMatrix())[childCity][
			self.currState.getIndex()]  # O(1)

		if cost > self.bssf._costOfRoute():  # O(1)
			self.prunedStates += 1  # O(1)
			return None  # state was pruned and not needed checking this now saves overall space in the priority queue
		# O(1) by setting the key as the lower bound cost and the
		# level it is on ensures that smaller costs and those closer to the bottom have higher priority
		key = (cost / currLevel)

		childMatrix = CostMatrix(self.cities[:])  # O(n) copies over cities
		childMatrix.setMatrix(
			[row[:] for row in (self.currState.getCostMatrix()).getMatrix()])  # time and space both O(n^2) copying 2d array
		childMatrix.setCost(cost)  # O(1)
		for i in range(0, len(self.cities)):  # O(n) sets things to inf that need to be inf
			childMatrix.setMatrixValue(i, self.currState.getIndex(), math.inf)  # O(1)
			childMatrix.setMatrixValue(childCity, i, math.inf)  # O(1)
		childMatrix.setMatrixValue(self.currState.getIndex(), childCity, math.inf)  # O(1)
		# Copies parent(childCity) state's stuff into childNode O(n^2) for copying 2d array time and space complexity
		childNode = copy.deepcopy((self.cities[:])[childCity])

		childNode.setCostMatrix(childMatrix)  #O(1)
		childNode.setIndexAndName(childCity, childName)  #O(1)
		childNode.setKey(key)  # O(1)
		childNode.setPath(self.currState.getPath())  # O(1)
		childNode.addToPath(childNode)  # O(1)
		return childNode