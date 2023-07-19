import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy.ma as ma
import math
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded, MapScatter

fn = '/Users/rohanchawla/Desktop/Ice Downloads/rtofs_glo.t00z.n00.cice_inst'
ds = nc.Dataset(fn) 
icedata = ds["aice"]
longitude = ds["TLON"]
latitude = ds["TLAT"]
adjacentEdges = {}
lons = []
lats = []
def main():
    # Create global map with no data using
    # PlateCarree projection and coastlines
    plot1 = CreatePlot()
    plot1.projection = 'plcarr'
    plot1.domain = 'global'
    plot1.add_map_features(['coastline', 'land', 'ocean'])
    plot1.add_xlabel(xlabel='longitude')
    plot1.add_ylabel(ylabel='latitude')
    

    fig = CreateFigure()
    fig.plot_list = [plot1]
    fig.create_figure()

    doDijkstras()

    plt.plot(lons, lats)

    plt.show()



def doDijkstras():
    calculateEdges()

    used = []
    shortest_root = {}
    shortest_cost = {}

    #Puts all of the points costsw at infinity and shortest roots to none
    for i in range (0, 3297):
        for j in range(0, 4500):
            shortest_root[(i, j)] = None
            shortest_cost[(i, j)] = float('inf')
    #Sets the cost of the starting point to 0
    shortest_cost[(0, 0)] = 0

    
    while (3296, 4499) not in used and len(used) < (3296*4499):
        min_val = float('inf')
        min_cost = None
        for min_vertex in shortest_cost:
            if min_vertex not in used and shortest_cost[min_vertex] < min_val:
                min_cost = min_vertex
                min_val = shortest_cost[min_vertex]
        if min_cost:
            used.append(min_cost)

        for vertex in adjacentEdges[min_cost]:
            if vertex not in used:
                old_cost = shortest_cost[vertex]
                new_cost = calculateCost(min_cost, vertex) + shortest_cost[min_cost]
                if old_cost > new_cost:
                    shortest_cost[vertex] = new_cost
                    shortest_root[vertex] = min_cost
        else:
            break
    
    reorder = []
    curr = (3296, 4499)
    reorder.append(curr)

    while shortest_root[curr] != None:
        curr = shortest_root[curr]
        reorder.append(curr)

    for i in range (0, len(reorder)):
        lons[i] = ma.getdata(longitude[0, reorder[i][1]])
        lats[i] = ma.getdata(latitude[reorder[i][0], 0])

    
#Calculates all the edge values for each point
def calculateEdges():
    #Nested for loop going through each of the points
    for i in range (3297-1):
        for j in range (4500-1):
            #If this is not the last row of cells, the above, diagonal, and right cells are the adjacent cells, else, only the above cell is adjacent
            #if(j < 4499 and i < 3296):
            adjacentEdges[(i, j)] = {(i+1, j) : calculateCost(i, j, i+1, j), (i+1, j+1) : calculateCost(i, j, i+1, j+1), (i, j+1) : calculateCost(i, j, i, j+1)}    
           # print("This is the i value:", i, "This is the j value:", j)
        if (i % 10 == 0):
            print(i)
            #elif(i == 3296 and j < 4499):
    for i in range(3295, 3296):
        for j in range(4499):
            dict = {(i, j+1) : calculateCost(i, j, i, j+1), (i-1, j) : calculateCost(i-1, j, i, j), (i-1, j+1) : calculateCost(i-1, j, i, j+1)}
            adjacentEdges[(i, j)] = dict
            #print("This is the i value:", i, "This is the j value:", j)
            #else:
    for i in range(3296):
        for j in range(4498, 4499):
            dict = {(i+1, j) : calculateCost(i, j, i+1, j)}
            adjacentEdges[(i, j)] = dict
           # print("This is the i value:", i, "This is the j value:", j)
#Calculates the cost of an edge between two points by the formula: distance + (distance * ice area)
def calculateCost(i1, j1, i2, j2):
    distance = calculate_distance(ma.getdata(latitude[i1, 0]), ma.getdata(longitude[0, j1]), ma.getdata(latitude[i2, 0]), ma.getdata(longitude[0, j2]))
    val = distance + distance * ma.getdata(icedata[0, i2, j2])
    print(val)
    return distance + distance * ma.getdata(icedata[0, i2, j2])


#Calculates the distance of two points based on the longitude and latitude points of each point
def calculate_distance(lat1, lon1, lat2, lon2):
    # Radius of the Earth in kilometers
    earth_radius = 6371

    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Calculate the distance
    distance = earth_radius * c

    return distance

if __name__ == '__main__':
    main()