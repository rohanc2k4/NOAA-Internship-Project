import numpy as np
import networkx as nx
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy.ma as ma
import math
from math import *
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded, MapScatter

fn = '/Users/rohanchawla/Desktop/Ice Downloads/rtofs_glo.t00z.n00.cice_inst'
ds = nc.Dataset(fn) 
icedata = ds["aice"]
longitude = ds["TLON"]
lons = ds.variables["TLON"][:,:]
lats = ds.variables["TLAT"][:,:]
latitude = ds["TLAT"]
aice = ds.variables["aice"][0,:,:]
hi = ds.variables["hi"][0,:,:]
nj = len(ds.dimensions["nj"])
ni = len(ds.dimensions["ni"])
#f = open('/Users/rohanchawla/Desktop/readme.txt', 'w')
#faster (~50 seconds) to used masked arrays than doubly nested loop (250 seconds)
lmask = ma.masked_array(lons > 2.*360.+180.)
lin = lmask.nonzero()
for k in range(0, len(lin[0])):
    i = lin[1][k]
    j = lin[0][k]
    print(ma.getdata(i),ma.getdata(j))
    lons[j,i] -= 3.*360.

lmask = ma.masked_array(lons > 1.*360.+180.)
lin = lmask.nonzero()
for k in range(0, len(lin[0])):
    i = lin[1][k]
    j = lin[0][k]
    print(ma.getdata(i),ma.getdata(j))
    lons[j,i] -= 2.*360.

#most (10.6 million of 14.7 million) rtofs points have lons > 180, so subtract 360 and 
# then correct the smaller number that are < -180 as a result
lons -= 360.

lmask = ma.masked_array(lons < -180.)
lin = lmask.nonzero()
for k in range(0, len(lin[0])):
    i = lin[1][k]
    j = lin[0][k]
    lons[j,i] += 1.*360.
def main():
    #Uses the user's input to find the polar class of the ship
    PC = int(input("What is the polar class of the ship vessel? (1-7)\n"))
    PossAnswers = [1, 2, 3, 4, 5, 6, 7]
    if(PC not in PossAnswers):
        raise Exception("Please select an answer between 1 and 7.") 

    # Create global map with no data using
    # PlateCarree projection and coastlines
    plot1 = CreatePlot()
    plot1.projection = 'plcarr'
    plot1.domain = 'global'
    plot1.add_map_features(['coastline', 'land', 'ocean'])
    plot1.add_xlabel(xlabel='longitude')
    plot1.add_ylabel(ylabel='latitude')
    
    #Sets the limit of the graph's longitude and latitude values
    plot1.set_xlim(-180, 0)
    plot1.set_ylim(30, 90)
    fig = CreateFigure()
    fig.plot_list = [plot1]
    fig.create_figure()
    
    lons1 = [180, 170, 170, 180]
    lats1 = [65, 65, 85, 85]

    lons2 = [-180, -50, -50, -180]
    lats2 = [85, 85, 65, 65]
    plt.plot(lons1, lats1, 'r')
    plt.plot(lons2, lats2, 'r')

    (beringnj, beringni) = find(lons, lats, -168.59, 65.68)
    (baffinnj, baffinni) = find(lons, lats, -74.0, 74.0)
    # Construct nodes -- limit area to keep run time manageable:
    latmin = 65.0
    latmax = 88.0
    #lonmin = 185.0-360.
    #lonmax = 290.0-360.
    lonmin = -175.0
    lonmax =  -70.0
    xmask = ma.masked_outside(lons, lonmin, lonmax)
    xin = xmask.nonzero()
    #print(len(xin), len(xin[0]))
    xmask = ma.logical_and(xmask, ma.masked_outside(lats, latmin, latmax))
    xin = xmask.nonzero()
    #print(len(xin), len(xin[0]))

    xmask = ma.logical_and(xmask, aice < 1000.)
    xin = xmask.nonzero()
    print(len(xin), len(xin[0]))



    #Not a directed graph
    G = nx.Graph()
    
    nodemap = np.full(shape = (nj, ni),fill_value=-1,dtype="int")
    for k in range(0, len(xin[0])):
        i = xin[1][k]
        j = xin[0][k]
        if (k%15000 == 0):
            print("adding nodes, k = ",k, flush=True)
  #debug print("node:",k,i,j,lats[j,i], lons[j,i], sst[j,i], flush=True)
        nodemap[j,i] = int(k)
        G.add_node(k, i = i, j =j, lat = lats[j,i], lon = lons[j,i], aice=aice[j,i], hi=hi[j,i])
    print("Done adding nodes, k=",k, flush=True)
    #print(G.nodes(data=True))
    distance = 0
    #For loop that goes through every cell value and sets an edge cost of each edge based on the calculateCost function
    for k in range(0, len(xin[0])):
        i = xin[1][k]
        j = xin[0][k]
        jp = j + 1
        jm = j - 1
        ip = i + 1
        im = i - 1
        n = nodemap[j,i]
        if (n == -1):
            continue

        if (im >= 0):
            if (nodemap[j,im] != -1):
               G.add_edge(n, nodemap[j,im], weight = calculateCost(G.nodes[n], G.nodes[nodemap[j, im]], PC))

        if (ip < ni):
            if (nodemap[j,ip] != -1):
                G.add_edge(n, nodemap[j,ip], weight = calculateCost(G.nodes[n], G.nodes[nodemap[j, ip]], PC))

        if (jp < nj ):
            if (nodemap[jp,i] != -1):
                G.add_edge(n, nodemap[jp,i], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jp, i]], PC))

            if (im >= 0):
                if (nodemap[jp,im] != -1):
                    G.add_edge(n, nodemap[jp,im], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jp, im]], PC))
                    
            if (ip < ni):
                if (nodemap[jp,ip] != -1):
                    G.add_edge(n, nodemap[jp,ip], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jp, ip]], PC))
                    
  #RG: a guess about the archipelago seam
        else:
            tmpi = i
            if (i < ni/2-1):
                tmpi = ni - 1 - i
            if (nodemap[j,tmpi] != -1):
                G.add_edge(n, nodemap[j,tmpi], weight = calculateCost(G.nodes[n], G.nodes[nodemap[j, tmpi]], PC))
                
        if (jm >= 0 ):
            if (nodemap[jm,i] != -1):
                G.add_edge(n, nodemap[jm,i], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jm, i]], PC))
                
            if (im >= 0):
                if (nodemap[jm,im] != -1):
                    G.add_edge(n, nodemap[jm,im], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jm, im]], PC))
                    
            if (ip < ni):
                if (nodemap[jm,ip] != -1):
                    G.add_edge(n, nodemap[jm,ip], weight = calculateCost(G.nodes[n], G.nodes[nodemap[jm, ip]], PC))
                    
    print("Have constructed graph, number of edges =",k, len(G.edges), flush=True)
    
    #Creates a dijkstra_path from the start node to the end node
    FoundPath = nx.dijkstra_path(G, nodemap[beringnj, beringni], nodemap[baffinnj, baffinni])
    lonVals = []
    latVals = []
    for i in range (0, len(FoundPath)):
        print(G.nodes[FoundPath[i]])
    
    #Finds the distance of all the nodes and adds them together, puts all the latitude values and longitude values in a list to be plotted
    for i in range (0, len(FoundPath)):
        if(i != len(FoundPath)-1):
            distance += calculate_distance(G.nodes[FoundPath[i]]['lat'], G.nodes[FoundPath[i]]['lon'], G.nodes[FoundPath[i+1]]['lat'], G.nodes[FoundPath[i+1]]['lon'])
        lonVals.append(G.nodes[FoundPath[i]]['lon'])
        latVals.append(G.nodes[FoundPath[i]]['lat'])

    print("Total Distance Travelled: ", str(round(distance, 2)), "KM")
    plt.plot(lonVals, latVals, marker = ".", markersize = 1, color = "purple")
    plt.show()

def find(lons, lats, lonin, latin):
  #debug print("lon, lat in:",lonin, latin, flush=True)
  tmpx = lons - lonin
  tmpy = lats - latin
  #debug print("x ",tmpx.max(), tmpx.min(), lons.max(), lons.min(), flush=True )

  xmask = ma.masked_outside(tmpx, -0.5, 0.5)
  xin = xmask.nonzero()
  wmask = ma.logical_and(xmask, ma.masked_outside(tmpy, +0.5, -0.5) )
  win = wmask.nonzero()

  imin = -1
  jmin = -1
  dxmin = 999.
  dymin = 999.
  dmin  = 999.
  for k in range(0, len(win[0]) ):
    i = win[1][k]
    j = win[0][k]
    #debug print(k,i,j,abs(tmpx[j,i]), abs(tmpy[j,i]), dxmin, dymin, dmin, flush=True)
    #if (abs(tmpx[j,i]) < dxmin and abs(tmpy[j,i]) < dymin):
    if (sqrt(tmpx[j,i]**2 + tmpy[j,i]**2) < dmin):
      imin = i
      jmin = j
      dxmin = abs(tmpx[j,i])
      dymin = abs(tmpy[j,i])
      dmin  = sqrt(tmpx[j,i]**2 + tmpy[j,i]**2)
  #print("dmin:",imin, jmin, dmin, dxmin, dymin)
  return (jmin,imin)

#Calculate Cost function that takes in the ice concentration, ice thickness, polar class, and distance
#and outputs a cost from the current node to the next node based off of inputs
def calculateCost(currNode, nextNode, PolarClass):
    #RIO = (aice*10)RV
    #If aice <= .1, return 0
    #If RIO < 0, return 99999
    iceCon = currNode["aice"]
    iceThick = nextNode["hi"] * 100
    cost = calculate_distance(currNode["lat"], currNode["lon"], nextNode["lat"], nextNode["lon"])
    #Considered Ice-Free
    if(iceCon <= .1):
        return cost
    #Costs are defined based off of the ice thickness that each ship can handle signified by its Polar Class
    if(PolarClass == 1 or PolarClass == 2 or PolarClass == 3 or PolarClass == 4):
        if(iceThick <= 70):
            cost = cost + (2 * (iceCon * 10))
        elif(iceThick <= 120):
            cost = cost + (4 * (iceCon * 10))
        else:
            cost = cost + (6 * (iceCon * 10))

    elif(PolarClass == 5 or PolarClass == 6):
        if(iceThick <= 70):
            cost = cost + (2 * (iceCon * 10))
        elif(iceThick <= 95):
            cost = cost + (4 * (iceCon * 10))
        elif(iceThick <= 120):
            cost = cost + (6 * iceCon * 10)
        else:
            return 999
        
    else:
        if(iceThick <= 30):
            cost = cost + (2 * (iceCon * 10))
        elif(iceThick <= 50):
            cost = cost + (4 * (iceCon * 10))
        elif(iceThick <= 70):
            cost = cost + (6 * (iceCon * 10))
        else:
            return 999
    return cost

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