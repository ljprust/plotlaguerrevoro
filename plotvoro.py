import yt
import matplotlib
matplotlib.rc("text", usetex=True)
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import sys
#sys.path.append('~/code/plotvoronoi')
from laguerre import *

nfiles = 1
first = 140
interval = 100
nref = 64
do_marks = False
filename = []
for i in range(0, nfiles) :
    num = 1000000 + first + i * interval
    numstr = str(num)
    cut = numstr[1:7]
    filename.append(cut)

Rsun = 7.0e10
Rpl  = 0.1*Rsun
RA   = 3.0*Rpl
#dPeriod = 3.5e14
#lim = dPeriod / 2. * 1.0001
#hbox = np.array([[-lim,lim],[-lim,lim],[-lim,lim]])

for file in filename :
    print 'starting dataset ' + file
    plt.clf()

    ds = yt.load('agb.'+file, n_ref=nref) #, bounding_box=hbox )
    ad = ds.all_data()
    cl = ds.arr(1.0, 'code_length')
    #time = ds.current_time
    #mycenter = [0., cmvel*time, 0.]
    #timestr = str(time.in_units('day'))[0:5]
    #plot.annotate_text( (0.02,0.02), 't = ' + timestr + ' ' + 'd', coord_system='axis' )

    #mirrorFile = 'wind.bc.' + file
    #mirrorData = np.loadtxt(mirrorFile)

    #mirrorRadius = mirrorData[4]
    #mirrorCenter = np.array([mirrorData[5], mirrorData[6], mirrorData[7]])
    #mirrorVel    = np.array([mirrorData[9], mirrorData[10], mirrorData[11]])

    mirrorCenter = np.array([-RA,0.,0.])
    posGas = (ad[('Gas','Coordinates')]/cl-mirrorCenter)/Rpl
    x = posGas[:,0]
    y = posGas[:,1]
    z = posGas[:,2]

    zShift = -1.9
    z = z-zShift
    circleRadius = np.sqrt(1.0-zShift*zShift)

    cutWidth = 5.0 #*Rsun # 10.0*Rsun
    cutThick = 0.5*cutWidth #*Rsun # 8.0*Rsun
    xBounds = [ -cutWidth/2.0, cutWidth/2.0 ]
    yBounds = [ -cutWidth/2.0, cutWidth/2.0 ]
    zBounds = [ -cutThick/2.0, cutThick/2.0 ]
    xbool  = np.logical_and(x>xBounds[0],x<xBounds[1])
    ybool  = np.logical_and(y>yBounds[0],y<yBounds[1])
    zbool  = np.logical_and(z>zBounds[0],z<zBounds[1])
    inbool = np.logical_and( np.logical_and(xbool,ybool), zbool )
    numIn = inbool.sum()
    print('number of particles = ',numIn)

    posGasIn = posGas[inbool,:]
    posGasPlane = posGasIn[:,0:2]
    #inBoundary = mirrorRadius*mirrorRadius > x[inbool]*x[inbool] + y[inbool]*y[inbool]
    zIn = z[inbool]

    #zIn = zIn/Rsun/10.0
    #posGasPlane = posGasPlane/Rsun/10.0

    #print(posGasPlane)
    #print(zIn)

    #posGasPlane = numpy.stack((x,y),axis=1)
    z = zIn

    zmax = numpy.absolute(z).max()
    #print('zmax = ',zmax)
    #print('x,y max = ',posGasPlane.max())
    K = zmax*zmax + 1.0
    S = posGasPlane
    R = numpy.sqrt( K-z*z )
    #print(R)

    print('finding power triangulation')
    tri_list, V = get_power_triangulation(S, R)
    print('getting voronoi cell map')
    voronoi_cell_map = get_voronoi_cells(S, V, tri_list)
    print('displaying results')
    display(S, R, tri_list, voronoi_cell_map, xBounds, yBounds, circleRadius, zShift)

    '''
    vor = Voronoi(posGasPlane)
    fig = voronoi_plot_2d(vor, show_vertices = False, show_points = False)

    colors = np.empty(numIn, dtype='string')
    for i in range(0,numIn) :
        colors[i] = 'w'
    colors[inBoundary] = 'r'

    for j in range(0,numIn):
        region = vor.regions[vor.point_region[j]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            xPoly = []
            yPoly = []
            for k in range(0,len(polygon)) :
                xPoly.append(polygon[k][0])
                yPoly.append(polygon[k][1])
            plt.fill( xPoly, yPoly, colors[j] )

    axes = plt.gca()
    axes.set_aspect(1)

    if do_marks:
        ad = ds.all_data()
        dm_pos = ad[("DarkMatter","Coordinates")]
        core = dm_pos[0][:]
        comp = mirrorCenter
        plot.annotate_marker( core, coord_system = 'data', plot_args={'color':'black'}, marker = '+')
        plot.annotate_marker( comp, coord_system = 'data', plot_args={'color':'green','s':25}, marker = 'x')

    saveas = 'laguerre' + file + '.pdf'
    plt.savefig(saveas)
    print('saved figure ' + saveas)
    plt.clf()
    '''
