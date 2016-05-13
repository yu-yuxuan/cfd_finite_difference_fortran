#!/usr/bin/python2.7

"""
This is codes to plot the results from cfd
"""
import sys
import optparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def two_d_plot(file,show):
    data=np.loadtxt(file)
    x=data[:,0]
    y=data[:,1]
    fig = plt.figure()
    plt.plot(x,y,'b-')  # connect points with a blue line
    plt.xlim(x.min()-0.1, x.max()+0.1)       # set limits in y for plot
    plt.ylim(y.min()-0.1, y.max()+0.1)       # set limits in y for plot
    plt.title(file)
    sep = '.'
    filename = file.split(sep, 1)[0]

    plt.savefig(filename+'.png')   # save figure as .png file
    plt.close(fig)

def three_d_plot(file,show):

    xloc=('two_d_initial_x_value' #0
            ,
    )
    yloc=('two_d_initial_y_value' #0
            ,
    )
    data=np.loadtxt(file)
    X=np.loadtxt(xloc[0]+".txt")
    Y=np.loadtxt(yloc[0]+".txt")
    u1,u2=data.shape
    X=range(u2)
    Y=range(u1)
    X,Y=np.meshgrid(X,Y)
    #    print X,Y
    fig = plt.figure()
    ax=fig.add_subplot(111, projection='3d')
#    ax.set_ylim(ax.get_ylim()[::-1])

    ax.plot_surface(X,Y,data, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False )
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.invert_yaxis()
    plt.title(file)
    sep = '.'
    filename = file.split(sep, 1)[0]
    plt.savefig(filename+'.png')   # save figure as .png file
    if show:
        # print('test')
        plt.show()
    plt.close(fig)

def three_d_contour(file,show):

    xloc = raw_input("xlocation(two_d_initial_x_value.txt):") or "two_d_initial_x_value.txt"
    yloc = raw_input("ylocation(two_d_initial_y_value.txt):") or "two_d_initial_y_value.txt"
    uloc = raw_input("ulocation(two_d_channel_u_results.txt):") or "two_d_channel_u_results.txt"
    vloc = raw_input("vlocation(two_d_channel_v_results.txt):") or "two_d_channel_v_results.txt"
    ploc = raw_input("plocation(two_d_channel_p_results.txt):") or "two_d_channel_p_results.txt"
    u=np.loadtxt(uloc)
    v=np.loadtxt(vloc)
    p=np.loadtxt(ploc)
    X=np.loadtxt(xloc)
    Y=np.loadtxt(yloc)
    X,Y=np.meshgrid(X,Y)
#    print X,Y
    fig = plt.figure()
    # plt.contourf(X,Y,p,alpha=0.5)    ###plnttong the pressure field as a contour
    # plt.colorbar()
    # plt.contour(X,Y,p)               ###plotting the pressure field outlines
    plt.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3]) ##plotting velocity
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(file)
    sep = '.'
    filename = file.split(sep, 1)[0]
    plt.savefig(filename+'.png')   # save figure as .png file
    if show:
        # print('test')
        plt.show()
    plt.close(fig)

def main():
    usage = "usage: %prog [options] txt_files"
    version="0.1"
    parser = optparse.OptionParser(version="%prog "+version, description="""\
    make a picture based on the txt files
                """)
    parser = optparse.OptionParser(usage, version="%prog 0.1.0", description="""plot with file""")
    parser.add_option("-d", "--dimension", type="int", default=2, help="plot 2d or 3d")
    parser.add_option("-t", "--type", default="", help="3d type contour=> c ")
    # parser.add_option("-t", "--drift", type="float", default=0, help="the drift number(time):for the true solution")
    parser.add_option("-v", action="store_true", dest="verbose",help="open show picture")
    parser.add_option("--out_dir", help="default OUT_DIR=.", default=".")
    parser.add_option("-f","--filename", default="", help="output filename")
    (options, args) = parser.parse_args()

    out_dir=options.out_dir
    dimension =options.dimension
    show=options.verbose
    print(show)
    type=options.type
    # drift=options.drift
    file=options.filename
    # sep = '.'
    # filename_no_ext = args[0].split(sep, 1)[0]
    if not file:
        file=args[0]
    if dimension==2:
        two_d_plot(file,show)
    elif type=="c":
        three_d_contour(file,show)
    else:
        three_d_plot(file,show)




if __name__ == "__main__":
    main()
