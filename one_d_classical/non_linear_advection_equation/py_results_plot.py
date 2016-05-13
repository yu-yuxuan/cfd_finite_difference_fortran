#!/usr/bin/python2.7

"""
This is codes to plot the results from cfd
"""
import sys
import optparse
import os
import numpy as np
import matplotlib.pyplot as plt

def two_d_plot(file):
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
def main():
    usage = "usage: %prog [options] txt_files"
    version="0.1"
    parser = optparse.OptionParser(version="%prog "+version, description="""\
    make a picture based on the txt files
                """)
    parser = optparse.OptionParser(usage, version="%prog 0.1.0", description="""plot with file""")
    parser.add_option("-d", "--dimension", type="int", default=2, help="plot 2d or 3d")
    # parser.add_option("-t", "--drift", type="float", default=0, help="the drift number(time):for the true solution")
    parser.add_option("--out_dir", help="default OUT_DIR=.", default=".")
    parser.add_option("--filename", default="", help="output filename")
    (options, args) = parser.parse_args()

    out_dir=options.out_dir
    dimension =options.dimension
    # drift=options.drift
    file=options.filename
    # sep = '.'
    # filename_no_ext = args[0].split(sep, 1)[0]
    if not file:
        file=args[0]
    if dimension==2:
        two_d_plot(file)



if __name__ == "__main__":
    main()
