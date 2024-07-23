import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append('/Users/xxxx/.local/pipx/venvs/ltspice/lib/python3.12/site-packages/')
import ltspice
from scipy.interpolate import griddata
from pylab import figure, cm
from matplotlib.colors import LogNorm

path_here = os.path.dirname(os.path.realpath(__file__))+'/'

params = {"Number of rows (n)": 15,
          "Number of columns (m)": 15,
          "Length distribution type": "constant", #can be constant or exponential
          "Length": 1.5e-3, #for constant doping
          "Length distribution address": path_here+"exp_length_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2.txt",
          "Doping distribution type": "poisson", #poisson or constant
          "ND": 1e16, #if constant
          "Doping distribution address": path_here+"poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2.txt", 
          "Applied voltage": 1,
          ".raw file": path_here+"Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3.raw",
          "Sim Dir": path_here,
          "Sim file": "Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2",
          "Voltage divider data?": True,
          }

def nearest_voltage(raw_file=params['.raw file'], SimDir=params["Sim Dir"], SimFile=params["Sim file"], n=params['Number of rows (n)'], m=params["Number of columns (m)"], VA=params['Applied voltage']):
    if raw_file == None:
        filepath = path_here+"\\"+SimDir+"\\"+SimFile+".raw"
    else:
        filepath = raw_file
    filesize = os.stat(filepath).st_size
    print(str(f'filesize = {filesize}'))
    LT = ltspice.Ltspice(filepath)
    LT.parse()
    V_source = LT.get_data('V1')
    
    indx = np.searchsorted(V_source, VA) #the voltage may not exactly match the desired applied voltage so have to find the closest one
    
    if abs(V_source[indx]-VA) < abs(V_source[indx-1]-VA):
        indx = indx
    else:
        indx = indx - 1
    VA_actual = V_source[indx] #this is the closest value from the simulation
    print(str(f'index = {indx}'))
    return indx, VA_actual

def load_files(n=params["Number of rows (n)"], m=params["Number of columns (m)"], SimDir=params["Sim Dir"], SimFile=params["Sim file"], L_file=params["Length distribution address"], L_dist=params["Length distribution type"], length=params["Length"], ND_const=params['ND'], doping_dist=params["Doping distribution type"], ND_file=params["Doping distribution address"]):
    Imap = np.loadtxt(path_here+'current_map_'+SimFile+'.txt')

    try:
        Vmap = np.loadtxt(path_here+'voltage_map_'+SimFile+'.txt')
    except:
        print("No voltage map found")
        Vmap = 0
    
    if L_file == None and L_dist=="exponential":
        length_matrix = np.loadtxt(path_here+"\\"+SimDir+"\\exp_length_dist_"+SimFile+".txt")
        l_display = length_matrix
    elif L_dist=="exponential":
        length_matrix = np.loadtxt(L_file)
        l_display = length_matrix
    elif L_dist=='constant' and L_file!=None: #constant
        length_matrix = np.full((n, 2*m-2), length)
        l_display = np.loadtxt(L_file)
    else:
        length_matrix = np.full((n, 2*m-2), length)
        l_display = length_matrix

    if ND_file == None and doping_dist=="poisson":
        doping_matrix = np.loadtxt(path_here+"\\"+SimDir+"\\poisson_doping_dist_"+SimFile+".txt")
    elif doping_dist == "poisson":
        doping_matrix = np.loadtxt(ND_file)
    else: #constant
        doping_matrix = np.full((n, 2*m-2), ND_const)
    
    print(f'length matrix = {length_matrix}')
    print(f'ND = {doping_matrix}')
    
    return Imap, Vmap, length_matrix, doping_matrix, l_display
    



def SimParser(length_matrix, doping_matrix, length_display, n=params["Number of rows (n)"], m=params["Number of columns (m)"]):
    x = np.zeros(4*m*n-3*n-m)
    y = np.zeros_like(x)
    ND = np.zeros_like(x)
    L = np.zeros_like(x)

    # if raw_file == None:
    #     filepath = path_here+"\\"+SimDir+"\\"+SimFile+".raw"
    # else:
    #     filepath = raw_file
    # LT = ltspice.Ltspice(filepath)
    # LT.parse()

    counter = 0 #to count through columns 
    grain_counter = 0
    grain_boundary_seq = np.arange(1,3*m-4, 3) 
    # prev_y = np.zeros(3*m-3)
    prev_y = np.zeros(3*m-3)
    print(f'Starting y value = {prev_y}')
    input_counter = []
    output_counter = []
    boundary_counter =[]
    # starting_grain = max(length_matrix[:,0]/2)
    for j in range(n):
        for i in range(3*m-3):
            print(f'horizontal i = {i}, j = {j}')
            if i == 0: #first grain
                x[counter]=length_matrix[j,grain_counter]/2
                y[counter]=length_matrix[j,grain_counter]/2 + prev_y[i]
                prev_length = x[counter] 
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2 #to account for the halving
                if j == 0:
                    input_counter.append(counter)
                ND[counter]=doping_matrix[j,grain_counter]
                L[counter]=length_display[j,grain_counter]
                grain_counter += 1
            elif i == 1:
                x[counter]=prev_length + length_matrix[j,grain_counter-1]/2 #-1 because the grain counter already increased to next boundary
                y[counter]=(length_matrix[j,grain_counter-1]/2+length_matrix[j,grain_counter]/2)/2 + prev_y[i] #take average height of two surrounding grains
                prev_length = x[counter]
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2
                boundary_counter.append(counter)
                print('Grain boundary:')
            elif i in grain_boundary_seq: #grain boundary
                x[counter]=prev_length + length_matrix[j,grain_counter-1]/4 #-1 because the grain counter already increased to next boundary
                y[counter]=(length_matrix[j,grain_counter-1]/2+length_matrix[j,grain_counter]/2)/2 + prev_y[i]
                prev_length = x[counter]
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2
                boundary_counter.append(counter)
                print('Grain boundary:')
            elif i == 3*m-4:# final grain in row
                x[counter]=prev_length + length_matrix[j,grain_counter]/2
                y[counter]=length_matrix[j,grain_counter]/2 + prev_y[i]
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2
                if j == n-1:
                    output_counter.append(counter)
                grain_counter = 0 #reset for next row
                ND[counter]=doping_matrix[j,grain_counter]
                L[counter]=length_display[j,grain_counter]
            elif grain_counter%2==0: #for grains next to the boundary
                x[counter]=prev_length + length_matrix[j,grain_counter]/2
                y[counter]=length_matrix[j,grain_counter]/2 + prev_y[i]
                prev_length = x[counter]
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2
                ND[counter]=doping_matrix[j,grain_counter]
                L[counter]=length_display[j,grain_counter]
                grain_counter += 1
            else: #for grains just after a boundary
                x[counter]=prev_length + length_matrix[j,grain_counter]/4
                y[counter]=length_matrix[j,grain_counter]/2 + prev_y[i]
                prev_length = x[counter]
                prev_y[i] = y[counter] + length_matrix[j,grain_counter]/2
                ND[counter]=doping_matrix[j,grain_counter]
                L[counter]=length_display[j,grain_counter]
                grain_counter += 1
            # print(f'i={i}, j={j}, x={x[counter]}, y={y[counter]}, L from matrix={length_matrix[j,grain_counter]}, previous_length={prev_length}')
            counter += 1

    # prev_yv = np.zeros(2*m-2)
    prev_yv = np.zeros(2*m-2)
    # x[counter] = 0
    # y[counter] = prev_yv[0] + length_matrix[0, 0]
    # prev_yv[0] = y[counter]
    # counter += 1
    
    input_counter.append(counter)

    for j in range(n-1):
        prev_length = 0
        x[counter] = length_matrix[j,0]/2
        y[counter] = prev_yv[0] + length_matrix[j, 0]
        prev_yv[0] = y[counter]
        boundary_counter.append(counter)
        counter += 1
        for i in range(0, 2*m-2):
            print(f'vertical i = {i}, j = {j}')
            if i in range(1, 2*m-2, 2):
                # if i == 2*m-3:
                #     x[counter] = prev_length + length_matrix[j,i]
                #     y[counter] = prev_yv[i] + length_matrix[j,i]
                #     prev_yv[i] = y[counter]
                #     counter += 1
                # else:
                x[counter] = prev_length + (length_matrix[j,i]/2+length_matrix[j+1,i]/2)/2
                y[counter] = prev_yv[i] + length_matrix[j,i]
                prev_yv[i] = y[counter]
                boundary_counter.append(counter)
                counter += 1
            else:
                prev_length += (length_matrix[j, i] + length_matrix[j+1, i])/2
    output_counter.append(counter-1)
    print(f'x = {x}\ny = {y}\nlength matrix = {length_matrix}')
    plt.figure()
    plt.plot(x, y, 'k.')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    # plt.show()
    # print(f'x shape = {x.shape}, y shape = {y.shape}, Imap shape = {Imap.shape}')
    return x, y, input_counter, output_counter, boundary_counter, ND, L

def griddata_plot(x, y, z, Vmap, input_counter, output_counter, boundary_counter, ND, Ldist, L=params['Length']):
    # target grid to interpolate to
    # xi = yi = np.arange(0,1.01,0.01)
    x = x*10**4
    y = y*10**4
    maxx = max(x)
    maxy = max(y)
    xi = np.linspace(0,maxx+L/2*10**4,1000)
    yi = np.linspace(0,maxy+L/2*10**4,1000)
    xi,yi = np.meshgrid(xi,yi)
    
    #current
    zi = griddata((x,y),z,(xi,yi),method='nearest')

    #voltage
    Vi = griddata((x,y),Vmap,(xi,yi),method='nearest')

    #doping
    zeros = np.where(ND==0)[0] 
    ND_reduced = np.delete(ND, zeros)
    L_reduced = np.delete(Ldist, zeros)
    x_reduced = np.delete(x, zeros)
    y_reduced = np.delete(y, zeros)


    # print(f'ND_reduced = {ND_reduced}')
    NDi_masked = griddata((x_reduced,y_reduced),ND_reduced,(xi,yi),method='nearest')
    Li = griddata((x_reduced,y_reduced),L_reduced*10**4,(xi,yi),method='nearest')
    # NDi = griddata((x,y),ND,(xi,yi),method='nearest') #this was NDi
    # NDi_masked = np.ma.masked_where(NDi == 0, NDi)

    # Exceeded Vlimit (ONLY works for divider plots)
    if params['Voltage divider data?']==True:
        VlimMap = Vmap > 1.235
        Xlim = x*VlimMap #includes grains
        Ylim = y*VlimMap
        Xlim = Xlim[zeros] #only grain boundaries
        Ylim = Ylim[zeros]
        toremove = np.where(Xlim==0)[0]
        Xlim=np.delete(Xlim,toremove)
        Ylim=np.delete(Ylim,toremove)
        np.savetxt(path_here+'x_XS_voltage'+params['Sim file']+'.txt', Xlim)
        np.savetxt(path_here+'y_XS_voltage'+params['Sim file']+'.txt', Ylim)
        print(f'VlimMap={VlimMap}')
    else:
        try:
            Xlim = np.loadtxt(path_here+'x_XS_voltage'+params['Sim file']+'.txt')
            Ylim = np.loadtxt(path_here+'y_XS_voltage'+params['Sim file']+'.txt')
        except:
            print("No voltage limit data available.")
            Xlim=np.nan
            Ylim=np.nan

    
    fig = plt.figure()
    plt.subplot(121)
    # ax = fig.add_subplot(111)
    # f = figure()
    # ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
    # axcolor = f.add_axes([0.90, 0.02, 0.03, 0.79])
    # plt.imshow(zi, origin='lower', cmap='viridis', alpha=0.5, label='Interpolated Grid')
    # plt.contourf(xi,yi,zi, cmap="Blues_r", norm=LogNorm())
    im1 = plt.imshow(zi, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar(im1).set_label(label="Current (A)", size=16)
    plt.plot(x,y,'k.', zorder=1)
    plt.scatter(x[input_counter], y[input_counter], s=175, marker='o', edgecolors='red', c='None', zorder=3)
    plt.scatter(x[output_counter], y[output_counter], s=175, marker='o', edgecolors='green', c='None', zorder=2)
    plt.scatter(x[boundary_counter], y[boundary_counter], marker='.', c='red', zorder=4)
    plt.scatter(Xlim, Ylim, marker='*', c='r', zorder=5, s=100)
    plt.xlabel('x ($\mu$m)',fontsize=16)
    plt.ylabel('y ($\mu$m)',fontsize=16)
    # plt.title("Current plot",fontsize=16)
    # fig = plt.figure()
    plt.subplot(122)
    # ax = fig.add_subplot(111)
    # plt.imshow(zi, origin='lower', cmap='viridis', alpha=0.5, label='Interpolated Grid')
    # plt.contourf(xi,yi,Vi, cmap="Blues_r")
    im2=plt.imshow(Vi, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar(im2).set_label(label="Voltage (V)", size=16)
    plt.scatter(x[boundary_counter], y[boundary_counter], marker='.', c='r', zorder=4, label='Grain boundaries')
    plt.plot(x,y,'k.', zorder=1, label='Grains')
    plt.scatter(x[input_counter], y[input_counter], s=175, marker='o', edgecolors='red', c='None', zorder=2, label='Input')
    plt.scatter(x[output_counter], y[output_counter], s=175, marker='o', edgecolors='green', c='None', zorder=3, label='Ground')
    plt.scatter(Xlim, Ylim, marker='*', c='r', zorder=5, label="Exceeded voltage limit", s=100)
    plt.xlabel('x ($\mu$m)',fontsize=16)
    # plt.ylabel('y ($\mu$m)',fontsize=16)
    # plt.title("Voltage plot",fontsize=16)
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, ncol=5)
    # plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    # plt.legend(bbox_to_anchor=(1, 0), loc="lower right", bbox_transform=fig.transFigure, ncol=3)
    plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 0.9), ncol=3, fontsize=12)
    plt.tight_layout()
    plt.figure()
    # plt.contourf(xi,yi,NDi_masked, cmap="Blues_r")
    plt.imshow(NDi_masked, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar().set_label(label="Doping concentration (cm$^{-3}$)", size=16)
    # plt.plot(x_reduced,y_reduced,'k.')
    # plt.scatter(x[input_counter], y[input_counter], s=100, marker='o', edgecolors='red', c='None')
    # plt.scatter(x[output_counter], y[output_counter], s=100, marker='o', edgecolors='green', c='None')
    # plt.scatter(x[boundary_counter], y[boundary_counter], s=100, marker='o', edgecolors='blue', c='None')
    plt.xlabel('x ($\mu$m)',fontsize=16)
    plt.ylabel('y ($\mu$m)',fontsize=16)
    # plt.title("Doping plot",fontsize=16)

    plt.figure()
    plt.imshow(Li, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar().set_label(label="Grain size ($\mu$m)", size=16)
    # plt.plot(x_reduced,y_reduced,'k.')
    # plt.scatter(x[input_counter], y[input_counter], s=100, marker='o', edgecolors='red', c='None')
    # plt.scatter(x[output_counter], y[output_counter], s=100, marker='o', edgecolors='green', c='None')
    # plt.scatter(x[boundary_counter], y[boundary_counter], s=100, marker='o', edgecolors='blue', c='None')
    plt.xlabel('x ($\mu$m)',fontsize=16)
    plt.ylabel('y ($\mu$m)',fontsize=16)
    # plt.title("Doping plot",fontsize=16)

    plt.figure()
    im1 = plt.imshow(zi, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar(im1).set_label(label="Current (A)", size=16)
    # plt.plot(x,y,'k.', zorder=1)
    # plt.scatter(x[input_counter], y[input_counter], s=175, marker='o', edgecolors='red', c='None', zorder=3)
    # plt.scatter(x[output_counter], y[output_counter], s=175, marker='o', edgecolors='green', c='None', zorder=2)
    # plt.scatter(x[boundary_counter], y[boundary_counter], marker='.', c='red', zorder=4)
    # plt.scatter(Xlim, Ylim, marker='*', c='r', zorder=5, s=100)
    plt.xlabel('x ($\mu$m)',fontsize=16)
    plt.ylabel('y ($\mu$m)',fontsize=16)

    plt.figure()
    im2=plt.imshow(Vi, extent=(xi.min(), xi.max(), yi.min(), yi.max()), origin='lower', aspect='auto', cmap='Blues_r')
    plt.colorbar(im2).set_label(label="Voltage (V)", size=16)
    # plt.scatter(x[boundary_counter], y[boundary_counter], marker='.', c='r', zorder=4, label='Grain boundaries')
    # plt.plot(x,y,'k.', zorder=1, label='Grains')
    # plt.scatter(x[input_counter], y[input_counter], s=175, marker='o', edgecolors='red', c='None', zorder=2, label='Input')
    # plt.scatter(x[output_counter], y[output_counter], s=175, marker='o', edgecolors='green', c='None', zorder=3, label='Ground')
    # plt.scatter(Xlim, Ylim, marker='*', c='r', zorder=5, label="Exceeded voltage limit", s=100)
    plt.xlabel('x ($\mu$m)',fontsize=16)
    plt.ylabel('y ($\mu$m)',fontsize=16)

    plt.show()

    # plt.savefig('interpolated.png',dpi=100)
    # plt.close(fig)
    return Vi, xi, yi


    
            


z, Vmap, length_matrix, doping_matrix, l_display = load_files()
x,y,in_counter,out_counter,boundary_counter,ND, Ldist=SimParser(length_matrix=length_matrix, doping_matrix=doping_matrix, length_display=l_display)
print(f'ND={ND}')
griddata_plot(x=x, y=y, z=z, Vmap=Vmap, input_counter=in_counter, output_counter=out_counter, boundary_counter=boundary_counter, ND=ND, Ldist=Ldist)