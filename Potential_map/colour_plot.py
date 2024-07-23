import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append('/Users/xxxx/.local/pipx/venvs/ltspice/lib/python3.12/site-packages/')
import ltspice

path_here = os.path.dirname(os.path.realpath(__file__))+'/'

params = {"Number of rows (n)": 15,
          "Number of columns (m)": 15,
          "Applied voltage": 1,
          ".raw file": path_here+"Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3.raw",
          "Sim Dir": "Report things 2 Random length and doping",
          "Sim file": "Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3",
          "Voltage names": ["NV1", "NV2", "NV3", "NV4", "NV5"],
          "Multiple sources locs": None,
        #   ([range(3), range(3,6), range(6,9), range(9,12), range(12,15)],["V1","V2","V3","V4","V5"]), #taken from array_model, not necessary if only one input voltage
          #([range(5), range(10,15)],["V1","V2"])
          }

def nearest_voltage(raw_file=params['.raw file'], n=params['Number of rows (n)'], m=params["Number of columns (m)"], VA=params['Applied voltage']):
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

def Vin(V_in=params["Multiple sources locs"], n=params['Number of rows (n)']):
    if V_in == None:
        V_in_map = 0
    else:
        num_sources = len(V_in[0])
        in_locs = V_in[0]
        in_names = V_in[1]
        V_in_map = np.zeros(n+1,dtype='object')
        for i in range(num_sources):
            if type(in_locs[i]) == range:
                for j in in_locs[i]:
                    V_in_map[j] = str(f'N{in_names[i]}')
            else:
                V_in_map[i] = str(f'N{in_names[i]}')
    return V_in_map

def SimParser(indx, Vinmap_defined, v_node=params['Voltage names'], n=params["Number of rows (n)"], m=params["Number of columns (m)"], raw_file=params[".raw file"], SimDir=params["Sim Dir"], SimFile=params["Sim file"]):
    row_map = np.zeros((n, 3*m-3))
    V_row_map = np.zeros_like(row_map)
    V_in_map = np.zeros(n+1, dtype='object')
    V_ground_map = np.zeros_like(V_in_map)
    Vnodemap= np.zeros_like(row_map, dtype='object')
    Imap_exp = []
    Vmap_exp = []

    filepath = raw_file
    LT = ltspice.Ltspice(filepath)
    LT.parse()
    # LT.get_data(str(f'I(Rx{j}y{i})'))

    #read all the rows (ignore vertical grain boundaries between rows)
    for i in range(n):
        for j in range(3*m-3):
            # print(str(f'I(Rx{j}y{i})'))
            if j == 3*m-4:
                # row_map[i, j] = np.log10(abs(LT.get_data(str(f'I(V1)'))[indx]))
                # Vnodemap[i,j] = str(f'V(x{j:02.0f}y{i:02.0f})')

                row_map[i, j] = LT.get_data(str(f'I(Rx{j}y{i})'))[indx]
                V_row_map[i, j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
                try:
                    LT.get_data(str(f'V(x{j+1:02.0f}y{i:02.0f})'))[indx]
                except: 
                    V_ground_map[i] = 1
                    print(f'NOTE: Ground connected to x{j:02.0f}y{i:02.0f}')
                
                #this is for final node but it isn't needed for the current
            elif j == 0:
                row_map[i, j] = LT.get_data(str(f'I(Rx{j}y{i})'))[indx]
                try:
                    V_row_map[i, j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
                except:
                    print(f"NOTE: Voltage source at node x{j:02.0f}y{i:02.0f}")
                    if params["Multiple sources locs"]==None:
                        for v in v_node:
                            try:
                                V_row_map[i, j] = LT.get_data(str(f'V({v})'))[indx]
                                V_in_map[i] = v
                                break
                            except:
                                V_row_map[i, j] = 0
                                print(f"WARNING: No voltage V({v}) nodes found between nodes x{j:02.0f}y{i:02.0f} and x{j+1:02.0f}y{i:02.0f}.") 
                    else:
                        V_row_map[i, j] = LT.get_data(str(f'V({Vinmap_defined[i]})'))[indx]
                        print(f"Voltage source {Vinmap_defined[i]} between nodes x{j:02.0f}y{i:02.0f} and x{j+1:02.0f}y{i:02.0f}.")
            else:
                # Vnodemap[i, j] = str(f'V(x{j:02.0f}y{i:02.0f})-V(x{j+1:02.0f}y{i:02.0f})')
                row_map[i, j] = LT.get_data(str(f'I(Rx{j}y{i})'))[indx]
                V_row_map[i,j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
            Imap_exp.append(row_map[i, j])
            Vmap_exp.append(V_row_map[i, j])
            # print(f'{v_node[0]}')
            # print(f"V={LT.get_data(str(f'V({v_node[0]})'))[indx]}")
    print(V_in_map)
    if params["Multiple sources locs"] != None:
        V_in_map=Vinmap_defined
    print(V_in_map)
    print(V_ground_map)
    col_map = np.zeros((n-1, 3*m-3))
    V_col_map = np.zeros_like(col_map)
    V_col_nodemap = np.zeros_like(col_map, dtype='object')
    for i in range(n-1): 
        for j in range(0,3*m-2,3):
            print(str(f'I(Rx{j}y{i}V), x{j:02.0f}y{i:02.0f}'))
            if j == 3*m-3:
                if V_ground_map[i] == 1 and V_ground_map[i+1] == 1:
                    V_col_map[i, j-1] = 0
                    print(f"WARNING: No resistor found between x{j:02.0f}y{i:02.0f} and x{j:02.0f}y{i+1:02.0f}")
                elif V_ground_map[i] == 1:
                    V_col_map[i, j-1] = 0
                    print(f"NOTE: Ground connected above x{j:02.0f}y{i:02.0f}")
                elif V_ground_map[i+1] == 1:
                    V_col_map[i, j-1] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
                    print(f"NOTE: Ground connected below x{j:02.0f}y{i:02.0f}")
                else:
                    V_col_map[i, j-1] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]

                # try:
                #    V_col_map[i, j-1] = np.log10(abs(LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx] - LT.get_data(str(f'V(x{j:02.0f}y{i+1:02.0f})'))[indx]))
                # except:
                #     try:
                #         V_col_map[i, j-1] = np.log10(abs(LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]))
                #         print(f"NOTE: Ground connected below x{j:02.0f}y{i:02.0f}")
                #     except:
                #         try:
                #            V_col_map[i, j-1] = np.log10(abs(-LT.get_data(str(f'V(x{j:02.0f}y{i+1:02.0f})'))[indx]))
                #            print(f"NOTE: Ground connected above x{j:02.0f}y{i:02.0f}") 
                #         except:
                #             V_col_map[i, j-1] = 0
                #             print(f"WARNING: No resistor found between x{j:02.0f}y{i:02.0f} and x{j:02.0f}y{i+1:02.0f}")

                # if i == n-2:
                #     V_col_map[i, j-1] = np.log10(abs(LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]))
                # else:
                #     V_col_map[i, j-1] = np.log10(abs(LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx] - LT.get_data(str(f'V(x{j:02.0f}y{i+1:02.0f})'))[indx]))
                #     # V_col_nodemap[i, j-1] = str(f'V(x{j:02.0f}y{i:02.0f})-V(x{j:02.0f}y{i+1:02.0f})')
                # #move the array index to the left
                try:
                    col_map[i, j-1] = LT.get_data(str(f'I(Rx{j}y{i}V)'))[indx]
                    Imap_exp.append(col_map[i, j-1])
                except:
                    print(str(f'No resistor found at I(Rx{j}y{i}V)'))
                    col_map[i, j-1] = 0
                    Imap_exp.append(col_map[i, j-1])
                # col_map[i, j-1] = str(f'I(Rx{j}y{i}V)')
                Vmap_exp.append(V_col_map[i, j-1])
                
            elif j == 0:
                try:
                    col_map[i, j] = LT.get_data(str(f'I(Rx{j}y{i}V)'))[indx]
                    Imap_exp.append(col_map[i, j])
                except:
                    print(str(f'No resistor found at I(Rx{j}y{i}V)'))
                    col_map[i, j] = 0
                    Imap_exp.append(col_map[i, j])
                # col_map[i, j] = str(f'I(Rx{j}y{i}V)')
                    
                if V_in_map[i] != 0 and V_in_map[i+1]!=0:
                    V_col_map[i, j]=0
                    print(f"WARNING: No resistor found between x{j:02.0f}y{i:02.0f} and x{j:02.0f}y{i+1:02.0f}")
                elif V_in_map[i] != 0:
                    V_col_map[i, j] = LT.get_data(str(f'V({V_in_map[i]})'))[indx]
                    print(f"NOTE: Voltage {V_in_map[i]} above node x{j:02.0f}y{i+1:02.0f}.")
                elif V_in_map[i+1] != 0:
                    V_col_map[i, j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
                    print(f"NOTE: Voltage {V_in_map[i+1]} below node x{j:02.0f}y{i:02.0f}.")
                else:
                    V_col_map[i, j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]

                # if i == 0:
                #     V_col_map[i, j] = np.log10(abs(LT.get_data(str(f'V(NV1)'))[indx] - LT.get_data(str(f'V(x{j:02.0f}y{i+1:02.0f})'))[indx]))
                # else:
                #     V_col_map[i, j] = np.log10(abs(LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx] - LT.get_data(str(f'V(x{j:02.0f}y{i+1:02.0f})'))[indx]))
                # # V_col_nodemap[i, j] = str(f'V(x{j:02.0f}y{i:02.0f})-V(x{j:02.0f}y{i+1:02.0f})')
                Vmap_exp.append(V_col_map[i, j])
            else:
                col_map[i, j] = LT.get_data(str(f'I(Rx{j}y{i}V)'))[indx]
                col_map[i, j-1] = col_map[i, j]
                Imap_exp.append(col_map[i, j])
                # col_map[i, j] = str(f'I(Rx{j}y{i}V)')
                # col_map[i, j-1] = str(f'I(Rx{j}y{i}V)')
                V_col_map[i, j] = LT.get_data(str(f'V(x{j:02.0f}y{i:02.0f})'))[indx]
                V_col_map[i, j-1] = V_col_map[i, j]
                Vmap_exp.append(V_col_map[i, j])
                # V_col_nodemap[i, j] = str(f'V(x{j:02.0f}y{i:02.0f})-V(x{j:02.0f}y{i+1:02.0f})')
                # V_col_nodemap[i, j-1] = V_col_nodemap[i, j]
    newImap = np.zeros((2*n-1, 3*m-3))
    Vmap = np.zeros_like(newImap)
    
    np.savetxt(path_here+'\\'+SimDir+'\\current_map_'+SimFile+'.txt', Imap_exp)
    np.savetxt(path_here+'\\'+SimDir+'\\voltage_map_'+SimFile+'.txt', Vmap_exp)

    for i in range(2*n-1):
        if i % 2 == 0:
            newImap[i, :] = row_map[int(i/2)]
            Vmap[i, :] = V_row_map[int(i/2)]
        else: 
            newImap[i, :] = col_map[int((i-1)/2)]
            Vmap[i, :] = V_col_map[int((i-1)/2)]

    print(str(f'Imap = {newImap}'))
    print(str(f'Vmap = {Vmap}'))
    print(str(f"Vmap exported = {Vmap_exp}"))
    masked_Imap = np.ma.masked_equal(newImap, 0)
    masked_Vmap = np.ma.masked_equal(Vmap, 0)
    plt.figure(1)
    I = plt.imshow(masked_Imap, origin="lower", extent=(0, 1e-3, 1e-3, 0), cmap='Blues_r')
    plt.colorbar(I)
    plt.figure(2)
    V = plt.imshow(masked_Vmap, origin="lower", extent=(0, 1e-3, 1e-3, 0), cmap='Blues_r')
    plt.colorbar(V)
    plt.show()
    return masked_Imap, masked_Vmap

indx, V_actual = nearest_voltage(VA=params["Applied voltage"])
V_in = Vin()
masked_Imap, masked_Vmap = SimParser(indx,Vinmap_defined=V_in)

# fig, ax = plt.subplots()
# I = plt.imshow(masked_Imap, extent=(0, 1e-3, 1e-3, 0), cmap='Blues_r')
# plt.colorbar(I)

# def update(frame):
#         indx, V_actual = nearest_voltage(frame)
#         masked_Imap= SimParser(indx)
#         return masked_Imap


# from matplotlib.animation import FuncAnimation

# ani = FuncAnimation(fig, update, frames=params["Applied voltage"], blit=True)
# plt.show()


#NEXT just test imshow on a length matrix just to see if it works

