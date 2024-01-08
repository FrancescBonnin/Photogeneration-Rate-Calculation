# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 10:00:35 2019

@author: Francesc Bonnín Ripoll

Aquest és el programa definitiu pel càlcul de la photogeneracio
"""
#-----------------------------------------------------------------------------#
#                              LLIBRERIES
#-----------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import sys


print("Data i hora : ", time.ctime())

t = time.time
start = t()
#-----------------------------------------------------------------------------#
#                              FUNCIONS
#-----------------------------------------------------------------------------#

class Parametres:
    step_z = 1.25e-6 #mm
    step_lambda = 0.5
    error = 0.001
    S = 1.0
    q = 1.602176634e-19      # Electron charge (C)  

def v_longitud(vector):
    return np.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])

def tabulated_function(xvalues, yvalues):
    return lambda x: np.interp(x, xvalues, yvalues)

def limiting_wavelengths(dataframe, wl_limit):
    
    dataframe = dataframe.drop(dataframe[dataframe.lambda_PV > wl_limit].index)
    
    return dataframe

def wavelength_selection(dataframe, wl_region):
    
    dataframe = dataframe.drop(dataframe[dataframe.lambda_PV != wl_region].index)
    dataframe.reset_index()
    return dataframe

    
def llegir_source(file_in):

    source = pd.read_csv(file_in, sep='\t', header=None, comment='#')
    S_source = source.iloc[0,0]
    S_PV = source.iloc[2,0]
    wl_min = source.iloc[3,0]
    wl_max = source.iloc[4,0]
    step_source = source.iloc[5,0]    
    rays_per_wl = source.iloc[6,0]
    # array_wl_source = source.iloc[7:]
    array_wl_source = np.arange(wl_min,wl_max + step_source, step_source)
    array_wl_source = pd.DataFrame(array_wl_source)
    factor = S_source/S_PV #Factor geomètric entre les superfícies

    return wl_min, wl_max, factor, S_source, S_PV, rays_per_wl, step_source, array_wl_source

def spectrum_to_constant_step(file_in, step, wl_min, wl_max):
#   The wavelenght has to be in nm

    data_array = pd.read_csv(file_in, sep='\t', header=None)
    wl_spectrum = data_array.iloc[:,0]
    I_spectrum = data_array.iloc[:,1]
    
    array_spectrum_inter = pd.DataFrame([[x,np.interp(x,wl_spectrum,I_spectrum)] for x in np.arange(wl_min,wl_max+step,step)])
    array_spectrum_inter.columns = ['wl','I']

    return array_spectrum_inter, array_spectrum_inter['wl'], array_spectrum_inter['I']


def agrupar_espectre(array_spectrum_inter, wl_min, wl_max):
#   The wavelenght has to be in nm
#   Agrupam els valors de l'espectre solar per poder tenir menys dades.
    n = 6
    
    array_espectre_i = array_spectrum_inter[n:-n]    
    array_espectre_mes1 = array_spectrum_inter[(n-1):-(n+1)]
    array_espectre_menys1 = array_spectrum_inter[(n+1):-(n-1)]
    array_espectre_mes2 = array_spectrum_inter[(n-2):-(n+2)]
    array_espectre_menys2 = array_spectrum_inter[(n+2):-(n-2)]
    
    array_espectre_menys2 = array_espectre_menys2.reset_index(drop=True)
    array_espectre_mes2 = array_espectre_mes2.reset_index(drop=True)
    array_espectre_menys1 = array_espectre_menys1.reset_index(drop=True)
    array_espectre_mes1 = array_espectre_mes1.reset_index(drop=True)
    array_espectre_i = array_espectre_i.reset_index(drop=True)   
    
    array_espectre_2nm = pd.DataFrame()
    array_espectre_2nm['wl'] = array_spectrum_inter['wl'][n:-n].reset_index(drop=True)
    array_espectre_2nm['I'] = array_espectre_i['I'] + array_espectre_mes1['I'] + array_espectre_menys1['I'] + 0.5*(array_espectre_mes2['I'] + array_espectre_menys2['I']).reset_index(drop=True)
    array_espectre_2nm = array_espectre_2nm[(array_espectre_2nm.wl % 2 == 0.) & (array_espectre_2nm.wl >= wl_min-4.0) & (array_espectre_2nm.wl <= wl_max+4.0)].reset_index(drop=True)
    
    return array_espectre_2nm, array_espectre_2nm['wl'], array_espectre_2nm['I']

def numero_fotons(wl_min, array_wl_spectrum, array_I_spectrum, array_rays_source, rays_per_wl, step, factor):
#La font envia X W/m2, sigui quin sigui el número de rajos. Augmentar el nombre de rajos és augmentar 
# la informació que enviam. S en m^2. 

#En aquesta funció estam fent grups de rajos amb longituds d'ona de 2 en 2 nm i feim una mitjana de la irradiancia d'aquestes longituds d'ona. 
    h = 6.62607E-34
    c = 299792458.0
    hc = h*c
    I0 = array_I_spectrum[0]
    fotons_rang_promig = 0.0
    energia_rang_promig = 0.0
    df=[]
    
    rays_per_wl = rays_per_wl/factor
    print('rays per wl', rays_per_wl)
    
    for index, row in array_rays_source.itertuples():    
        lambda_f = row
        irradiance = tabulated_function(array_wl_spectrum,array_I_spectrum)(lambda_f) #irradiancia = W/m^2
        IM = (irradiance + I0)/2.0 * step
        en = hc/(lambda_f * 1e-9) # longitud d'ona en m
        num_fot = IM/en   #numero de fotons per m2     

        energia_rang_promig = (energia_rang_promig + num_fot * en)/rays_per_wl
        fotons_rang_promig = (fotons_rang_promig + num_fot)/rays_per_wl #fotons/ray        
        df.append((lambda_f, energia_rang_promig, fotons_rang_promig))

#Dividim pel nombre de rajos perque després li assignarem a cada raig un nombre de fotons.
        
        I0 = irradiance
        energia_rang_promig = 0.0
        fotons_rang_promig = 0.0 
        num_fot = 0.0
        
       
    array_num_fotons = pd.DataFrame(df)
    array_num_fotons.columns = ['lambda','energia_raig','fotons_raig']


    return array_num_fotons

def llegir_PV_values(file_in):
#   Adaptat pels fitxers source_wavelengths de 20/11/2019

    data_array = pd.read_csv(file_in, sep=' ', header=None, usecols = [0,1,2,3,4,5,6,7,8,9,10], skiprows=1, float_precision = 'round_trip')
    data_array.columns = ['x1','y1','z1','x2','y2','z2','en1','en2','lambda_PV','alfa','theta']
    
#    data_array['weigth'] = 1.0
    data_array['x1'] = np.round(data_array['x1'] , decimals = 10)
    data_array['y1'] = np.round(data_array['y1'] , decimals = 10)
    data_array['z1'] = np.round(data_array['z1'] , decimals = 10)
    data_array['x2'] = np.round(data_array['x2'] , decimals = 10)
    data_array['y2'] = np.round(data_array['y2'] , decimals = 10)
    data_array['z2'] = np.round(data_array['z2'] , decimals = 10)
#------------------------------------------------------------------------------  
    
    data_array['v_x'] = data_array['x2']-data_array['x1']
    data_array['v_y'] = data_array['y2']-data_array['y1']
    data_array['v_z'] = data_array['z2']-data_array['z1']
    
    del data_array['x1']
    del data_array['y1']
    del data_array['x2']
    del data_array['y2']
    
    vector = data_array[['v_x', 'v_y','v_z']].values
    data_array['vector'] = vector.tolist()
    data_array['longitud_vector'] = data_array['vector'].map(v_longitud)
    data_array['cos_theta_calculat'] = np.abs(data_array['v_z'])/data_array['longitud_vector'] 
    
    del data_array['v_x']
    del data_array['v_y']
    del data_array['vector']
    del data_array['longitud_vector']
    
    data_array['material'] = pd.read_csv(file_in, sep=' ', header=None, usecols = [11], skiprows=1)

    
    return data_array

def llegir_PV_values_1mat(file_in):
#   Adaptat pels fitxers source_wavelengths de 20/11/2019

    data_array = pd.read_csv(file_in, sep=' ', header=None, usecols = [0,1,2,3,4,5,6,7,8,9,10], skiprows=1, float_precision = 'round_trip')
    
    data_array.columns = ['x1','y1','z1','x2','y2','z2','en1','en2','lambda_PV','alfa','theta']
    

    data_array['x1'] = np.round(data_array['x1'] , decimals = 10)
    data_array['y1'] = np.round(data_array['y1'] , decimals = 10)
    data_array['z1'] = np.round(data_array['z1'] , decimals = 10)
    data_array['x2'] = np.round(data_array['x2'] , decimals = 10)
    data_array['y2'] = np.round(data_array['y2'] , decimals = 10)
    data_array['z2'] = np.round(data_array['z2'] , decimals = 10)
#------------------------------------------------------------------------------  
    
    data_array['v_x'] = data_array['x2']-data_array['x1']
    data_array['v_y'] = data_array['y2']-data_array['y1']
    data_array['v_z'] = data_array['z2']-data_array['z1']
    
    del data_array['x1']
    del data_array['y1']
    del data_array['x2']
    del data_array['y2']
    
    vector = data_array[['v_x', 'v_y','v_z']].values
    data_array['vector'] = vector.tolist()
    data_array['longitud_vector'] = data_array['vector'].map(v_longitud)
    data_array['cos_theta_calculat'] = np.abs(data_array['v_z'])/data_array['longitud_vector'] 
    
    del data_array['v_x']
    del data_array['v_y']
    del data_array['vector']
    del data_array['longitud_vector']

    
    return data_array
    
def filtrar_PV_values_3materials(data_array, HTM_width, ETM_width):  
    
    z_initial = data_array['z1'][0]
    PV_width = np.round(np.abs(data_array['z2'][1]-data_array['z1'][1]), decimals = 10)
    z_final = np.round(z_initial - PV_width - HTM_width - ETM_width, decimals = 10) #decimals de l'ordre de nanometres
    
    numrays_i = data_array['z1'].size
    
    data_array = data_array.drop(data_array[data_array.v_z == 0.0].index)
    data_array = data_array.drop(data_array[data_array.z1 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z2 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z1 < z_final].index) 
    data_array = data_array.drop(data_array[data_array.z2 < z_final].index)
    
    del data_array['v_z']
    
    data_array = data_array.reset_index(drop=True)
    
    numrays_f = data_array['z1'].size
    
    # print('rajos inicials = ', numrays_i)
    # print('rajos finals = ',numrays_f)
    # print('rajos borrats = ', numrays_i-numrays_f)
    # print('.......................................................')
    
    return data_array, z_final, z_initial, PV_width
    
def filtrar_PV_values_2materials(data_array, HTM_width):  
    
    z_initial = data_array['z1'][0]
    PV_width = np.round(np.abs(data_array['z2'][0]-data_array['z1'][0]), decimals = 10)
    z_final = np.round(z_initial - PV_width - HTM_width, decimals = 10) #decimals de l'ordre de nanometres
    
    numrays_i = data_array['z1'].size
    
    data_array = data_array.drop(data_array[data_array.v_z == 0.0].index)
    data_array = data_array.drop(data_array[data_array.z1 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z2 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z1 < z_final].index) 
    data_array = data_array.drop(data_array[data_array.z2 < z_final].index)
    
    del data_array['v_z']
    
    data_array = data_array.reset_index(drop=True)
    
    numrays_f = data_array['z1'].size
    
    # print('rajos inicials = ', numrays_i)
    # print('rajos finals = ',numrays_f)
    # print('rajos borrats = ', numrays_i-numrays_f)
    # print('.......................................................')
    
    return data_array, z_final, z_initial, PV_width
    
def filtrar_PV_values_HTM(data_array): 
    
    z_initial = np.around(data_array['z1'][0],decimals = 10)
    PV_width = np.around(np.abs(data_array['z2'][0]-data_array['z1'][0]),decimals = 10)
    z_final = np.around(z_initial-PV_width, decimals = 10) #decimals de l'ordre de nanometres
    
    print('z final =', z_final,'z inicial =', z_initial,'total width =', PV_width    )
    
    numrays_i = data_array['z1'].size
    
    data_array = data_array.drop(data_array[data_array.v_z == 0.0].index)
    print(data_array['z1'].size)
    data_array = data_array.drop(data_array[data_array.z1 > z_initial].index)
    print(data_array['z1'].size)
    data_array = data_array.drop(data_array[data_array.z2 > z_initial].index)
    print(data_array['z1'].size)
    data_array = data_array.drop(data_array[data_array.z1 < z_final].index) 
    print(data_array['z1'].size)
    data_array = data_array.drop(data_array[data_array.z2 < z_final].index)
    print(data_array['z1'].size)
    
    del data_array['v_z']
    
    data_array = data_array.reset_index(drop=True)
    
    numrays_f = data_array['z1'].size
    
    print('rajos inicials = ', numrays_i)
    print('rajos finals = ',numrays_f)
    print('rajos borrats = ', numrays_i-numrays_f)
    print('.......................................................')
    
    return data_array, z_final, z_initial, PV_width 

def filtrar_PV_values_1material(data_array): 
    
    z_initial = np.around(data_array['z1'][0],decimals = 10)
    PV_width = np.around(np.abs(data_array['z2'][0]-data_array['z1'][0]),decimals = 10)
    z_final = np.around(z_initial-PV_width, decimals = 10) #decimals de l'ordre de nanometres
    
    numrays_i = data_array['z1'].size
    
    data_array = data_array.drop(data_array[data_array.v_z == 0.0].index)
    data_array = data_array.drop(data_array[data_array.z1 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z2 > z_initial].index)
    data_array = data_array.drop(data_array[data_array.z1 < z_final].index) 
    data_array = data_array.drop(data_array[data_array.z2 < z_final].index)
    
    del data_array['v_z']
    
    data_array = data_array.reset_index(drop=True)
    
    numrays_f = data_array['z1'].size
    
    # print('rajos inicials = ', numrays_i)
    # print('rajos finals = ',numrays_f)
    # print('rajos borrats = ', numrays_i-numrays_f)
    # print('.......................................................')
    
    return data_array, z_final, z_initial, PV_width 
    

#-----Funcions per fer el tractament de les parassitic losses-----#
    
def llegir_reflectancia(file_in):
    
    data_array_R = pd.read_csv(file_in, sep=' ', header=None,skiprows=1)
    data_array_R.columns = ['lambda_source','theta','R_s','R_p','T_s','T_p']
    data_array_R['T_promig'] = 0.5*(data_array_R['T_s']+data_array_R['T_p'])
  
    return data_array_R
    
    
def separacio_PV_values(data_array, width):
    
    array_in = data_array[data_array['theta'] == 0.0]
    array_trap = data_array[data_array['theta'] != 0.0]
    
    return array_in, array_trap
#-----------------------------------------------------------------#
    
def Calcul_Photogeneration_sense_mirall(data_array, array_z, S):
    
    
    array_z_reflectit = array_z[::-1]
    z_initial = np.abs(data_array['z1'][0])

    i = 0
    j = 0
    k = 0
    l = 0
    
    T = np.zeros(np.size(array_z))
    T_inicial = 0.0
    T_final = 0.0
    
    for index, row in data_array.iterrows():
        z1_0 = row['z1']
        delta_z = np.abs(row['v_z'])
        alfa = row['alfa']
        cos = row['cos_theta_calculat']
        n0 = row['n0']
        

        if np.abs(z1_0) == z_initial:
            
            if cos == 1.0:
#Aqui entren els rajos amb angle d'incidencia 0.
                T_inicial = T_inicial + (alfa * 10) * n0  *  1e-4
                T = T + (alfa * 10) * (n0 * 1e-4) * np.exp(-alfa*array_z)  #z esta en mm, alfa està en mm^-1 , n0 en m^-2 i T en cm^-3
                i = i + 1
            else:
#Aqui entren els rajos amb angle d'incidència diferent a 0
                T_inicial = T_inicial + (alfa * 10) * n0  * 1e-4 
#                T = T + (alfa/cos * 10) * n0 * 1e-4 * np.exp(-alfa*array_z/cos)  #z esta en mm, alfa està en mm^-1, n0 en m^-2 i T en cm^-3
                T = T + (alfa/cos * 10) * n0 * 1e-4 * np.exp(-alfa*array_z/cos)  #z esta en mm, alfa està en mm^-1, n0 en m^-2 i T en cm^-3
                k = k + 1
        else:
#Aqui entren les reflexions
            if cos == 0.0: 
                l = l + 1  
                T_final = T_final + (alfa * 10) * n0 * 1e-4
            else:    
                T_final = T_final + (alfa * 10) * n0 * 1e-4
                T = T + (alfa/cos * 10) * n0 *1e-4 * np.exp(-alfa*array_z_reflectit/cos)  #alfa està en mm^-1 i n0 en m^-2            
                j = j + 1
#    inclourem un filtre pel càlcul de la població de fotons. Quan la quantitat d'aquests sigui lo suficientment petita, deixarem de realitzar el càlcul, que mai serà 0.
#
#    print('rajos entrants amb incidencia 0 = ',i,', rajos entrants amb incidencia no 0 = ',k)
#    print('rajos rebotats = ',j,', rajos amb cos nul =',l)
#    print('numero total rajos = ', i+j+k+l)
#    print('nombre de línies en el fitxer = ', data_array.shape[0])


    return data_array, T, T_inicial, T_final
    
def Calcul_Photogeneration_BC(data_array, array_z, step_z, z_initial, z_final):
#    Calculam la photogeneracio G(z) tenint en compte miralls ideals en les condicions de contorn.
    
    G = np.zeros(np.size(array_z))

    index_total_z = array_z.size
    ordre_step = np.int(-np.log10(step_z))+1 
    
    i = 0
    j = 0
    k = 0
    l = 0
    num = 0
    
    for index, row in data_array.iterrows():
        z1_0 = row['z1']
        z2_0 = row['z2']
        en1 = row['en1']
        alfa = row['alfa']
        cos = row['cos_theta_calculat']
        n0 = row['n0']   
        
        if np.abs(z1_0) == z_initial:
                                                                                                            #incidence from the top
            i = i + 1            
            delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)   
            BC_index = delta_z.size
            
            G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos) 

        else:
            if np.abs(z1_0) == np.abs(z_final):                                                                     #reflection from the bottom
                j = j + 1
                delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)
                delta_z = delta_z[::-1]
                BC_index = delta_z.size
                
                G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

            else:                                                                                            
                if z2_0 == z_initial:                                                                       #lateral reflection to the top
                    k = k + 1             
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    delta_z = delta_z[::-1]
                    BC_index = delta_z.size
                    
                    G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

                else:                                                                                       #lateral reflection to the bottom
                    l = l + 1
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    BC_index = delta_z.size
                    
                    G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  

    # print('.......................................................')
    # print('rajos desde z=0 = ', i )
    # print('rajos laterals up= ', k, ', rajos laterals down = ', l)
    # print('rajos desde z=width = ', j )
    # print('suma rajos = ', i+j+k+l)
    # print('nombre de línies en el fitxer = ', data_array.shape[0])
    # print('rajos que no compleixen la relacio = ', num )


    return data_array, G

def Calcul_Photogeneration_BC_2mat(data_array, array_z, step_z, z_initial, z_final):
#    Calculam la photogeneracio G(z) tenint en compte miralls ideals en les condicions de contorn.
    

    G = np.zeros(np.size(array_z))

    index_total_z = array_z.size
    ordre_step = np.int(-np.log10(step_z))+1 
    
    i = 0
    j = 0
    k = 0
    l = 0
    num = 0
    
    for index, row in data_array.iterrows():
        z1_0 = row['z1']
        z2_0 = row['z2']
        en1 = row['en1']
        alfa = row['alfa']
        cos = row['cos_theta_calculat']
        n0 = row['n0']   
        
        if np.abs(z1_0) == z_initial:
                                                                                                            #incidence from the top
            i = i + 1            
            delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)   
            BC_index = delta_z.size
            
            G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos) 

        else:
            if np.abs(z1_0) == np.abs(z_final):                                                                     #reflection from the bottom
                j = j + 1
                delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)
                delta_z = delta_z[::-1]
                BC_index = delta_z.size
                
                G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

            else:                                                                                            
                if z2_0 == z_initial:                                                                       #lateral reflection to the top
                    k = k + 1             
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    delta_z = delta_z[::-1]
                    BC_index = delta_z.size
                    
                    G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

                else:                                                                                       #lateral reflection to the bottom
                    l = l + 1
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    BC_index = delta_z.size
                    
                    G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  

    # print('.......................................................')
    # print('rajos desde z=0 = ', i )
    # print('rajos laterals up= ', k, ', rajos laterals down = ', l)
    # print('rajos desde z=width = ', j )
    # print('suma rajos = ', i+j+k+l)
    # print('nombre de línies en el fitxer = ', data_array.shape[0])
    # print('rajos que no compleixen la relacio = ', num )


    return data_array, G
    
def Calcul_Photogeneration_BC_3mat(data_array, array_z, step_z, z_initial, z_final, array_z_ETM):
#    Calculam la photogeneracio G(z) tenint en compte miralls ideals en les condicions de contorn.
    
    G = np.zeros(np.size(array_z))

    index_total_z = array_z.size
    index_ETM = array_z_ETM.size
    ordre_step = np.int(-np.log10(step_z))+1 
    
    i = 0
    j = 0
    k = 0
    l = 0
    num = 0
    
    for index, row in data_array.iterrows():
        z1_0 = row['z1']
        z2_0 = row['z2']
        en1 = row['en1']
        alfa = row['alfa']
        cos = row['cos_theta_calculat']
        n0 = row['n0']   

        if np.abs(z1_0) == np.abs(z_initial):
                                                                                                            #incidence from the top
            i = i + 1            
            delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)   
            BC_index = delta_z.size
            
            G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos) 

        else:
            if np.abs(z1_0) == np.abs(z_final):                                                                     #reflection from the bottom
                j = j + 1
                delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)
                delta_z = delta_z[::-1]
                BC_index = delta_z.size
                
                G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

            else:                                                                                            
                if z2_0 == z_initial:                                                                       #lateral reflection to the top
                    k = k + 1             
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    delta_z = delta_z[::-1]
                    BC_index = delta_z.size
                    
                    G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)

                else:                                                                                       #lateral reflection to the bottom
                    if cos == 1.0:
                        delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                        BC_index = delta_z.size
                        # print (z1_0,z2_0, BC_index, G[(index_total_z-BC_index):].size)
                        G[index_ETM:BC_index+index_ETM] = G[(index_ETM):BC_index+index_ETM] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  
                        # print(G.size,G)
                        # sys.exit()
                    else:
                        l = l + 1
                        delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                        BC_index = delta_z.size
                        
                        G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  
                    
    print('.......................................................')
    print('rajos desde z=0 = ', i )
    print('rajos laterals up= ', k, ', rajos laterals down = ', l)
    print('rajos desde z=width = ', j )
    print('suma rajos = ', i+j+k+l)
    print('nombre de línies en el fitxer = ', data_array.shape[0])
    print('rajos que no compleixen la relacio = ', num )
    

    return data_array, G
    
    
def Calcul_Photogeneration_BC_4G(data_array, array_z, step_z):
    
    z_initial = np.abs(data_array['z1'][0])
    z_final = np.abs(data_array['z2'][0])
    
    G = np.zeros(np.size(array_z))  
    G1 = np.zeros(np.size(array_z))
    G2 = np.zeros(np.size(array_z))
    G3 = np.zeros(np.size(array_z))
    G4 = np.zeros(np.size(array_z))
    
    sum_pvpath_1 = 0.0
    sum_pvpath_2 = 0.0
    sum_pvpath_3 = 0.0
    sum_pvpath_4 = 0.0
    
    index_total_z = array_z.size
    ordre_step = np.int(-np.log10(step_z))+1 
    i = 0
    j = 0
    k = 0
    l = 0
    num = 0
    
    for index, row in data_array.iterrows():
        z1_0 = row['z1']
        z2_0 = row['z2']       
        en1 = row['en1'] 
        en2 = row['en2'] 
        alfa = row['alfa']
        cos = row['cos_theta_calculat']
        n0 = row['n0']
        
        if np.abs(z1_0) == z_initial:
                                                                                                            #incidence from the top
            i = i + 1            
            delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)   
            BC_index = delta_z.size
            
            G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos) 
            G1[:BC_index] = G1[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos) 
            sum_pvpath_1 = sum_pvpath_1 + (n0*(en1-en2))

        else:
            if np.abs(z1_0) == z_final:                                                                     #reflection from the bottom
                j = j + 1
                delta_z = np.arange(0.0, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z)
                delta_z = delta_z[::-1]
                BC_index = delta_z.size
                
                G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)
                G2[(index_total_z-BC_index):] = G2[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)
                sum_pvpath_2 = sum_pvpath_2 + (n0*(en1-en2))
                
            else:                                                                                            
                if z2_0 == z_initial:                                                                       #lateral reflection to the top
                    k = k + 1             
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    delta_z = delta_z[::-1]
                    BC_index = delta_z.size
                    
                    G[:BC_index] = G[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)
                    G3[:BC_index] = G3[:BC_index] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)
                    sum_pvpath_3 = sum_pvpath_3 + (n0*(en1-en2))

                else:                                                                                       #lateral reflection to the bottom
                    l = l + 1
                    delta_z = np.arange(step_z, np.around(np.abs(z1_0-z2_0), decimals = ordre_step) + step_z*0.1, step_z) 
                    BC_index = delta_z.size
                    
                    G[(index_total_z-BC_index):] = G[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  
                    G4[(index_total_z-BC_index):] = G4[(index_total_z-BC_index):] + 1.0/cos * (alfa * 10.0) * (en1*n0 * 1e-4) * np.exp(-alfa*delta_z/cos)  
                    sum_pvpath_4 = sum_pvpath_4 + (n0*(en1-en2))

    # print('...............................')
    integral_1 = (1e4)*np.trapz(G1,array_z*0.1)
    integral_2 = (1e4)*np.trapz(G2,array_z*0.1)
    integral_3 = (1e4)*np.trapz(G3,array_z*0.1)
    integral_4 = (1e4)*np.trapz(G4,array_z*0.1)    
        
    # print('rajos desde z=0 = ', i )
    # print('rajos laterals up= ', k, ', rajos laterals down = ', l)
    # print('rajos desde z=width = ', j )
    # print('suma rajos = ', i+j+k+l)
    # print('nombre de línies en el fitxer = ', data_array.shape[0])
    # print('rajos que no compleixen la relacio = ', num )
    # print('-------------------------------------------')
    # print('TOP',integral_1/sum_pvpath_1*100.0)
    # print('BOTTOM',integral_2/sum_pvpath_2*100.0)
    # print('UP',integral_3/sum_pvpath_3*100.0)
    # print('DOWN',integral_4/sum_pvpath_4*100.0)


    return data_array, G, G1, G2, G3, G4
    
def return_step(file_in, error):
#   alfa està en mm^-1
#   les distàncies també estan en mm

    alfa_file = pd.read_csv(file_in, sep=' ', header=None, skiprows=1)
    alfa_file.columns = ['x1','y1','z1','x2','y2','z2','en1','en2','lambda_PV','alfa','theta'] 
    alfa_file['width'] = np.abs(alfa_file['z2'] - alfa_file['z1'])
    width = alfa_file['width'].iloc[0] #en mm
    df = []
    df2 = []

    for index, row in alfa_file.iterrows():
        alfa = row['alfa']
        lambda_PV = row['lambda_PV'] 
        if alfa != 0.:
            delta_x = -np.log(1.-error)/alfa
        else:
            delta_x = 0.0
        df.append((alfa,delta_x))
        df2.append((alfa,lambda_PV))
    df_alfa = pd.DataFrame(df2)
    df_alfa.columns = ['alfa','lambda_PV']
    mask = df_alfa.duplicated()
    df_alfa = df_alfa[~mask]
    df_alfa = df_alfa.reset_index(drop=True)
    df_delta_x = pd.DataFrame(df)
    df_delta_x.columns = ['alfa','step']

    return df_delta_x, df_alfa, width

def discretize_z_1mat(step_z, width):
#Aquesta funció fa una discretització de la dimensió z segons cada step definit per cada
#coeficient d'absorció. 

    zmin = 0.0
    zmax = width  #mm
    array_z = np.arange(zmin,zmax + step_z*0.1, step_z)
    
    
    print('Material Width =', width*1e6,'nm')

    return array_z
    
def discretize_z_2mat(step_z, width, htm_width):
#Aquesta funció fa una discretització de la dimensió z segons cada step definit per cada
#coeficient d'absorció. 

    zmin = 0.0
    zmax = width + htm_width #mm
    array_z = np.arange(zmin,zmax + step_z*0.1, step_z)
    
    n_split = np.int(np.round(width/step_z)+1)
    
    array_z_PV, array_z_HTM = np.split(array_z, [n_split])
    
    
    print('Perovskite Width =', width*1e6,'nm')
    print('HTM Width =', htm_width*1e6,'nm')
    
    return array_z, array_z_PV, array_z_HTM
    
def discretize_z_3mat(step_z, width, htm_width, etm_width):
#Aquesta funció fa una discretització de la dimensió z segons cada step definit per cada
#coeficient d'absorció. 

    zmin = 0.0
    zmax = width + htm_width + etm_width #mm
    array_z = np.arange(zmin,zmax + step_z*0.2, step_z)
    array_z_ETM = np.arange(zmin,etm_width + step_z*0.2, step_z)
    
    print('ETM Width =', etm_width*1e6,'nm')
    print('Perovskite Width =', width*1e6,'nm')
    print('HTM Width =', htm_width*1e6,'nm')
    
    return array_z, array_z_ETM
    
def separar_PV_path(data_array):
    
    data_array_PV = data_array[data_array.material == 'PV']
    data_array_HTM = data_array[data_array.material == 'Spiro']
    data_array_PV = data_array_PV.reset_index()
    data_array_HTM = data_array_HTM.reset_index()
    
    print('Separació per material')
    print('suma = ', data_array_HTM.shape[0]+data_array_PV.shape[0])
    print('initial data array = ',data_array.shape[0])
    print("diferencia = ", -data_array_HTM.shape[0]-data_array_PV.shape[0]+data_array.shape[0])
    
    return data_array_PV, data_array_HTM

def separar_PV_path_z_2mat(data_array, HTM_width):
    
    z_initial = np.around(data_array['z1'][0],decimals = 10)
    PV_width = np.around(np.abs(data_array['z2'][0]-data_array['z1'][0]),decimals = 10)
    z_final_PV = np.around(z_initial-PV_width, decimals = 10)
    z_final_HTM = np.around(z_initial-PV_width-HTM_width, decimals = 10)
    
    
    data_array_PV = data_array[(data_array['z1'] <= z_initial) & (data_array['z1'] >=z_final_PV) & (data_array['z2'] >= z_final_PV) & (data_array['z2'] <= z_initial)]
    data_array_HTM = data_array[(data_array['z1'] <= z_final_PV) & (data_array['z1'] >= z_final_HTM) & (data_array['z2'] >= z_final_HTM) & (data_array['z2'] <= z_final_PV)]
    data_array_PV = data_array_PV.reset_index()
    data_array_HTM = data_array_HTM.reset_index()
        
    return data_array_PV, data_array_HTM

def separar_PV_path_z_3mat(data_array, HTM_width, ETM_width):
    
    z_initial = np.around(data_array['z1'][0],decimals = 10)
    z_initial_PV = np.around(data_array['z1'][1],decimals = 10)
    PV_width = np.around(np.abs(data_array['z2'][1]-data_array['z1'][1]),decimals = 10)
    z_final_PV = np.around(z_initial_PV-PV_width, decimals = 10)
    z_final_HTM = np.around(z_initial-PV_width-HTM_width-ETM_width, decimals = 10)
    
    data_array_ETM = data_array[(data_array['z1'] <= z_initial) & (data_array['z1'] >=z_initial_PV) & (data_array['z2'] >= z_initial_PV) & (data_array['z2'] <= z_initial)]
    data_array_PV = data_array[(data_array['z1'] <= z_initial_PV) & (data_array['z1'] >=z_final_PV) & (data_array['z2'] >= z_final_PV) & (data_array['z2'] <= z_initial_PV)]
    data_array_HTM = data_array[(data_array['z1'] <= z_final_PV) & (data_array['z1'] >= z_final_HTM) & (data_array['z2'] >= z_final_HTM) & (data_array['z2'] <= z_final_PV)]
    
    data_array_ETM = data_array_ETM.reset_index()
    data_array_PV = data_array_PV.reset_index()
    data_array_HTM = data_array_HTM.reset_index()
    
    print('Separació per z')
    print('suma = ', data_array_ETM.shape[0]+data_array_HTM.shape[0]+data_array_PV.shape[0])
    print('initial data array = ',data_array.shape[0])   
    print("diferencia = ", -data_array_ETM.shape[0]-data_array_HTM.shape[0]-data_array_PV.shape[0]+data_array.shape[0])
    
    return data_array_PV, data_array_HTM, data_array_ETM
    
    
def resultats(data_array, array_G, array_z, S, q):
    fotons_PV_path = (data_array['n0']*(-data_array['en2']+data_array['en1'])).sum() #n0 esta en m-2
    
    #Càlcul de la integral de G per a la comprovació del nombre de fotons.
    integral = (S*1e4)*np.trapz(array_G,array_z*0.1) #tot ha d'estar en cm. Això dona el fotos absorbits per m2
    
    J_ph = 1e3*q*np.trapz(array_G,array_z*0.1) # mA/cm2

    
    return fotons_PV_path, integral, J_ph

def resultats_separats(array_G_PV, array_G_HTM, array_z_HTM, array_z_PV, S, q):
    
    #Càlcul de la integral de G per a la comprovació del nombre de fotons.
    integral_PV = (S*1e4)*np.trapz(array_G_PV,array_z_PV*0.1) #tot ha d'estar en cm. Això dona el fotos absorbits per m2
    integral_HTM = (S*1e4)*np.trapz(array_G_HTM,array_z_HTM*0.1)
    
    J_ph_PV = 1e3*q*np.trapz(array_G_PV,array_z_PV*0.1) # mA/cm2
    J_ph_HTM = 1e3*q*np.trapz(array_G_HTM,array_z_HTM*0.1) # mA/cm2
    
    return integral_PV, integral_HTM, J_ph_PV, J_ph_HTM

def resultats_separats_3mat(array_G_ETM, array_G_PV, array_G_HTM, array_z_HTM, array_z_PV, array_z_ETM, S, q):
    
    #Càlcul de la integral de G per a la comprovació del nombre de fotons.
    integral_ETM = (S*1e4)*np.trapz(array_G_ETM,array_z_ETM*0.1)
    integral_PV = (S*1e4)*np.trapz(array_G_PV,array_z_PV*0.1) #tot ha d'estar en cm. Això dona el fotos absorbits per m2
    integral_HTM = (S*1e4)*np.trapz(array_G_HTM,array_z_HTM*0.1)
    
    J_ph_ETM = 1e3*q*np.trapz(array_G_ETM,array_z_ETM*0.1) # mA/cm2
    J_ph_PV = 1e3*q*np.trapz(array_G_PV,array_z_PV*0.1) # mA/cm2
    J_ph_HTM = 1e3*q*np.trapz(array_G_HTM,array_z_HTM*0.1) # mA/cm2
    
    return integral_ETM, integral_PV, integral_HTM, J_ph_ETM, J_ph_PV, J_ph_HTM

def fotons_data(data_array, lambda_source, fotons_rang_source, z_initial):
    
    data_array['n0'] =  tabulated_function(lambda_source, fotons_rang_source)(data_array['lambda_PV']) #la usam si 'weight' és 1.0 sempre.
    nz0 = data_array['n0'][data_array.z1 == z_initial][data_array['theta'] == 0.0].sum()
    
    return nz0, data_array

def fotons_data_ETM(data_array, lambda_source, fotons_rang_source, z_initial):
    
    #la usam si 'weight' és 1.0 sempre.
    data_array['n0'] =  tabulated_function(lambda_source, fotons_rang_source)(data_array['lambda_PV'])
    nz0 = data_array['n0'][data_array.z1 == z_initial][data_array['theta'] == 180.0].sum() #   n0 són els fotons per àrea (m^-2). theta = 180 per casuistica del pvpath
    fotons_PV_path = (data_array['n0']*(-data_array['en2']+data_array['en1'])).sum()
    
    return nz0, fotons_PV_path, data_array

def fotons_data_2mat(data_array, lambda_source, fotons_rang_source, z_initial):
    
    #la usam si 'weight' és 1.0 sempre.
    data_array['n0'] =  tabulated_function(lambda_source, fotons_rang_source)(data_array['lambda_PV'])
    
    nz0 = data_array['n0'][data_array.z1 == z_initial][data_array['theta'] == 0.0].sum() #   n0 són els fotons per àrea (m^-2)
    
    fotons_PV_path = (data_array['n0']*(-data_array['en2']+data_array['en1'])).sum()
    
    return nz0, fotons_PV_path, data_array

    
def llegir_th_spectral(file_in):
#   Adaptat pels fitxers source_wavelengths de 20/11/2019

    th_spectral = pd.read_csv(file_in, sep='\s+', header=None, skiprows=1, float_precision = 'round_trip')
    th_spectral.columns = ['lambda_PV','efficiency_Th_absorbed']

    return th_spectral

def llegir_th_integral(file_in):
#   Adaptat pels fitxers source_wavelengths de 20/11/2019

    th_integral = pd.read_csv(file_in, sep='\s+', header=None, skiprows=1, float_precision = 'round_trip')
    th_integral.columns = ['power_absorbed','irradiance_emitted','efficiency']
    return th_integral

