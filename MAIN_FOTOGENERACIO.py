# -*- coding: utf-8 -*-
"""
Created on Thu May 07 10:17:57 2020

@author: Francesc Bonnín Ripoll
"""
#-----------------------------------------------------------------------------#
#                              LLIBRERIES
#-----------------------------------------------------------------------------#
from Funcions_Photogeneration_BC_mirall_many_material_DEFINITIU import *
from Photogeneration_plots import G_Plot
import pathlib as Path
from os import listdir
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import isfile, join


#-----------------------------------------------------------------------------#
#                              PARÀMETRES
#-----------------------------------------------------------------------------#
    
step_z = 1.25e-6 #mm
S = 1.0
q = 1.602176634e-19 
step_lambda = 0.5
wl_min_spectrum = 295.0
wl_max_spectrum = 805.0

array_wl_astm = np.arange(300.0, 800.0 + step_lambda , step_lambda*2.0)
array_wl_astm = pd.DataFrame(array_wl_astm)

one_material = True
two_material = False
two_material_separats = False
three_materials = False
three_materials_separats = False

prints = False

limiting_wavelength = False
wl_limit = 564.0 #nm

HTM = 'organic'
# HTM = 'inorganic'

HTM_mat = 'Spiro'
ETM_mat = 'TiO2'

if HTM == 'organic': 
    HTM_width = 225e-6 #nm
if HTM == 'inorganic':
    HTM_width = 85e-6 #nm
if one_material:
    HTM_width = 0.0 #nm
    
# ETM_width = 56.5e-6 #nm
    
# ETM_width = 41e-6 #nm



#-----------------------------------------------------------------------------#
#                                Fitxers
#-----------------------------------------------------------------------------#  
 
cwd = Path.Path().resolve()
path = str(cwd.parents[0])

carpeta_output = 'OUTPUTS' 

carpeta_input = 'INPUTS'

carpeta_spectrum = 'Spectrum_Files'

file_ASTM = path + '\Photogeneration_Rate_Calculation\\'+ carpeta_spectrum + '\ASTMG173-total.txt'

array_spectrum_inter, wl_array_inter, I_array_inter = spectrum_to_constant_step(file_ASTM, step_lambda*2., wl_min_spectrum, wl_max_spectrum)

dir_PVpath = path + '\Photogeneration_Rate_Calculation\\' + carpeta_input

dir_outputs = path + '\Photogeneration_Rate_Calculation\\' + carpeta_output

fitxers_carpeta = [dir_PVpath + '\\' + f for f in listdir(dir_PVpath) if isfile(join(dir_PVpath, f))]

num_fitxers = len(fitxers_carpeta)

print("Casos a fer: ", num_fitxers/2)

fitxers_PV_path = fitxers_carpeta[:int(num_fitxers/2)]
fitxers_source = fitxers_carpeta[int(num_fitxers/2):]

for pv_path, source in zip(fitxers_PV_path, fitxers_source):
    parts = pv_path.split('_values')
    
    file_source = source
    file_pvpath = pv_path  
    
    if not two_material_separats:
        output_G = open(dir_outputs + '\\' + "G" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data" + parts[1] , "w+") 
    
    print('-------------------------------------------------------------------')
    print('fitxer entrada: ', parts[1])
    print('-------------------------------------------------------------------')

    
    wl_min, wl_max, factor, S_source, S_PV, rays_per_wl, step_source, array_wl_source = llegir_source(file_source)
    
    array_espectre_2nm, wl_array, I_array = agrupar_espectre(array_spectrum_inter, wl_min, wl_max)
    
    data_array = llegir_PV_values(file_pvpath)
    
    step_source = step_source*0.25 #Estam sumant 4 vegades el mateix nombre

    array_num_fotons = numero_fotons(wl_min, wl_array, I_array, array_wl_source, rays_per_wl, step_source, factor)
    
    #Per obtenir l'espectre cada 1 nm.
    array_num_fotons_astm = numero_fotons(wl_min, wl_array_inter, I_array_inter, array_wl_astm, rays_per_wl, step_source*2., factor)
    array_num_fotons_astm = array_num_fotons_astm [::-1]
    
    #data_array_reduced = reduccio_data(data_array,rays_per_wl)

    if three_materials_separats:
        # data_array_PV, data_array_HTM = separar_PV_path(data_array)
        data_array_PV, data_array_HTM, data_array_ETM = separar_PV_path_z_3mat(data_array, HTM_width, ETM_width)
        
        data_array_ETM, z_final_ETM, z_initial_ETM, ETM_width  = filtrar_PV_values_1material(data_array_ETM)
        data_array_PV, z_final_PV, z_initial_PV, PV_width  = filtrar_PV_values_1material(data_array_PV)
        data_array_HTM, z_final_HTM, z_initial_HTM, HTM_width  = filtrar_PV_values_1material(data_array_HTM)
        
        array_z_ETM = discretize_z_1mat(step_z, ETM_width)
        array_z_PV = discretize_z_1mat(step_z, PV_width)
        array_z_HTM = discretize_z_1mat(step_z, HTM_width)
        
        nz0_ETM, fotons_PV_path_ETM, data_array_ETM = fotons_data_ETM(data_array_ETM, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial_ETM)
        nz0_PV, fotons_PV_path_PV, data_array_PV = fotons_data_2mat(data_array_PV, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial_PV)
        nz0_HTM, fotons_PV_path_HTM, data_array_HTM = fotons_data_2mat(data_array_HTM, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial_HTM) 
        
        
        data_array_ETM, array_G_ETM = Calcul_Photogeneration_BC(data_array_ETM, array_z_ETM, step_z, z_initial_ETM, z_final_ETM)
        data_array_PV, array_G_PV = Calcul_Photogeneration_BC(data_array_PV, array_z_PV, step_z, z_initial_PV, z_final_PV)
        data_array_HTM, array_G_HTM = Calcul_Photogeneration_BC(data_array_HTM, array_z_HTM, step_z, z_initial_HTM, z_final_HTM)
        
        integral_ETM, integral_PV, integral_HTM, J_ph_ETM, J_ph_PV, J_ph_HTM = resultats_separats_3mat(array_G_ETM, array_G_PV, array_G_HTM, array_z_HTM, array_z_PV, array_z_ETM, S, q)
        
        print('-------------------------Results per Separat ---------------------------------')
        print('Photocurrent PV (mA/cm2) = ', J_ph_ETM)
        print('Photocurrent PV (mA/cm2) = ', J_ph_PV)
        print('Photocurrent HTM (mA/cm2) = ', J_ph_HTM)  
        print('Photocurrent TOTAL (mA/cm2) = ', J_ph_ETM + J_ph_HTM + J_ph_PV)  
        print('Number of absorbed photons (ETM) (1/m2) = ',integral_ETM)
        print('Number of absorbed photons (PV) (1/m2) = ',integral_PV)
        print('Number of absorbed photons (HTM) (1/m2) = ',integral_HTM)
        print('Number of absorbed photons (TOTAL) (1/m2) = ',integral_ETM + integral_PV + integral_HTM)
        print('------------------------- RATIOS ------------------------------------')
        print('RATIO integral and pvpath (ETM)= ', integral_ETM/fotons_PV_path_ETM*100.0)
        print('RATIO integral and pvpath (PV)= ', integral_PV/fotons_PV_path_PV*100.0)
        print('RATIO integral and pvpath (HTM)= ', integral_HTM/fotons_PV_path_HTM*100.0)
        print('RATIO integral and pvpath (TOTAL)= ', (integral_ETM + integral_HTM + integral_PV)/(fotons_PV_path_ETM + fotons_PV_path_HTM + fotons_PV_path_PV)*100.0)
        print('RATIO absorbed/sent photons ETM= ', integral_ETM/nz0_ETM*100.0)
        print('RATIO absorbed/sent photons PV= ', integral_PV/nz0_PV*100.0)
        print('RATIO absorbed/sent photons HTM= ', integral_HTM/nz0_HTM*100.0)
        print('RATIO absorbed/sent photons TOTAL= ', (integral_ETM + integral_HTM + integral_PV)/(nz0_PV)*100.0)
        
        df_plot_ETM = pd.DataFrame({'z($\mu m$)':array_z_ETM,'Generation Rate G(n/cm3)':array_G_ETM})
        df_plot_PV = pd.DataFrame({'z($\mu m$)':array_z_PV,'Generation Rate G(n/cm3)':array_G_PV})
        df_plot_HTM = pd.DataFrame({'z($\mu m$)':array_z_HTM,'Generation Rate G(n/cm3)':array_G_HTM})
        df_data_ETM = pd.DataFrame({'Width(nm)':ETM_width,'Absorbed Photons':integral_ETM, 'Ratio absorbed':integral_ETM/nz0_ETM*100.0, 'Check':integral_ETM/fotons_PV_path_ETM*100.0,'Photocurrent_PV(mA/cm2)': J_ph_ETM}, index=[0])
        df_data_PV = pd.DataFrame({'Width(nm)':PV_width,'Absorbed Photons':integral_PV, 'Ratio absorbed':integral_PV/nz0_PV*100.0, 'Check':integral_PV/fotons_PV_path_PV*100.0,'Photocurrent_PV(mA/cm2)': J_ph_PV}, index=[0])
        df_data_HTM = pd.DataFrame({'Width(nm)':HTM_width,'Absorbed Photons':integral_HTM, 'Ratio absorbed':integral_HTM/nz0_HTM*100.0, 'Check':integral_HTM/fotons_PV_path_HTM*100.0,'Photocurrent_HTM(mA/cm2)': J_ph_HTM}, index=[0])
     
        output_G = open(dir_outputs + '\\' + "G_PV" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data_PV" + parts[1] , "w+") 
        
        df_plot_PV.to_csv(output_G, header=True, index=False, sep='\t') 
        df_data_PV.to_csv(output_data, header=True, index=False, sep='\t')
        
        output_G = open(dir_outputs + '\\' + "G_HTM" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data_HTM" + parts[1] , "w+") 
        
        df_plot_HTM.to_csv(output_G, header=True, index=False, sep='\t') 
        df_data_HTM.to_csv(output_data, header=True, index=False, sep='\t')

        output_G = open(dir_outputs + '\\' + "G_ETM" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data_ETM" + parts[1] , "w+") 
        
        df_plot_ETM.to_csv(output_G, header=True, index=False, sep='\t') 
        df_data_ETM.to_csv(output_data, header=True, index=False, sep='\t')
        
        output_G.close()
        output_data.close()
        
        array_z_HTM = array_z_HTM*1e6
        array_z_PV = array_z_PV*1e6
        array_z_ETM = array_z_ETM*1e6
        G_Plot(array_z_ETM, array_G_ETM, ETM_mat, parts[1].split('_PV.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        G_Plot(array_z_PV, array_G_PV, HTM_mat, parts[1].split('_PV.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        G_Plot(array_z_HTM, array_G_HTM, HTM_mat, parts[1].split('_HTM.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        
        array_G_tot = np.append(np.append(array_G_ETM, array_G_PV), array_G_HTM)
        array_z_tot = np.append(array_z_ETM, np.append(array_z_PV, array_z_HTM + PV_width*1e6) + ETM_width*1e6)
        G_Plot(array_z_tot, array_G_tot, HTM_mat+'_'+ETM_mat, parts[1].split('_total')[0], (ETM_width + PV_width + HTM_width)*1e6, dir_outputs)
        # sys.exit()
    if two_material_separats:
            
        data_array_PV, data_array_HTM = separar_PV_path_z_2mat(data_array, HTM_width)
        
        data_array_PV, z_final_PV, z_initial_PV, PV_width  = filtrar_PV_values_1material(data_array_PV)
        data_array_HTM, z_final_HTM, z_initial_HTM, HTM_width  = filtrar_PV_values_1material(data_array_HTM)
        # sys.exit()
        array_z_PV = discretize_z_1mat(step_z, PV_width)
        array_z_HTM = discretize_z_1mat(step_z, HTM_width)
        
        nz0_PV, fotons_PV_path_PV, data_array_PV = fotons_data_2mat(data_array_PV, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial_PV)
        nz0_HTM, fotons_PV_path_HTM, data_array_HTM = fotons_data_2mat(data_array_HTM, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial_HTM) 
        
        if limiting_wavelength: 
            data_array_HTM = limiting_wavelengths(data_array_HTM, wl_limit)
        
        data_array_PV, array_G_PV = Calcul_Photogeneration_BC(data_array_PV, array_z_PV, step_z, z_initial_PV, z_final_PV)
        data_array_HTM, array_G_HTM = Calcul_Photogeneration_BC(data_array_HTM, array_z_HTM, step_z, z_initial_HTM, z_final_HTM)
        
        integral_PV, integral_HTM, J_ph_PV, J_ph_HTM = resultats_separats(array_G_PV, array_G_HTM, array_z_HTM, array_z_PV, S, q)
        if prints:
            print('-------------------------Results per Separat ---------------------------------')
            print('Photocurrent PV (mA/cm2) = ', J_ph_PV)
            print('Photocurrent HTM (mA/cm2) = ', J_ph_HTM)  
            print('Photocurrent TOTAL (mA/cm2) = ', J_ph_HTM + J_ph_PV)  
            print('Number of absorbed photons (PV) (1/m2) = ',integral_PV)
            print('Number of absorbed photons (HTM) (1/m2) = ',integral_HTM)
            print('Number of absorbed photons (TOTAL) (1/m2) = ',integral_PV + integral_HTM)
            print('------------------------- RATIOS ------------------------------------')
            print('RATIO integral i pvpath (PV)= ', integral_PV/fotons_PV_path_PV*100.0)
            print('RATIO integral i pvpath (HTM)= ', integral_HTM/fotons_PV_path_HTM*100.0)
            print('RATIO integral i pvpath (TOTAL)= ', (integral_HTM + integral_PV)/(fotons_PV_path_HTM + fotons_PV_path_PV)*100.0)
            print('RATIO absorbed/sent photons PV= ', integral_PV/nz0_PV*100.0)
            print('RATIO absorbed/sent photons HTM= ', integral_HTM/nz0_HTM*100.0)
            print('RATIO absorbed/sent photons TOTAL= ', (integral_HTM + integral_PV)/(nz0_PV)*100.0)
        
        df_plot_PV = pd.DataFrame({'z($\mu m$)':array_z_PV,'Generation Rate G(n/cm3)':array_G_PV})
        df_plot_HTM = pd.DataFrame({'z($\mu m$)':array_z_HTM,'Generation Rate G(n/cm3)':array_G_HTM})
        df_data_PV = pd.DataFrame({'Width(nm)':PV_width,'Absorbed Photons':integral_PV, 'Ratio absorbed':integral_PV/nz0_PV*100.0, 'Check':integral_PV/fotons_PV_path_PV*100.0,'Photocurrent_PV(mA/cm2)': J_ph_PV}, index=[0])
        df_data_HTM = pd.DataFrame({'Width(nm)':HTM_width,'Absorbed Photons':integral_HTM, 'Ratio absorbed':integral_HTM/nz0_HTM*100.0, 'Check':integral_HTM/fotons_PV_path_HTM*100.0,'Photocurrent_HTM(mA/cm2)': J_ph_HTM}, index=[0])
     
        output_G = open(dir_outputs + '\\' + "G_PV" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data_PV" + parts[1] , "w+") 
        output_G = dir_outputs + '\\' + "G_PV" + parts[1]
        output_data = dir_outputs + '\\' + "Data_PV" + parts[1]
        
        df_plot_PV.to_csv(output_G, header=True, index=False, sep='\t') 
        df_data_PV.to_csv(output_data, header=True, index=False, sep='\t')
        
        output_G = open(dir_outputs + '\\' + "G_HTM" + parts[1] , "w+")
        output_data = open(dir_outputs + '\\' + "Data_HTM" + parts[1] , "w+")
        output_G = dir_outputs + '\\' + "G_HTM" + parts[1]
        output_data = dir_outputs + '\\' + "Data_HTM" + parts[1]
        
        df_plot_HTM.to_csv(output_G, header=True, index=False, sep='\t') 
        df_data_HTM.to_csv(output_data, header=True, index=False, sep='\t')
        
        # output_G.close()
        # output_data.close()
        
        print('-------RESULTS NUMBER PHOTONS--------', integral_HTM)
        
        array_z_HTM = array_z_HTM*1e6
        array_z_PV = array_z_PV*1e6
        G_Plot(array_z_PV, array_G_PV, HTM_mat, parts[1].split('_PV.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        G_Plot(array_z_HTM, array_G_HTM, HTM_mat, parts[1].split('_HTM.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        
        array_G_tot = np.append(array_G_PV, array_G_HTM)
        array_z_tot = np.append(array_z_PV, array_z_HTM + PV_width*1e6)
        G_Plot(array_z_tot, array_G_tot, HTM_mat, parts[1].split('_total.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        
    else:
        if one_material:
            print('One Material')
            data_array, z_final, z_initial, PV_width  = filtrar_PV_values_1material(data_array)   
            array_z = discretize_z_1mat(step_z, PV_width)
            nz0, data_array = fotons_data(data_array, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial)
            data_array, array_G = Calcul_Photogeneration_BC(data_array, array_z, step_z, z_initial, z_final)

        if two_material:
            data_array, z_final, z_initial, PV_width = filtrar_PV_values_2materials(data_array, HTM_width)
            array_z, array_z_PV, array_z_HTM = discretize_z_2mat(step_z, PV_width, HTM_width)
            nz0, data_array = fotons_data(data_array, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial)
            data_array, array_G = Calcul_Photogeneration_BC_2mat(data_array, array_z, step_z, z_initial, z_final)
        
        if three_materials:
            data_array, z_final, z_initial, PV_width  = filtrar_PV_values_3materials(data_array, HTM_width, ETM_width)
            array_z,array_z_ETM = discretize_z_3mat(step_z, PV_width, HTM_width, ETM_width)
            nz0, data_array = fotons_data(data_array, array_num_fotons['lambda'], array_num_fotons['fotons_raig'], z_initial)
            data_array, array_G = Calcul_Photogeneration_BC_3mat(data_array, array_z, step_z, z_initial, z_final, array_z_ETM)
                    
        fotons_PV_path, integral, J_ph = resultats(data_array, array_G, array_z, S, q)
        
        df_plot = pd.DataFrame({'z($\mu m$)':array_z,'Generation Rate G(n/cm3)':array_G})
        df_data = pd.DataFrame({'Width(nm)':PV_width,'Absorbed Photons':integral, 'Ratio absorbed':integral/nz0*100.0, 'Check':integral/fotons_PV_path*100.0,'Photocurrent(mA/cm2)': J_ph}, index=[0])
    
        df_plot.to_csv(output_G, header=True, index=False, sep='\t', line_terminator='\n') 
        df_data.to_csv(output_data, header=True, index=False, sep='\t')

        output_G.close()
        output_data.close()
        
        print('-------------------------Results---------------------------------')
        print('Incident number of photons z=0 (1/m2): ', nz0)
        print('Photocurrent (mA/cm2) = ', J_ph)
        print('Number of absorbed photons (TOTAL) (1/m2) = ',integral)
        print('Number of photons of PV_path file(1/m2) = ', fotons_PV_path)
        
        print('-------------------------RATIOS------------------------------------')
        print('RATIO integral and pvpath (TOTAL)= ', integral/fotons_PV_path*100.0)
        print('RATIO absorbed/sent photons = ', integral/nz0*100.0)
        
        print('-------------------------PLOTS---------------------------------')
        array_z = array_z*1e6
        
        G_Plot(array_z, array_G, HTM_mat, parts[1].split('.txt')[0], (PV_width + HTM_width)*1e6, dir_outputs)
        
    plt.show()
