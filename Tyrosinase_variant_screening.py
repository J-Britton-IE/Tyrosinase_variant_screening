# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:19:04 2021

@author: UCD-Jim
"""
## Standalone software for tyrosinase colony screening

import PySimpleGUI as sg
import matplotlib.pyplot as plt
import pandas as pd
# Note the matplot tk canvas import
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import screen_function as screen
import seaborn as sns


# VARS CONSTS:
_VARS = {'window': False}


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def delete_figure(figure):
    figure.get_tk_widget().forget()
    plt.close()
# \\  -------- PYSIMPLEGUI -------- //


AppFont = 'Any 16'
sg.theme('LightGrey')

font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 5}

plt.rc('font', **font)

tab1 = [[sg.Canvas(key ='figCanvas')]]
tab2 = [[sg.Canvas(key = 'figCanvas2')]]
tab3 = [[sg.Canvas(key = 'figCanvas3')]]
tab4 = [[sg.Canvas(key = 'figCanvas4')]]
tab5 = [[sg.Canvas(key = 'figCanvas5')]]
tab6 = [[sg.Canvas(key = 'figCanvas6')]]

buttons = [[sg.Frame('', [
           [sg.Button('Clear', use_ttk_buttons=True), 
           sg.Button('Help', use_ttk_buttons= True), 
           sg.Button('Exit', use_ttk_buttons= True)]], border_width = 5, element_justification='right')]]


data_outs  = [[sg.Frame('Data',
                        [[sg.TabGroup([[#sg.Tab('', tab0),
                                        sg.Tab('BCA', tab1), 
                                        sg.Tab('Quinone Endpoint', tab2), 
                                        sg.Tab('Control Quinone', tab3),
                                        sg.Tab('Quinone timecourse', tab4), 
                                        sg.Tab('Sample rate bar', tab5),
                                        sg.Tab('Sample rate heatmap', tab6)
                  ]],)],
                [sg.Text('High activity wells:'), sg.Input('', size = (50,2), key ='max_result')],
                [sg.Text('Control min:'), sg.Input('', size =(15, 2), key = 'control_min')],
                [sg.Text('Control mean:'), sg.Input('', size =(15, 2), key = 'control_mean')],
                [sg.Text('Control max:'), sg.Input('', size =(15, 2), key = 'control_max')]
                ], border_width = 5)]]

analyse_plot = [[sg.Frame('', [[sg.Button('Analyze', use_ttk_buttons= True), sg.Button('Plot', use_ttk_buttons= True)]], border_width= 5)]]

savings   = [[sg.Frame('',[[sg.Text('Results file name: '), sg.Input('', size=(55, 2), key = 'resultfilename')],
          
          [sg.FolderBrowse('Save folder', size=(10, 1)), sg.Input(size=(58, 1), key='resultfolder')], 
          
          [sg.Button('Save')]], border_width = 5)]]


layout = [[sg.FileBrowse('Activity File', size=(10, 1)), sg.Input('', size=(60, 2), key ='myFile')],
          
          [sg.FileBrowse('BCA File', size =(10, 1)), sg.Input('', size = (60, 2), key = 'BCA_file')],
          
          [sg.Column(analyse_plot, element_justification= 'left')],
          
          [ sg.Column(data_outs, element_justification='left', vertical_alignment='top')],
          
          [sg.Column(savings, element_justification='left', vertical_alignment = 'top')],
          
          
          [sg.Column(buttons, element_justification='right', vertical_alignment= 'bottom')]]
         
                                                                    
_VARS['window'] = sg.Window('Tyrosinase variant activity screen version 1',
                            layout,
                            finalize=True,
                            resizable=True)
# MAIN LOOP

while True:
    event, values = _VARS['window'].read(timeout=200)
    
    results = []
    
    
    
    if event == 'Analyze':
        
        if values['Activity File'] == '' or values['BCA File'] == '':
            sg.Popup('Please ensure an activity datafile and BCA datafile have been selected',
                 title='Warning - No datafile',
                 keep_on_top= True)
        
        elif values['Activity File'] != '' and values['BCA File'] != '':
            
            myActivityFile  = values['Activity File'];
            myActivityFile  = pd.read_excel(myActivityFile)
            myBCAFile  = values['BCA File'];
            myBCAFile  = pd.read_excel(myBCAFile)
            
            try:
                results           = screen.screen(myBCAFile, myActivityFile);
            except:
                sg.Popup('Unable to analyze.\n'
                         'Check input file and input parameters.',
                         title = 'Warning!',
                         keep_on_top= True)
                continue
                
            time                                       = results[0]
            BCA_ControlProtConc_activity_assay         = results[1]
            BCA_samplesProtConc                        = results[2]
            FinalQuinoneData_DF                        = results[3]
            controlQuinoneDta                          = results[4]
            sampleQuinoneData                          = results[5]
            ProtAdj_Qunione_DF                         = results[6]
            above_control_mean                         = results[7]
            above_control_max                          = results[8]
            above_control_min                          = results[9]
            MeanControlActivity                        = results[10]
            ControlActivityMax                         = results[11]
            ControlActivityMin                         = results[12]
        
    elif event == 'Plot':
                
        fig = plt.figure()
        BCA_heatmap = sns.heatmap(BCA_samplesProtConc, cmap= "RdPu", annot = True, annot_kws = {"fontsize": 8}, linewidth = 0.3);
        plt.title('Mg protein loaded per test well (BCA)')
        draw_figure(_VARS['window']['figCanvas'].TKCanvas, fig)
        
        fig2 = plt.figure()
        Quinone_heatmap = sns.heatmap(FinalQuinoneData_DF, cmap= "hot_r", annot = True, annot_kws = {"fontsize": 8}, linewidth = 0.3);
        plt.title('Final Quinone formation per well (mM)')
        draw_figure(_VARS['window']['figCanvas2'].TKCanvas, fig2)
    
        fig3 = plt.figure()
        plt.plot(time, controlQuinoneDta);
        plt.xlabel('Time (mins)');
        plt.ylabel('Quinone conc. mM');
        plt.title('Screening assay control wells Quinone formation');
        draw_figure(_VARS['window']['figCanvas3'].TKCanvas, fig3)
        
        
        fig4 = plt.figure()
        plt.plot(time, sampleQuinoneData);
        plt.xlabel('Time (mins)');
        plt.ylabel('Quinone conc. mM');
        plt.title('Screening assay sample wells Quinone formation');
        draw_figure(_VARS['window']['figCanvas4'].TKCanvas, fig4);
        
        
        def create_column_label(row):
            return f"{row['level_0']}_{row['level_1']}";
        
        melted_ProtAdjQuinone_DF = ProtAdj_Qunione_DF.stack().reset_index()
        melted_ProtAdjQuinone_DF['Column'] = melted_ProtAdjQuinone_DF.apply(create_column_label, axis = 1);
        
        fig5 = plt.figure()
        plt.bar(melted_ProtAdjQuinone_DF['Column'], melted_ProtAdjQuinone_DF[0])
        
        
        plt.ylabel('Quinone mM/min/mg protein');
        plt.title("Protein adjusted Quinone formation per well - Control limits")
        plt.plot([-10, 100],[ControlActivityMin, ControlActivityMin], "k--");
        plt.plot([-10, 100], [ControlActivityMax, ControlActivityMax], "k--");
        plt.plot([-10, 100], [MeanControlActivity, MeanControlActivity], "k--");
        draw_figure(_VARS['window']['figCanvas5'].TKCanvas, fig5)
        
        fig6 = plt.figure()
        ProtAdj_Quinone_heatmap = sns.heatmap(ProtAdj_Qunione_DF, cmap= "hot_r", annot = True, annot_kws = {"fontsize": 8}, linewidth = 0.3);
        plt.title('Protein adjusted Quinone formation per min (mM)')
        draw_figure(_VARS['window']['figCanvas6'].TKCanvas, fig6)
        
       
            
        _VARS['window']['max_result'](above_control_min);
        _VARS['window']['control_min'](ControlActivityMin);
        _VARS['window']['control_mean'](MeanControlActivity);
        _VARS['window']['control_max'](ControlActivityMax);
       
        
    elif event == 'Help':
        sg.Popup('Microtiter plate ctivity analysis tool.\n'
                 'Ensure plate set up as follows: \n'
                 ' - Row A = Blanks ONLY. \n'
                 ' - Every row should contain 1 protein concentration. \n'
                 ' - Columns 1 - 12 should contain replicates. \n',
                 title='Help',

                 keep_on_top= True)
    
    elif event == 'Clear':
        _VARS['window']['myFile']('')
        _VARS['window']['epsilon']('')
        _VARS['window']['replicateNumber']('')
        _VARS['window']['stepWidth']('')
        _VARS['window']['Row_B']('')
        _VARS['window']['Row_C']('')
        _VARS['window']['Row_D']('')
        _VARS['window']['Row_E']('')
        _VARS['window']['Row_F']('')
        _VARS['window']['Row_G']('')
        _VARS['window']['Row_H']('')
        _VARS['window']['max_result']('')
        _VARS['window']['resultfilename']('')
        _VARS['window']['resultfolder']('')
        _VARS['window']['r2_lim']('')
               
    elif event == 'Save':
        
        results_location = values['Save folder']
    
        _VARS['window']['resultfolder'](results_location)
      
        #save_spot = str(results_location) + "/" + str(values['resultfilename'] + '.xlsx')
        
        save_folder = str(values['resultfilename'])
        
        save_file_path = str(results_location) + "/" + str(save_folder)
        
        save_spot = str(results_location) + "/" + str(save_folder) + "/" + str(values['resultfilename'] + '.xlsx')
        
        os.mkdir(save_file_path)
        
        resultsFileName = str(values['resultfilename']) + '.xlsx'
        
        writer = pd.ExcelWriter(save_spot, engine='xlsxwriter')
        
        myActivityFile.to_excel(writer, sheet_name = 'Activity_Raw', header = False, index = False)
        
        myBCAFile.to_excel(writer, sheet_name = 'BCA_raw', header = False, index = False)
        
        BCA_samplesProtConc.to_excel(writer, sheet_name = 'BCA_result', header = True, index = True)
        
        FinalQuinoneData_DF.to_excel(writer, sheet_name = 'Activity_end', header = True, index = True)
        
        control_Q_with_time = pd.concat([time, controlQuinoneDta], axis = 1)
        
        control_Q_with_time.to_excel(writer, sheet_name = 'Control_quinone', header = True, index = False)
        
        sample_Q_with_time = pd.concat([time, sampleQuinoneData], axis = 1)
        
        sample_Q_with_time.to_excel(writer, sheet_name = 'Sample_quinone', header = True, index = False)
        
        ProtAdj_Qunione_DF.to_excel(writer, sheet_name = 'Sample_rates', header = True, index = True)
              
        writer.save() 
       
        fig.savefig(str(save_file_path) + '/' + '1_BCA_Heat_plot.png')
        fig2.savefig(str(save_file_path) + '/' + '2_Quinone_End_heat_plot.png')
        fig3.savefig(str(save_file_path) + '/' + '3_control_Quinone_plot.png')
        fig4.savefig(str(save_file_path) + '/' + '4_sample_quinone_plot.png')
        fig5.savefig(str(save_file_path) + '/' + '5_sample_rate_plot.png')
        fig6.savefig(str(save_file_path) + '/' + '6_sample_rate_heat_plot.png')
        
    elif event == sg.WIN_CLOSED or event == 'Exit':
        break
    
_VARS['window'].close()
