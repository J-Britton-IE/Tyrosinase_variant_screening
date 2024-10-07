# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:24:48 2023

@author: UCD
"""

## Colony activity screen function ##

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def screen(BCA_file, activity_file):
    
    epsilon               = 716; # Molar extinction coefficient for desired product. 716
    pathlength            = 0.622;

## From here is copied BCA_analysis.py #####################################

    proteinCoefficient = 0.9898

    dilutionFactor = 5

    BCA_File =  BCA_file;

    BCA_File     = BCA_File.drop(BCA_File.index[0:14]);

    BCA_File.dropna(how = 'all', axis = 1, inplace = True);

    BCA_File.fillna(0, inplace = True);

    BCA_rawData = BCA_File;

##Step 1 gather blanks

    BCA_blanks = BCA_rawData.drop(BCA_rawData.columns[2:13], axis =1);

    BCA_blanks = BCA_blanks.drop(BCA_blanks.columns[0], axis = 1);

    BCA_blanks = BCA_blanks.drop(BCA_blanks.index[0:5]);

    BCA_blankMean = BCA_blanks.mean(axis = 0);

    BCA_blankMean = BCA_blankMean.values[0];

##Step 2 gather controls

    BCA_controls = BCA_rawData.drop(BCA_rawData.columns[2:13], axis = 1);

    BCA_controls = BCA_controls.drop(BCA_controls.columns[0], axis = 1);

    BCA_controls = BCA_controls.drop(BCA_controls.index[0]);

    BCA_controls = BCA_controls.drop(BCA_controls.index[4:9]);

    BCA_controlMean = BCA_controls.mean(axis = 0);

##Step  3 adjust control values by blank values

    BCA_controlBlkAdj = BCA_controls.sub(BCA_blankMean);

##Step 4 divide by protein coefficient (from standard curve)

    BCA_controlProtConc = BCA_controlBlkAdj.div(proteinCoefficient);

##Step 5 multiply by dilution factor

    BCA_controlProtConc = BCA_controlProtConc*dilutionFactor;

    BCA_ControlProtConc_activity_assay = BCA_controlProtConc/10;

    BCA_ControlProtConc_activity_assay = BCA_ControlProtConc_activity_assay.to_numpy();

    BCA_samples = BCA_rawData.drop(BCA_rawData.columns[0:2], axis = 1);

    BCA_samples = BCA_samples.drop(BCA_samples.index[0]);

    BCA_samplesBlkAdj = BCA_samples.sub(BCA_blankMean);

    BCA_samplesProtConc = BCA_samplesBlkAdj.div(proteinCoefficient);

    BCA_samplesProtConc = BCA_samplesProtConc*dilutionFactor;

    BCA_ProtConc_activity_assay = BCA_samplesProtConc/10; 

    BCA_ProtConc_activity_assay = BCA_ProtConc_activity_assay.set_axis(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], axis = 0);

    BCA_ProtConc_activity_assay = BCA_ProtConc_activity_assay.set_axis(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'], axis = 'columns');

##Make a nice heatmap
    
    BCA_samplesProtConc[BCA_samplesProtConc < 0] =0;
    
    BCA_samplesProtConc = BCA_samplesProtConc.set_axis(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], axis = 0);

    BCA_samplesProtConc = BCA_samplesProtConc.set_axis(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'], axis = 'columns');

    plt.figure(num = "A");

    BCA_heatmap = sns.heatmap(BCA_samplesProtConc, cmap= "RdPu", annot = True, annot_kws = {"fontsize": 8}, linewidth = 0.3);

    plt.title('Mg protein loaded per test well (BCA)')

### Final steps --> Divide activity per min by Protein conc. loaded.

## 1 get file and sort data
    myFile     = activity_file;

    myFile     = myFile.drop(myFile.index[0:12]);

    myFile.dropna(how = 'all', axis = 1, inplace = True);

    myFile.fillna(0, inplace = True);

    rawData = myFile;

    rawData = rawData.drop(rawData.columns[0], axis = 1);

##Organize timepoints
    time = myFile.drop(myFile.columns[0:1], axis = 1); # Get new array with time values.

    time = time.drop(time.columns[1::], axis = 1);

    time = time.drop(time.index[0]);

    time.columns = ['Times'];

    time_numeric = len(time);

    time_numeric_list = []

    for i in range(time_numeric):
        tee = (i * 1);
        time_numeric_list.append(tee)

    time_numeric_list = pd.DataFrame(time_numeric_list);

    time_numeric_list.columns = ['Times']

    time = time_numeric_list;

##Organize blanks
    blanks = rawData.drop(rawData.columns[3::], axis =1);

    blanks = blanks.drop(blanks.columns[0], axis = 1);

    blanks = blanks.drop(blanks.index[0]);

    blank_means = blanks.mean(axis = 'columns');


##Organize controls and subtract blank
    controls = rawData.drop(rawData.columns[0:3], axis = 1);

    controls = controls.drop(controls.index[0]);

    controls = controls.drop(controls.columns[4::], axis = 1);

    controls.drop(controls.index[0]);

#controls = controls.drop(controls.columns[0], axis = 1);

    controlBlkAdj = controls.subtract(blank_means, axis = 'index');

    controlBlkAdj.reset_index(inplace = True);

##Convert control wells OD400 to quinone formed.

    controlQuinone_mM = controlBlkAdj/(epsilon * pathlength); #convert OD value to uM

    controlQuinone_mM = controlQuinone_mM * 1000; # convert to mM

    controlQuinoneDta = controlQuinone_mM;

    controlQuinoneDta = controlQuinoneDta.drop(controlQuinoneDta.columns[0], axis = 1);

    #controlQuinoneDta.rename(columns = {0: "Time (mins)", 2: "1", 3: "2", 4: "3"}, inplace = True);

##Organize samples

    samples = rawData.drop(rawData.columns[0:7], axis = 1);

    samples.columns = samples.iloc[0]

    samples = samples[1:]

    sampleBlkAdj = samples.subtract(blank_means, axis = 'index');

    sampleBlkAdj.reset_index(inplace = True);

    sampleQuinone_mM = sampleBlkAdj/(epsilon * pathlength);

    sampleQuinone_mM = sampleQuinone_mM * 1000;

    #sampleQuinoneData = pd.concat([time, sampleQuinone_mM], axis = 1, ignore_index = True);
    sampleQuinoneData = sampleQuinone_mM;

    sampleQuinoneData = sampleQuinoneData.drop(sampleQuinoneData.columns[0], axis = 1);

    #sampleQuinoneData.rename(columns = {0: "Time (mins)"}, inplace = True);

    samplesFirstRow = sampleQuinoneData.iloc[[0]].values[0];

    sampleQuinoneData = sampleQuinoneData.apply(lambda row: row - samplesFirstRow, axis =1);
    
    sampleQuinoneData[sampleQuinoneData < 0] = 0;
## 
    sampleQuinonePerMin = sampleQuinoneData/60;

    #sampleQuinonePerMin = sampleQuinonePerMin.drop(sampleQuinonePerMin.columns[0], axis = 1);

##Making heat map of final quinone values

    sampleValsFinal = sampleQuinonePerMin.values[-1].tolist();

    #sampleQuinoneData2 = sampleQuinoneData.drop(sampleQuinoneData.columns[0], axis = 1);

    sampleValsFinal_array = sampleQuinoneData.values[-1];

    FinalQuinoneData_DF = pd.DataFrame(np.reshape(sampleValsFinal_array, (8, 11)), columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'], index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

    FinalQuinoneData_DF[FinalQuinoneData_DF < 0] =0; #Turn negative values into zer0


## Plot controls quinone formation

    control_OD_plot = plt.figure();
    plt.plot(time, controlQuinoneDta);
    plt.xlabel('Time (mins)');
    plt.ylabel('Quinone conc. mM');
    plt.title('Screening assay control wells Quinone formation');

## Plot sample Quinone formation

    sample_OD_plot = plt.figure()
    plt.plot(time, sampleQuinoneData);
    plt.xlabel('Time (mins)');
    plt.ylabel('Quinone conc. mM');
    plt.title('Screening assay sample wells Quinone formation');
##Ok - Blanks and controls are plotted. Samples are plotted but need to re-associate the sample number with plate location. 
##Before this, need to calculate Quinone formation rates per well. Then rearrange data into same format as BCA. Then divide Quinone formation by protein loaded.
##Then select for better than controls.

## Easy way - adjust all values by the T0 for that well, Then take endpoint sample and divide by 60 to get mM/min?
## Hard way - for each well, find a place where the rate is linearly increasing at highest rate. Divide by selected time frame. Voila.


## Subtract first row from all rows in control data

    controlFirstRow = controlQuinoneDta.iloc[[0]].values[0];

    controlQuinoneDta = controlQuinoneDta.apply(lambda row: row - controlFirstRow, axis =1);

    controlQuinonePerMin = controlQuinoneDta/60;

    #controlQuinonePerMin = controlQuinonePerMin.drop(controlQuinonePerMin.columns[0], axis = 1);

    controlValsFinal = controlQuinonePerMin.values[-1].tolist();

    ##Controls sorted - blanked by T0 for each well. Divided last timepoint by 60 to get rough estimate of mM/min Quinone and saved as list.

## Ok - so have a list of all samples, all controls and their rough Quninone formation per min. Now adjust quinone per protein loaded.

    arr = np.array(sampleValsFinal);
    
    sampleValsFinal_DF = pd.DataFrame(np.reshape(arr, (8, 11)), columns = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'], index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

## From here is copied BCA_analysis.py #####################################

### Final steps --> Divide activity per min by Protein conc. loaded.

    ProtAdj_Qunione_DF = sampleValsFinal_DF / BCA_ProtConc_activity_assay; 

    controlValsFinal = pd.DataFrame(controlValsFinal);

    BCA_ControlProtConc_activity_assay = pd.DataFrame(BCA_ControlProtConc_activity_assay);

    controlProtAdj_Quinone = controlValsFinal / BCA_ControlProtConc_activity_assay;

    controlProtAdj_Quinone = pd.DataFrame(controlProtAdj_Quinone);

## Ok - so far, the values for each test well (Protein adjusted Quinone formation per min) and the control wells are saved. 
## Next step - find out if any of the test wells have a higher value than the controls.

    MeanControlActivity = controlProtAdj_Quinone.mean();

    MeanControlActivity = MeanControlActivity.values[0];

    ControlActivityMax = controlProtAdj_Quinone.max();

    ControlActivityMax = ControlActivityMax.values[0];

    ControlActivityMin = controlProtAdj_Quinone.min();

    ControlActivityMin = ControlActivityMin.values[0];

    ProtAdj_Qunione_DF[ProtAdj_Qunione_DF < 0] =0; #Turn negative values into zero

    ProtAdjQuinone_plot = ProtAdj_Qunione_DF.plot.bar(rot = 0, legend = False);

    ProtAdjQuinone_plot.set_ylabel('Quinone mM/min/mg protein');

    ProtAdjQuinone_plot.set_title("Protein adjusted Quinone formation per well - Control limits")

    ProtAdjQuinone_plot.plot([-10, 100],[ControlActivityMin, ControlActivityMin], "k--");

    ProtAdjQuinone_plot.plot([-10, 100], [ControlActivityMax, ControlActivityMax], "k--");

    ProtAdjQuinone_plot.plot([-10, 100], [MeanControlActivity, MeanControlActivity], "k--");

    plt.figure(num= "C");

    ProtAdj_Quinone_heatmap = sns.heatmap(ProtAdj_Qunione_DF, cmap= "hot_r", annot = True, annot_kws = {"fontsize": 8}, linewidth = 0.3);

    plt.title('Protein adjusted Quinone formation per min (mM)')

##Come back to this - trying to add in some text to identify control lines
#ProtAdjQuinone_plot.text(ControlActivityMax, ControlActivityMax, "Control max", ha = 'center', bbox = dict(facecolor = 'w', alpha=0.5));

#ProtAdjQuinone_plot.text(ControlActivityMin, ControlActivityMin, "Control min", ha = 'center', bbox = dict(facecolor = 'w', alpha=0.5));

#ProtAdjQuinone_plot.text(MeanControlActivity, MeanControlActivity, "Control mean", ha = 'center', bbox = dict(facecolor = 'w', alpha=0.5));
    
##Ok so now have all my values, just need to select good sample wells. 
##Step 1 - find any values above the mean of controls.
##Step 2 - Make the graph nice, want to illustrate the controls on the bar char with the samples
#        - Make it so that there is a dotted line going from min and max control activities across the bar.


##ChatGPT function to return indexes of wells with activity higher than X.
    def select_values_above_control(df, control_value):
        selected_positions = []
        
        for column in df.columns:
            for index, value in df[column].items():
                if value > control_value:
                    selected_positions.append((index, column))
                
        return selected_positions


    above_control_mean = select_values_above_control(ProtAdj_Qunione_DF, MeanControlActivity);

    above_control_max = select_values_above_control(ProtAdj_Qunione_DF, ControlActivityMax);

    above_control_min = select_values_above_control(ProtAdj_Qunione_DF, ControlActivityMin);
    
    above_control_min = above_control_min;


    return time, BCA_ControlProtConc_activity_assay, BCA_samplesProtConc, FinalQuinoneData_DF, controlQuinoneDta, sampleQuinoneData, ProtAdj_Qunione_DF,  above_control_mean,  above_control_max,  above_control_min, MeanControlActivity, ControlActivityMax, ControlActivityMin 
