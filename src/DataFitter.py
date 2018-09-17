"""
DataFitter
Created by Jerry LaRue, larue@chapman.edu, 1/2016
Last modified by Jerry LaRue, larue@chapman.edu, 9/2018
"""

import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import GaussFitter

class SpectraFitter ( object ) :
    
    def __init__ ( self ) :
        
        pass
    
    def Go ( self, Parameters ) :
        
        self.Folder_Input = Parameters['Folder_Analysis']
        self.File_Input = Parameters['File_Input']
        self.Folder_Output = Parameters['Folder_Analysis']
        self.File_Output = Parameters['File_Output']
        
        print('----------------------------------------')
        print('Fitting spectra from: ')
        print(self.Folder_Input + self.File_Input)
        
        ##### Check for files #####
        
        self.Success = True
        if os.path.isfile(self.Folder_Input + self.File_Input) :
            f = h5py.File(self.Folder_Input + self.File_Input, 'r')
        else :
            self.Success = False
            print(self.File_Input + ' file missing')
        
        if os.path.isfile(self.Folder_Input + self.File_Input) :
            
            # General
            if not 'runs' in f :
                self.Success = False
                print('Run list missing')
            
            # Energy
            if not 'BinnedData/E_bin_centers' in f :
                self.Success = False
                print('Energy values missing')
            
            # Delay
            if not 'BinnedData/delays_fs' in f :
                self.Success = False
                print('Delay Values missing')
            
            # Spectra
            if not 'BinnedData/XAS_2dmatrix' in f :
                self.Success = False
                print('Spectral data missing')
            if Parameters['Fit_Spectra_NumGauss'] == 0 :
                self.Success = False
                print('Need to have at least 1 gauss to fit spectra')
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                if not par.Spectra_Reference in f :
                    self.Success = False
                    print('Reference spectra missing')
        
        if self.Success :
            
            ##### Load data #####
            
            # General
            Runs = f['runs'][...]
            
            # Energy
            Energy_Values = f['BinnedData/E_bin_centers'][...]
            
            # Delay
            Delay_Values = f['BinnedData/delays_fs'][...]
            
            # Spectra
            Spectra_Values = f['/BinnedData/XAS_2dmatrix'][...]
            Spectra_Values = np.transpose(Spectra_Values)
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                Spectra_Reference = f[par.Spectra_Reference][...]
            
            f.close()
            
        if self.Success :
            
            ###### Trim data #####
            
            # Data Range
            X_Index_Min = (np.abs(Energy_Values - Parameters['Fit_ROI_Min'])).argmin()
            X_Index_Max = (np.abs(Energy_Values - Parameters['Fit_ROI_Max'])).argmin()
            
            # Data
            Spectra_Values = Spectra_Values[X_Index_Min:X_Index_Max]
            Spectra_Values = np.transpose(Spectra_Values)
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                Spectra_Reference = Spectra_Reference[X_Index_Min:X_Index_Max]
                Spectra_Reference = np.transpose(Spectra_Reference)
            Energy_Values = Energy_Values[X_Index_Min:X_Index_Max]
            
            ##### Fit reference spectra #####
            
            Parameters_Fit = np.zeros((0))
            Parameters_Fit_Reference = np.zeros((0))
            Parameters_Fixed = np.zeros((0))
            
            # Fit parameters for reference data
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                if Parameters['Fit_Reference_NumGauss'] > 0 :
                    Parameters_Fit = np.append(Parameters_Fit, [Parameters['Fit_Reference_Amplitude1'],Parameters['Fit_Reference_Position1'],Parameters['Fit_Reference_Width1']])
                    Parameters_Fixed = np.append(Parameters_Fixed, [True, True, True])
                if Parameters['Fit_Reference_NumGauss'] > 1 :
                    Parameters_Fit = np.append(Parameters_Fit, [Parameters['Fit_Reference_Amplitude2'],Parameters['Fit_Reference_Position2'],Parameters['Fit_Reference_Width2']])
                    Parameters_Fixed = np.append(Parameters_Fixed, [True, True, True])
                if Parameters['Fit_Reference_NumGauss'] > 2 :
                    Parameters_Fit = np.append(Parameters_Fit, [Parameters['Fit_Reference_Amplitude3'],Parameters['Fit_Reference_Position3'],Parameters['Fit_Reference_Width3']])
                    Parameters_Fixed = np.append(Parameters_Fixed, [True, True, True])
                Parameters_Fit_Found = GaussFitter.multigaussfit(Energy_Values, Spectra_Reference, ngauss=Parameters['Fit_Reference_NumGauss'], params=Parameters_Fit)
                Parameters_Fit_Reference = Parameters_Fit_Found[0]
                
            # Fit parameters for spectral data
            if Parameters['Fit_Spectra_NumGauss'] > 0 :
                Parameters_Fit = np.append(Parameters_Fit_Reference, [Parameters['Fit_Spectra_Amplitude1'],Parameters['Fit_Spectra_Position1'],Parameters['Fit_Spectra_Width1']])
                Parameters_Fixed = np.append(Parameters_Fixed, [Parameters['Fit_Spectra_Amplitude1_Fix'], Parameters['Fit_Spectra_Position1_Fix'], Parameters['Fit_Spectra_Width1_Fix']])
            
            if Parameters['Fit_Spectra_NumGauss'] > 1 :
                Parameters_Fit = np.append(Parameters_Fit, [Parameters['Fit_Spectra_Amplitude2'],Parameters['Fit_Spectra_Position2'],Parameters['Fit_Spectra_Width2']])
                Parameters_Fixed = np.append(Parameters_Fixed, [Parameters['Fit_Spectra_Amplitude2_Fix'], Parameters['Fit_Spectra_Position2_Fix'], Parameters['Fit_Spectra_Width2_Fix']])
            
            if Parameters['Fit_Spectra_NumGauss'] > 2 :
                Parameters_Fit = np.append(Parameters_Fit, [Parameters['Fit_Spectra_Amplitude3'],Parameters['Fit_Spectra_Position3'],Parameters['Fit_Spectra_Width3']])
                Parameters_Fixed = np.append(Parameters_Fixed, [Parameters['Fit_Spectra_Amplitude3_Fix'], Parameters['Fit_Spectra_Position3_Fix'], Parameters['Fit_Spectra_Width3_Fix']])
            
            ##### Fit signal spectra #####
            
            Fit_Energy_Values = np.zeros((0))
            j = 0
            while Parameters['Fit_ROI_Min'] + j * Parameters['Fit_X_Delta'] <= Parameters['Fit_ROI_Max'] :
                Fit_Energy_Values = np.append(Fit_Energy_Values, Parameters['Fit_ROI_Min'] + j * Parameters['Fit_X_Delta'])
                j = j + 1
            Fit_Spectra_Values = np.zeros((len(Delay_Values),len(Fit_Energy_Values)))
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                Fit_Spectra_Reference = np.zeros((len(Fit_Energy_Values)))
                Fit_Data_Difference = np.zeros((len(Delay_Values),len(Fit_Energy_Values)))
                Fit_Data_Sum = np.zeros((len(Delay_Values),len(Fit_Energy_Values)))
                Fit_Data_Contrast = np.zeros((len(Delay_Values),len(Fit_Energy_Values)))
                j = 0
                while j < len(Fit_Energy_Values) :
                    k = 0
                    while k < Parameters['Fit_Reference_NumGauss'] :
                        Fit_Spectra_Reference[j] = Fit_Spectra_Reference[j] + Parameters_Fit[3 * k] * np.exp( -(Fit_Energy_Values[j] - Parameters_Fit[3 * k + 1])**2 / (2.0*Parameters_Fit[3 * k + 2]**2) )
                        k = k + 1
                    j = j + 1
            
            i = 0
            while i < len(Delay_Values) :
                
                # Fit the Data
                Intensities = np.array([0.1])
                ChiSquares = np.zeros((0))
                Streak = 0
                Go = True
                Counter = 0
                while Go :
                    if Streak > 10 :
                        Go = False
                        Index = (np.abs(ChiSquares - min(ChiSquares))).argmin()
                    else :
                        Index = -1
                    Parameters_Fit_Intensities = np.zeros((0))
                    Parameters_Fit_Intensities = np.append(Parameters_Fit_Intensities,Parameters_Fit)
                    if Parameters['Fit_Reference_NumGauss'] > 0 :
                        j = 0
                        while j < Parameters['Fit_Reference_NumGauss'] :
                            Parameters_Fit_Intensities[j * 3] = Parameters_Fit[j * 3] * Intensities[Index]
                            j = j + 1
                    else :
                        Go = False
                    Parameters_Fit_Found = GaussFitter.multigaussfit(Energy_Values, Spectra_Values[i], ngauss=Parameters['Fit_Reference_NumGauss']+Parameters['Fit_Spectra_NumGauss'], params=Parameters_Fit_Intensities, err=True, fixed=Parameters_Fixed)
                        #limitedmin=[True,False,True,True,False,True,True,False,True,True,False,True])
                    Parameters_Fit_Spectra = Parameters_Fit_Found[0]
                    Counter_Max = 500
                    if Counter > Counter_Max :
                        Go = False
                        print('Warning: Fit failed after ' + str(Counter_Max) + ' iterations.')
                    if Go :
                        if len(ChiSquares) > 0 :
                            if Parameters_Fit_Found[-1] > ChiSquares[-1] :
                                Streak = Streak + 1
                        ChiSquares = np.append(ChiSquares,Parameters_Fit_Found[-1])
                        Intensities = np.append(Intensities, Intensities[-1] + 0.01)
                    Counter = Counter + 1
                if self.Success :
                    
                    ##### Make fits #####
                    
                    j = 0
                    while j < len (Fit_Energy_Values) :
                        k = 0
                        while k < Parameters['Fit_Reference_NumGauss'] + Parameters['Fit_Spectra_NumGauss'] :
                            Fit_Spectra_Values[i][j] = Fit_Spectra_Values[i][j] + Parameters_Fit_Spectra[3 * k] * np.exp( -(Fit_Energy_Values[j] - Parameters_Fit_Spectra[3 * k + 1])**2 / (2.0*Parameters_Fit_Spectra[3 * k + 2]**2) )
                            k = k + 1
                        j = j + 1
                    if Parameters['Fit_Reference_NumGauss'] > 0 :
                        Fit_Data_Difference[i,:] = Fit_Spectra_Values[i,:] - Fit_Spectra_Reference
                        Fit_Data_Sum[i,:] = Fit_Spectra_Values[i,:] + Fit_Spectra_Reference
                        Fit_Data_Contrast[i,:] = Fit_Data_Difference[i,:] / Fit_Data_Sum[i,:] * 100
                    
                    ##### Plot data #####
                    
                    # Figure Data
                    Fit_Spectra_Reference_Gauss = np.zeros((Parameters['Fit_Reference_NumGauss'],len(Fit_Energy_Values)))
                    Fit_Spectra_Values_Gauss = np.zeros((Parameters['Fit_Reference_NumGauss'] + Parameters['Fit_Spectra_NumGauss'], len(Fit_Energy_Values)))
                    j = 0
                    while j < len(Fit_Energy_Values) :
                        k = 0
                        while k < Parameters['Fit_Reference_NumGauss'] :
                            Fit_Spectra_Reference_Gauss[k][j] = Fit_Spectra_Reference_Gauss[k][j] + Parameters_Fit_Reference[3 * k] * np.exp( -(Fit_Energy_Values[j] - Parameters_Fit_Reference[3 * k + 1])**2 / (2.0*Parameters_Fit_Reference[3 * k + 2]**2) )
                            k = k + 1
                        k = 0
                        while k < Parameters['Fit_Reference_NumGauss'] + Parameters['Fit_Spectra_NumGauss'] :
                            Fit_Spectra_Values_Gauss[k][j] = Fit_Spectra_Values_Gauss[k][j] + Parameters_Fit_Spectra[3 * k] * np.exp( -(Fit_Energy_Values[j] - Parameters_Fit_Spectra[3 * k + 1])**2 / (2.0*Parameters_Fit_Spectra[3 * k + 2]**2) )
                            k = k + 1
                        j = j + 1
                    
                    # Reference & Signal
                    plt.figure(figsize = [4,4])
                    plt.plot(Energy_Values, Spectra_Values[i])
                    plt.plot(Fit_Energy_Values, Fit_Spectra_Values[i])
                    if Parameters['Fit_Reference_NumGauss'] > 0 :
                        plt.plot(Energy_Values, Spectra_Reference)
                        plt.plot(Fit_Energy_Values, Fit_Spectra_Reference)
                        plt.plot(Fit_Energy_Values, Fit_Data_Difference[i])
                    plt.xlabel('Energy')
                    plt.ylabel('Intensity [au]')
                    plt.title('Data ' + str(int(Delay_Values[i])) + ' fs')
                    plt.show()
                    
                    # Total intensities
                    Intensity_Reference = 0
                    Intensity_Signal = 0
                    j = 0
                    while j < Parameters['Fit_Reference_NumGauss'] + Parameters['Fit_Spectra_NumGauss'] :
                        if j < Parameters['Fit_Reference_NumGauss'] :
                            Intensity_Reference = Intensity_Reference + Parameters_Fit_Reference[j * 3] * Parameters_Fit_Reference[j * 3 + 2]
                        Intensity_Signal = Intensity_Signal + Parameters_Fit_Spectra[j * 3] * Parameters_Fit_Spectra[j * 3 + 2]
                        j = j + 1
                    
                    print('********************')
                    print('Delay = ' + str(Delay_Values[i]) + ' fs')
                    if Parameters['Fit_Reference_NumGauss'] > 0 :
                        print('Reference Intensity = ' + str(Intensity_Reference))
                    print('Signal Intensity = ' + str(Intensity_Signal))
                i = i + 1
            
            ##### Save to file #####
            
            print('********************')
            
            # Output File
            if self.Folder_Output + self.File_Output == self.Folder_Input + self.File_Input :
                f = h5py.File(self.Folder_Output + self.File_Output, 'a')
            else :
                if not os.path.exists(self.Folder_Output):
                    os.makedirs(self.Folder_Output)
                f = h5py.File(self.Folder_Output + self.File_Output, 'w')
            dt = h5py.special_dtype(vlen=bytes)
            
            # Group name            
            if self.Folder_Output + self.File_Output == self.Folder_Input + self.File_Input :
                print('Appending fits to:')
                print(self.Folder_Output + self.File_Output)
            else :
                print('Saving fits to:')
                print(self.Folder_Output + self.File_Output)
                
                # General
                dataSet = f.create_dataset('Run_List', data = Runs, dtype = 'int32')
                
                # Energy
                dataSet = f.create_dataset('Energy/Values', data = Energy_Values, dtype = np.dtype('float64'))
                
                # Delay
                dataSetText = f.create_dataset('Delay/Values', data = Delay_Values, dtype = np.dtype('float64'))
            
            # Clear data for overwriting
            if 'Fit/' in f :
                del f['Fit/']
            
            # Energy
            dataSet = f.create_dataset('Fit/Energy_Values', data = Fit_Energy_Values, dtype = np.dtype('float64'))
            
            # Fits
            dataSet = f.create_dataset('Fit/Signal', data = Fit_Spectra_Values, dtype = np.dtype('float64'))
            if Parameters['Fit_Reference_NumGauss'] > 0 :
                dataSet = f.create_dataset('Fit/Reference', data = Fit_Spectra_Reference, dtype = np.dtype('float64'))
            
            f.close()
            
            print('Done\n')
            
        else :
            print('Fitting cancelled\n')
    
    def Parameters ( self ) :
    
        if self.Success :
            return Parameters
        else :
            return ''
    
    def Folder_Output ( self ) :
        
        if self.Success :
            return self.Folder_Output
        else :
            return ''
    
    def File_Output ( self ) :
        
        if self.Success :
            return self.File_Output
        else :
            return ''
    
    def Success ( self ) :
        
        return self.Success