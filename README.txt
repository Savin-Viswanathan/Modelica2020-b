This download contains the OceanEngineering library with component models to simulate the hydrodynamic response of free floating and catenary moored cylindrical objects as discussed in the paper titled 'Modelica Component Models for Non-diffracting Floating Objects and Quasi-static Catenary Moorings' submitted to the American Modelica Conference 2020.

Instruction to generate the plots discussed in the above paper follows:
********************************************************************************************************************************************************************************
a. Load the current OceanEngineering library in OMEdit.

****Figure 5-a,b
a. Under SampleSimulations>FreeFloatingCylinder, simulate FreeFloatingCylinderInRegularWaves.
b. Plot variable 'z' under the cylindricalBuoy_Free1 model to get the plot for heave shown in Figure 5a of the paper.
c. Plot variable 'x' under the cylindricalBuoy_Free1 model to get the plot for surge shown in Figure 5b of the paper.
d. For orcaflex results, plots are generated from the orcaflex simulation file named 'FreeFloatingCylinderInRegularWaves'.

****Figure 6-a,b
a. Simulate  SampleSimulations>FreeFloatingCylinder>FreeFloatingCylinderInCurrent and plot 'x' under the cylindricalBuoy_Free1 for the surge in current.
b. Simulate SampleSimulations>FreeFloatingCylinder>FreeFloatingCylinderInWavesAndCurrent and plot 'x' under the cylindricalBuoy_Free1 for surge in waves and current.
d. Orcaflex plots are generated from the orcaflex simulation files named 'FreeFloatingCylinderInCurrent' and 'FreeFloatingCylinderInWavesAndCurrent'.

****Figure 7-a,b,c
a. Simulate SampleSimulations>FixedCylinder>FixedCylinderInCurrent and plot 'MF_x' under the cylindricalBuoy_Fixed1 for Morision loads in current.
b. Simulate SampleSimulations>FixedCylinder>FixedCylinderInRegularWaves and plot 'MF_x' under the cylindricalBuoy_Fixed1 for Morision loads in  waves.
c. Simulate SampleSimulations>FixedCylinder>FixedCylinderInRegularWavesAndCurrent and plot 'MF_x' under the cylindricalBuoy_Fixed1 for Morision loads in  waves and current.
d. Orcaflex plots are generated from the orcaflex simulations file named 'FixedCylinderInCurrent', 'FixedCylinderInWaves','FixedCylinderInWavesAndCurrent'. Since Orcaflex plots the connection X-force, which is the reaction to the morison force and hence opposite in magnitude, the direction of waves and current is specified at 180 to facilitate direct plotting.

****Figure 8-a,b,c
a. Simulate SampleSimulations>MooredCylinder>ThCalculator and plot X Vs Th in an array parametric plot to get the plots for horizontal tensions.
b. Specifying different values for the parameter spm_chain will generate plots for mooring chains of different specific masses.
c. Corresponding Orcaflex plots are generated from the orcaflex simulation files named 'HorizontalTensionsSpm10','HorizontalTensionsSpm16' and 'HorizontalTensionsSpm43', by placing the top end of the catenary at different x locations manually with z=0 in all cases.

****Figure 9
a. Set x=60 in the equation section of SampleSimulations>MooredCylinder>CylindricalBuoy_Fixed.
b. Simulate SampleSimulations>MooredCylinder>CatenaryShape.
c. In the results, under the catenaryMooring_MfC1 model, plot the variables z_lnk_plot Vs x_lnk_plot to get the catenary shape. This shape does not change depending on the specified current since this is a theoretical shape calculated based on the horizontal tension look-up table.
d. Re-simulate with x=80 for corresponding catenary shape.
e. Orcaflex plots are generated from the simulation file named 'HorizontalTensionsSpm10' by specifying the buoy x location accordingly and manually noting down the values of the nodal coordinates for both the cases x=60 and x=80, in the presence or absence of a current.

****Figure 10-a,b,c,d
a. Simulate SampleSimulations>MooredCylinder>MooredCylinderInCurrentMfC and plot variable 'x' under the cylindricalBuoy1 model to get the x-coordinate of the buoy in current and variable 'z3' go get the position of the keel as the buoy drifts away.
b. Simulate SampleSimulations>MooredCylinder>MooredCylinderInWavesAndCurrentMfC and repeat plotting steps above for corresponding plots.
c. Orcaflex plots are generated from the orcaflex simulation files named 'MooredBuoyInCurrent' and 'MooredBuoyInWavesAndCurrent' 

****Figure11-a,b,c
a. Simulate ampleSimulations>MooredCylinder>MooredCylinderInWavesAndCurrentMfCW and plot 'x' and 'z3' as above to get the surge and heave responses. Corresponding orcaflex results are plotted from the orcaflex simulation file named 'MooredBuoyInWaveAndCurrentReduced'.
b. Plot azl[2] under the catenaryMooring_MfCW1 component model to plot the vertical acceleartion of the second segment from the fairlead.
*************************************************************************************************************************************************
for clarification, send mail to savinvis@gmail.com
*************************************************************************************************************************************************