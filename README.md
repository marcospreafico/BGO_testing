# Setup 
Before doing anything run the `setup.sh` script to create the required folders to run 

# Getting data
To get data, copy them from the machine in SGM. If connected to INFN web or using a VPN, you can copy them using the following command: 
`scp streamdaq@193.206.147.141:/home/streamdaq/Documents/BGO_test/DAQ/run_name/RAW/SDataR_run_name.root ./data/` 
If not connected to the INFN web or using a VPN, ask to an adult to get the data for you

I strongly suggest to rename the files getting rid of the CoMPASS naming: `SDataR_260824.root` -> `260824.root`

# Running the reconstruction
In order to run the reconstruction follow the next steps. 

Note that in order to change the run to be used you need to change the file name in all the `.C` files. 

### 1. Phe calibration
Run the `calib_phe.C` script on a `phe_xxxxxx.root` file. This creates a file in the `calib` folder with the phe charge for each channel and a pdf file in the `pdf` folder with the fit results. 
**Important** always check that all the fits are good otherwise some channels may be miscalibrated. 
### 2. Data preprocessing 
Run the `preprocess.C` script on a data `xxxxxx.root` file. This creates a `xxxxxx_preprocessed.root` file in the `out` folder and a pdf with all the monitoring plots in the `pdf` folder.
**Note** The code requires a calibration file with the same name as the run. Since not all the time a calibration has been done before the run start, a workaround is to copy a calibration file for the same batch with the same name of the run you want to analyse. 
### 3. Event building
Run the `event_building.C` script. I strongly suggest to run in batch mode with `root -l -q -b event_building.C` as this script creates hundreds of canvases. The output is a `xxxxxx_evt.root` file containing the summarised information in the `out` folder and a file with all the canvases produced in the `pdf` folder. 
### 4. Data fitting 
If multiples files are available, merge them with `hadd` to create a single file with all the data available. Then use the `ana_fit.C` script to perform the fit to get the attenuation length. Again run this code in batch mode. The output in `pdf` will contain the results of all the fits, as well as the estimate of the LY trend. 

