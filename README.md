# Tic timing supporting data and files

  * Reads text file of timing of tics, 
  * Produces tic modulated random walk trajectory
  * Determines fractal dimension of the trajectory using Fourier analysis

## Inputs

  * Max time - the maximum time at which a tic can be recorded, in seconds. 
   	 * Example: If tics are monitored for ten minutes, the max time is 600 seconds.
  * dt - time resolution of tic detection, in seconds. All tic detections will be rounded to the same decimal place as dt. Additionally, any consecutive tics which are recorded as occuring at the same time will be separated by dt.
  	 *  Example: Tics can only be distinguished if they occur more than 0.2 seconds apart (dt = 0.2). All recorded tic times will be rounded to the first decimal place. If two consecutive tics are recorded at 10 seconds, the first tic time will remain at 10 seconds and the second tic time will be changed to 10.2 seconds. 
  *  Filename - full path to input data file. Each row reflects the tic timing of an individual patient. The first column is the patient identifier, and subsequent columns contain times when a tic is detected, in seconds. Must be space-delimited or tab-delimited.
  	 *  Example: See [example_input.txt](https://github.com/beelerpayton/tic_timing_supporting_data_and_files/blob/main/example_input.txt)
  *  Output directory - directiry where output data will be written. Must not end in '/'.
  	 *  Good Example: /Users/username/Desktop
  	 *  Bad Example: /Users/username/Desktop/

## Outputs

  * Df.csv - Shows patient ID and corresponding fractal dimension.
  * trajectories.csv - First column is time (in seconds), subsequent columns are the position of walkers at the time given by the first column. First row is header with patient ID.
  * temp.txt - temporary file that reads input data line-by-line. Can be discarded after code has finished running.

## Quality assurance/troubleshooting

  * Values for fractal dimension should be roughly between 1 and 2. If values lie significantly out of this range, verify that trajectory is being generated properly by examining the trajectory (stored in trajectory.csv output file).
  * Output to screen should be the number of trajectories completed followed by the patient ID.
   	 * Example screen output: <br />
   	   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1 patient_ID_#1 <br />
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 patient_ID_#2 <br />
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.... <br />
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;.... <br />
    If the patient ID does not match what is expected, the input file is not being read correctly, and the input file format should be checked. 

## Installation

This code is written in C++, and can be downloaded and compiled using Xcode (Mac), Microsoft Visual Studio (Windows), or another appropriate compiler.

## Publications

  * [Beeler P., Jensen N.O., Kim S., Robichaux-Viehoever A., Schlaggar B., Greene D., Black K., and Chakrabarty R.K. (2022) Fractality of tics as a quantitative assessment tool for Tourette syndrome. Journal of the Royal Society Interface.](https://doi.org/10.1098/rsif.2021.0742)
