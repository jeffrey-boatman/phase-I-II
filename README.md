=======
# phase-I-II
Estimating toxcitity and efficacy from a phase I-II clinical trial and a follow up cohort for more efficient estimation

The code is based on Kusha's plan B MS project.

***Prior to running the code, the user must create the following directories in the same directory where the code files are located:
  'output'
  'RData'
  'tables'
  'plots'
  
 ***note: the R files were designed to run on departmental servers in parallel. The user must update these to run locally or remotely as appropriate.
 
To run the code:

  1 - create the directories as described above
  
  2 - run 'simulate-data.R' to simulate the datasets and save them to the 'RData' directory
  
  3 - run 'apply-estimators.R' to apply the estimators to the saved data sets and store the results in the 'output' directory
  
  4 - run 'plots-and-tables.R' to create the tables and write them to the 'tables' directory and to generate plots and save them to the 'plots' directory.


