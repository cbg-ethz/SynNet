# SynNet
Hi.

This package includes two main components: *01_Optimize_Paramters* and *02_SeekNet*.


## 01_Optimize_Paramters
This package can be used to optimize circuit biochemical paramters using on-off ratio for a set of given circuits. 

You can run it by enterin the folder and running ```Launch_Global_Optimization.m```. This is tested on Matlab 2015. All parameters are set in ```Launch_Global_Optimization.m``` file, you can edit the file to change your choice of circuits, optimization parameters ranges, binarization point, input expression distributions, and other things. Advanced users can change the circuit evaluation functions in ```Auxilary_Functions``` folder to apply this software to other in-situ classification technologies with a different biochemical transfer function.

## 02_SeekNet
This package can be used to learn classifier circuit for a given dataset for a specific technology given fixed biochemical paramters (continuous mode), or using the general optimization strategy for any boolean-like circuit technology (boolean mode).


You can run an example case provided by running: ```Seek_Net('Example_Data/D04_EX_BreastCancer_Constraints.txt')```

Input: The example data includes three mandatory (Expresion data, Annotation metadata, and analysis Constraints) and two optional input files (mature miRNA sequences, and blacklisted miRNAs to be excluded). 

NOTE: There's one MEX file, fastAUC, used in this package which is precompiled for 64bit mac, you might need to recompile it for your computer, check ```02_SeekNet/Auxilary_Functions/fastAUC/install.m``` to learn how to do so.

Output: For each analysis the output inludes a set of figures with self displanatory names in .fig and .pdf format. Output also inludes the following text files:

1. ```[Analysis name]-Used_Constants.log```
A carbon-copy from the Constraint file used for the analaysis.
 
2. ```[Analysis name].txt```
This is the main file you want to look into. For each analysis case a file be produced like this that you can import to excel and has all the details of the analaysis, the performance, the circuit details the genes included and excluded. 

3. ```[Analysis name]-Prune.txt```
This is the same file as the pervious one, just made for the pruned circuit.

4. ```Summary_Report_[Constraint file name].txt```
The metadata file provided in the Constraints file might have more than one analysis, a summary of all of the analysis will be given in this file. This is usefull when you have several analyses done in a batch.

5. ```Summary_Report_[Constraint file name]-CV.txt```
This is similar to the file above except that it's only generated if cross-validation is on, and it will convey the cross-validation resutls.


I hope this helps to get you started. For problems, try again, if it fails please feel free to contact: pejman.m@gmail.com


best
Pej

