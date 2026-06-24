## Brief description
This code:  

* Automatically fits time-dependent correlations using the ***GNUPLOT*** software [\[1\]](https://gnuplot.sourceforge.net/). Correlations must be provided by the user.
* Uses the fitted parameters to derive physical quantities using models that are related to the type of input data under consideration.   
* Estimate errors for the fitted parameters and the physical quantities (only if the input correlations contain error bars).  

The implemented types of input data that can be fitted are:   

* Orientational Correlation Functions (OCFs) of reactive [\[2\]](https://www.sciencedirect.com/science/article/pii/S0006349580851095) and non-reactive species [\[3\]](https://pubs.acs.org/doi/10.1021/jacs.1c08552).  
* Transfer Correlation Functions (TCFs) of reactive species [\[4\]](https://pubs.acs.org/doi/10.1021/jp0554036).  
* Special Pair Correlation Functions (SPCFs) involving the closest donor/acceptor pair in reactive systems.  
* Mean Square Displacements (MSD).  

For a detailed description of the implemented directives for data analysis, the user is provided with a [manual](./manual.pdf). 

## Disclaimer
The ALC does not fully guarantee the code is free of errors and assumes no legal responsibility for any incorrect outcome or loss of data.

## Author 
Ivan Scivetti (ALC, SCD, STFC)

## Structure of files and folders
ALC_FACT contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains the tests files (in .tar format) needed for CI purposes. The user should execute the available scripts of the [***tools***](./tools) folder to run the test automatically and verify the code has been installed properly (see the [build_code.md](./build_code.md) file for instructions).  
* [***examples***](./examples): folder with example cases to help familiarising with the code.  
* [***scripts***](./scripts): contains scripts for data processing.  
* [***source***](./source): contains the source code. Files have the *.F90* extension  
* [***tools***](./tools): shell files for building, compiling and testing the code automatically.  
* [.gitignore](./.gitignore): instructs Git which files to ignore.  
* [CMakeLists.txt](./CMakeLists.txt): sets the framework for code building and testing with CMake. This file must ONLY be modified to add test cases.  
* [LICENSE](./LICENSE): BSD 3-Clause License for ALC_FACT.
* README.md: this file.  
* [build_code.md](./build_code.md): steps to build, compile and run tests using the CMake platform.    
* [manual.pdf](./manual.pdf): ALC_FACT manual.  

## Dependencies
The user must have access to the following software (locally):

* GNU-Fortran (11.2.0) or Intel-Fortran (ifx 2023.1.0)
* Gnuplot (5.4)
* Cmake (3.16)
* Make (4.2.1)
* Git (2.34.1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work.

## Getting started

### Obtaining the code
The user can clone the code locally by executing the following command with the SSH protocol
```sh
$ git clone git@github.com:stfc/alc_fact.git
```
Instead, if the user wants to use the HTTPS protocol it must execute
```sh
$ git clone https://github.com/stfc/alc_fact.git
```
Both ways generate the ***alc_fact*** folder as the root directory. Alternatively, the code can be downloaded from any of the available assets.


### Building and testing the code
Details can be found in file [build_code.md](./build_code.md)

### Using the code
To execute the code, the user needs the SETTINGS file (with directives to set the fitting procedure) and the input file that contains correlation data to be fitted. The name of such file must be specified in the SETTINGS file with the ***filename*** directive. Relevant information of the defined settings as well as the fitted parameters and the computed physical quantities are printed in the OUTPUT file. The user is referred to the [manual](./manual.pdf) of the code for a detailed description of the implemented functionalities. In the [***examples***](./examples) folder the user will find several examples to help familiarising with the code.

## Acknowledgements  
* Ada Lovelace Centre for funding.  
* Paul Donaldson (CLF-ULTRA) and Gilberto Teobaldi (SCD) for scientific discussions and support.  
* Sunita Jones, Marion Samler and Elizabeth Bain for assistance with the licensing process.  
* Lesley Mansfield for project management support.  
