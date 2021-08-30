# NLPAR
NLPAR (Non-Local Pattern Averaging and Reindexing) is a method for enhancing noisy EBSD patterns for increased indexing success. This repository does not include an indexing method, only the pattern enhancement method.  The advantage of NLPAR over other methods is that that a large number of patterns can be averaged together without significnatly degrading the original signal.  

## Installation 

### MacOS:

The library by default is linked against OpenMP, which is not included with MacOS. 

#### Option 1: Install OpenMP with Homebrew
Go to https://brew.sh and follow the instructions to install Homebrew.  Then:

    $> brew install libomp

Then go into the "makefile" and change the lines that look like:

    #OPENMP=-Xpreprocessor -fopenmp -lomp
    #CCOMP= clang 
    #CPPCOMP= clang++
    OPENMP = -fopenmp
    CCOMP = /usr/local/bin/gcc
    CPPCOMP = /usr/local/bin/g++

To look like this: 
    
    OPENMP=-Xpreprocessor -fopenmp -lomp
    CCOMP= clang 
    CPPCOMP= clang++
    #OPENMP = -fopenmp
    #CCOMP = /usr/local/bin/gcc
    #CPPCOMP = /usr/local/bin/g++

Then:

    $>make clean
    $>make -f makefile

#### Option 2: Install GCC with OpenMP 
The version I have found for this is at http://hpc.sourceforge.net/
Follow the directions for installing the latest gcc complier.  Once gcc with OpenMP is installed you should be able to (from the NLPAR Source Directory).

Go into the "makefile" and change the lines that look like:

    OPENMP=-Xpreprocessor -fopenmp -lomp
    CCOMP= clang 
    CPPCOMP= clang++
    #OPENMP = -fopenmp
    #CCOMP = /usr/local/bin/gcc
    #CPPCOMP = /usr/local/bin/g++

To look like this: 
    
    #OPENMP=-Xpreprocessor -fopenmp -lomp
    #CCOMP= clang 
    #CPPCOMP= clang++
    OPENMP = -fopenmp
    CCOMP = /usr/local/bin/gcc
    CPPCOMP = /usr/local/bin/g++

    $>make clean
    $>make -f makefile

#### Option 3: Run as a single threaded process  
Open the make file in an editor and ...
Comment out the line:

    #OPENMP=-Xpreprocessor -fopenmp -lomp
    CCOMP= clang 
    CPPCOMP= clang++
    #OPENMP = -fopenmp
    #CCOMP = /usr/local/bin/gcc
    #CPPCOMP = /usr/local/bin/g++

Then save changes and at the command prompt:

    $>make clean
    $>make -f makefile





## Running NLPAR

The assumption is that you have a running version of IDL.  I have tested this with IDL 8.5 - 8.7.  I do not see any reason why it would not work with earlier versions, nor later, in this code I do not use any advanced IDL data structures (HASH, LIST, DICTIONARY).  

_Not fully supported:_ If you do not want to purchase a IDL license, the open source GDL has worked for me with these purposes: https://github.com/gnudatalanguage/gdl To install GDL on mac systems, you can use Homebrew https://brew.sh 
Once Homebrew is installed at a command prompt:

        $>brew tap brewsci/science
        $>brew install gnudatalanguage

Then from the command prompt:

        $>gdl
And follow the instuctions below for IDL. 


File compatability: The codes expects an UP1/2 version 3 file as an input. Version 3 files contain the scan grid size, which is required for processing.  If the up1/2 file is version 1, then the user needs to input the number of columns and rows using the nColScan and nRowScan keywords. 

The scan grid must be a square grid.  NLPAR does not currently support Hex Grids (a tool should be released soon that converts a up1/2 file from hex to square for processing, then back to a hex).

From the IDL prompt get an estimate of lambda:

    IDL> nlpar_lambda_opt, './scan8.up2', /up2

    Approximate RAM usage: 1.53687GB
    Block Number:  0

    ****************
    Guess at lambda
    Low:      0.645804  High:       1.17495  Mean:      0.885553

Then use this in a call to NLPAR:

    IDL> nlparv5, 'scan8.up2' , 'scan8nlpar.up2', lambda=1.17, searchWindow=4, /up2

    Approximate RAM usage: 1.71886GB
    Block Number:  0


This will generate a new up2 file with the processed patterns.  To reindex the patterns, copy the associated .osc or .ang file, change the name so that the file root is identical to the nlpar processed pattern file and reindex using the EDAX software.  

Further details of nlpar program options can be found in the 'nlparv5.pro' source.

