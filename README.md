# SCANNET: A software tool for identifying communities in protein networks

## Project Description

 O SCANNET is a program that integrates all algorithms of the [Protein Similarity Networks (PSN)] method into a unique tool to perform all the steps needed to identify communities in protein networks as a single workflow.

## Technical information 
 
 SCANNET was developed in C language using the NetBeans IDE 7.4 environment.
 
## Installation

##### Generating the Executable File
To compile  and generate the executable file you must have a C compiler installed on your computer. 
* Windows OS  
We suggest installing the compiler MinGW (Minimalist GNU for Windows) , available at: http://www.mingw.org/

* Linux OS  
GCC (GNU Compiler Collection) for Linux or Mac OS, available at: https://gcc.gnu.org/.  

After installation compiler, is required run the command line below to compile the scannet.c and generate the executable file.

			gcc -o <EXECUTABLE NAME> scannet.c
	
<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig1.png" /> </p>
    
<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig2.png" /> </p>
			
For enjoy the full potential of SCANNET, we suggest you install on your computer the GNUPLOT software, available at: http://www.gnuplot.info/.	
			
## How Use 

##### Running Program	
For start the program you need only double-click on the generated executable file.

##### Setting Input File  
    
SCANNET use as input, a file (.txt or .dat) containing the similarity matrix between the proteins that will be analyzed. So you must provide the name of the input file. The used input file must be in the same SCANNET run directory.

<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig3.png" /> </p>

##### Setting Protein Order File
During the communities identification process, the SCANNET reassembles the proteins of the similarity matrix, changing the initial order these proteins. To facilitate the proteins identification in the generated output (Dendrogram, Complex Network and Color Matrix), the SCANNET generates an output file with the final proteins order. If the user wanted has this feature, it's necessary set the file (.csv or .dat or .txt) with the initial protein order. The protein order file is optional, but if used must be in the same SCANNET run directory.

<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig4.png" /> </p>

##### Select Similarity Value
The PSN method assumes that a set of 101 networks is generated from the dataset, each one for a specific threshold of similarity S in the range [0,100]. Afterwards, it calculates the distance between all pairs of subsequent networks in the network set, in order to select the similarity value that best represents the communities in the dataset. Such network is generally the one with the largest distance from its immediately subsequent network. In some cases, more than one similarity value could be selected. To facilitate the work of selecting the best similarity value, SCANNET displays a chart with all distances between pairs of neighboring networks and prompts the user to select the similarity value on which he will continue the analytical process.

<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig5.png" /> </p>
	
The chart is characterized by a series of peaks, where each peak represents a point at which the obtained networks has greater chances of suffer important alterations topological, or to be divided in new communities. After the analysis of the graph the user can select the similarity value to be analyzed.		

<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig6.png" /> </p>

> **Note:** For the SCANNET displays the chart with all distances between pairs of neighboring networks is necessary you have installed on your computer GNUPLOT software, available at: http://www.gnuplot.info/.
 
##### End of Execution
At the end of the run, the SCANNET generates the critical network, the dendrogram, the color representation of the neighborhood matrix and the ordered proteins, each one in a different output file. It also generate a log file with all the information of program execution. The output files can be read by several graphical analysis softwares such as Origin Lab®, Matlab®, GNU Plot, etc.

<p align="center"> <img src="https://dl.dropboxusercontent.com/u/20113863/Fig7.png" /> </p>

## Next Steps
	
Develop a graphical user interface more intuitive and powerful that present graphically on the screen the complex network, color matrix and the dendrogram generated by SCANNET.


[Protein Similarity Networks (PSN)]:http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001131
