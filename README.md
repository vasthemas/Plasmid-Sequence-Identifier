# Plasmid Sequence Identifier

Plasmid Sequence Identifier parses .seq, .dna, .fa, .gbk files and identifies insert between given primer regions. 

- Translates insert into Amino Acid sequence 
- Converts insert into reverse compelement
- Stores all data into an excel file
    
# Installation 
### PyPi 

Install the following packages with pip 

xlwt: http://xlwt.readthedocs.io/en/latest/
```sh
$ pip install xlwt
```
Biopython:https://biopython.org/wiki/Documentation
```sh
$ pip install biopython:
```
snapgene_reader: https://github.com/IsaacLuo/SnapGeneFileReader
```sh
$ pip install snapgene_reader
```


# How to use 
Go to main/python/scripts 

```sh
$ python Plasmid_Sequence_Identifier.py
```
![](https://i.imgur.com/IdGI2Dy.png)

![](https://i.imgur.com/BoSEaav.png)


When prompted;Enter the path name of the folder you like to analyze. If the path is found, a list of all the files in the directory will be shown. 

You will then asked to enter the upstream and downstream primer regions. Make sure to enter a space between two primers, to indicate distinct seqeuences. If the regions are mixed up, the program will resolve the issue by finding out which region is really upstream. 

Finally, you will be prompted for a file name of the excel file. The file will be saved in the folder being analzyed. 

For each file, you will be told if the region was found and display the grabbed seqence. NaN is what is shown when nothing is found. 

All the data from the program will be stored in a excel file, including the translation, reverse compelement and its translation. 

![](https://i.imgur.com/sDx4Z0z.png)


You can use the test_sequence directory to test the script. 


If you are more comfortable with jupyter-notebook, one is provided with the exact same code in the main/python/notebooks directory. 



