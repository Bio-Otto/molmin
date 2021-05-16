# molmin ![molmim v0.1 ff69b4](https://img.shields.io/badge/<molmim>-<v0.1>-<ff69b4>)
 Macro-Molecule mutater and minimizer

[![Powered by |Ozbek' Lab](https://github.com/Bio-Otto/Example_MD_Scripts/blob/master/PoweredByOzbekLab.png)](http://compbio.bioe.eng.marmara.edu.tr/)


A Package for mutating and minimizing

For terminal usage just use; 

```sh
$ python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4
```

* __-p__  -->  Give the absolute path of your pdb file. 
* __-pff__  -->  Protein Forcefield (The program defaultly will use ```amber96``` forcefield)
* __-wff__  -->  Water Forcefield (The program defaultly will use ```tip3p``` forcefield)
* __-mut__  -->  The program will apply the mutation described as "target_residue-residue_number-desired_residue" (For example: ```-mut ASN-120-ASP```).
* __-mut_ch__  -->  The program will mutate selected residue according to indicated chain ID.



## Then run it.

### Dependencies

molmin uses a number of open source projects to work properly:

* __OpenMM__ - A high performance toolkit for molecular simulation. 
* __pdbfixer__ 

And of course molmim v0.1 itself is open source with a [public repository][molmim] on GitHub.

### Also you can check full functional parameters with typing 

```sh
$ python minimizer.py -h
```

For production just type...

```sh
$ python minimizer.py -p <pdb file> -pff <amber96> -wff <tip3p> -mut <ASP-121-ASN> -mut_ch <A>
```

The Program applying minimize and mutation according to yur choises using powerfull OpenMM Molecular Dynamic Toolkit, which also supports the Cuda platform. 

# New Features!

  - Mutation
  - Fixing non-standart residues
  - Minimizing with heterogens or ions

### Installation [Not Supporting Yet!]

Open your favorite Terminal and run these commands.

```sh
$ cd molmim
$ python setup.py --install
```


License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


[MDPERTOOL]: <https://github.com/Bio-Otto/molmim>
