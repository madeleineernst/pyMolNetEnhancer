# pyMolNetEnhancer

### Table of contents

* [Install](#install)
* [Map MS2LDA substructural information to mass spectral molecular networks (classical)](#Mass2Motifs_to_Network_Classical)
* [Map MS2LDA substructural information to mass spectral molecular networks (feature based)](#Mass2Motifs_to_Network_FeatureBased)
* [Map chemical class information to mass spectral molecular networks](#ChemicalClasses_to_Network)
* [Dependencies](#dependencies)
* [License](#licence)

## Install

Install with:

 `pip install pyMolNetEnhancer`
 
 
## Map MS2LDA substructural information to mass spectral molecular networks (classical) <a name="Mass2Motifs_to_Network_Classical"></a>
 
In order to map substructural information to a mass spectral molecular network you need to:
 
* [Create a molecular network](https://ccms-ucsd.github.io/GNPSDocumentation/quickstart/) through the Global Natural Products Social Molecular Networking (GNPS) platform
* Create an LDA experiment on [http://ms2lda.org/](http://ms2lda.org/) using the MGF clustered spectra downloaded from GNPS:

<img src="IMG/DownloadCulsteredMGF.png"/>

Then execute the code in [Example_notebooks/Mass2Motifs_2_Network_Classicaly.ipynb](https://github.com/madeleineernst/pyMolNetEnhancer/blob/master/Example_notebooks/Mass2Motifs_2_Network_Classical.ipynb) line by line.
The only things you need to specify are:

<ol>
  <li>Your GNPS job ID 
  <img src="IMG/GNPSJobID.png"/></li>
  <li>Your MS2LDA job ID
  <img src="IMG/MS2LDAJobID.png"/>
  <b>Note</b>: Depending on the size of this file, a server connection timeout may occur. Alternatively, you may download the file manually at http://ms2lda.org/: <br>
  <img src="IMG/Export_MS2LDA_Summary.jpg"/></li>
  <li>User-defined parameters for mapping the Mass2Motifs onto the network
  <img src="IMG/Mass2Motif_2_Network_Parameters.png"/>
  <tt>prob</tt>: minimal probability score for a Mass2Motif to be included. Default is 0.01. <br>
  <tt>overlap</tt>: minimal overlap score for a Mass2Motif to be included. Default is 0.3. <br>
  <b>Important</b>: The probability and overlap thresholds can be set within the ms2lda.org app as well under the Experimental Options tab. It is recommendable to do so when inspecting results in the web app. Importantly, the summary table contains filtered motif-document relations using the set thresholds in the web app. <br>
  <tt>top</tt>: This parameter specifies how many most shared motifs per molecular family (network component index) should be shown. Default is 5.
</li>
</ol>

## Map MS2LDA substructural information to mass spectral molecular networks (feature based) <a name="Mass2Motifs_to_Network_FeatureBased"></a>

In order to map substructural information to a mass spectral molecular network created through the feature based workflow you need to:

* [Create a feature based molecular network](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/) through the Global Natural Products Social Molecular Networking (GNPS) platform
* Create an LDA experiment on [http://ms2lda.org/](http://ms2lda.org/) using the MGF file created within MZmine (see [GNPS documentation](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/))

Then execute the code in [Example_notebooks/Mass2Motifs_2_Network_FeatureBased.ipynb](https://github.com/madeleineernst/pyMolNetEnhancer/blob/master/Example_notebooks/Mass2Motifs_2_Network_FeatureBased.ipynb) line by line.
The only things you need to specify are:

<ol>
  <li>Your GNPS job ID 
  <img src="IMG/GNPSJobID.png"/></li>
  <li>Your MS2LDA job ID
  <img src="IMG/MS2LDAJobID.png"/>
  <b>Note</b>: Depending on the size of this file, a server connection timeout may occur. Alternatively, you may download the file manually at http://ms2lda.org/: <br>
  <img src="IMG/Export_MS2LDA_Summary.jpg"/></li>
  <li>User-defined parameters for mapping the Mass2Motifs onto the network
  <img src="IMG/Mass2Motif_2_Network_Parameters.png"/>
  <tt>prob</tt>: minimal probability score for a Mass2Motif to be included. Default is 0.01. <br>
  <tt>overlap</tt>: minimal overlap score for a Mass2Motif to be included. Default is 0.3. <br>
  <b>Important</b>: The probability and overlap thresholds can be set within the ms2lda.org app as well under the Experimental Options tab. It is recommendable to do so when inspecting results in the web app. Importantly, the summary table contains filtered motif-document relations using the set thresholds in the web app. <br>
  <tt>top</tt>: This parameter specifies how many most shared motifs per molecular family (network component index) should be shown. Default is 5.
</li>
</ol>

## Map chemical class information to mass spectral molecular networks <a name="ChemicalClasses_to_Network"></a>

In order to map chemical class information to a mass spectral molecular network you need to:

* Create a molecular network using the [classical](https://ccms-ucsd.github.io/GNPSDocumentation/quickstart/) or [feature based](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/) workflow through the Global Natural Products Social Molecular Networking (GNPS) platform
* Perform <i>in silico</i> structure annotation using [Network Annotation Propagation](https://ccms-ucsd.github.io/GNPSDocumentation/nap/) (NAP), [DEREPLICATOR](https://ccms-ucsd.github.io/GNPSDocumentation/dereplicator/) or another tool of preference for <i>in silico</i> structure annotation 

Then execute the code in [Example_notebooks/ChemicalClasses_2_Network.ipynb](https://github.com/madeleineernst/pyMolNetEnhancer/blob/master/Example_notebooks/ChemicalClasses_2_Network.ipynb) line by line.
The only things you need to specify are:

<ol>
  <li>Your GNPS job ID 
  <img src="IMG/GNPSJobID.png"/></li>
  <li>Your DEREPLICATOR job ID(s)
  <img src="IMG/DereplicatorJobID.png"/></li>
  <li>Your NAP job ID(s)
  <img src="IMG/NAPJobID.png"/>
  </li>
</ol>

You can specify as many <i>in silico</i> annotation outputs as you wish. If you import results from applications different than NAP or DEREPLICATOR make sure that your input file is tab separated and includes a column named 'Scan' containing numeric identifiers matching the numeric node identifiers in the GNPS network and a column named 'SMILES' containing SMILES structures.
Make sure that you include all results as dataframe list items in the 'matches' object. The object 'gnpslib' contains all GNPS library hits:

 `matches = [gnpslib, nap, derep]`
 
In this notebook we use [ChemAxon's molconvert](https://docs.chemaxon.com/display/docs/Molconvert) to convert SMILES to InChIKeys. Make sure to have molconvert installed and add the path to the environment:
 
```
path = '/Applications/MarvinSuite/bin/'
os.environ['PATH'] += ':'+path
```
 



## Dependencies

python 3.6.5, collections 0.6.1, csv 1.0, functools, joblib 0.13.0, json 2.0.9, multiprocessing, networkx 2.1, operator, os, pandas 0.22.0, rdkit, re 2.2.1, requests 2.18.4, sys, time
 
## Licence
This repository is available under the following licence https://github.com/madeleineernst/pyMolNetEnhancer/blob/master/LICENSE
