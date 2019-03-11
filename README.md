# pyMolNetEnhancer

### Table of contents

* [Install](#install)
* [Map MS2LDA substructural information to mass spectral molecular networks (classical)](#Mass2Motifs_to_Network_Classical)
* [Help and Troubleshooting](#help-and-troubleshooting)
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
  <li>The MS2LDA job ID
  <img src="IMG/MS2LDAJobID.png"/></li>
  <li>User-defined parameters for mapping the Mass2Motifs onto the network
  <img src="IMG/Mass2Motif_2_Network_Parameters.png"/>
  where <br>
  <tt>prob</tt>: minimal probability score for a Mass2Motif to be included. Default is 0.01. <br>
  <tt>overlap</tt>: minimal overlap score for a Mass2Motif to be included. Default is 0.3. <br>
   <b>Important</b>: The probability and overlap thresholds can be set within the ms2lda.org app as well under the Experimental Options tab. It is recommendable to do so when inspecting results in the web app. Importantly, the summary table contains filtered motif-document relations using the set thresholds in the web app. <br>
  <tt>top</tt>: This parameter specifies how many most shared motifs per molecular family (network component index) should be shown. Default is 5.
</li>
</ol>
 
## Licence
This repository is available under the following licence https://github.com/madeleineernst/pyMolNetEnhancer/blob/master/LICENSE
