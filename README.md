# Bioinformatic Learning Lab
Script for BILL project analysis


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#Acknowledgement">Acknowledgement</a></li>
        <li><a href="#prerequisites">Prerequisites</a></li>       
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- GETTING STARTED -->
## Getting Started


### Acknowledgement

This script was concieved in the first place to only sort and filter VCF files in order to understand the informations they store for the needs of the BiLL Project. The script was later on upgraded to answer further questions such as mutation selection and pertinence, and also to create simple visualisations such as mutation depth/coverage or mutation type distribution.

### Prerequisites

In order to run optimally our SamReader.py, this are configuration and packages our script uses:

* Python 3
  ```sh
  # If your Python version is under 3.0, the script could not work, please install Python 3.0:
  
  sudo apt install python3
  ```
  
* Packages used: re, os, sys (basic Python packages) pandas, numpy (dataframes) ,metaplotlib.pyplot (graphs)


## Usage

### Input

Two valid VCF Files with the '.VCF' extension.

### Launch the script

```sh
./VCFReader.py  file1.VCF file2.VCF 


```

### Output

* Two summary files of the two VCF files in input, sorted and filtered.
* Graphic visualisation as well as a comparative dataframe of the input files :

The script consists of two functions, the first one allows the user to extract every mutation that is located in both input files. Second function takes as an input the position of a mutation and outputs all the necessary informations ( chromosome, name, tag, localisation, sequence, allelic frequency, depth...) if this mutation exists at least once in any of the two input files.




<!-- LICENSE -->
## License

This program is free software: you can **redistribute** it and/or **modify** it under the terms of the **GNU General Public License** as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

See the GNU General Public License for more details ([GNU website](https://www.gnu.org/licenses/)).



<!-- CONTACT -->
## Contact

MECHKOURI Nazim - [sabri-nazim.mechkouri@etu.umontpellier.fr](mailto:sabri-nazim.mechkouri@etu.umontpellier.fr) 

<p align="right">(<a href="#top">back to top</a>)</p>
