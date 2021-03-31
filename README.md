# Pre-trained classifiers
Pre-trained classifiers can be found here:
https://drive.google.com/drive/folders/14llv-YqOWYKUykXUPozA0BdNXkBCxxcv?usp=sharing

# Supplementary data
Supplementary data containing DNM calls can be found here:
https://docs.google.com/spreadsheets/d/1OSAlX-XnBzWZQrpzt_tfHd58Jhz6FqO05BBgrJUx4gI/edit?usp=sharing

# Usage
```
usage: 
 _______  __   __  __    _  _______  __   __  ______   __    _  __   __ 
|       ||  | |  ||  |  | ||       ||  | |  ||      | |  |  | ||  |_|  |
|  _____||  |_|  ||   |_| ||_     _||  |_|  ||  _    ||   |_| ||       |
| |_____ |       ||       |  |   |  |       || | |   ||       ||       |
|_____  ||_     _||  _    |  |   |  |       || |_|   ||  _    ||       |
 _____| |  |   |  | | |   |  |   |  |   _   ||       || | |   || ||_|| |
|_______|  |___|  |_|  |__|  |___|  |__| |__||______| |_|  |__||_|   |_|

Version 0.1.0.1    Authors: Danny Antaki, Aojie Lian, James Guevara    
                   Contact: j3guevar@ucsd.health.edu
---------------------------------------------------------------------------------
    synthdnm  -f <in.fam>  -v  <in.vcf.gz>  [-oLgkVh]
    
necessary arguments:
  
  -v, -vcf    PATH    VCF file
  -f, -fam    PATH    PLINK pedigree (.fam/.ped) file
  
optional arguments:

  -g  --gen    STR     human reference genome version [default: hg38]
  
  -h, --help           show this message and exit
     
synthdnm.py: error: the following arguments are required: -v/--vcf, -f/--fam

```
