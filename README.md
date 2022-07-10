# InSAR_atmo_corrections
This project is made to correct the atmospheric noise in InSAR interferogramms

## requirements : 
> **OS**: Macintosh 12.4
> 
> **Python version** : 3.9
> 
> **Packages** : see `requirements.txt`

## Installation :

`pip install -r requirements.txt`

## Organisation :
![project flowchart](/Volumes/Pierre_2TO/Stage_M1_Albino/Rapport/Images/methodology/STAGE_M1_methodologie_V4_full.png)

### Input data
To run this programm, the required input data are : 
> 1. [LiCSAR](https://comet.nerc.ac.uk/comet-lics-portal/) database: 
>    - Interferogramm(s) (date1_date2.geo.unw.tiff)
>    - DEM (XXXX.geo.hgt.tiff)
>    - Incidence angle array (XXXX.geo.U.tiff)
> 
> 2. [GACOS](http://www.gacos.net) database :
>    - Gacos atmospheric model at date1 
>    - Gacos atmospheric model at date2 

### Functiuns 
To correct the input interferogramm from atmospheric noise with GACOS data, run the function call
`Gacos_correction_tiff.py` wich is the main programm who call different other functiuns.
You will need to modify the different paths.
