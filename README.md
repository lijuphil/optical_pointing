# List of sources used for COMAP optical pointing
The script sorts the sources from a star catalogue into specific elevation- and azimuth-wide bins, using the ranges specified by the user.\
\
To run the script, for e.g. elevation window of 30--40 deg, and 20-deg wide azimuth chunks:
``` 
python pointing_sources.py -sl star_catalog.cat -pl pathfinder -el 30 40 -az 20
```
The basic command line plotting option can be activated by setting -p to True as follows:
``` 
python pointing_sources.py -sl star_catalog.cat -pl pathfinder -el 30 40 -az 20 -p True
```
