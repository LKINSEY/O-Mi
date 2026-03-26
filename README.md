# O-Mi
Scripts for Working With Various Open Microscope Datasets

# Tutorial
- requires a windows machine and anaconda (or miniconda3)
- working directory should be a location with enough storage capacity to hold data you want to work with (example URLs will take up several Gb of data)
In conda bash (change C:/temp to be dir of choice):
```
conda create -n openMicroEnv python=3.9 conda-forge::zeroc-ice==3.6.5 omero-py
conda activate openMicroEnv
cd C:/path/to/PlaceRepo
git clone https://github.com/LKINSEY/O-Mi.git
cd C:/path/to/PlaceRepo/O-Mi
pip install -r requirements.txt
```

# Run function to acquire imaging data
- imageURLs.txt is contained in repo with example URLs contained
- can also make imageURLs.txt a remote location with better capacities (to be tested in future)

simple 1 line command to store data at chosen destination
```
python OMiReqs.py dir/containing/imageURLs.txt
```

# Note:
- line 313 should be edited to args = [(url, urlDir, (1000,1000,1000)) for url in urls] for now***
