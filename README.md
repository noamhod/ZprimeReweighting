# ZprimeReweighting

Put all flies in the pythia8xxx/examples directory
Prequisits:
- pythia8215 or later
- ROOT v5 (not yet 6), copmpatile with the available Python version

Put all files found here in the pythia82xx/examples directory and have fun :)

To run,
- make main02
- edit the number of events to generate in main02.cc
- execute `./main02 DY`
- execute `./main02 ZP`
- etc
- execute `python makeMatrix.py` to make all libraries for this configuration
- execute `python run2.py` to reweight the SM sample into this configuration and store the weights in the tree
- execute `python template.py` to plot the tempaltes
- execute `python run3.py` to validate a single template
