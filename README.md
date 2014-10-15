# Tools for Resonant HbbHgg analysis


## How to install
  * First setup a CMSSW 611 release
  * Clone the project, setup the environment
    
        git clone git@github.com:ResonantHbbHgg/Selection.git
        cd Selection
        source setup.sh
    
  * Compile everything

        make

## How to make the plot-trees
  * Edit do_selection.sh with the necessary information, then:
    
        sh do_selection.sh

## How to make the limit-trees
  * After making the plottrees, edit do_quicktrees.sh with the necessary information, then:
    
        sh do_quicktrees.sh

