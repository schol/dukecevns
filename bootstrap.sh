#!/bin/sh
cd ..
if [ ! -d COHERENTProposal2018 ]; then
  git clone https://code.ornl.gov/COHERENT/COHERENTProposal2018.git
fi
cd COHERENTProposal2018
git pull
cd assumptions
rsync -rtphv eff gs jsonfiles qf *.cc *.root ../../dukecevns

cd ../../dukecevns
# fetch json submodule if json/ is empty
if [ ! "$(ls -A json)" ]; then git submodule update --init; fi

make sns_rates
