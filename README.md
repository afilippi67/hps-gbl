hps-gbl
=======

GBL implementation for HPS


Setup user name: git config --global user.name "name" git config --global user.email "name@email.com"

Clone the remote repository: (Recommended) git clone https://github.com/perhansson/hps-gbl.git 

Checkout the GBL repository: svn checkout http://svnsrv.desy.de/public/GeneralBrokenLines/tags/V01-15-02 GeneralBrokenLines

Other requirements:
- numpy, pyROOT


Run GBL over 100 events:

python gbltst-hps.py gbloutputfile.gbl -n 100



Basic help for git usage:

Fetch latest from remote (local and remote is already configured with clone): git pull

Check status of local repository: git status

Add files to be committed to the local repository: git add file(s)

Commit files to be committed to the local repository (add -a to commit all staged files): git commit -m "Message" file(s)

Push your local changes to the remote (local and remote is configured by clone): git push

