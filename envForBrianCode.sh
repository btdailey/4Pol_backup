export ROOTSYS=/usr/local/root-5.34.18/bin/root
. thisroot.sh
export COAST_DIR=/rh5stuff/32bit/src/coast-v3r2-test/install
export COAST_USER_LIB=/rh5stuff/32bit/src/coast-v3r2-test/THRadio
source /data/anita/btdailey/analysis_info/anita_new_root/eventReaderRoot/eventReaderSetup.sh
export PATH=${ROOTSYS}/bin:/rh5stuff/64bit/bin:/usr/local/bin:/bin:/usr/bin:$PATH:
export LD_LIBRARY_PATH=/home/dailey.110/analysis:/data/anita/btdailey/Healpix_stuff/cfitsio:${ROOTSYS}/lib:/home/dailey.110/Healpix_3.30/src/cxx/Healpix_cxx:/home/dailey.110/Healpix_3.30/lib:$LD_LIBRARY_PATH

