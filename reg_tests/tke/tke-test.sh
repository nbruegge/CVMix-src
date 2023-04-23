#!/bin/bash

path_output='/Users/nbruegge/work/src/CVMix-src/reg_tests/tke/out/'

#cat > ${path_output}/input.nl << EOF
cat > input.nl << EOF
&cvmix_nml
mix_type = 'tke'
nlev     = 50
max_nlev = 50
/
! TKE parameteris
&tke_nml
dtime = 60.
nt = 1000
nt_output = 20
path_out = "${path_output}"
/
EOF

# (1) Load required routines
. ../common/environ.sh
. ../common/usage.sh
. ../common/parse_inputs.sh
. ../common/build.sh
. ../common/run.sh

rm -rf path_output
mkdir path_output

parse_inputs $@
build
run

## (2) Look at output
#if [ "${USE_NETCDF}" == "netcdf" ]; then
#  ncdump -v Tdiff data_LMD.nc
#  ncdump -v Mdiff data_PP1d.nc
#else
#  cat data_LMD.out
#  cat data_PP1d.out
#fi
