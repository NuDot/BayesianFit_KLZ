#!/bin/bash -l

source /project/snoplus/datacleaning/snoing/install/env_rat-dev.sh
# source /project/snoplus/datacleaning/rat-6.16.8/env.sh

#$$ -P snoplus

#$$ -l h_rt=${TIME}

#$$ -j y

#$$ -pe omp 1

SimulationSpectrumDir=/projectnb/snoplus/KLZanalysis/SimulationSpectrum-tmp-20190815
SimulationSpectrumFile_XeLS=$SimulationSpectrumDir/PhotonSimulationSpectrum_run015404_0.073_13.1_1.107_1.00_XeLS.root
SimulationSpectrumFile_KamLS=$SimulationSpectrumDir/PhotonSimulationSpectrum_run015404_0.073_13.1_1.107_1.00_KamLS.root
SimulationSpectrumFile_film=$SimulationSpectrumDir/PhotonSimulationSpectrum_run015404_0.073_13.1_1.107_1.00_film.root

DataDir=/projectnb/snoplus/KLZ_TAUP/data
EnergySpectrumDir=/projectnb/snoplus/KLZanalysis/EnergySpectrum-temporary-190215
LiveTimeDir=/projectnb/snoplus/KLZ_TAUP/livetime
BiPoDir=/projectnb/snoplus/KLZanalysis/BiPo_zen800
RunListDir=/projectnb/snoplus/KLZ_TAUP/data
ParameterDir=/projectnb/snoplus/KLZ_TAUP/vector-file-XeLS-tmp-v3/target_all_free_v3
Eth=0.5
UppEth=4.8

export DOUBLEBETA_ANALYSIS_DATA_DIR=$DataDir
export DOUBLEBETA_ANALYSIS_ENERGY_SPECTRUM_DIR=$EnergySpectrumDir
export DOUBLEBETA_ANALYSIS_LIVETIME_DIR=$LiveTimeDir
export DOUBLEBETA_ANALYSIS_BIPO_DIR=$BiPoDir
export DOUBLEBETA_ANALYSIS_KLG4SIM_SPECTRUM_DIR=$Klg4simSpectrumDir
export DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_XELS=$SimulationSpectrumFile_XeLS
export DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_KAMLS=$SimulationSpectrumFile_KamLS
export DOUBLEBETA_ANALYSIS_SIMULATION_SPECTRUM_FILE_FILM=$SimulationSpectrumFile_film
export DOUBLEBETA_ANALYSIS_RUNLIST_DIR=$RunListDir
export DOUBLEBETA_ANALYSIS_PARAMETER_DIR=$ParameterDir
export DOUBLEBETA_ANALYSIS_ETH=$Eth
export DOUBLEBETA_ANALYSIS_UPP_ETH=$UppEth

source /project/snoplus/KLZ/bat_new/bat_install/bat/bat_parallel.sh
make clean
make
./plot_parameter
