1. Download DeepMD source code from github

2. Untar/unzip the downloaded file

3. Go to the source directory 
cd deepmd-kit-[version-number]/source

4. For Expanse, install python 3.9 through Conda
module load gpu/0.17.3b anaconda3/2021.05/kfluefz
conda create -n py39_env python=3.9
conda activate py39_env

5. Install necessary packages:
conda install -c conda-forge python-devel
pip install -U cmake
pip install --upgrade tensorflow-cpu

6. Create directory to install DeepMD:
mkdir ~/opt/DeepMD/ -p
deepmd_root=~/opt/DeepMD

7. Run cmake:
cmake -DENABLE_TENSORFLOW=TRUE -DUSE_TF_PYTHON_LIBS=TRUE -DPython_EXECUTABLE=$(which python) -DCMAKE_INSTALL_PREFIX=$deepmd_root ..


