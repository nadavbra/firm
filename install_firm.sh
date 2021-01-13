#! /bin/tcsh -e


# Determine the Python version to use.

# If the environment variable PYTHON_CMD exists, then just use it.
if ($?PYTHON_CMD) then
    set python_cmd = $PYTHON_CMD
    goto determine_final_python_version
endif

ask_default_python_version:
printf "Do you want to install geneffect in your default Python version (`python --version |& cat`)? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set python_cmd = "python"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative version of Python: "
    set python_cmd = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_python_version
endif

determine_final_python_version:
echo "Will install geneffect for the Python version `$python_cmd --version |& cat` ($python_cmd)"


# Determine a temp working directory.

# If a the environment variable PYTHON_CMD exists, then just use it.
if ($?TEMP_DIR) then
    set temp_dir = $TEMP_DIR
    goto determine_final_temp_dir
endif

ask_default_temp_dir:
printf "Do you want to use /tmp as a temporary working directory? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set temp_dir = "/tmp"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set temp_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_temp_dir
endif

determine_final_temp_dir:
set temp_dir = `echo $temp_dir | sed 's:/*$::'`
echo "Will use $temp_dir as a temporary working directory."


# Determine the data directory.

# If a the environment variable DATA_DIR exists, then just use it.
if ($?DATA_DIR) then
    set data_dir = $DATA_DIR
    goto determine_final_data_dir
endif

ask_data_dir:
printf "Do you want to use ~/data as the directory for all the data files required by geneffect? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set data_dir = "~/data"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set data_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_data_dir
endif

determine_final_data_dir:
set data_dir = `echo $data_dir | sed 's:/*$::'`
echo "Will use $data_dir as the directory for data files."


# Install general dependencies.

echo "Installing dependencies..."
$python_cmd -m pip install numpy pandas biopython sklearn cython interval_tree 


# Install geneffect.

echo "Installing geneffect..."
# TODO replace with wget
wget https://raw.githubusercontent.com/nadavbra/geneffect/master/install_geneffect.sh -O ${temp_dir}/install_geneffect.sh
chmod a+x ${temp_dir}/install_geneffect.sh
setenv PYTHON_CMD $python_cmd
setenv TEMP_DIR $temp_dir
setenv DATA_DIR $data_dir
${temp_dir}/install_geneffect.sh
rm -f ${temp_dir}/install_geneffect.sh


# Install firm.

mkdir -p ${temp_dir}/firm
git clone https://github.com/nadavbra/firm.git ${temp_dir}/firm
cd ${temp_dir}/firm
$python_cmd ./setup.py install
cd -
rm -fr ${temp_dir}/firm


# Set the Pfam HMM profiles.

echo "Installing hmmer3..."
mkdir -p ${temp_dir}/hmmer
cd ${temp_dir}/hmmer
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd -
cd ${temp_dir}/hmmer/hmmer-3.1b2-linux-intel-x86_64
./configure
make
sudo make install
rehash
cd -
rm -fr ${temp_dir}/hmmer

echo "Preparing Pfam HMM profiles..."
mkdir -p ~/data/pfam/hmm/profiles
cd ~/data/pfam/hmm
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.dat.gz
gzip -d Pfam-A.hmm.dat.gz
hmmpress Pfam-A.hmm
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/active_site.dat.gz
gzip -d active_site.dat.gz

echo "Succefully installed firm."