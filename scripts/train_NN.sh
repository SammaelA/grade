export LD_LIBRARY_PATH=/usr/local/cuda-11.0/lib64 ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export PATH=/usr/local/cuda-11.0/bin${PATH:+:${PATH}}
python3 neural_selection.py