# Tesi_GPU
CUDA implementation of Tesi (https://github.com/FRSYH/Tesi)

Compiled with CUDA 9.2 nvcc for NVIDIA GT 1050 compute capability 6.1
> nvcc gm.cu -o gm -w -Xcompiler " -openmp" -gencode arch=compute_61,code=sm_61 -lcudadevrt -rdc=true

Execute
> ./gm --input input_file --output output_file

Performance delle varie tipologie dei kernel su input contenuto nel file input.txt
1 solo thread per tutto il kernel -> 2.1 sec
1 thread per riga da ridurre -> 1.1 sec
1 thread per cella della matrice da ridurre -> TODO
