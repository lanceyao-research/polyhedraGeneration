# polyhedraGeneration
MATLAB codes for generating projection images of polyhedra used in article: Lehan Yao, Zhiheng Lyu, Jiahui Li, Qian Chen, "No Ground Truth Needed: Unsupervised Sinogram Inpainting for Nanoparticle Electron Tomography (UsiNet) to Correct Missing Wedges" npj Computation Materials 2024 10, 28  
<br/><br/>
The codes requires the MATLAB version of astra-toolbox, which can be downloaded and compiled from https://github.com/astra-toolbox/astra-toolbox. The version I used is 1.9.0. Newer versions might work but are not tested.  
Seven script files used for generating different polyhedra are included: cube, coreshell cube, octahedron, triangular prism, tetrahedron, tetrahexahedron, and random polydron. Those shapes account for the most common nanoparticles encountered in synthesis and self-assembly. The generation codes considers randomness indlucing random orientation, size distribution and vertex roundess.  
<br/><br/>
The output includes projection images and their corresponding volumetric images. The total number of generated image pairs could be adjusted in the script. The output folders should first be created before running the codes.  
