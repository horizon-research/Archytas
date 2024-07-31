# Archytas

This repo contains the code of Archytas: A Framework for Synthesizing and Dynamically Optimizing Accelerators for Robotic Localization.

To run the demo, use python to run `localization_demo.py`.



## Build the project
```
git clone https://github.com/horizon-research/Archytas.git
cd Archytas
mkdir build
cd build
cmake ..
make
```

This demo contains a large-scale localization problem.

## Citation

If you think this work is useful in your research, please consider cite our paper:
```
@inproceedings{liu2021archytas,
  title={Archytas: A framework for synthesizing and dynamically optimizing accelerators for robotic localization},
  author={Liu, Weizhuang and Yu, Bo and Gan, Yiming and Liu, Qiang and Tang, Jie and Liu, Shaoshan and Zhu, Yuhao},
  booktitle={MICRO-54: 54th Annual IEEE/ACM International Symposium on Microarchitecture},
  pages={479--493},
  year={2021}
}
```