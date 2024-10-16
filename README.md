# AdaHRBF: Gradient-Adaptive Hermite-Birkhoff Radial Basis Function Interpolants for Three-dimensional Stratigraphic Implicit Modeling
Three-dimensional Stratigraphic Implicit Modeling Using Gradient-Adaptive Hermite-Birkhoff Radial Basis Function Interpolants

This is a [Matlab](https://ww2.mathworks.cn/en/products/matlab.html) version of implicit modeling method using an Hermite-Birkhoff radial basis function (HRBF) formulation with an adaptive gradient magnitude (AdaHRBF) for continuous three-dimesnional stratigrahic potential field (3D SPF) modeling of multiple stratigraphic interfaces.  

As described in **AdaHRBF: Gradient-Adaptive Hermite-Birkhoff Radial Basis Function Interpolants for Three-dimensional Stratigraphic Implicit Modeling** by Baoyi Zhang<sup>1</sup>, 
Linze Du<sup>1</sup>,
Umair Khan<sup>2</sup>, 
Yongqiang Tong<sup>1,3</sup>,
Lifang Wang<sup>3,4</sup> and
Hao Deng<sup>1*</sup>.   
<sup>1</sup>Key laboratory of Metallogenic Prediction of Nonferrous Metals & Geological Environment Monitoring (Ministry of Education) / School of Geosciences & Info-Physics, Central South University; <sup>2</sup>Institute of Deep-Sea Science and Engineering, Chinese Academy of Sciences; <sup>3</sup>Wuhan ZGIS Science & Technology Co. Ltd.; <sup>4</sup>School of Geomatics and Geography, Hunan Vocational College of Engineering

## Requirments

```
matlab2021a
```

## Dataset

**According to the comprehensive stratigraphic column, the burial depth of each stratigraphic interface relative to the top surface of the Quaternary was used as the attribute value of the SPF.We extracted the geological boundary points with 3D coordinates from 2D cross-sections. The attribute points and attitude points of each stratigraphic interface and fault plane extracted from the geological map and cross-sections were used as the original dataset for 3D SPF modeling. .** 

**The proposed stratigraphic potential field modeling by AdaHRBF interpolants can provide a suitable basic model for subsequent geosciences numerical simulation.**

The attribute points and attitude points dataset is available in Data file.

## Modeling

Run Experiment_AdaHRBF.m to start building a model with multiple stratigraphic interfaces.

## License

This extension  is released under a creative commons license which allows for personal and research use only. 
For a commercial license please contact the authors.
