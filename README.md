#SLR: 

![SLR](README_TOP.png)  
IBL image: [sIBL Archive](http://www.hdrlabs.com/sibl/archive.html)  

SLRはモンテカルロレイトレーシングに基づいたレンダラーです。  
SLR is a renderer based on Monte Carlo ray tracing.  

SLRは次の要素から構成されています。  
SLR consists of the following components.  
* libSLR - レンダリングコア / rendering core
* libSLRSceneGraph - シーン管理・読み込み機能 / scene managing & loading
* HostProgram

##特徴 / Features
* Full Spectral Rendering (Monte Carlo Spectral Sampling)  
  (For RGB resources, RGB->Spectrum conversion is performed using Meng-Simon's method \[Meng2015\].)
* RGB Rendering
* Various BSDF Types
    * Ideal Diffuse (Lambert) BRDF
    * Ideal Specular BRDF and BSDF
    * Oren-Nayer BRDF \[Oren1994\]
    * Improved Ward-Dür BRDF \[Moroder2010\]
    * Ashikhmin-Shirley BRDF \[Ashikhmin2000\]
    * GGX Microfacet BRDF and BSDF \[Walter2007\] with visible normal sampling \[Heitz2014\]
    * Mixed BSDF
    * ~~Layered BSDF~~ (TODO)
* Procedural Textures
    * Checkerboard (color & normal & float)
    * Cell Noise \[Worley1996\] like Voronoi diagram (color & normal & float)
* Bump Mapping (Normal Map)
* Correct handling of non-symmetric scattering due to shading normals \[Veach1996, 1997\]
* ~~Volume Rendering~~ (TODO)
* Image Based Environmental Light
* Camera
    * Perspective Camera with Depth of Field (thin-lens model)
    * Environment Camera
* Motion Blur
    * Camera Motion Blur
    * Object Motion Blur
    * ~~Deformation Blur~~ (TODO)
* Geometry Instancing
* Acceleration Structure
    * Binned SAH BVH
    * Split BVH \[Stich2009\]
    * QBVH \[Dammertz2008\] constructed by collapsing SBVH
* Light Transport Algorithms
    * Unidirectional Path Tracing \[Kajiya1986\] with MIS
    * Bidirectional Path Tracing \[Veach1994, 1997\]
    * ~~Adaptive MCMC Progressive Photon Mapping~~ \[Hachisuka2011\]  
(has been dropped from current SLR implementation.)
* SLR Custom Language (C/Python-like syntax) for flexible scene description

##動作環境 / Confirmed Environment
現状以下の環境で動作を確認しています。  
I've confirmed that the program runs correctly on the following environment.

* OS X 10.11
* Windows 8.1
* MacBook Pro Retina Late 2013

動作させるにあたっては以下のライブラリが必要です。  
It requires the following libraries.

* OpenEXR 2.2
* libpng 1.6
* assimp 3.2

##注意 / Note
モデルデータやテクスチャーを読み込むシーンファイルがありますが、それらアセットはリポジトリには含まれていません。

There are some scene files loading a model data and textures, but those assets are NOT included in this repository.

##参考文献 / References
[Ashikhmin2000] "An Anisotropic Phong BRDF Model"  
[Dammertz2008] "Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of Incoherent Rays"  
[Hachisuka2011] "Robust Adaptive Photon Tracing Using Photon Path Visibility"  
[Heitz2014] "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"    
[Kajiya1986] "THE RENDERING EQUATION"  
[Meng2015] "Physically Meaningful Rendering using Tristimulus Colours"  
[Moroder2010] "A New Ward BRDF Model with Bounded Albedo"  
[Oren1994] "Generalization of Lambert’s Reflectance Model"  
[Stich2009] "Spatial Splits in Bounding Volume Hierarchies"  
[Veach1994] "Bidirectional Estimators for Light Transport"  
[Veach1996] "Non-symmetric Scattering in Light Transport Algorithms"  
[Veach1997] "ROBUST MONTE CARLO METHODS FOR LIGHT TRANSPORT SIMULATION"  
[Walter2007] "Microfacet Models for Refraction through Rough Surfaces"  
[Worley1996] "A Cellular Texture Basis Function"  

----
2016 [@Shocker_0x15](https://twitter.com/Shocker_0x15)
