# SLR: 

![SLR](README_TOP.png)  
IBL image: [sIBL Archive](http://www.hdrlabs.com/sibl/archive.html)  

SLRはモンテカルロレイトレーシングに基づいたレンダラーです。  
SLR is a renderer based on Monte Carlo ray tracing.  

SLRは次の要素から構成されています。  
SLR consists of the following components.  
* libSLR - レンダリングコア / rendering core
* libSLRSceneGraph - シーン管理・読み込み機能 / scene managing & loading
* HostProgram

## 特徴 / Features
* Full Spectral Rendering (Monte Carlo Spectral Sampling)  
  (For RGB resources, RGB->Spectrum conversion is performed using Meng-Simon's method \[Meng2015\].)
* RGB Rendering
* BSDFs
    * Ideal Diffuse (Lambert) BRDF
    * Ideal Specular BRDF/BSDF
    * Oren-Nayer BRDF \[Oren1994\]
    * Improved Ward-Dür BRDF \[Moroder2010\]
    * Ashikhmin-Shirley BRDF \[Ashikhmin2000\]
    * GGX Microfacet BRDF/BSDF \[Walter2007\] with visible normal sampling \[Heitz2014\]
    * Mixed BSDF
    * ~~Layered BSDF~~ (TODO)
* Volume Rendering \[Pauly1999, Raab2008, Novák2014\]
    * Participating Media
        * Homogeneous Media
        * Heterogeneous Media (trilinear-interpolated grid)
        * Media enclosed by surfaces
    * Phase Functions
        * Isotropic
        * Henyey-Greenstein
        * Schlick
* Texture Types
    * Constant
    * Image
    * Procedural
        * Checkerboard (color & normal & float)
        * Cell Noise \[Worley1996\] like Voronoi diagram (color & normal & float)
        * Improved Perlin Noise \[Perlin2002\] (float)
        * Analytic Sky Texture \[Hosek2012\] (color)
* Bump Mapping (Normal Map)
* Light Types
    * Area (Polygonal) Light
    * Infinitely Distant Image Based Environmental Light
    * ~~Infinitely Distant Directional Light~~ (TODO)
* Light Distributions
    * Diffuse
    * Ideal Directional
* Camera Types
    * Perspective Camera with Depth of Field (thin-lens model)
    * Environment (Equirectangular) Camera
* Motion Blur
    * Camera Motion Blur
    * Object Motion Blur
    * ~~Deformation Blur~~ (TODO)
* Geometry Instancing
* Acceleration Structure Types
    * Standard BVH (median, mid-point, binned SAH)
    * Split BVH \[Stich2009\]
    * QBVH \[Dammertz2008\] constructed by collapsing SBVH
* Light Transport Algorithms
    * Path Tracing \[Kajiya1986\] with MIS (+ volumetric variant)
    * Bidirectional Path Tracing \[Veach1994, 1997\] (+ volumetric variant)
    * ~~Adaptive MCMC Progressive Photon Mapping~~ \[Hachisuka2011\]  
(has been dropped from current SLR implementation.)
* Correct handling of non-symmetric scattering due to shading normals \[Veach1996, 1997\]
* SLR Custom Language (C/Python-like syntax) for flexible scene description

## 動作環境 / Confirmed Environment
現状以下の環境で動作を確認しています。  
I've confirmed that the program runs correctly on the following environment.

* OS X 10.12.3 & Xcode 8.2.1 (Apple LLVM-Clang 8.0.0)
* Windows 8.1 & Visual Studio 2015 (MSVC 14.0, character-encoding conversion required)
* MacBook Pro Retina Late 2013

動作させるにあたっては以下のライブラリが必要です。  
It requires the following libraries.

* OpenEXR 2.2
* libpng 1.6
* assimp 3.2

## 注意 / Note
モデルデータやテクスチャーを読み込むシーンファイルがありますが、それらアセットはリポジトリには含まれていません。

There are some scene files loading model data and textures, but those assets are NOT included in this repository.

## 参考文献 / References
[Ashikhmin2000] "An Anisotropic Phong BRDF Model"  
[Dammertz2008] "Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of Incoherent Rays"  
[Hachisuka2011] "Robust Adaptive Photon Tracing Using Photon Path Visibility"  
[Heitz2014] "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"  
[Hosek2012] "An Analytic Model for Full Spectral Sky-Dome Radiance"  
[Kajiya1986] "THE RENDERING EQUATION"  
[Meng2015] "Physically Meaningful Rendering using Tristimulus Colours"  
[Moroder2010] "A New Ward BRDF Model with Bounded Albedo"  
[Novák2014] "Residual ratio tracking for estimating attenuation in participating media"  
[Oren1994] "Generalization of Lambert’s Reflectance Model"  
[Pauly1999] "Robust Monte Carlo Methods for Photorealistic Rendering of Volumetric Effects"  
[Perlin2002] "Improving Noise"  
[Raab2008] "Unbiased Global Illumination with Participating Media"  
[Stich2009] "Spatial Splits in Bounding Volume Hierarchies"  
[Veach1994] "Bidirectional Estimators for Light Transport"  
[Veach1996] "Non-symmetric Scattering in Light Transport Algorithms"  
[Veach1997] "ROBUST MONTE CARLO METHODS FOR LIGHT TRANSPORT SIMULATION"  
[Walter2007] "Microfacet Models for Refraction through Rough Surfaces"  
[Worley1996] "A Cellular Texture Basis Function"  

----
2017 [@Shocker_0x15](https://twitter.com/Shocker_0x15)
