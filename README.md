#SLR: 

SLRはモンテカルロレイトレーシングに基づいたレンダラーです。

SLR is a renderer based on Monte Carlo ray tracing.

##特徴 / Features
* Unidirectional Path Tracing (with MIS)
* Adaptive MCMC Progressive Photon Mapping
* Multiple BSDF Layers
* Image Based Environmental Light
* Bump Mapping (Normal Map)
* Depth of Field
* Camera / Object Motion Blur
* Geometry Instancing
* Binned SAH BVH

##動作環境 / Confirmed Environment
現状以下の環境で動作を確認しています。  
I've confirmed that the program runs correctly on the following environment.

* OS X 10.9.5
* MacBook Pro Retina Late 2013

動作させるにあたっては以下のライブラリが必要です。  
It requires the following libraries.

* OpenEXR 2.10
* libpng 1.6.14
* assimp (currently obtained from GitHub for supporting the .assbin format)

##注意 / Note
...

----
2015 [@Shocker_0x15](https://twitter.com/Shocker_0x15)
