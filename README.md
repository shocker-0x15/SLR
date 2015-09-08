#SLR: 

SLRはモンテカルロレイトレーシングに基づいたレンダラーです。

SLR is a renderer based on Monte Carlo ray tracing.

##特徴 / Features
* ~~Full Spectral Rendering~~ (under Dev.)  
  (For RGB resources, RGB->Spectrum conversion is performed using Meng-Simon's method [1].)
* Various BSDF Types (including Mixed and ~~Layered~~ (under Dev.) BSDF)
* Image Based Environmental Light
* Bump Mapping (Normal Map)
* Depth of Field
* Camera / Object Motion Blur
* Geometry Instancing
* Binned SAH BVH
* Unidirectional Path Tracing \[2\] with MIS
* Adaptive MCMC Progressive Photon Mapping \[3\]  
  (current implementation is buggy and does not support IBL.)

[1] "Physically Meaningful Rendering using Tristimulus Colours"  
[2] "THE RENDERING EQUATION"  
[3] "Robust Adaptive Photon Tracing Using Photon Path Visibility"

##動作環境 / Confirmed Environment
現状以下の環境で動作を確認しています。  
I've confirmed that the program runs correctly on the following environment.

* OS X 10.9.5
* MacBook Pro Retina Late 2013

動作させるにあたっては以下のライブラリが必要です。  
It requires the following libraries.

* OpenEXR 2.2
* libpng 1.6
* assimp (currently obtained from GitHub for supporting the .assbin format)

##注意 / Note
モデルデータやテクスチャーを読み込むコードが書かれていますが、それらアセットはリポジトリには含まれていません。

There are code for loading a model data and textures, but those assets are NOT included in this repository.

----
2015 [@Shocker_0x15](https://twitter.com/Shocker_0x15)
