Online Multiple Camera Multiple Target Tracking (MCMTT)
=============

This repository is a Visual C++ implementation of  [Online Scheme for Multiple Camera Multiple Target Tracking Based on Multiple Hypothesis Tracking](http://ieeexplore.ieee.org/document/7517399/), but with some legacy codes like frame grabbing with a commercial camera. The code runs best under a multi-core system.

This repository also contains two datasets names PETS2009.S2.L1 and PILSNU. The formal one is widely used in an MCMTT field and the later one is made by ourselves. All dataset contain frame images, calibration information, pedestrian detection results and ground truths of the tracking. You can test the code after unzip files under 'OnlineMCMTT/dataset/' folder and set a target dataset in 'OnlineMCMTT/data/settings.xml' file.

The code is guaranteed to run with Microsoft Visual Studio and OpenCV 3.1. You can run the code with other versionㄴ of OpenCV by modifying linking options in the project preference. But I strongly recommend you to use 3.1 or later versions because the previous versions have different function calls for feature point tracking.

Requirements
------------

The following softwares and libraries are required:

- (Operating System) Microsoft's Windows
- Microsoft's Visual Studio 2015
- OpenCV 3.1


()()  
(' ')HAANJU.YOO