# Online Multiple Camera Multiple Target Tracking (MCMTT)


This repository is a Visual C++ implementation of  [Online Scheme for Multiple Camera Multiple Target Tracking Based on Multiple Hypothesis Tracking](http://ieeexplore.ieee.org/document/7517399/), but with some legacy codes like frame grabbing with a commercial camera. The code runs best under a multi-core framework.

This repository also contains two datasets names PETS2009.S2.L1 and PILSNU. The formal one is widely used in a MCMTT field and the later one is made by ourselves. All dataset contains frame images, calibration informations, pedestrian detection results and ground truthes of the tracking. You can test the code after unzip files under 'OnlineMCMTT/dataset/' folder and set a target dataset in 'OnlineMCMTT/data/settings.xml' file.

The code is guaranteed to run with Microsoft Visual Studio and OpenCV 3.1. You can run the code with other version of OpenCV by modifying linking options in the project preference. But I strongly recommand you to use 3.1 or later versions because the previous versions have different function calls for feature point tracking.

Requirements
------------

The following softwares and libraries are required:

- Microsoft Visual Studio 2015
- OpenCV 3.1

Notice
------------
- 2017.03.21: There is something wrong at the evaluation code. So, don't trust the evaluation result.

()()

('')HAANJU.YOO