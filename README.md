# GH_Wind
Wind Simulation (FFD) plugin for Rhinoceros Grasshopper.

Please refer to this publication for citation: [Waibel et al. (2017)](http://www.ibpsa.org/proceedings/BS2017/BS2017_582.pdf). It also includes validation and comparison studies with OpenFOAM and wind tunnel experiments.

<br><br>

## Installation 
The current build can be found in the Tutorials folder. Place [FastFluidSolverMT.dll](https://github.com/christophwaibel/GH_Wind/tree/master/Tutorials/FastFluidSolverMT.dll) and [GHWind.gha](https://github.com/christophwaibel/GH_Wind/tree/master/Tutorials/GHWind.gha) into your Rhino Grasshopper components folder and you are good to go.

## Examples
Try out the `*.gh` and `*.3dm` files in the [Tutorials](https://github.com/christophwaibel/GH_Wind/tree/master/Tutorials) folder as examples. The folder also contains two optimization examples (`200528_Highrise.gh`and `200528_Breathability.gh`) from [Waibel et al. (2019)](http://www.ibpsa.org/proceedings/BS2019/BS2019_210621.pdf). See the paper for a description of the problems.

![alt text](https://github.com/christophwaibel/GH_Wind/blob/master/Documentation/wind_optimization.png "Airflow optimization. Left: Surface pressure minimization of a high rise building. Right: Improving natural ventilation of 4 buildings. From: [Waibel et al. (2019)](http://www.ibpsa.org/proceedings/BS2019/BS2019_210621.pdf)")
*Airflow optimization. Left: Surface pressure minimization of a high rise building. Right: Improving natural ventilation of 4 buildings. From: [Waibel et al. (2019)](http://www.ibpsa.org/proceedings/BS2019/BS2019_210621.pdf)*

<br><br>

## Some GIFs

![alt text](https://github.com/christophwaibel/GH_Wind/blob/master/Documentation/slide0005_image017.gif "Image from Rhino")

*Capture from Rhino. Serves only as an illustrative example, don't make such a small domain in your simulations!*


![alt text](https://github.com/christophwaibel/GH_Wind/blob/master/Documentation/image23.gif "Image from Rhino")

*Capture from Rhino, 50x50x40 mesh, 10 times speed up*
