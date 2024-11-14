# crater-id
Read a crater data file, create triads of craters with associated invariants, search database for matching triads based on images, then estimate the camera position using only the image for navigation. This is based on the work of Christian "Lunar Crater Identification in Digital Images" (2020) and my own work "Optical Navigation using Quadrics and Indirect Invariant Indexing" (2022).

Starting with a set (3+) of identified craters on the lunar surface, we backproject the image ellipses to the disk quadrics which represent the crater rims. This backprojection gives the distance from the camera to the quadric, albeit with two possible orientations for the normal. Using the relationship between the know quadrics lets us choose the correct normals. With the quadrics known in both the camera and world frame, we can estimate (in the least squares sense) the camera's pose (position and orientation).

## Installing
This program requires OpenCV and Eigen3 libraries, which manage the 2D visualizations and matrix operations, respectively. If you want to have VTK for 3D visualization, this is useful to see where things are in 3D space. I have the Python scripts that use PyVista, which is a wrapper for VTK and is very useful. 

I use Meson as the build tool which is good for my purposes here. To download and install,

```bash
git clone https://github.com/drenshaw/crater-id.git
cd crater-id
meson setup <builddir>
cd <builddir>
```
where ```<builddir>``` can be named whatever you want (I just use "build" or "builddir"). Once this is done, we need to build:

```bash
cd <builddir>
meson compile
```

## Running tests
I have a set of over 100 unit tests testing everything from basic matrix math to homography operations. You can enable or disable the test suites (e.g., math, navigation, conics) as you see fit, which are included in the `meson.build` file in the `test` directory. Ensure that you run `meson compile` again if you make changes in that file. Then you can run the tests by running:

```bash
./test/test_cid
```
