
DIRECTORIES
===========
./curvature
380 files containing the per-vertex curvature and curvature derivatives
information for our mesh data set. Line i has 6 numbers for vertex i: the
first two are the curvatures along the major and minor principle axes u and v;
the next four are the curvature derivatives dcurv(u,u), dcurv(u,v), dcurv(v,u)
and dcurv(v,v). The curvature related information is computed with "trimesh2",
see "www.cs.princeton.edu/gfx/proj/trimesh2" for more details


./seg
Segmentation results from human(Benchmark) and 7 algorithms (5 taking number
of segments as input: FitPrim, KMeans, NormCuts, RandCuts, RandWalks; 2 not
taking number of segments as input: CoreExtra and ShapeDiam). You can add new
segmentation results here as ".seg" files. 

In a ".seg" file, there are nFaces number of lines. Line i contains the
segment number of face i. Segment numbers should be within 0 to (numSegments-1).


./off
380 mesh models in "off" format (sequence number:1-260 and 281-400) from 
Watertight Track of SHREC 2007 (http://watertight.ge.imati.cnr.it/)

They belong to 19 categories:
     1 -  20  Human
    21 -  40  Cup
    41 -  60  Glasses
    61 -  80  Airplane
    81 - 100  Ant
   101 - 120  Chair
   121 - 140  Octopus
   141 - 160  Table
   161 - 180  Teddy
   181 - 200  Hand
   201 - 220  Plier
   221 - 240  Fish
   241 - 260  Bird
  (261 - 280) Spring (excluded from our study)
   281 - 300  Armadillo
   301 - 320  Bust
   321 - 340  Mech
   341 - 360  Bearing
   361 - 380  Vase
   381 - 400  Fourleg

