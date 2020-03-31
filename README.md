# NWTMasterSample

R Code to accompany a Northwest Territories Freshwater Master Sample built using Halton Iterative Partitioning.

## HIP Master Sample:

A HIP Master Sample works by first discretizing space using a Halton frame (Robertson et al. 2017). 
For Northwest Territories these are roughly 100 m x 60 m, or a football field. The frame's box is defined as a bounding box
of the entire NWT, and then given a random start and single permutations introduced in Robertson et al. (2018) each Halton box
is given an ordering which we refer to as the Master Sample Index (MSI). 

When there are large gaps in the resource this leads to many skips in Halton boxes for sample selection resulting in
poor spatial balance, roughly simple random sampling. As a result, the MSI is not necessarily an effective Master Sample on its own.
Instead, we use it to stitch together Halton Iterative Partitioning (HIP) samples (Robertson et al. 2018). For any individual sample,
we partition the space using base 2 along the x-axis and base 3 along the y-axis into equal amounts of resource into $B = 2^{J_1} \times 3^{J_2}$ boxes
that are determined by the user. We then choose the HIP box that has the minimum MSI within it as the start of the sample. Then, 
based on random permutations, we select the remaining samples using HIP methodology with the minimum MSI in each HIP, or primary, box. 
For the individual user the sample is a well balanced HIP sample, equi-probable and well spatially balanced. They are also flexible for
any sample size and configuration required by the user that HIP provides.

The MSI offers a coordinating ordering insuring that the samples match and spatially scale between users. It adds an underlying determinism
to HIP which makes this a Master Sample. If B is very large for a small region, then it may not spatially scale very well but will have excellent
spatial balance. On the contrary, if B = 1, this is equivalent to simple random sampling and the minimum set of MSI are chosen.

We are currently working on a manuscript to accompany this work so please be patient for it to make more sense in the future.

## R Package Functionality

At this stage this is just geared for lakes in NWT. However, in the background we have existing code to select random points on
stream networks as well. This works by using one dimensional BAS on the stream network within a Halton box to select any number
of points randomly based on the MSI. We have implemented this for a New Zealand freshwater streams Master Sample but not wrapped it in this
particular version. The code is included but the user will require a bit more work to do this which we currently don't provide documentation for.
More to come but feel free to contact me if you want some base code that will do this.

To adapt this version of the R package for another Master Sample, it only needs the random seeds and bounding box shifted as well
as the base projection hard coded in the functions changed.

In the background we also provide the code to create the Halton Iterative Partitioning boxes. All this will be documentated in the future
as this package is generalized and grows.

For now you can use the function below to create the HIP boxes.
hip must be a data.table with columns X, Y, index = MasterSampleIndex.
Easiest way to create it is to call getIndividualBoxIndices() first. With that output you can then call 
getHIPBoxes().

To check spatial balance use. Best to test this against multiple different partitions to see the impact on balance.
getBalance()

To sample a river in the style of the NZ FW Master Sample you can use:
getPolyRast()

This code essentially is driven by:
getHipSample()
You can pass any random seed, bounding box etc and return a HIP sample from this.

Please see the attached vignettes for some example use with lakes in the Upper Coppermine River Basin.

To install:

library(devtools)

install_github("paul-vdb/NWTMasterSample", build_vignettes = TRUE)

library(NWTMasterSample)

browseVignettes("NWTMasterSample")

