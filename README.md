The function takes 3D line in XYZ format and decimation. The function finds circles passing through 3 points (triplet). Triplets are formed by decimated points. 
Decimation stands for points down sampling. I.e. decimation = 1 means that 3 adjacent points are taken, while decimation = 2 means 3 second points are taken. Obviously greater decimation means lower resolution.

In UI select start/end points of the line and try various decimations. If decimation steps are too high due to large amount of points try to reduce points range first. Then select decimation and select required points range again.
