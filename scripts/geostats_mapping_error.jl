using GeoStats
using DataFrames

grid = CartesianGrid(10,10)
# grid = CartesianGrid((10,10), (1,1), (1,1)) # works!
# grid = CartesianGrid((10,10), (0,0), (1,1)) # Does not work!

coordinates = [
  7  6;
  4  4]
 # 4  4  4  4  1  3;
 # 6  7  3  5  5  7]
# table = DataFrame(value=rand(size(coordinates,2)))
table = DataFrame(value=[11000, -100])
cdomain = PointSet(coordinates)
sdata = georef(table, cdomain)
mapping = map(sdata, grid, (:value,), NearestMapping())[:value]

@info coordinates'
@info mapping
@info length(mapping) == size(coordinates,2)

dlocs = Int64[m[1] for m in mapping]
D = [centroid(grid, i) for i in dlocs]
