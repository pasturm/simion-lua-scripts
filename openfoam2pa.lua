--[[
openfoam2pa.lua
Converts OpenFOAM data file to SIMION PA files.

OpenFOAM is an open source software for computational fluid dynamics (CFD) 
https://openfoam.org/

This will create one SIMION PA for each scalar component (T, U_x, U_y, U_z, p).

How to export OpenFOAM data:
In ParaView apply a MergeBlocks filter and then a PointVolumeInterpolator 
filter with the box size and the resolution as desired (does not need to 
be the same as in the SIMION geometry PA). Note that a resolution of 
e.g. 10x10x10 gives 11x11x11 grid points.
Then save the data as csv file (Choose T, U and p arrays and select
Point Data).

How to run this script:
Adjust the filename, number of grid points, domain size and symmetry
and run the script.

Patrick Sturm
(c) 2018-2021 TOFWERK
--]]

-------------------------------------------------------------------------------
local filename = "C:/Users/pst/OpenFOAM/data.csv"

-- number of grid points (same as in csv input file)
local nx = 11
local ny = 11
local nz = 11

-- domain size, in mm (same as in csv input file)
local xmin = -5
local xmax = 5
local ymin = -5
local ymax = 5
local zmin = -5
local zmax = 5

-- set symmetry and mirroring type, e.g. '3dplanar[xy]' 
local symmetry = '3dplanar'  
-------------------------------------------------------------------------------


local dx = (xmax - xmin) / (nx-1)
local dy = (ymax - ymin) / (ny-1)
local dz = (zmax - zmin) / (nz-1)

-- string split http://lua-users.org/wiki/SplitJoin
function split(s, pat)
  local st, g = 1, s:gmatch("()("..pat..")")
  local function getter(self, segs, seps, sep, cap1, ...)
    st = sep and seps + #sep
    return self:sub(segs, (seps or 0) - 1), cap1 or sep, ...
  end
  local function splitter(self)
    if st then return getter(self, st, g()) end
  end
  return splitter, s
end

local function parse_numbers(line)
  local nums = {}
  for word in split(line, ",") do
    local val = tonumber(word)
    table.insert(nums, val)
  end
  return nums
end


--[[
Convert OpenFOAM file (filename) to SIMION PA (pa_filename_prefix).
--]]
local function convert_file(filename, pa_filename_prefix)

  -- Open input file.
  local fh = assert(io.open(filename, 'rb'))

  -- Read header data.
  local line = fh:read'*l'
  local header = {}
  for word in split(line, ",") do
    table.insert(header, word)
  end

  -- Load columns of values.
  local cols = {}
  for i=1, #header do
    table.insert(cols, {})
  end
  for line in fh:lines() do
    local nums = parse_numbers(line)
    for i=1,#nums do
      table.insert(cols[i], nums[i])
    end
  end

  -- Create PA object.
  local pa = simion.pas:open()
  pa:size(nx,ny,nz)  -- Note for example that nx = 11 means ten grid points that define the boundary of 10 cells, giving a physical size of 10 grid units.
  pa.symmetry = symmetry  -- set symmetry and mirroring type
  pa.refinable = false  -- prevent SIMION prompting to refine this PA.
  pa.dx_mm = dx
  pa.dy_mm = dy
  pa.dz_mm = dz

  -- Copy values from input array to SIMION PA.
  local names = {"t", "vx", "vy", "vz", "p", "rho"}
  for j=1,5 do
    local col = cols[j]
    for i=0,#col-1 do
      local xi = i % nx
      local yi = math.floor(i/nx) % ny
      local zi = math.floor(i/(nx*ny))
      pa:potential(xi,yi,zi, col[i+1])
    end
    pa:save(pa_filename_prefix .. '_' .. names[j] .. '.pa')
  end

  -- Close input file.
  fh:close()

end

-- Convert files.
convert_file(filename, 'openfoam')
