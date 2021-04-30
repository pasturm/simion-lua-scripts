--[[
openfoam2pa.lua
Converts OpenFOAM data file to SIMION PA files.

OpenFOAM is an open source software for computational fluid dynamics (CFD) 
https://openfoam.org/

This will create one SIMION PA for each scalar component.

In ParaView apply a PointVolumeInterpolator filter with the box size 
the same as the PA size and the resolution and cell size as desired.
Then save the data as csv file (Choose T, U and p arrays and select
Point Data).

Some adjustment to this program may be necessary for your own purposes.

Patrick Sturm
(c) 2018-2021 TOFWERK
--]]

local filename = "C:/Users/pst/OpenFOAM/data.csv"

-- grid size (same as in CSV input file)
local nx = 52
local ny = 52
local nz = 104

-- domain size, in mm (same as in CSV input file)
local xmin = -13
local xmax = 13
local ymin = -13
local ymax = 13
local zmin = -12
local zmax = 40

local dx = (xmax - xmin) / nx
local dy = (ymax - ymin) / ny
local dz = (zmax - zmin) / nz


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
  pa:size(nx+1,ny+1,nz+1)  -- Note for example that nx = 10 means ten grid points that define the boundary of 9 cells, giving a physical size of 9 grid units.
  pa.refinable = false  -- prevent SIMION prompting to refine this PA.
  pa.dx_mm = dx
  pa.dy_mm = dy
  pa.dz_mm = dz

  -- Copy values from input array to SIMION PA.
  local names = {"t", "vx", "vy", "vz", "p", "rho"}
  for j=1,5 do
    local col = cols[j]
    for i=0,#col-1 do
      local xi = i % (nx+1)
      local yi = math.floor(i/(nx+1)) % (ny+1)
      local zi = math.floor(i/((nx+1)*(ny+1)))
      pa:potential(xi,yi,zi, col[i+1])
    end
    pa:save(pa_filename_prefix .. '_' .. names[j] .. '.pa')
  end

  -- Close input file.
  fh:close()

end

-- Convert files.
convert_file(filename, 'openfoam')
