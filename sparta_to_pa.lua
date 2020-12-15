--[[
sparta_to_pa.lua
Converts SPARTA data file to SIMION PA files.

SPARTA is a parallel DSMC code for performing simulations of low-density 
gases in 2D or 3D.
http://sparta.sandia.gov/

This will create one SIMION PA for each scalar component.

Some adjustment to this program may be necessary for your own purposes.
Here it is assumed that the SPARTA data files have been produced by the
following SPARTA commands:
compute   1 grid all air u v w nrho temp
fix       1 ave/grid all 1 10000 10000 c_1[*]
dump      1 grid all 10000 data.grid.* id xc yc zc f_1[*]

Assuming 3D simulation and SI units.

The SPARTA output then needs to be pre-processed in ParaView:
  * MergeBlocks
  * CellDataToPointData
  * Save Data... (CSV File, default settings)

The number of grid cells and box size needs to be the same as in the
SPARTA input configuration file.

v2018-01-23

Patrick Sturm
(c) 2018 TOFWERK
--]]


-- grid size and box dimensions (same as in SPARTA input file)
local nx = 240 
local ny = 72
local nz = 72

local xmin = 0
local xmax = 0.12
local ymin = 0
local ymax = 0.018
local zmin = 0
local zmax = 0.018

local dx = (xmax - xmin) / nx * 1000  -- mm
local dy = (ymax - ymin) / ny * 1000  -- mm
local dz = (zmax - zmin) / nz * 1000  -- mm


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

-- Rounds to nearest integer.
local function round(x) return math.floor(x + 0.5) end


--[[
Convert SPARTA file (filename) to SIMION PA (pa_filename_prefix).
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
  local xs = cols[10]
  local ys = cols[11]
  local zs = cols[12]


  -- Create PA object.
  local pa = simion.pas:open()
  pa:size(nx+1,ny+1,nz+1)  -- Note for example that nx = 10 means ten grid points that define the boundary of 9 cells, giving a physical size of 9 grid units.
  pa.refinable = false  -- prevent SIMION prompting to refine this PA.
  pa.mirror_x = false
  pa.mirror_y = true
  pa.mirror_z = true
  pa.dx_mm = dx
  pa.dy_mm = dy
  pa.dz_mm = dz

  -- Copy values from input array to SIMION PA.
  local names = {"vx", "vy", "vz", "p", "t"}
  -- velocity, m/s
  for j=1,3 do
    local col = cols[j]
    for i=1,#col do
      local xi = round((xs[i]*1000 - xmin) / dx)
      local yi = round((ys[i]*1000 - ymin) / dy)
      local zi = round((zs[i]*1000 - zmin) / dz)
      pa:potential(xi,yi,zi, col[i])
    end
    pa:save(pa_filename_prefix .. '_' .. names[j] .. '.pa')
  end
  -- temperature, K
  local col = cols[5]
  for i=1,#col do
    local xi = round((xs[i]*1000 - xmin) / dx)
    local yi = round((ys[i]*1000 - ymin) / dy)
    local zi = round((zs[i]*1000 - zmin) / dz)
    pa:potential(xi,yi,zi, col[i])
  end
  pa:save(pa_filename_prefix .. '_' .. names[5] .. '.pa')
  -- pressure, Pa
  local col_rho = cols[4]
  for i=1,#col do
    local xi = round((xs[i]*1000 - xmin) / dx)
    local yi = round((ys[i]*1000 - ymin) / dy)
    local zi = round((zs[i]*1000 - zmin) / dz)
    pa:potential(xi,yi,zi, col_rho[i]*8.31446*col[i]/6.022e23)
  end
  pa:save(pa_filename_prefix .. '_' .. names[4] .. '.pa')

  -- Close input file.
  fh:close()

end

-- Convert files.
convert_file('C:/Users/pst/SPARTA/Skimmer_BSQ_sparta_12/data.grid.40000.csv', 'sparta')
