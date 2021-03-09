--[[
sparta_to_pa.lua
Converts SPARTA grid data to SIMION PAs (potential arrays).

SPARTA is a parallel DSMC code for performing simulations of low-density 
gases in 2D or 3D.
http://sparta.sandia.gov/

This will create one SIMION PA for each scalar component.

The SPARTA grid data first need to be converted to ParaView .pvd format
with grid2paraview.py (from SPARTA tools). Then in ParaView the cell 
data are interpolated to the grid points using the CellDataToPointData 
filter. If the SPARTA grid is not a regular grid, then the data need to be 
resampled: Import a regular grid into ParaView and resample the 
interpolated point data to the new grid using the ResampleWithDataset filter 
(with the regular grid as the destination mesh). Finally the data can be 
saved as csv file (Save data... with default settings).

Some adjustment to this script may be necessary for your own purposes.

Here, it is assumed that it is a 3D simulation with y and z mirroring 
and that all values are in SI units.

Patrick Sturm
(c) 2021 TOFWERK
--]]

local filename = "C:/Users/pst/SPARTA/paraview_export_test/data.csv"

-- box size and number of grid cells (same as in SPARTA csv file)
local xmin = 0  -- m
local xmax = 1  -- m
local ymin = 0  -- m
local ymax = 1  -- m
local zmin = 0  -- m
local zmax = 1  -- m

local nx = 4 
local ny = 4
local nz = 4


-- data columns in csv file
local col_x = 13  -- x coordinate, m
local col_y = 14  -- y coordinate, m
local col_z = 15  -- z coordinate, m
local col_vx = 2  -- x velocity, m/s
local col_vy = 3  -- y velocity, m/s
local col_vz = 4  -- z velocity, m/s
local col_rho = 5  -- number density, particles/m^3
local col_t = 6  -- temperature, K

-- cell size
local dx = (xmax-xmin)/nx*1000  -- mm
local dy = (ymax-ymin)/ny*1000  -- mm
local dz = (zmax-zmin)/nz*1000  -- mm


-- string split function, http://lua-users.org/wiki/SplitJoin
local function split(s, pat)
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


-- parse_numbers function
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
	local fh = assert(io.open(filename, "rb"))

	-- Read header data.
	local line = fh:read"*l"
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
	local x = cols[col_x]
	local y = cols[col_y]
	local z = cols[col_z]
	local vx = cols[col_vx]
	local vy = cols[col_vy]
	local vz = cols[col_vz]
	local rho = cols[col_rho]
	local t = cols[col_t]

	-- Convert the coordinates to PA integer points.
	for i=1,#x do
		x[i] = round((x[i]*1000-xmin)/dx)
		y[i] = round((y[i]*1000-ymin)/dy)
		z[i] = round((z[i]*1000-zmin)/dz)
	end


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
	-- x-velocity, m/s
	for i=1,#x do
		pa:potential(x[i],y[i],z[i], vx[i])
	end
	pa:save(pa_filename_prefix.."_vx.pa")
	-- y-velocity, m/s
	for i=1,#x do
		pa:potential(x[i],y[i],z[i], vy[i])
	end
	pa:save(pa_filename_prefix.."_vy.pa")
	-- z-velocity, m/s
	for i=1,#x do
		pa:potential(x[i],y[i],z[i], vz[i])
	end
	pa:save(pa_filename_prefix.."_vz.pa")
	-- temperature, K
	for i=1,#x do
		pa:potential(x[i],y[i],z[i], t[i])
	end
	pa:save(pa_filename_prefix.."_t.pa")
	-- pressure, Pa
	for i=1,#x do
		pa:potential(x[i],y[i],z[i], rho[i]*8.31446*t[i]/6.022e23)
	end
	pa:save(pa_filename_prefix.."_p.pa")

	-- Close input file.
	fh:close()

end

-- Convert files.
convert_file(filename, "sparta")
