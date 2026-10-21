--[[
COMSOL -> SIMION magnetic-field PA converter

Input: COMSOL text file
Output: Bx.pa, By.pa, Bz.pa

Patrick Sturm
(c) TOFWERK, 2026
--]]

local filename = "20261009 Export magnetic field.txt"
local output_bxpa = "Bx.pa"
local output_bypa = "By.pa"
local output_bzpa = "Bz.pa"

-- comsol grid
local dx = 0.5
local dy = 0.5
local dz = 0.5

local xmin = -10.0
local ymin = -20.0
local zmin = -150.0

local xmax = 10.0
local ymax = 180.0
local zmax = 30.0

local nx = (xmax - xmin)/dx + 1
local ny = (ymax - ymin)/dy + 1
local nz = (zmax - zmin)/dz + 1

local tesla_to_gauss = 10000

-- create PAs
local function make_pa(filename)
    local pa = simion.pas:open()
    pa:size(nx, ny, nz)
    pa.symmetry = '3dplanar'
    pa.dx_mm = dx
    pa.dy_mm = dy
    pa.dz_mm = dz
    pa.potential_type = 'magnetic'
    pa.refinable = false
    return pa
end

local bxpa = make_pa(output_bxpa)
local bypa = make_pa(output_bypa)
local bzpa = make_pa(output_bzpa)

-- open COMSOL file
local fh = assert(io.open(filename, "r"))

--parse line: x y z Bx By Bz
local function parse_line(line)
    local values = {}
    for token in line:gmatch("%S+") do
        local value = tonumber(token)
        table.insert(values, value)
    end
    return values
end

function round(num)
    return num >= 0 and math.floor(num + 0.5) or math.ceil(num - 0.5)
end

-- read file and set potential in PA
for line in fh:lines() do

    -- Ignore COMSOL comments/header
    if not line:match("^%s*%%") and line:match("%S") then

        -- Parse line
        local vals = parse_line(line)
        local x, y, z, bx, by, bz = unpack(vals)

        local xi = round((x*1000 - xmin) / dx)  -- m to mm conversion
        local yi = round((y*1000 - ymin) / dy)  -- m to mm conversion
        local zi = round((z*1000 - zmin) / dz)  -- m to mm conversion

        bxpa:point(xi, yi, zi, bx*tesla_to_gauss, false)
        bypa:point(xi, yi, zi, by*tesla_to_gauss, false)
        bzpa:point(xi, yi, zi, bz*tesla_to_gauss, false)
    end
end

fh:close()

bxpa:save(output_bxpa)
bypa:save(output_bypa)
bzpa:save(output_bzpa)

bxpa:close()
bypa:close()
bzpa:close()
