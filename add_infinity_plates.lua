--[[
 Script to add "infinity plates" to a gridless planar mirror PA.

 Patrick Sturm
 (c) 2025 TOFWERK
 --]]

local path = "E:/Patrick/SIMION/InfinityPlates/"
local filename = "GLPM3_v2"
local pa = simion.pas:open(path..filename..".pa#")
local pa0 = simion.pas:open(path..filename..".pa0")

-- number of infinity electrodes
local n1 = 9  -- number of plates between V5 and Drift/V4
local n2 = 8  -- number of plates between Drift/V4 and V3
local n3 = 5  -- number of plates between V3 and V2
local n4 = 6  -- number of plates between V2 and V1

-- Electrode voltages
local V1 =  1580
local V2 =  1115
local V3 =  927.5
local V4 =  0
local V5 =  -1000
pa0:fast_adjust{[1] = V1, [2] = V2, [3] = V3, [4] = V4, [5] = V5}

-- calculate midpoints (i.e. voltages of infinity electrodes)
local midpoints = {}
for i=1,n1 do
	midpoints[i] = V5 + (i-1)*(V4-V5)/(n1-1)
end
for i=1,n2-1 do
  midpoints[#midpoints+1] = V4 + i*(V3-V4)/(n2-1)
end
for i=1,n3-1 do
	midpoints[#midpoints+1] = V3 + i*(V2-V3)/(n3-1)
end
for i=1,n4-1 do
	midpoints[#midpoints+1] = V2 + i*(V1-V2)/(n4-1)
end

-- break points (i.e. gaps between infinity electrodes) 
local breaks = {}
breaks[1] = V5
breaks[#midpoints+1] = V1
for i=2,#midpoints do
	breaks[i] = (midpoints[i] + midpoints[i-1])/2
end

-- function to bin the value x into the intervals given by "breaks".
local function bin(x, breaks)
	local lo = 1
	local hi = #breaks
	while ((hi - lo) >= 2) do
		local new = math.floor((hi + lo)/2)
		if (x >= breaks[new]) then
			lo = new
		else
			hi = new
		end
	end
	return lo
end

-- replace edge of mirror with infinity plates
local z_edge = 220/pa.dz_mm  -- z position of edge
local H = 17.5/pa.dy_mm  -- half-height of GLPM
local d = 10/pa.dy_mm  -- height of electrode edge
local x_pb = 525/pa.dx_mm  -- x position of back plate
pa:fill { 
	function(x,y,z)
		if (y<=H+d and x<(x_pb+d) and z>=z_edge and z<=(z_edge+d)) then
			local v = pa0:potential(x,y,0)
			local el = bin(v, breaks) + 5  -- assign electrode numbers 6,7,.. to infinity plates
			return el, true
		end
	end, surface='none', units = 'gu'
}

-- circle function to generate larger gaps
local function circle(x, y, r, array)
  local dx, dy, err = r, 0, 1-r
  while dx >= dy do
    for i=-dx+1,dx do
      table.insert(array, {x+i, y+dy})
      table.insert(array, {x+i, y-dy})
      table.insert(array, {x+dy, y+i})
      table.insert(array, {x-dy, y+i})
    end
    dy = dy + 1
    if err < 0 then
      err = err + 2 * dy + 1
    else
      dx = dx-1
      err = err + 2 * (dy - dx) + 1
    end
  end
end

-- make 2D array of gap points
local gap = 1  -- gap width radius in grid units, gives 2*gap+1 wide gap.
local gap_array = {}
for x=0,x_pb do
	for y=0,H+d do
		local diff1 = pa:potential(x,y,z_edge) - pa:potential(x+1,y,z_edge)
		local diff2 = pa:potential(x,y,z_edge) - pa:potential(x,y+1,z_edge)
		if (diff1~=0 or diff2~=0) then
			circle(x,y,gap,gap_array)
		end
	end
end

-- add gaps
for i,v in pairs(gap_array) do
	if (v[1]>=0 and v[2]>=0) then
		for z=z_edge,z_edge+d do
			pa:electrode(v[1],v[2],z, false)
		end
	end
end

pa:save(path..filename.."_"..tostring(n1)..tostring(n2)..tostring(n3)..tostring(n4).."_infinity.pa#")

-- make pa+ file.
local paplus = "local n1 = "..tostring(n1).."\nlocal n2 = "..tostring(n2)..
               "\nlocal n3 = "..tostring(n3).."\nlocal n4 = "..tostring(n4)..[[

local map1, map2, map3, map4, map5 = {}, {}, {}, {}, {}
for i=1,(5+n1+n2+n3+n4-3) do
  map1[i] = 0
  map2[i] = 0
  map3[i] = 0
  map4[i] = 0
  map5[i] = 0
end
map1[1] = 1
map2[2] = 1
map3[3] = 1
map4[4] = 1
map5[5] = 1

local from = 6
local to = from + n1 - 1
for i=from,to do
  map4[i] = (i-from)/(n1-1)
  map5[i] = 1 - (i-from)/(n1-1)
end

from = to
to = from + n2 - 1
for i=from,to do
  map3[i] = (i-from)/(n2-1)
  map4[i] = 1 - (i-from)/(n2-1)
end

from = to
to = from + n3 - 1
for i=from,to do
  map2[i] = (i-from)/(n3-1)
  map3[i] = 1 - (i-from)/(n3-1)
end

from = to
to = from + n4 - 1
for i=from,to do
  map1[i] = (i-from)/(n4-1)
  map2[i] = 1 - (i-from)/(n4-1)
end

potential_array {
  scalable_electrodes = {
    [1] = map1,
    [2] = map2,
    [3] = map3,
    [4] = map4,
    [5] = map5
  }
}
]]
local fh = assert(io.open(path..filename.."_"..tostring(n1)..tostring(n2)..tostring(n3)..tostring(n4).."_infinity.pa+", "wb"))
fh:write(paplus)
fh:close()

pa:refine {convergence=1e-5}
simion.pas:close()
