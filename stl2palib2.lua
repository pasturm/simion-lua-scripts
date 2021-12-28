--[[
stl2palib2.lua

Library to convert an electrode geometry from a STL file to a SIMION potential array file.

The main advantages compared to the stl2pa command from SIMION SL Tools are:
1. The STL surface is filled with solid points (avoids the solid points 
   strategy issues of SL Tools < 8.2.0.11).
2. The electrode surface enhancement feature can be used (leads to higher 
   field accuracies).
3. Multiple electrodes can be generated from a single STL file.
4. Array mirroring can be used.
5. Bounding boxes can easily be added.
6. Ideal and real wire grids can be added.

Note: SIMION 8.2.0.11-20210626 SL Tools now contains similar functionality. But there are still bugs
and it is slower than this lua script.

Example:
	-- load STL2PA library
	local STL2PA = simion.import("./stl2palib2.lua")

	-- STL input file
	local stl_filename = "my_geometry.stl"

	-- region min and max points, mm
	local xmin = -50
	local xmax = 130
	local ymin = -105
	local ymax = 130
	local zmin = -30
	local zmax = 30

	-- grid cell size, mm
	local dx_mm = 1
	local dy_mm = 1
	local dz_mm = 1

	-- surface enhancement
	local surface = "fractional"

	-- bounding box
	local boundingbox = "-x+x-y+y-z+z"

	-- array symmetry and mirroring
	local symmetry = "3dplanar[yz]"

	-- assign voltages to fully connected surfaces closest to the given point
	-- {voltage, {x, y, z}}
	local electrodes = {
	 	{1, {0,-75,-70}},
	 	{2, {0,-75,-66}}
	}

	-- run stl2pa conversion
	STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, surface, symmetry, electrodes)

	-- modify or add one electrode in .pa#
	-- STL2PA.modify(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, surface, electrodes)

	-- add bounding box in .pa# file
	-- STL2PA.boundingBox(stl_filename, boundingbox, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm)

	-- remove an electrode from .pa#
	-- removeElectrode(stl_filename, electrode_no)


Patrick Sturm
(c) 2020-2022 TOFWERK
--]]

local STL2PA = {}


-- Convert string of bytes (4 bytes) to IEEE-754 floating point ---------------
-- http://lua-users.org/lists/lua-l/2010-03/msg00911.html
local function bytes2float(x)
	local sign = 1
	local mantissa = string.byte(x, 3) % 128
	for i = 2, 1, -1 do mantissa = mantissa * 256 + string.byte(x, i) end
	if string.byte(x, 4) > 127 then sign = -1 end
	local exponent = (string.byte(x, 4) % 128) * 2 + math.floor(string.byte(x, 3) / 128)
	if exponent == 0 then return 0 end
	mantissa = (math.ldexp(mantissa, -23) + 1) * sign
	return math.ldexp(mantissa, exponent - 127)
end


-- Convert string of bytes (4 bytes) to 32-bit integer ------------------------
local function bytes2int(x)
	local i1,i2,i3,i4 = string.byte(x, 1, 4)
	-- return i1 + 2^8*i2 + 2^16*i3 + 2^24*i4
	-- return i1 + 256*i2 + 65536*i3 + 16777216*i4
	return ((i4*256 + i3)*256 + i2)*256 + i1
end


-- Check if STL file is in ASCII or binary format -----------------------------
local function isascii(stl_filename)
	local fh = io.open(stl_filename, "r")
	local line = fh:read("*l")
	fh:close()
	if string.find(line, "^solid ") then
		return true
	else
		return false
	end
end


-- Read a STL file ------------------------------------------------------------
local function readSTL(stl_filename)
	local t_faces = {}
	local current_face = {}
	if isascii(stl_filename) then
		-- parse stl file (from http://lua-users.org/wiki/StlToObj)
		local fh = io.open(stl_filename, "r")
		local count = 0
		for line in fh:lines() do
			if string.find(line, "^%s*facet normal") then
				local x,y,z = string.match(line, "(%S+)%s(%S+)%s(%S+)$")
				table.insert(current_face, {tonumber(x),tonumber(y),tonumber(z)})
			elseif string.find(line, "^%s*vertex") then
				local x,y,z = string.match(line, "(%S+)%s(%S+)%s(%S+)$")
				table.insert(current_face, {tonumber(x),tonumber(y),tonumber(z)})
				count = count + 1
				if count == 3 then
					table.insert(t_faces,current_face)
					current_face = {}
					count = 0
				end
			end
		end
		fh:close()
	else
		local fh = io.open(stl_filename, "rb")
		fh:seek("cur", 80)  -- ignore header
		local triangles = bytes2int(fh:read(4))  -- number of triangles
		for i=1,triangles do
			local x = bytes2float(fh:read(4))  -- normal vector
			local y = bytes2float(fh:read(4))
			local z = bytes2float(fh:read(4))
			table.insert(current_face, {x,y,z})
			local x = bytes2float(fh:read(4))  -- vertex 1
			local y = bytes2float(fh:read(4))
			local z = bytes2float(fh:read(4))
			table.insert(current_face, {x,y,z})
			local x = bytes2float(fh:read(4))  -- vertex 2
			local y = bytes2float(fh:read(4))
			local z = bytes2float(fh:read(4))
			table.insert(current_face, {x,y,z})
			local x = bytes2float(fh:read(4))  -- vertex 3
			local y = bytes2float(fh:read(4))
			local z = bytes2float(fh:read(4))
			table.insert(current_face, {x,y,z})
			fh:read(2)  -- ignore attribute byte count
			table.insert(t_faces,current_face)
			current_face = {}
		end
		fh:close()
	end
	return t_faces
end


-- Calculate Axis Aligned Bounding Boxes of the triangles ---------------------
-- Note: only 2D projection is needed.
local function boundingBoxes(t_faces)
	local t_bb = {}
	for _,v in ipairs(t_faces) do
		local min_x = math.min(v[2][1], v[3][1], v[4][1])
		local max_x = math.max(v[2][1], v[3][1], v[4][1])
		local min_y = math.min(v[2][2], v[3][2], v[4][2])
		local max_y = math.max(v[2][2], v[3][2], v[4][2])
		table.insert(t_bb, {min_x,max_x,min_y,max_y})
	end
	return t_bb
end


-- Calculate STL size ---------------------------------------------------------
local function stlSize(t_faces)
	local min_x = math.huge
	local max_x = -math.huge
	local min_y = math.huge
	local max_y = -math.huge
	local min_z = math.huge
	local max_z = -math.huge
	for _,v in ipairs(t_faces) do
		min_x = math.min(min_x, v[2][1], v[3][1], v[4][1])
		max_x = math.max(max_x, v[2][1], v[3][1], v[4][1])
		min_y = math.min(min_y, v[2][2], v[3][2], v[4][2])
		max_y = math.max(max_y, v[2][2], v[3][2], v[4][2])
		min_z = math.min(min_z, v[2][3], v[3][3], v[4][3])
		max_z = math.max(max_z, v[2][3], v[3][3], v[4][3])
	end
	return {min_x, max_x, min_y, max_y, min_z, max_z}
end


-- Check if two rectangles overlap (or touch each other) ----------------------
local function doRectaglesOverlap(xmin1,xmax1,ymin1,ymax1,xmin2,xmax2,ymin2,ymax2)
	local tol = 0.01  -- add some tolerance so no triangles are missed when ray is slightly tilted from z axis
	-- if one rectangle is on left side of other 
	if (xmin1 > (xmax2+tol) or xmin2 > (xmax1+tol)) then
	    return false
	-- if one rectangle is above other 
	elseif (ymin1 > (ymax2+tol) or ymin2 > (ymax1+tol)) then
	    return false
	else
		return true
	end
end


-- Hash the triangle bounding boxes to a 2D grid ------------------------------
-- The output is a 3D array, where the first two dimensions are
-- the grid point indices and the third dimension is a vector
-- where the n-th element is 0 if the n-th bounding box is inside the grid cell
-- and nil if it is completely outside. 
local function map2Grid(t_faces, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
	local x_min = math.max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
	local x_max = math.min(math.floor(t_size[2]/dx_mm)*dx_mm, xmax)
	local y_min = math.max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
	local y_max = math.min(math.floor(t_size[4]/dy_mm)*dy_mm, ymax)
	local t_bb = boundingBoxes(t_faces)
	local t_hash = {}
	for i=1,(x_max-x_min)/dx_mm+1 do  -- loop over x axis
		t_hash[i] = {}
		-- convert the grid point indices to min and max
		-- values of the grid cell boundary
		local x0 = (i-1)*dx_mm + x_min
		local x1 = (i)*dx_mm + x_min
		for j=1,(y_max-y_min)/dy_mm+1 do  -- loop over y axis
			t_hash[i][j] = {}
			local y0 = (j-1)*dy_mm + y_min
			local y1 = (j)*dy_mm + y_min
			for k,v in ipairs(t_bb) do  -- loop over bounding boxes
				-- check if the bounding box touches the grid cell
				if (doRectaglesOverlap(x0,x1,y0,y1,v[1],v[2],v[3],v[4])) then
					t_hash[i][j][k] = 0
				end
			end
		end
	end
	return t_hash
end

local eps = 1e-10  -- floating point tolerance (machine epsilon is about 2e-16)

-- Get the unique elements of a table -----------------------------------------
-- https://stackoverflow.com/questions/20066835/lua-remove-duplicate-elements
local function getUniqueElements(t_in)
	local hash = {}
	local res = {}
	for _,v in ipairs(t_in) do
		v = math.floor(v/eps+0.5)*eps  -- account for floating point tolerance
	   if (not hash[v]) then
	       res[#res+1] = v
	       hash[v] = true
	   end
	end
	return res
end

-- check whether a grid point is inside the STL -------------------------------
-- adapted from http://geomalgorithms.com/a06-_intersect-2.html
-- all points on STL surface are also considered to be inside
local function isInsideSTL(x,y,z, t_faces, xmin, xmax, ymin, ymax, zmin, zmax, t_hash, dx_mm, dy_mm)
	if (x<xmin or x>xmax or y<ymin or y>ymax or z<zmin or z>zmax) then -- outside of STL
		return false
	end
	local i = math.floor((x-xmin)/dx_mm) + 1  -- hash table grid index
	local j = math.floor((y-ymin)/dy_mm) + 1
	local r_array = {}  -- use this to check whether ray hits edge of two connected triangles
	for l,_ in pairs(t_hash[i][j]) do  -- loop through all relevant triangles
		-- Note: ipairs does not work here, it stops at the first nil -> use pairs 
		-- ray-triangle intersection testing:
		-- triangle vertices
		local v0_1 = t_faces[l][2][1]
		local v0_2 = t_faces[l][2][2]
		local v0_3 = t_faces[l][2][3]
		local v1_1 = t_faces[l][3][1]
		local v1_2 = t_faces[l][3][2]
		local v1_3 = t_faces[l][3][3]
		local v2_1 = t_faces[l][4][1]
		local v2_2 = t_faces[l][4][2]
		local v2_3 = t_faces[l][4][3]
		-- trianlge vectors
		local u1 = v1_1-v0_1
		local u2 = v1_2-v0_2
		local u3 = v1_3-v0_3
		local v1 = v2_1-v0_1
		local v2 = v2_2-v0_2
		local v3 = v2_3-v0_3
		local n1 = t_faces[l][1][1]
		local n2 = t_faces[l][1][2]
		local n3 = t_faces[l][1][3]
		-- ray direction
		local l1 = 1e-5
		local l2 = 1e-5
		local l3 = math.sqrt(1-2e-10)
		-- distance to plane (https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection)
		local r = (n1*(v0_1-x)+n2*(v0_2-y)+n3*(v0_3-z))/(l1*n1+l2*n2+l3*n3)
		if (r>=-eps) then  -- ray goes towards triangle
			-- intersect point of ray and plane
			local I1 = x+l1*r
			local I2 = y+l2*r
			local I3 = z+l3*r
		 	local uu = (u1*u1 + u2*u2 + u3*u3)
			local uv = (u1*v1 + u2*v2 + u3*v3)
			local vv = (v1*v1 + v2*v2 + v3*v3)
			local w1 = I1-v0_1
			local w2 = I2-v0_2
			local w3 = I3-v0_3
			local wu = (w1*u1 + w2*u2 + w3*u3)
			local wv = (w1*v1 + w2*v2 + w3*v3)
			local D = uv*uv - uu*vv
			local s = (uv*wv - vv*wu)/D
			local t = (uv*wu - uu*wv)/D
			-- I is inside or on edge or on corner of T
			if (s>=-eps and s<=1+eps and t>=-eps and (s+t)<=1+eps) then
				if (math.abs(r)<eps) then
					return true
				end
				table.insert(r_array, r)
			end
		end
	end
	if #r_array==0 then return false end  -- no intersection
	local res = getUniqueElements(r_array)
	local count = #res
	if (count%2==1) then
		return true
	else
		return false
	end
end


-- Split path/filename into path, filename and extension ----------------------
local function splitPath(path)
	return string.match(path, "(.-)([^\\/]-)%.?([^%.\\/]*)$")  -- path, filename without extension, extension
end


-- write header ---------------------------------------------------------------
local function writeHeader()
	io.write("*********************************************\n")
	io.write("STL2PA CONVERSION\n")
	io.write("Copyright (c) 2020-2022 TOFWERK\n")
	io.write("Author: Patrick Sturm\n")
	io.write("*********************************************\n\n")
	io.flush()
end


-- add bounding box in .pa# file -----------------------------------
function STL2PA.addBoundingBox(stl_filename, bounding, xmin, xmax, ymin, ymax, 
									 zmin, zmax, dx_mm, dy_mm, dz_mm)
	local path,name = splitPath(stl_filename)
	local file = path..name..".pa#"

	io.write("Adding bounding box "..bounding.."\n")
	io.flush()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(file)  -- open pa# file

	if (string.match(bounding, "%-x")) then
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(0,yi,zi, 0,true)
			end
		end
	end

	local ximax = (xmax-xmin)/dx_mm
	if (string.match(bounding, "%+x")) then
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(ximax,yi,zi, 0,true)
			end
		end
	end

	if (string.match(bounding, "%-y")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,0,zi, 0,true)
			end
		end
	end

	local yimax = (ymax-ymin)/dy_mm
	if (string.match(bounding, "%+y")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,yimax,zi, 0,true)
			end
		end
	end

	if (string.match(bounding, "%-z")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,0, 0,true)
			end
		end
	end

	local zimax = (zmax-zmin)/dz_mm
	if (string.match(bounding, "%+z")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,zimax, 0,true)
			end
		end
	end

	pa:save(file)
	simion.pas:close()  -- remove all PAs from RAM.
end


-- remove bounding box in .pa# file -----------------------------------
function STL2PA.removeBoundingBox(stl_filename, bounding, xmin, xmax, ymin, ymax, 
									 zmin, zmax, dx_mm, dy_mm, dz_mm)
	local path,name = splitPath(stl_filename)
	local file = path..name..".pa#"

	io.write("Removing bounding box "..bounding.."\n")
	io.flush()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(file)  -- open pa# file

	if (string.match(bounding, "%-x")) then
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(0,yi,zi, 0,false)
			end
		end
	end

	local ximax = (xmax-xmin)/dx_mm
	if (string.match(bounding, "%+x")) then
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(ximax,yi,zi, 0,false)
			end
		end
	end

	if (string.match(bounding, "%-y")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,0,zi, 0,false)
			end
		end
	end

	local yimax = (ymax-ymin)/dy_mm
	if (string.match(bounding, "%+y")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,yimax,zi, 0,false)
			end
		end
	end

	if (string.match(bounding, "%-z")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,0, 0,false)
			end
		end
	end

	local zimax = (zmax-zmin)/dz_mm
	if (string.match(bounding, "%+z")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,zimax, 0,false)
			end
		end
	end

	pa:save(file)
	simion.pas:close()  -- remove all PAs from RAM.
end


-- remove one electrode of an existing PA# array ------------------------------
function STL2PA.removeElectrode(stl_filename, electrode_no)
	local path,name = splitPath(stl_filename)

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(path..name..".pa#")  -- open pa# file

	-- delete electrode
	for xi,yi,zi in pa:points() do
		if pa:potential(xi,yi,zi) == electrode_no then
			pa:point(xi,yi,zi, 0, false)
		end
	end

	pa:save(path..name..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


-- find closest triangle vertex ------------------------------------------------------
local function find_closest_triangle(point, t_faces)
	dist = math.huge
	k = 0 
	for i=1,#t_faces do
		for j=2,4 do
			local x = t_faces[i][j][1]
			local y = t_faces[i][j][2]
			local z = t_faces[i][j][3]
			local d = (point[1]-x)^2 + (point[2]-y)^2 + (point[3]-z)^2
			if (d<dist) then
				-- if (d==0) then 
				-- 	return i 
				-- end
				dist = d
				k = i
			end
		end
	end
	return k
end


-- make array of vertices keys (as string) and triangle indices ---------------
local function make_vertices_string_array(t_faces)
	local vertices = {}
	for i=1,#t_faces do
		for j=2,4 do
			local x = t_faces[i][j][1]
			local y = t_faces[i][j][2]
			local z = t_faces[i][j][3]
			local s = string.format("%.6f,%.6f,%.6f", x,y,z)  -- round to 6 decimal digits
			table.insert(vertices, {s,i})
		end
	end
	return vertices
end


-- map vertices to triangles --------------------------------------------------
-- returns a table where the i-th element contains a table with all triangle 
-- indices connected to the i-th vertex
local function map_vertices2triangles(vertices)
	local t_hash = {}
	for i=1,#vertices do  -- loop through vertices
		if (not t_hash[vertices[i][1]]) then  -- vertex does not exist yet
			t_hash[vertices[i][1]] = {}  -- generate it
		end
		table.insert(t_hash[vertices[i][1]], vertices[i][2])  -- add triangle index to corresponding vertex
	end
	return t_hash 
end


-- make t_neighbour[i] --------------------------------------------------------
-- indices of triangles that are neighbours of i-th triangle
local function get_neighbours(t_faces)
	local vertices = make_vertices_string_array(t_faces)
	local map = map_vertices2triangles(vertices)
	local t_neighbour = {}
	for i=1,#t_faces do
		t_neighbour[i] = {}
		local lookup = {}
		for j=1,3 do
			local test_vertex = vertices[3*(i-1)+j][1]  -- loop through all vertices
			for k=1,#map[test_vertex] do  -- loop through all triangles that belong to one vertex
				local tmp = map[test_vertex][k]  -- triangle index
				if i~=tmp and lookup[tmp]==nil then  -- if triangle index ~= i, it is a neighbour triangle. But only add it once.
					table.insert(t_neighbour[i], tmp)
					lookup[tmp] = true
				end
			end
		end
	end
	return t_neighbour
end


-- get list of all triangles connected to the i-th triangle -------------------
local function get_island(t_faces, i)
	t_neighbour = get_neighbours(t_faces)
	isAdded = {}
	for i=1,#t_faces do
		isAdded[i] = false
	end
	isAdded[i] = true
	local island = t_neighbour[i]  -- initialize with all neighbours of i-th triangle
	for k,v in pairs(island) do
		isAdded[v] = true
	end	
	local j = 1
	local n = #island
	while j <= n do  -- loop through all triangles in island 
		for k,v in pairs(t_neighbour[island[j]]) do  -- loop through all neighbours of island[j]-th triangle
			if not isAdded[v] then  -- if v-th trianlgle has not been added yet
				island[n+1] = t_neighbour[island[j]][k]  -- append to island
				n = n + 1
				isAdded[v] = true
			end
		end
		j = j + 1
	end
	return island
end


-- STL to PA conversion -------------------------------------------------------
-- conversion of all connected triangle that are closest to point
function STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, 
	surface, symmetry, electrodes)
	local surfenhance = surface or "none"  -- surface enhancement
	local e = electrodes or {0, {0,0,0}}
	local path,name = splitPath(stl_filename)
	local start_time, end_time1, end_time2 = 0, 0, 0

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	local pa = simion.pas:open()
	pa:size((xmax-xmin)/dx_mm+1, (ymax-ymin)/dy_mm+1, (zmax-zmin)/dz_mm+1)
	pa.dx_mm = dx_mm
	pa.dy_mm = dy_mm
	pa.dz_mm = dz_mm
	pa.symmetry = symmetry or "3dplanar"

	local t_faces = readSTL(stl_filename)
	
	for i=1,#e do
		io.write(i.."/"..#e.." Finding connected STL surfaces... ")
		io.flush()
		start_time = os.clock()
		local ii = find_closest_triangle(e[i][2], t_faces)  -- index of closest triangle
		local island = get_island(t_faces, ii)  -- indices of triangles belonging to group of ii-th triangle
		local t_faces2 = {}
		for j,v in pairs(island) do
			table.insert(t_faces2, t_faces[v])
		end
		local t_size = stlSize(t_faces2)
		end_time1 = os.clock()
		io.write(" ["..string.format("%.3f", end_time1-start_time).." s]\n")
		local t_hash = map2Grid(t_faces2, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
		io.write(i.."/"..#e.." Filling PA points... ")
		io.flush()
		local xminstl = math.max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
		local xmaxstl = math.min(math.floor(t_size[2]/dx_mm)*dx_mm, xmax)
		local yminstl = math.max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
		local ymaxstl = math.min(math.floor(t_size[4]/dy_mm)*dy_mm, ymax)
		local zminstl = t_size[5]
		local zmaxstl = t_size[6]
		pa:fill { 
			function(x,y,z)
				if isInsideSTL(x+xmin,y+ymin,z+zmin, t_faces2, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl, t_hash, dx_mm, dy_mm) then
					return e[i][1], true
				else
					if i==1 then  -- do not overwrite other electrodes
						return 0, false
					end
				end
			end, 
			surface=surfenhance
		}
		end_time2 = os.clock()
		io.write(" ["..string.format("%.3f", end_time2-end_time1).." s]\n")
	end

	pa:save(path..name..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


-- modify (or add) one or more electrodes of an existing PA# array ---------------------
function STL2PA.modify(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, surface, electrodes)
	local surfenhance = surface or "none"  -- surface enhancement
	local e = electrodes or {0, {0,0,0}}
	local path,name = splitPath(stl_filename)
	local start_time, end_time1, end_time2 = 0, 0, 0

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	local pa = simion.pas:open(path..name..".pa#")  -- open pa# file
	local dx_mm = pa.dx_mm
	local dy_mm = pa.dy_mm
	local dz_mm = pa.dz_mm

	local t_faces = readSTL(stl_filename)

	for i=1,#e do
		-- delete electrode
		for xi,yi,zi in pa:points() do
			if pa:potential(xi,yi,zi) == e[i][1] then
				pa:point(xi,yi,zi, 0, false)
			end
		end
		local ii = find_closest_triangle(e[i][2], t_faces)  -- index of closest triangle
		local island_list = get_island(t_faces, ii)  -- indices of triangles belonging to group of ii-th triangle
		local t_faces2 = {}
		for j,v in pairs(island_list[1]) do
			if t_faces[v][1][3]~=0 then  -- ignore verticle faces
				table.insert(t_faces2, t_faces[v])
			end
		end

		local t_size = stlSize(t_faces2)
		io.write(i.."/"..#e.." Generating STL hash table... ")
		io.flush()
		start_time = os.clock()
		local t_hash = map2Grid(t_faces2, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
		end_time1 = os.clock()
		io.write(" [Finished in "..string.format("%.3f", end_time1-start_time).." s]\n")
		io.write(i.."/"..#e.." Building PA... ")
		io.flush()
		local xminstl = math.max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
		local xmaxstl = math.min(math.floor(t_size[2]/dx_mm)*dx_mm, xmax)
		local yminstl = math.max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
		local ymaxstl = math.min(math.floor(t_size[4]/dy_mm)*dy_mm, ymax)
		local zminstl = t_size[5]
		local zmaxstl = t_size[6]
		pa:fill { 
			function(x,y,z)
				if isInsideSTL(x+xmin,y+ymin,z+zmin, t_faces2, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl, t_hash, dx_mm, dy_mm) then
					return e[i][1], true
				end
			end, 
			surface=surfenhance
		}
		end_time2 = os.clock()
		io.write(" [Finished in "..string.format("%.3f", end_time2-end_time1).." s]\n")
	end

	pa:save(path..name..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


-- add ideal grid in .pa# file -----------------------------------
function STL2PA.addIdealGrid(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, voltage)

	local path,name = splitPath(stl_filename)
	local file = path..name..".pa#"
	local v = voltage or 0

	io.write("Adding ideal grid...\n")
	io.flush()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(file)  -- open pa# file

	pa:fill { 
		function(x,y,z)
			if ((x>=xmin and x<=xmax) and (y>=ymin and y<=ymax) and (z>=zmin and z<=zmax)) then
				return v, true
			end
		end, surface='none' 
	}

	pa:save(file)
	simion.pas:close()  -- remove all PAs from RAM.
end


-- add real grid in .pa# file -----------------------------------
function STL2PA.addRealGrid(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax,
	voltage, pitch, radius, dir_wire, dir_plane)

	local path,name = splitPath(stl_filename)
	local file = path..name..".pa#"
	local v = voltage or 0  -- electrode voltage
	local p = pitch or 0  -- grid pitch (mm)
	local r = radius or 0  -- wire radius (mm)
	local dw = dir_wire or "z"  -- wire direction, "x" or "y" or "z" 
	local dp = dir_plane or "y"  -- grid direction, orthogonal to wire direction

	io.write("Adding real grid...\n")
	io.flush()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(file)  -- open pa# file

	if dw == "x" and dp == "y" then
		if zmin~=zmax then
			print("WARNING: zmin~=zmax. zmin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((x>=xmin and x<=xmax) and (y>=ymin and y<=ymax)) then
					if (math.sqrt((z-zmin)^2 + (y%p)^2) < r or math.sqrt((z-zmin)^2 + (-y%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	elseif dw == "y" and dp == "x" then
		if zmin~=zmax then
			print("WARNING: zmin~=zmax. zmin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((x>=xmin and x<=xmax) and (y>=ymin and y<=ymax)) then
					if (math.sqrt((z-zmin)^2 + (x%p)^2) < r or math.sqrt((z-zmin)^2 + (-x%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	elseif dw == "z" and dp == "y" then
		if xmin~=xmax then
			print("WARNING: xmin~=xmax. xmin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((y>=ymin and y<=ymax) and (z>=zmin and z<=zmax)) then
					if (math.sqrt((x-xmin)^2 + (y%p)^2) < r or math.sqrt((x-xmin)^2 + (-y%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	elseif dw == "y" and dp == "z" then
		if xmin~=xmax then
			print("WARNING: xmin~=xmax. xmin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((y>=ymin and y<=ymax) and (z>=zmin and z<=zmax)) then
					if (math.sqrt((x-xmin)^2 + (z%p)^2) < r or math.sqrt((x-xmin)^2 + (-z%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	elseif dw == "x" and dp == "z" then
		if ymin~=ymax then
			print("WARNING: ymin~=ymax. ymin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((x>=xmin and x<=xmax) and (z>=zmin and z<=zmax)) then
					if (math.sqrt((y-ymin)^2 + (z%p)^2) < r or math.sqrt((y-ymin)^2 + (-z%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	elseif dw == "z" and dp == "x" then
		if ymin~=ymax then
			print("WARNING: ymin~=ymax. ymin taken for grid plane.")
		end 
		pa:fill { 
			function(x,y,z)
				if ((x>=xmin and x<=xmax) and (z>=zmin and z<=zmax)) then
					if (math.sqrt((y-ymin)^2 + (x%p)^2) < r or math.sqrt((y-ymin)^2 + (-x%p)^2) < r) then
						return v, true
					end
				end
			end, surface='fractional' 
		}
	end

	pa:save(file)
	simion.pas:close()  -- remove all PAs from RAM.
end


return STL2PA
