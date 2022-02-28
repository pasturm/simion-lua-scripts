--[[
stl2palib.lua

Library to convert an electrode geometry from STL files to a SIMION potential array file.

The main advantages compared to the stl2pa command from SIMION SL Tools (< 8.2.0.11) are:
1. The STL surface is filled with solid points (avoids the solid points 
   strategy issues of SL Tools).
2. The electrode surface enhancement feature can be used (leads to higher 
   field accuracies).
3. Array mirroring can be used.
4. Bounding boxes can easily be added.

Multi-electrode systems are handled as in SL Tools with the "%" character
as a placeholder in the STL file name.

The input STL needs to be watertight and must not contain joined faces.

Note: SIMION 8.2.0.11-20210626 SL Tools now contains similar functionality. But there are still bugs
and it is slower than this lua script.

Example:
	-- load STL2PA library
	local STL2PA = simion.import("./stl2palib.lua")

	-- STL input file
	local stl_filename = "my_geometry-%.stl"

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

	-- run stl2pa conversion
	STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, surface, boundingbox, symmetry)

	-- -- modify or add one electrode in my_geometry.pa#
	-- local electrode_no = 2  -- my_geometry-2.stl
	-- STL2PA.modify(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, electrode_no, surface, boundingbox)


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
function STL2PA.readSTL(stl_filename)
	local t_faces = {}
	local current_face = {}
	local z0
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
function STL2PA.stlSize(t_faces)
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
	-- if one rectangle is on left side of other 
	if (xmin1 > xmax2 or xmin2 > xmax1) then
	    return false
	-- if one rectangle is above other 
	elseif (ymin1 > ymax2 or ymin2 > ymax1) then
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
	local x_max = math.min(math.ceil(t_size[2]/dx_mm)*dx_mm, xmax)
	local y_min = math.max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
	local y_max = math.min(math.ceil(t_size[4]/dy_mm)*dy_mm, ymax)
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


-- Get the unique elements of a table -----------------------------------------
-- https://stackoverflow.com/questions/20066835/lua-remove-duplicate-elements
local function getUniqueElements(t_in)
	local hash = {}
	local res = {}
	for _,v in ipairs(t_in) do
	   if (not hash[v]) then
	       res[#res+1] = v
	       hash[v] = true
	   end
	end
	return res
end

-- check whether a grid point is inside the STL -------------------------------
-- from http://geomalgorithms.com/a06-_intersect-2.html
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
		if (r>=0) then  -- ray goes towards triangle
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
			if (s>=0 and s<=1 and t>=0 and (s+t)<=1) then
				if (r==0) then
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


-- Find all STL files
-- Note: To convert multiple overlapping STL files into single PA file, have your CAD software create one STL file for each electrode. 
-- Each file must be numbered by the potential it will be assigned in the PA file (e.g. myfile-1.stl, myfile-2.stl, and myfile-3.stl). 
-- Then set stl_filename to a template for the input file name using the "%" character as a placeholder (e.g. "my-file-%.stl"). 
-- The SL Tools will search and read all matching STL files and assign these voltages of 1V, 2V, and 3V respectively.
-- Adapted from https://stackoverflow.com/questions/5303174/how-to-get-list-of-directories-in-lua
local function getSTLfiles(stl_filename)
	local path,name,ext = splitPath(stl_filename)
	local basename =  string.match(stl_filename, "^.+/(.+)%-%%")  -- filename up to %
    local i, t, popen = 0, {}, io.popen
    local pfile = popen('dir "'..path..'" /b')
    for filename in pfile:lines() do
    	if string.match(filename, basename.."%-") and (string.match(filename, ".stl") or string.match(filename, ".ast")) then
	        i = i + 1
	        t[i] = filename
	    end
    end
    pfile:close()
    return t
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


-- add bounding box -----------------------------------------------------------
local function addBoundingBox(bounding, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm)
	if (string.match(bounding, "%-x")) then
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(0,yi,zi, 0,true)
			end
		end
	end
	if (string.match(bounding, "%+x")) then
		local ximax = (xmax-xmin)/dx_mm
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
	if (string.match(bounding, "%+y")) then
		local yimax = (ymax-ymin)/dy_mm
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
	if (string.match(bounding, "%+z")) then
		local zimax = (zmax-zmin)/dz_mm
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,zimax, 0,true)
			end
		end
	end
end

-- STL to PA conversion -------------------------------------------------------
function STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, surface, boundingbox, symmetry)
	local surfenhance = surface or "none"  -- surface enhancement
	local bounding = boundingbox or ""  -- grounded bounding box
	local path,name = splitPath(stl_filename)
	local start_time, end_time1, end_time2 = 0, 0, 0

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open()
	pa:size((xmax-xmin)/dx_mm+1, (ymax-ymin)/dy_mm+1, (zmax-zmin)/dz_mm+1)
	pa.dx_mm = dx_mm
	pa.dy_mm = dy_mm
	pa.dz_mm = dz_mm
	pa.symmetry = symmetry or "3dplanar"

	if string.match(stl_filename, "%-%%%.") then  -- check if multiple STL files have to be converted
		-- find all STL files
		local files = getSTLfiles(stl_filename)
		for i,v in ipairs(files) do
			
			local potential = tonumber(string.match(v, "^.+%-(%d+)"))
			local t_faces = STL2PA.readSTL(path..v)
			local t_size = STL2PA.stlSize(t_faces)
			io.write(i.."/"..#files.." Generating STL hash table for "..v.."... ")
			io.flush()
			start_time = os.clock()
			local t_hash = map2Grid(t_faces, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
			end_time1 = os.clock()
			io.write(" [Finished in "..string.format("%.3f", end_time1-start_time).." s]\n")
			io.write(i.."/"..#files.." Building PA... ")
			io.flush()
			local xminstl = max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
			local xmaxstl = min(math.ceil(t_size[2]/dx_mm)*dx_mm, xmax)
			local yminstl = max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
			local ymaxstl = min(math.ceil(t_size[4]/dy_mm)*dy_mm, ymax)
			local zminstl = t_size[5]
			local zmaxstl = t_size[6]
			pa:fill { 
				function(x,y,z)  -- in grid units
					if isInsideSTL(x+xmin,y+ymin,z+zmin, t_faces, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl, t_hash, dx_mm, dy_mm) then
						return potential, true
					else
						if i==1 then  -- do not overwrite other electrodes
							return 0, false
						end
					end
				end,
				surface=surfenhance
			}
			end_time2 = os.clock()
			io.write(" [Finished in "..string.format("%.3f", end_time2-end_time1).." s]\n")
		end
	else
		local t_faces = STL2PA.readSTL(stl_filename)
		local t_size = STL2PA.stlSize(t_faces)
		io.write("Generating STL hash table for "..stl_filename.."... ")
		io.flush()
		start_time = os.clock()
		local t_hash = map2Grid(t_faces, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
		end_time1 = os.clock()
		io.write(" [Finished in "..string.format("%.3f", end_time1-start_time).." s]\n")
		io.write("Building PA...\n")
		io.flush()
		local xminstl = math.max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
		local xmaxstl = math.min(math.ceil(t_size[2]/dx_mm)*dx_mm, xmax)
		local yminstl = math.max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
		local ymaxstl = math.min(math.ceil(t_size[4]/dy_mm)*dy_mm, ymax)
		local zminstl = t_size[5]
		local zmaxstl = t_size[6]
		pa:fill { 
			function(x,y,z)
				if isInsideSTL(x+xmin,y+ymin,z+zmin, t_faces, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl, t_hash, dx_mm, dy_mm) then
					return 1, true
				else
					return 0, false
				end
			end, 
			surface=surfenhance
		}
		io.write("\n")
	end

	addBoundingBox(bounding, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm)

	pa:save(path..string.gsub(name, "%-%%", "")..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


-- modify (or add) one electrode of an existing PA# array ---------------------
function STL2PA.modify(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, electrode_no, surface, boundingbox)
	local surfenhance = surface or "none"  -- surface enhancement
	local bounding = boundingbox or ""  -- grounded bounding box
	local path,name = splitPath(stl_filename)

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(path..string.gsub(name, "%-%%", "")..".pa#")  -- open pa# file

	-- read new STL file
	local files = getSTLfiles(stl_filename)
	local t_faces = {}
	for i,v in ipairs(files) do
		local potential = tonumber(string.match(v, "^.+%-(%d+)"))
		if (potential==electrode_no) then
			t_faces = STL2PA.readSTL(path..v)
			break
		end
	end
	if next(t_faces)==nil then 
		print("Error: no STL file with electrode_no "..electrode_no.." found.")
		return
	end
	local t_size = STL2PA.stlSize(t_faces)

	-- delete electrode
	for xi,yi,zi in pa:points() do
		if pa:potential(xi,yi,zi) == electrode_no then
			pa:point(xi,yi,zi, 0, false)
		end
	end

	local dx_mm = pa.dx_mm
	local dy_mm = pa.dy_mm
	local dz_mm = pa.dz_mm
	io.write("Generating STL hash table for "..path..string.gsub(name, "%-%%", "-"..electrode_no)..".stl".."... ")
	io.flush()
	local t_hash = map2Grid(t_faces, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
	start_time = os.clock()
	local t_hash = map2Grid(t_faces, t_size, dx_mm, dy_mm, xmin, xmax, ymin, ymax)
	end_time1 = os.clock()
	io.write(" [Finished in "..string.format("%.3f", end_time1-start_time).." s]\n")
	io.write("Building PA...\n")
	io.flush()
	local xminstl = max(math.floor(t_size[1]/dx_mm)*dx_mm, xmin)
	local xmaxstl = min(math.ceil(t_size[2]/dx_mm)*dx_mm, xmax)
	local yminstl = max(math.floor(t_size[3]/dy_mm)*dy_mm, ymin)
	local ymaxstl = min(math.ceil(t_size[4]/dy_mm)*dy_mm, ymax)
	local zminstl = t_size[5]
	local zmaxstl = t_size[6]
	pa:fill { 
		function(x,y,z)
			if isInsideSTL(x+xmin,y+ymin,z+zmin, t_faces, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl, t_hash, dx_mm, dy_mm) then
				return electrode_no, true
			end
		end, 
		surface=surfenhance
	}
	io.write("\n")

	if (electrode_no==0) then
		addBoundingBox(bounding, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm)
	end

	pa:save(path..string.gsub(name, "%-%%", "")..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


-- add or modifiy bounding box in .pa# file -----------------------------------------------------------
function STL2PA.boundingBox(stl_filename, bounding, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm)
	local path,name = splitPath(stl_filename)
	local file = path..string.gsub(name, "%-%%", "")..".pa#"

	writeHeader()
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
	else
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
				pa:point(ximax,yi,zi, 0,true)
			end
		end
	else
		for yi=0,(ymax-ymin)/dy_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(ximax,yi,zi, 0,false)
			end
		end
	end

	if (string.match(bounding, "%-y")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,0,zi, 0,true)
			end
		end
	else
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
				pa:point(xi,yimax,zi, 0,true)
			end
		end
	else
		for xi=0,(xmax-xmin)/dx_mm do
			for zi=0,(zmax-zmin)/dz_mm do
				pa:point(xi,yimax,zi, 0,false)
			end
		end
	end

	if (string.match(bounding, "%-z")) then
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,0, 0,true)
			end
		end
	else
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
				pa:point(xi,yi,zimax, 0,true)
			end
		end
	else
		for xi=0,(xmax-xmin)/dx_mm do
			for yi=0,(ymax-ymin)/dy_mm do
				pa:point(xi,yi,zimax, 0,false)
			end
		end
	end

	pa:save(file)
	simion.pas:close()  -- remove all PAs from RAM.
end


-- remove one electrode of an existing PA# array ---------------------
function STL2PA.removeElectrode(stl_filename, electrode_no)
	local path,name = splitPath(stl_filename)

	writeHeader()

	simion.pas:close()  -- remove all PAs from RAM.
	pa = simion.pas:open(path..string.gsub(name, "%-%%", "")..".pa#")  -- open pa# file

	-- delete electrode
	for xi,yi,zi in pa:points() do
		if pa:potential(xi,yi,zi) == electrode_no then
			pa:point(xi,yi,zi, 0, false)
		end
	end

	pa:save(path..string.gsub(name, "%-%%", "")..".pa#")
	simion.pas:close()  -- remove all PAs from RAM.
end


return STL2PA
