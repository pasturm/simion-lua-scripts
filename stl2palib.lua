--[[
stl2palib.lua

STL to SIMION potential array conversion library.

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
and it is much slower than this lua script.

Example:
  -- load STL2PA library
  local STL2PA = simion.import("./stl2palib.lua")

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

  -- array symmetry and mirroring
  local symmetry = "3dplanar[yz]"

  -- assign voltages to fully connected surfaces closest to the given point
  -- {voltage, {x, y, z}}
  local electrodes = {
    {1, {0,-75,-70}},
    {2, {0,-75,-66}}
  }

  -- run stl2pa conversion
  STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, 
            surface, symmetry, electrodes, "z")

  -- add electrode to existing .pa#
  local add = true
  STL2PA.convert(stl_filename, xmin, xmax, ymin, ymax, zmin, zmax, dx_mm, dy_mm, dz_mm, 
            surface, symmetry, electrodes, "z", add)

  -- add bounding box in .pa# file
  -- local boundingbox = "-x+x-y+y-z+z"
  -- local add = true  -- false to remove boundary
  -- STL2PA.boundingBox(stl_filename, boundingbox, add)

  -- remove an electrode from .pa#
  -- STL2PA.removeElectrode(stl_filename, electrode_no)

Author: Patrick Sturm
Copyright (c) 2020-2026 TOFWERK
--]]
 
local STL2PA = {}

-- Split path/filename into path, filename and extension ----------------------
local function splitPath(path)
  return string.match(path, "(.-)([^\\/]-)%.?([^%.\\/]*)$")
end

-- Convert 4-byte binary string to unsigned 32-bit integer ------------------------
local function read_uint32(str, offset)
  local i1, i2, i3, i4 = string.byte(str, offset, offset+3)
  return ((i4*256 + i3)*256 + i2)*256 + i1
end

-- Convert 4-byte binary string to little-endian IEEE-754 single-precision float ----
-- b4 = [sign:1][exp high 7] | b3 = [exp low 1][mant high 7] | b2,b1 = mantissa
local function read_float32(str, offset)
  offset = offset or 1
  local b1, b2, b3, b4 = string.byte(str, offset, offset + 3)
  local sign = 1
  if b4 > 127 then sign = -1 end
  local mantissa = ((b3%128)*256 + b2)*256 + b1
  local exponent = (b4%128) * 2 + math.floor(b3/128)
  if exponent == 0 then return 0 end
  return sign * (1+mantissa/2^23) * 2^(exponent-127)
end

-- Read a STL file ------------------------------------------------------------
local function readSTL(filename)
  local fh = assert(io.open(filename,"rb"))
  local s = assert(fh:read("*a"))
  fh:close()
 
  local triangles = {}
  local binary = #s >= 84 and #s == 84+50*read_uint32(s,81)
 
  if binary then
    local n = read_uint32(s,81)
    local p = 85
    for _=1,n do
      triangles[#triangles+1] = {
        read_float32(s,p+12),read_float32(s,p+16),read_float32(s,p+20),
        read_float32(s,p+24),read_float32(s,p+28),read_float32(s,p+32),
        read_float32(s,p+36),read_float32(s,p+40),read_float32(s,p+44)
      }
      p = p+50
    end
  else
    local v = {}
    for a,b,c in s:gmatch("[vV]ertex%s+(%S+)%s+(%S+)%s+(%S+)") do
      v[#v+1] = tonumber(a)
      v[#v+1] = tonumber(b) 
      v[#v+1] = tonumber(c)
    end
    for i=0,#v-9,9 do
      triangles[#triangles+1] = {v[i+1],v[i+2],v[i+3],v[i+4],v[i+5],v[i+6],v[i+7],v[i+8],v[i+9]}
    end
  end
  return triangles
end
 
local function qround(x, scale)
  if x >= 0 then return math.floor(x*scale+0.5) end
  return -math.floor(-x*scale+0.5)
end
 
-- Build vertex-connected triangle components using a union-find data structure.
local function build_components(triangles)
  local vertex_scale = 1e6
  local keyid,parent,rank = {},{},{}
  local vx,vy,vz = {},{},{}
  local ta = {}
 
  local function vertex(x,y,z)
    local key = qround(x,vertex_scale)..","..qround(y,vertex_scale)..","..qround(z,vertex_scale)
    local id = keyid[key]
    if not id then
      id = #parent+1
      keyid[key] = id
      parent[id] = id
      rank[id] = 0
      vx[id],vy[id],vz[id] = x,y,z
    end
    return id
  end
 
  -- find the root of the tree and compress the path
  local function find(a)
    local root = a
    while parent[root] ~= root do
      root = parent[root]
    end
    while parent[a] ~= root do
      local next_a = parent[a]
      parent[a] = root
      a = next_a
    end
    return root
  end
 
  -- connect two vertices
  local function union(a,b)
    a,b = find(a),find(b)
    if a == b then return end
    if rank[a] < rank[b] then a,b = b,a end
    parent[b] = a
    if rank[a] == rank[b] then rank[a] = rank[a]+1 end
  end
 
  for i=1,#triangles do
    local t = triangles[i]
    local a = vertex(t[1],t[2],t[3])
    local b = vertex(t[4],t[5],t[6])
    local c = vertex(t[7],t[8],t[9])
    union(a,b)
    union(a,c)
    ta[i] = a
  end
 
  local rootcid = {}
  local C = {}
 
  for i=1,#triangles do
    local root = find(ta[i])
    local cid = rootcid[root]
 
    if not cid then
      cid = #C+1
      rootcid[root] = cid
      C[cid] = {
        triangle = {},
        xmin = math.huge,xmax = -math.huge,
        ymin = math.huge,ymax = -math.huge,
        zmin = math.huge,zmax = -math.huge
      }
    end
 
    local c = C[cid]
    local t = triangles[i]
    c.triangle[#c.triangle+1] = i
 
    for p=1,7,3 do
      local x,y,z = t[p],t[p+1],t[p+2]
      c.xmin = math.min(c.xmin,x)
      c.xmax = math.max(c.xmax,x)
      c.ymin = math.min(c.ymin,y)
      c.ymax = math.max(c.ymax,y)
      c.zmin = math.min(c.zmin,z)
      c.zmax = math.max(c.zmax,z)
    end
  end
 
  local vc = {}
  for i=1,#vx do vc[i] = rootcid[find(i)] end
 
  return C,vx,vy,vz,vc
end
 
-- Map each electrode seed to the component containing its nearest stored vertex.
-- If several seeds select one component, the last seed retains precedence.
local function mapSeeds(E,vx,vy,vz,vc)
  local mapped,last = {},{}
 
  for e=1,#E do
    local p = E[e][2]
    local best,bi = math.huge,1
 
    for i=1,#vx do
      local x = p[1]-vx[i]
      local y = p[2]-vy[i]
      local z = p[3]-vz[i]
      local d = x*x+y*y+z*z
 
      if d < best then
        best,bi = d,i
      end
    end
 
    mapped[e] = vc[bi]
    last[mapped[e]] = e
  end
 
  local S = {}
 
  for e=1,#E do
    local cid = mapped[e]
 
    if last[cid] == e then
      S[#S+1] = {cid=cid, V=E[e][1], seed=e}
    end
  end
 
  return S,#E-#S
end
 
-- Precompute the spatial lookup for one component.
--
-- Ordinary nonparallel triangles:
--   R = projected A/B/C plus affine wi(u,v)
--   B = projected cell -> R record IDs
-- Projected-degenerate/ray-parallel triangles:
--   D  = complete 3D triangle coordinates
--   DB = projected cell -> D record IDs
local function preprocess(c,triangles,iu,iv,iw,u0,v0,du,dv,nu,nv)
  local R,B = {},{}
  local D,DB = {},{}
  local ulast = u0+(nu-1)*du
  local vlast = v0+(nv-1)*dv
 
  local function clip(value,lo,hi)
    if value < lo then return lo end
    if value > hi then return hi end
    return value
  end
 
  local function insert(cells,id,umin,umax,vmin,vmax)
    if umax < u0 or umin > ulast or
       vmax < v0 or vmin > vlast
    then
      return
    end
    local i0 = clip(math.floor((umin-u0)/du),0,nu-1)
    local i1 = clip(math.floor((umax-u0)/du),0,nu-1)
    local j0 = clip(math.floor((vmin-v0)/dv),0,nv-1)
    local j1 = clip(math.floor((vmax-v0)/dv),0,nv-1)
    for i=i0,i1 do
      local row = i*nv
      for j=j0,j1 do
        local key = row+j+1
        local cell = cells[key]
        if not cell then
          cell = {}
          cells[key] = cell
        end
        cell[#cell+1] = id
      end
    end
  end
 
  for _,ti in ipairs(c.triangle) do
    local t = triangles[ti]
    local au,av,aw = t[iu],t[iv],t[iw]
    local bu,bv,bw = t[iu+3],t[iv+3],t[iw+3]
    local cu,cv,cw = t[iu+6],t[iv+6],t[iw+6]
    local d = (bu-au)*(cv-av) - (bv-av)*(cu-au)
 
    if d ~= 0 then
      -- Normalize the projected triangle to counter-clockwise winding.
      if d < 0 then
        bu,cu = cu,bu
        bv,cv = cv,bv
        bw,cw = cw,bw
        d = -d
      end
 
      local a1 = (bv-cv)/d
      local b1 = (cu-bu)/d
      local c1 = -a1*cu-b1*cv
      local a2 = (cv-av)/d
      local b2 = (au-cu)/d
      local c2 = -a2*cu-b2*cv
      local acw = aw-cw
      local bcw = bw-cw
      local wu = a1*acw+a2*bcw
      local wv = b1*acw+b2*bcw
      local wc = cw+c1*acw+c2*bcw
 
      local base = #R
      local rid = base/9+1
 
      R[base+1],R[base+2] = au,av
      R[base+3],R[base+4] = bu,bv
      R[base+5],R[base+6] = cu,cv
      R[base+7],R[base+8],R[base+9] = wu,wv,wc
 
      insert(B, rid, math.min(au,bu,cu), math.max(au,bu,cu), math.min(av,bv,cv), math.max(av,bv,cv))
    else
      local base = #D
      local did = base/9+1
 
      D[base+1],D[base+2],D[base+3] = au,av,aw
      D[base+4],D[base+5],D[base+6] = bu,bv,bw
      D[base+7],D[base+8],D[base+9] = cu,cv,cw
 
      insert(DB,did, math.min(au,bu,cu), math.max(au,bu,cu), math.min(av,bv,cv), math.max(av,bv,cv))
    end
  end
 
  c.R,c.B = R,B
  c.D,c.DB = D,DB
end
 
-- Inclusive 3D surface test for a triangle whose u-v projection is degenerate.
local function on_parallel_triangle(u,v,w,D,o,eps)
  local au,av,aw = D[o+1],D[o+2],D[o+3]
  local bu,bv,bw = D[o+4],D[o+5],D[o+6]
  local cu,cv,cw = D[o+7],D[o+8],D[o+9]
  local ux,uy,uz = bu-au,bv-av,bw-aw
  local vx,vy,vz = cu-au,cv-av,cw-aw
  local px,py,pz = u-au,v-av,w-aw
  local nx = uy*vz-uz*vy
  local ny = uz*vx-ux*vz
  local nz = ux*vy-uy*vx
  local n2 = nx*nx+ny*ny+nz*nz
  if n2 == 0 then return false end
  local plane = px*nx+py*ny+pz*nz
  if plane*plane > eps*eps*n2 then return false end
  local uu = ux*ux+uy*uy+uz*uz
  local uv = ux*vx+uy*vy+uz*vz
  local vv = vx*vx+vy*vy+vz*vz
  local wu = px*ux+py*uy+pz*uz
  local wv = px*vx+py*vy+pz*vz
  local den = uu*vv-uv*uv
  if den <= 0 then return false end
  local s = (vv*wu-uv*wv)/den
  local t = (uu*wv-uv*wu)/den
  local beps = eps*math.sqrt(math.max(uu,vv)/den)
  return s >= -beps and t >= -beps and s+t <= 1+beps
end
 
function STL2PA.convert(filename, xmin, xmax, ymin, ymax, zmin, zmax,
  dx, dy, dz, surface, symmetry, E, direction, add)
  E = E or {{0,{0,0,0}}}
  direction = direction or "z"
  surface = surface or "none"
  add = add or false
 
  assert(direction == "x" or direction == "y" or direction == "z")
 
  local path,name = splitPath(filename)
  simion.pas:close()
 
  local pa
  local nx = math.floor((xmax-xmin)/dx+0.5)+1
  local ny = math.floor((ymax-ymin)/dy+0.5)+1
  local nz = math.floor((zmax-zmin)/dz+0.5)+1
 
  if add then
    pa = simion.pas:open(path..name..".pa#")
    dx,dy,dz = pa.dx_mm,pa.dy_mm,pa.dz_mm
    nx,ny,nz = pa.nx,pa.ny,pa.nz
  else
    pa = simion.pas:open()
    pa:size(nx,ny,nz)
    pa.dx_mm,pa.dy_mm,pa.dz_mm = dx,dy,dz
    pa.symmetry = symmetry or "3dplanar"
  end
 
  local start = os.clock()
  io.write("Reading STL... ")
  io.flush()
 
  local triangles = readSTL(filename)
 
  io.write(#triangles.." triangles\nFinding connected components...\n")
  io.flush()
 
  local C,vx,vy,vz,vc = build_components(triangles)
  local S,reused = mapSeeds(E,vx,vy,vz,vc)
 
  io.write(#E.." electrode seed(s)\n")
  io.write(#S.." unique component(s) selected\n")
  io.write(reused.." duplicate component selection(s) ignored\n")
 
  local iu,iv,iw,u0,v0,du,dv,nu,nv
  if direction == "z" then
    iu,iv,iw,u0,v0,du,dv,nu,nv = 1,2,3,xmin,ymin,dx,dy,nx,ny
  elseif direction == "y" then
    iu,iv,iw,u0,v0,du,dv,nu,nv = 1,3,2,xmin,zmin,dx,dz,nx,nz
  else
    iu,iv,iw,u0,v0,du,dv,nu,nv = 3,2,1,zmin,ymin,dz,dy,nz,ny
  end
 
  local G = {}
  local ulast = u0+(nu-1)*du
  local vlast = v0+(nv-1)*dv
 
  for _,q in ipairs(S) do
    local c = C[q.cid]
    q.c = c
    if direction == "z" then
      c.umin,c.umax = c.xmin,c.xmax
      c.vmin,c.vmax = c.ymin,c.ymax
      c.wmin,c.wmax = c.zmin,c.zmax
    elseif direction == "y" then
      c.umin,c.umax = c.xmin,c.xmax
      c.vmin,c.vmax = c.zmin,c.zmax
      c.wmin,c.wmax = c.ymin,c.ymax
    else
      c.umin,c.umax = c.zmin,c.zmax
      c.vmin,c.vmax = c.ymin,c.ymax
      c.wmin,c.wmax = c.xmin,c.xmax
    end
 
    c.surface_eps = 1e-9 * math.max(c.umax-c.umin, c.vmax-c.vmin, c.wmax-c.wmin, du, dv, 1)
 
    preprocess(c,triangles,iu,iv,iw,u0,v0,du,dv,nu,nv)
 
    if not (c.umax < u0 or c.umin > ulast or c.vmax < v0 or c.vmin > vlast) then
      local i0 = math.max(0,math.floor((c.umin-u0)/du))
      local i1 = math.min(nu-1,math.floor((c.umax-u0)/du))
      local j0 = math.max(0,math.floor((c.vmin-v0)/dv))
      local j1 = math.min(nv-1,math.floor((c.vmax-v0)/dv))
 
      for i=i0,i1 do
        local row = i*nv
        for j=j0,j1 do
          local key = row+j+1
          local cell = G[key]
          if not cell then
            cell = {}
            G[key] = cell
          end
          cell[#cell+1] = q
        end
      end
    end
  end
 
  -- Determine whether an exact point (u,v,w) lies in a selected component.
  local function classify(u,v,w)
    if u < u0 or u > ulast or
       v < v0 or v > vlast
    then
      return nil
    end
 
    local i = math.floor((u-u0)/du)
    local j = math.floor((v-v0)/dv)
 
    -- Protect the exact global maximum from a rounded quotient of nu/nv.
    if i >= nu then i = nu-1 end
    if j >= nv then j = nv-1 end
 
    local key = i*nv+j+1
    local candidate = G[key]
    if not candidate then return nil end
 
    -- Reverse traversal preserves the existing last-selected precedence.
    for n=#candidate,1,-1 do
      local q = candidate[n]
      local c = q.c
      local eps = c.surface_eps
 
      if u >= c.umin and u <= c.umax and
         v >= c.vmin and v <= c.vmax and
         w >= c.wmin-eps and w <= c.wmax+eps
      then
        -- Projected-degenerate triangles are surface-only exceptions.
        local db = c.DB[key]
        if db then
          local D = c.D
          for m=1,#db do
            local o = (db[m]-1)*9
            if on_parallel_triangle(u,v,w,D,o,eps) then
              return q.V
            end
          end
        end
 
        local b = c.B[key]
        if b then
          local R = c.R
          local count = 0
          for m=1,#b do
            local o = (b[m]-1)*9
            local au,av = R[o+1],R[o+2]
            local bu,bv = R[o+3],R[o+4]
            local cu,cv = R[o+5],R[o+6]
 
            -- Three independent query-relative edge functions. Since the
            -- triangle is CCW, nonnegative values are inclusive coverage.
            local eA = (bu-u)*(cv-v) - (bv-v)*(cu-u) -- edge B->C
            local eB = (cu-u)*(av-v) - (cv-v)*(au-u) -- edge C->A
            local eC = (au-u)*(bv-v) - (av-v)*(bu-u) -- edge A->B
 
            if eA >= 0 and eB >= 0 and eC >= 0 then
              local wi = R[o+7]*u+R[o+8]*v+R[o+9]
              local distance = wi-w
 
              -- Surface membership is inclusive and independent of parity.
              if distance >= -eps and distance <= eps then
                return q.V
              end
 
              if distance > eps then
                if eA > 0 and eB > 0 and eC > 0 then
                  count = count+1
                else
                  -- Half-open top-left rule. Directed edge P->Q owns an exact
                  -- zero when dv>0, or when horizontal and du<0. At a vertex,
                  -- all incident zero-edge conditions for this triangle pass.
                  local ownA = eA > 0 or (eA == 0 and (cv > bv or (cv == bv and cu < bu)))
                  local ownB = eB > 0 or (eB == 0 and (av > cv or (av == cv and au < cu)))
                  local ownC = eC > 0 or (eC == 0 and (bv > av or (bv == av and bu < au)))
 
                  if ownA and ownB and ownC then
                    count = count+1
                  end
                end
              end
            end
          end
          if count%2 == 1 then return q.V end
        end
      end
    end
 
    return nil
  end
 
  io.write("Preparation: "..string.format("%.3f",os.clock()-start).." s\n")
  io.write("Filling PA... ")
  io.flush()
 
  local fill_start = os.clock()
 
  pa:fill{
    function(x,y,z)
      local xm,ym,zm = xmin+x,ymin+y,zmin+z
      local V
 
      if direction == "z" then
        V = classify(xm,ym,zm)
      elseif direction == "y" then
        V = classify(xm,zm,ym)
      else
        V = classify(zm,ym,xm)
      end
 
      if V ~= nil then return V,true end
      if not add then return 0,false end
    end,
    surface=surface
  }
 
  io.write("["..string.format("%.3f",os.clock()-fill_start).." s]\n")
  pa:save(path..name..".pa#")
  simion.pas:close()
end
 
-- Add or remove a bounding box in a PA# file. ---------------------------
function STL2PA.boundingBox(filename,bounding,add,V)
  local path,name = splitPath(filename)
  local file = path..name..".pa#"
 
  if add == nil then add = true end
  V = add and (V or 0) or 0
 
  simion.pas:close()
  local pa = simion.pas:open(file)
  local X,Y,Z = pa.nx-1,pa.ny-1,pa.nz-1
 
  if bounding:match("%-x") then
    for y=0,Y do for z=0,Z do pa:point(0,y,z,V,add) end end
  end
  if bounding:match("%+x") then
    for y=0,Y do for z=0,Z do pa:point(X,y,z,V,add) end end
  end
  if bounding:match("%-y") then
    for x=0,X do for z=0,Z do pa:point(x,0,z,V,add) end end
  end
  if bounding:match("%+y") then
    for x=0,X do for z=0,Z do pa:point(x,Y,z,V,add) end end
  end
  if bounding:match("%-z") then
    for x=0,X do for y=0,Y do pa:point(x,y,0,V,add) end end
  end
  if bounding:match("%+z") then
    for x=0,X do for y=0,Y do pa:point(x,y,Z,V,add) end end
  end
 
  pa:save(file)
  simion.pas:close()
end
 
-- Remove one electrode from an existing PA# array. ----------------------------
function STL2PA.removeElectrode(filename,electrode_no)
  local path,name = splitPath(filename)
  local file = path..name..".pa#"
 
  simion.pas:close()
  local pa = simion.pas:open(file)
 
  for xi,yi,zi in pa:points() do
    if pa:potential(xi,yi,zi) == electrode_no then
      pa:point(xi,yi,zi,0,false)
    end
  end
 
  pa:save(file)
  simion.pas:close()
end
 
-- Add an ideal zero-thickness grid region. --------------------------------
function STL2PA.addIdealGrid(filename,xmin,xmax,ymin,ymax,zmin,zmax,voltage)
  local path,name = splitPath(filename)
  local file = path..name..".pa#"
  local v = voltage or 0
 
  io.write("Adding ideal grid...\n")
  io.flush()
  simion.pas:close()
 
  local pa = simion.pas:open(file)
 
  pa:fill{
    function(x,y,z)
      if x >= xmin and x <= xmax and
         y >= ymin and y <= ymax and
         z >= zmin and z <= zmax
      then
        return v,true
      end
    end,
    surface="none"
  }
 
  pa:save(file)
  simion.pas:close()
end

-- add real grid in .pa# file -----------------------------------
local function wire_distance(value,pitch)
  return math.min(value%pitch,(-value)%pitch)
end
 
function STL2PA.addRealGrid(filename,xmin,xmax,ymin,ymax,zmin,zmax,
  voltage,pitch,radius,dir_wire,dir_plane)
  local path,name = splitPath(filename)
  local file = path..name..".pa#"
  local v = voltage or 0
  local p = pitch or 0
  local r = radius or 0
  local dw = dir_wire or "z"
  local dp = dir_plane or "y"
 
  assert(p > 0,"grid pitch must be > 0")
  assert(r >= 0,"grid radius must be >= 0")
 
  local valid = {xy=true, yx=true, zy=true, yz=true, xz=true, zx=true}
  assert(valid[dw..dp],"invalid wire/grid direction combination "..dw.."/"..dp)
 
  io.write("Adding real grid...\n")
  io.flush()
  simion.pas:close()
 
  local pa = simion.pas:open(file)
  local r2 = r*r
 
  pa:fill{
    function(x,y,z)
      local inside,dplane,dwire
 
      if dw == "x" and dp == "y" then
        inside = x >= xmin and x <= xmax and y >= ymin and y <= ymax
        dplane = z-zmin
        dwire = wire_distance(y,p)
      elseif dw == "y" and dp == "x" then
        inside = x >= xmin and x <= xmax and y >= ymin and y <= ymax
        dplane = z-zmin
        dwire = wire_distance(x,p)
      elseif dw == "z" and dp == "y" then
        inside = y >= ymin and y <= ymax and z >= zmin and z <= zmax
        dplane = x-xmin
        dwire = wire_distance(y,p)
      elseif dw == "y" and dp == "z" then
        inside = y >= ymin and y <= ymax and z >= zmin and z <= zmax
        dplane = x-xmin
        dwire = wire_distance(z,p)
      elseif dw == "x" and dp == "z" then
        inside = x >= xmin and x <= xmax and z >= zmin and z <= zmax
        dplane = y-ymin
        dwire = wire_distance(z,p)
      else
        inside = x >= xmin and x <= xmax and z >= zmin and z <= zmax
        dplane = y-ymin
        dwire = wire_distance(x,p)
      end
 
      if inside and dplane*dplane+dwire*dwire < r2 then
        return v,true
      end
    end,
    surface="fractional"
  }
 
  pa:save(file)
  simion.pas:close()
end
 
return STL2PA
