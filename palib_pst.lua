--[[
 Some utility functions for dealing with PA's.

 D.Manura, 2011-12-12, 2008-02.
 (c) 2011 Scientific Instrument Services, Inc. (Licensed under SIMION 8.1)

 adapted by
 Patrick Sturm
 (c) 2020 TOFWERK
--]]

local PAL = {}

--[[
 Gets the base file name (i.e. without extension) of the given PA.
--]]
function get_pa_basename(pa)  -- pst: make this function global
  local name = (pa.filename or 'blank.pa'):gsub('%.[^%.]+$', '')
  return name
end

--[[
 Creates a new empty PA of the size/symmetry required by the
 given coordinate definition.
--]]
function PAL.create_pa(coords)
  assert(coords.xmin, "missing xmin field")
  assert(coords.ymin, "missing ymin field")
  assert(coords.zmin, "missing zmin field")
  assert(coords.xmax, "missing xmax field")
  assert(coords.ymax, "missing ymax field")
  assert(coords.zmax, "missing zmax field")
  assert(coords.symmetry, "missing symmetry field")
  assert(coords.dx_mm, "missing dx_mm field")
  assert(coords.dy_mm, "missing dy_mm field")
  assert(coords.dz_mm, "missing dz_mm field")

  local dx = coords.xmax - coords.xmin
  local dy = coords.ymax - coords.ymin
  local dz = coords.zmax - coords.zmin
  
  local nx = math.floor(dx / coords.dx_mm + 1 + 0.5)
  local ny = math.floor(dy / coords.dy_mm + 1 + 0.5)
  local nz = math.floor(dz / coords.dz_mm + 1 + 0.5)

  local pa = simion.pas:open()
  pa:size(nx, ny, nz)
  pa.symmetry = coords.symmetry
  pa.dx_mm = coords.dx_mm
  pa.dy_mm = coords.dy_mm
  pa.dz_mm = coords.dz_mm
  
  return pa  
end

--[[
 Creates a new empty PA of the same size/symmetry as the given PA.
 If `off` is not `nil`, then the number grid points in each dimension
 of `destpa` will be `off` plus that of `pa` (typically `off` is -1, 0, or 1).
--]]
function PAL.create_similar_pa(pa, off)
  off = off or 0
  local expect_nx = pa.nx + off
  local expect_ny = pa.ny + off
  local expect_nz = pa.nz == 1 and 1 or pa.nz + off
  local newpa = simion.pas:open()
  newpa:size(expect_nx, expect_ny, expect_nz)
  newpa.symmetry = pa.symmetry
  newpa.dx_mm, newpa.dy_mm, newpa.dz_mm = pa.dx_mm, pa.dy_mm, pa.dz_mm
  newpa.potential_type = pa.potential_type
  if pa.potential_type == 'magnetic' then newpa.ng = pa.ng end
  newpa.refinable = pa.refinable
  return newpa
end

--[[
 Ensures that PA `destpa` has the same format (e.g. size, symmetry,
 potential type, etc.) as the given PA `pa'.  `pa` will be
 changed if necessary.
 If `off` has the same meaning as in PAL.create_similar_pa.
--]]
function PAL.make_similar_pa(destpa, pa, off)
  off = off or 0
  local expect_nx = pa.nx + off
  local expect_ny = pa.ny + off
  local expect_nz = pa.nz == 1 and 1 or pa.nz + off
  if destpa.nx~=expect_nx or destpa.ny~=expect_ny or destpa.nz~=expect_nz then
    destpa:size(expect_nx, expect_ny, expect_nz)
  end
  if destpa.symmetry ~= pa.symmetry then
    destpa.symmetry = pa.symmetry
  end
  if destpa.dx_mm~=pa.dx_mm or destpa.dy_mm~=pa.dy_mm or destpa.dz_mm~=pa.dz_mm then
    destpa.dx_mm, destpa.dy_mm, destpa.dz_mm = pa.dx_mm, pa.dy_mm, pa.dz_mm
  end
  if destpa.potential_type ~= pa.potential_type then
    destpa.potential_type = pa.potential_type
  end
  if pa.potential_type == 'magnetic' and destpa.ng ~= pa.ng then
    destpa.ng = pa.ng
  end
  if destpa.refinable ~= pa.refinable then
    destpa.refinable = pa.refinable
  end
  return destpa
end

--[[
 Creates a new empty PA suitable to hold space-charge for the given PA.
 Note: the space-charge array has one less grid point in each dimension.
--]]
function PAL.create_charge_pa(pa)
  local newpa = simion.pas:open()
  newpa:size(pa.nx-1, pa.ny-1, math.max(pa.nz-1, 1))
  newpa.symmetry = pa.symmetry
  newpa.dx_mm, newpa.dy_mm, newpa.dz_mm = pa.dx_mm, pa.dy_mm, pa.dz_mm
  local basename = get_pa_basename(pa)
     -- file name without extension
  newpa.filename = basename..'-charge.pa'
  newpa.refinable = false
  return newpa
end

--[[
 Returns a PA with filename `name`.  If the PA is already in memory,
 it is returned.  Directory paths are ignored in matching.
 If no such PA is found, it will be created in memory and returned.
 If the PA is created, it will have the same format as the PA `template_pa`
 in terms of size, symmetry, type, etc..  `off` has the same meaning
 as in PAL.create_similar_pa.
--]]
function PAL.cache_pa(name, template_pa, off)
  for i=1,#simion.pas do
    if (simion.pas[i].filename or ''):gsub('^.*[/\\]', '') == name then
      PAL.make_similar_pa(simion.pas[i], template_pa, off)
      return simion.pas[i]
    end
  end
  local newpa = PAL.create_similar_pa(template_pa, off)
  newpa.filename = name
  return newpa
end

--[[
 Fills PA `pa` with values returned by given function `f`.
 `f` is passed coordinates (x,y,z) in the system units (mm).
 `f` is expected to return (potential, electrode), where `potential` is
 the potential at that point and `electrode` is true or false to indicate
 whether the point should be an electrode point or not.
 The table `coords` is used to determine how to convert between system (mm) units
 and PA (gu) units.
 
 This will proper handle the case where `pa` represents space-charge (i.e.
 has one less grid point in each dimension), in which case values are sampled
 from the centers of cells not the vertices.
--]]
function PAL.fill_array_from_function(pa, coords, f)
  local nx_expect = math.floor((coords.xmax - coords.xmin) / coords.dx_mm + 1 + 0.5)
  local hx,hy,hz
  if pa.nx == nx_expect then -- Potential arrays sample on grid vertices
    hx,hy,hz = 0,0,0
  elseif pa.nx == nx_expect - 1 then  -- Space-charge arrays sample in center of grid cells
    hx,hy,hz = 0.5,0.5,0.5
  else
    error('PA and mesh size mismatch')
  end
  if pa.nz == 1 then hz = 0 end

  for xu,yu,zu in pa:points() do
    local x = (xu+hx-coords.xwo)*coords.dx_mm   -- Convert grid units (gu) to system units (mm).
    local y = (yu+hy-coords.ywo)*coords.dy_mm
    local z = (zu+hz-coords.zwo)*coords.dz_mm
    local potential, is_electrode = f(x,y,z)
    if potential    ~= nil then pa:potential(xu,yu,zu, potential) end
    if is_electrode ~= nil then pa:electrode(xu,yu,zu, is_electrode) end
  end
end

--[[
 This is a combination of `create_pa` and `fill_array_from_function` functions.
--]]
function PAL.build_pa_from_function(coords, f)
  local pa = PAL.create_pa(coords)
  PAL.fill_array_from_function(pa, coords, f)
  return pa
end

--[[
 Converts from system (mm) coordinates to PA (gu) coordinates.
--]]
function PAL.system_to_pa_units(coords, x,y,z)
  local xu = coords.xwo + x / coords.dx_mm
  local yu = coords.ywo + y / coords.dy_mm
  local zu = coords.zwo + z / coords.dz_mm
  return xu, yu, zu
end

--[[
 Creates a new PA whose potentials are the difference of the two given PAs.
 Summary statistics (maximum differences in elecrtode and non-elecrtode points)
 are also reported.
 This can be useful for visualizing how two PA's differ.
--]]
function PAL.diff_pas(apa, bpa)
  for _,name in ipairs{'nx','ny','nz','symmetry','dx_mm','dy_mm','dz_mm'} do
    if apa[name] ~= bpa[name] then
      error("PA's differ in "..name.." "..apa[name].." "..bpa[name])
    end
  end
  local diffpa = PAL.create_similar_pa(apa)
  local electrodes_same = true
  local max_diff_e = 0
  local max_diff_n = 0
  for x,y,z in diffpa:points() do
    local v1,e1 = apa:point(x,y,z)
    local v2,e2 = bpa:point(x,y,z)
    local dv = v2-v1
    diffpa:point(x,y,z, dv, e1 or e2)
    if e1 ~= e2 then electrodes_same = false end
    if e1 or e1 then
      max_diff_e = math.max(max_diff_e, math.abs(dv))
    else
      max_diff_n = math.max(max_diff_n, math.abs(dv))
    end
  end
  diffpa.refinable = false  -- This PA not intended to be refined.
  diffpa.filename = get_pa_basename(apa)..'-diff.pa'
  print('max |Delta V| non-electrode=', max_diff_n)
  print('max |Delta V| electrode=', max_diff_e)
  print('electrodes_same=', electrodes_same)
  return diffpa
end

--[[
 Creates an empty potential array object with size and symmetry
 defined by the grid specification and system size.

 `is_charge` is a Boolean indicating whether the array is will be a
 space-charge arrays.  Space-charge arrays are sized slightly
 differently than potential arrays.
 
 `gridspec` is a table containing these fields:
 
    `symmetry`
       symmetry for the array (e.g. '2dplanar')
    `d_mm`
       The length of a grid unit (mm per grid unit).
       to use in potential arrays.
    `dx_mm`
    `dy_mm`
    `dz_mm`
       If any of these specified, it overrides any value
       of d_mm in the given direction.  If dz_mm still
       unspecified, defaults to dy_mm.
    `N`
       A non-negative integer.  The grid density of any
       space-charge arrays is further reduced by a factor
       of 2^N compared to d_mm value given above.
 
 `size` is a table `{wx=wx, wy=wy, wz=wz}` describing the lengths of
   the solution space in mm in the x, y, and z dimensions.
--]]
function PAL.make_array(gridspec, size, is_charge)
  local pa = simion.pas:open()
  local d_mm  = gridspec.d_mm
  local dx_mm = gridspec.dx_mm or d_mm
  local dy_mm = gridspec.dy_mm or d_mm
  local dz_mm = gridspec.dz_mm or d_mm or dy_mm
  local N = is_charge and gridspec.N or 0
  local wx = size.wx        / dx_mm
  local wy = size.wy        / dy_mm
  local wz = (size.wz or 0) / dz_mm
  if is_charge then
    pa:size(wx/2^N, wy/2^N, math.max(1, wz/2^N))
  else
    pa:size(wx+1, wy+1, wz+1)
  end
  pa.symmetry = gridspec.symmetry
  pa.dx_mm, pa.dy_mm, pa.dz_mm = dx_mm, dy_mm, dz_mm
  return pa
end

--[[
 Builds a potential array that represents only the boundary
 conditions (electrode points) defined in the given
 system.  All non-electrode points are set to 0 V.
 The array is sized according to gridspec.
 
 Typically, one would later pass this array to the SIMION
 Refine function.
--]]
function PAL.build_boundary_potential_array(gridspec, system)
  local pa = PAL.make_array(gridspec, system.size, false)
  local d_mm  = gridspec.d_mm
  local dx_mm = gridspec.dx_mm or d_mm
  local dy_mm = gridspec.dy_mm or d_mm
  local dz_mm = gridspec.dz_mm or d_mm or dy_mm
  for xi,yi,zi in pa:points() do
    local x,y,z = xi*dx_mm, yi*dy_mm, zi*dz_mm
    local v            = system.potential(x,y,z)
    local is_electrode = system.boundary (x,y,z)
    pa:point(xi,yi,zi, is_electrode and v or 0, is_electrode)
  end
  return pa
end


--[[
 Builds a potential array containing the charge densities (in units
 of C/mm^3) in the given system.  The array is sized
 according to gridspec.
 
 Typically, one would later pass this array as a parameter in the
 SIMION Refine function for Poisson solving.
--]]
function PAL.build_charge_array(gridspec, system)
  local pa = PAL.make_array(gridspec, system.size, true)
  local N = gridspec.N
  local d_mm  = gridspec.d_mm
  local dx_mm = 2^N * (gridspec.dx_mm or d_mm)
  local dy_mm = 2^N * (gridspec.dy_mm or d_mm)
  local dz_mm = 2^N * (gridspec.dz_mm or d_mm or dy_mm)
  for xi,yi,zi in pa:points() do
    local x,y,z = (xi+0.5)*dx_mm,
                  (yi+0.5)*dy_mm,
                  (zi+0.5)*dz_mm
    local rho = system.charge(x,y,z)  -- (C/mm^3)
    pa:potential(xi,yi,zi, rho)
  end
  return pa
end

--[[
 Builds a potential array containing the potentials throughout all
 space in the given system `system`.  The array is sized according to
 `gridspec`.
 
 The potentials in here are theoretical potentials.  It can be
 instructive to compare potentials generated by the SIMION Refine
 function with these potentials.
--]]
function PAL.build_theoretical_potential_array(gridspec, system)
  local pa = PAL.make_array(gridspec, system.size, false)
  local d_mm  = gridspec.d_mm
  local dx_mm = gridspec.dx_mm or d_mm
  local dy_mm = gridspec.dy_mm or d_mm
  local dz_mm = gridspec.dz_mm or d_mm or dy_mm
  for xi,yi,zi in pa:points() do
    local x,y,z = xi*dx_mm, yi*dy_mm, zi*dz_mm
    local vt           = system.potential(x,y,z)
    local is_electrode = system.boundary (x,y,z)
    assert(vt == vt and vt+1 ~= vt, 'value is NaN or +-infinity')
    pa:point(xi,yi,zi, vt, is_electrode)
  end  
  return pa
end

--[[
 Returns an array (`pad`) that is the difference of arrays `pa1` and
 `pa2`.  That is, for each point `(x,y,z)` in the arrays,
 `pad(x,y,z) = pa1(x,y,z) - pa2(x,y,z)`.  It is assumed that the
 arrays have the same dimensions.
 
 Also returns the maximum absolute and relative error between the
 two arrays.
--]]
function PAL.build_difference_array(pa1, pa2)
  assert(pa1.nx == pa2.nx and pa1.ny == pa2.ny and pa1.nz == pa2.nz,
         "array dimensions don't match")

  -- Create difference array (pad = pa1 - pa2)
  local pad = simion.pas:open()
  pad:size(pa1:size())
  pad.symmetry = pa1.symmetry
  pad.dx_mm, pad.dy_mm, pad.dz_mm = pa1.dx_mm, pa1.dy_mm, pa1.dz_mm
  for xi,yi,zi in pa1:points() do
    local dvalue = pa1:potential(xi,yi,zi) - pa2:potential(xi,yi,zi)
    pad:potential(xi,yi,zi, dvalue)
  end

  -- Compute magnitude of |pad|.
  local max_abs_err = 0
  local maxv, minv = -math.huge, math.huge
  local abs = math.abs
  local min = math.min
  local max = math.max
  for xi,yi,zi in pa1:points() do
    local vt = pa1:potential(xi,yi,zi)
    local va = pa2:potential(xi,yi,zi)
    assert(va == va, 'NaN')
    assert(va+1 ~= va, 'inf')
    local abs_err = abs(va - vt)
    --local rel_err = abs_err / abs(vt)
    max_abs_err = max(max_abs_err, abs_err)
    maxv = max(maxv, vt)
    minv = min(minv, vt)
  end
  local max_rel_err = max_abs_err / (maxv - minv)

  return pad, max_abs_err, max_rel_err
end

--[[
 Computes the total charge in the space charge PA `chargepa`.
 Symmetry and mirroring is taken into account.
 For 2D planar arrays, this return C/mm.  For all others arrays, this returns C.
--]]
function PAL.total_charge(chargepa)
  local Q_sum = 0
  if chargepa.symmetry_type == '2dcylindrical' then
    for x,r in chargepa:points() do
      local cell_volume = ((r+1)^2-r^2)*math.pi*chargepa.dy_mm^2*chargepa.dx_mm  -- mm^3
      Q_sum = Q_sum + chargepa:potential(x,r,0) * cell_volume  -- (C/mm^3)*mm^3 = C
    end
  else -- 2D or 3D planar
    local cell_volume = chargepa.dy_mm*chargepa.dx_mm
                        *(chargepa.symmetry_type == '3dplanar' and chargepa.dy_mm or 1)
                        -- mm^2 (2D) or mm^3 (3D)
    for x,y,z in chargepa:points() do
      Q_sum = Q_sum + chargepa:potential(x,y,z) * cell_volume  -- C/mm (2D planar) or C (otherwise)
    end
  end
  if chargepa.mirror_x then Q_sum = Q_sum * 2 end
  if chargepa.symmetry_type ~= '2dcylindrical' then
    if chargepa.mirror_y then Q_sum = Q_sum * 2 end
    if chargepa.mirror_z then Q_sum = Q_sum * 2 end
  end
  return Q_sum
end

return PAL
