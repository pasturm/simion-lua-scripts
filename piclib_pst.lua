--[[
 piclib.lua
 SIMION utility functions for space-charge particle
 tracing (particle-in-cell -- PIC).

 Note: This code assumes that you properly set the ion_cwf variable
 yourself.  See code below and examples.

 D.Manura
 (c) 2008-2014 Scientific Instrument Services, Inc. (Licensed under SIMION 8.1)

 with some enhancements by
 Patrick Sturm
 (c) 2020 TOFWERK
--]]

local M = {_VERSION='0.1.20140721'}

-- check SIMION version.
assert(simion.VERSION >= simion.VERSION'8.1.0.43', 'SIMION 8.1.0.43 or above required.')

-- Whether particles represent 'charge' or 'current'.
local mode = 'current'

-- iff instance number n is a PA instance to Poisson solve,
-- then pa_info[n] is that charge array.
-- iff instance number n is a charge density array, pa_info[n] = 'charge'
-- In all other cases, pa_info[n] = nil.
local pa_info = {}
M.pa_info = pa_info

-- charge_skip[n] is number of grid cells in each dimension of PA instance to Poisson solve
-- for each corresponding cell in the charge density array.
local charge_skip = {}

-- insts[n] == simion.wb.instances[n]
local insts = {}

-- volume[n] is volume in mm^3 of one grid cell of PA instance number n.
-- For 2D cylindrical arrays, this refers to cells touching the axis.
-- For 2D planar arrays, this is an area instead.
local volumes = {}

-- Whether space-charge effects shall be on (1=yes,0=no).  Setting to
-- 0 disables this module, which is useful for quickly comparing the
-- system with and without space-charge effects.
adjustable PIC_enable = 1

-- Re-refine array every this number of time steps.  This is a
-- positive integer.  If your trajectory quality factor (TQual) is
-- zero, then in each time step particles move at most approximately 1
-- grid unit.  It's probably wastefull to re-refine prior to particles
-- advancing at least approximately charge_skip[i] grid units since the
-- charge density array won't change much.
-- This is only used when mode == 'charge'.
adjustable PIC_refine_period = 1

-- PE and contour display shall be updated every this number of
-- re-refines.  This is a positive integer.  Normally leave this at 1
-- to update display on each refine.  You might want to increase this
-- if display updates are slow and you don't care to see frequent
-- updates.
-- This is only used when mode == 'charge'.
adjustable PIC_display_period = 1

-- Radius (Gaussian standard deviation) to smooth (blur) over in the
-- space-charge array.  This is in units of grid units in the
-- space-charge array (i.e. multiples of charge_skip[i]).  Set to 0 for no
-- smoothing.  Smoothing makes the charge density more continuous
-- since tracing a relatively small number of representative particles
-- and a high grid density for the space-charge array can cause the
-- charge density distribution to be bumpy (e.g. like a poorly
-- bucketed histogram).  You can use the PIC_display_scaled_charge
-- option below to see if your charge density array is smooth.
adjustable PIC_smooth_radius = 5

-- Refining convergence objective (volts).  Default is 0.001 V in
-- SIMION. Range is 1e-1 to 1e-7.  Reduce value for greater accuracy
-- (and longer run-time).  Increase for shorter run-time (lower
-- accuracy).
adjustable PIC_refine_convergence = 1e-5

-- Number of reruns.
-- This is only used when mode == 'current'.
adjustable PIC_iterations = 10

-- These variables are only used when mode == 'current'.  They are
-- intended to improve the ability of the iterations to converge to
-- solution by slowly ramping up the current from some fraction
-- `PIC_current_ratio_start` of the final current, in steps of
-- size `PIC_current_ratio_step`.  For example, if `PIC_iterations`
-- is 6, `PIC_current_ratio_start` is 0.7, and `PIC_current_ratio_step`
-- is 0.1, then the currents in the 6 iterations would be 0.7, 0.8,
-- 0.9, 1.0, 1.0, and 1.0 times the final current.  Be sure that at
-- least two iterations are done at the final current.
-- Set `PIC_current_ratio_start` to 1 (the default) if you don't want
-- to ramp up the currents.
adjustable PIC_current_ratio_start = 1
adjustable PIC_current_ratio_step  = 0.1


-- These options can be used to scale the space-charge array values so
-- that its values are within a range PIC_display_scaled_charge_1 to
-- PIC_display_scaled_charge_2 volts so that they can be seen in the
-- PE and contour views.  Set both to the same value (e.g. 0) to
-- disable option.  Normally this option is disabled.  This option is
-- only for display purposes to help the user visualize the charge
-- density distribution using PE and contour plots.  It does not
-- affect the calculation since the scaling occurs after the charge
-- density array has been used in the refine procedure.
adjustable PIC_display_scaled_charge_1 = 0
adjustable PIC_display_scaled_charge_2 = 0

-- pst:
-- Perform initial refine (1=yes, 0=no) with no space charges 
-- (clears any previous refine that may have had charge data).
adjustable PIC_clear_refine = 1

-- Whether any SIMION initialize segment that needs to be called has
-- been called.  This can be used to avoid errors.  You may set this
-- to false at your top level and then set this to true inside your
-- SIMION initialize segment.  If the PIC module's other_actions
-- segment detects this value to be false, then it knows you wanted
-- initialize to be called but was not, and it raises an error.
M.is_init = true

M.charge_segment = {}   -- segments for mode == 'charge'
M.current_segment = {}  -- segments for mode == 'current'

-- For 'current' mode, this is the fraction of the final current
-- applies for the present iteration.
-- See PIC_current_ratio_start/PIC_current_ratio_step.
M.current_ratio = 1

-- Whether to use skipped point option in Refine.
-- Use false if the charge distribution does not change much between refines.
M.refine_skipped_point = true

-- Variables for triggering periodic updates.
local last_tof = 0     -- ion_time_of_flight in last time step
local rebuild = {}     -- rebuild[i] == true iff charge array or PA instance i is being rebuilt
local ndisplay = 0     -- display update counter
-- local nstep = 0        -- time step counter
adjustable nstep = 0 -- time step counter, [pst]

-- Cache for quick access.
local math  = math
local floor = math.floor


--[[
 Gets the base file name (i.e. without extension) of the given PA.
--]]
local function get_pa_basename(pa)
  local name = (pa.filename or 'blank.pa'):gsub('%.[^%.]+$', '')
  return name
end


--[[
 Re-scales value in charge array (pa) to be between values minv and
 maxv.  (This is used for display purposes to support
 PIC_display_scaled_charge.)
--]]
local function scale_charge_for_display(pa, minv, maxv)
  -- Find min and max values in original array.
  local minv_old,maxv_old = pa:potentials_minmax()

  -- Recale values in array.
  local span = maxv_old - minv_old
  local b = (span ~= 0) and (maxv - minv)/span or 0
  local a = minv - minv_old * b
  pa:potentials_scale(a, b)
end
M.scale_charge_for_display = scale_charge_for_display


--[[
 Applies charge density rho (unknown units) to cell (xi,yi,zi) and
 nearby cells in array pa_charge.  Charge is smoothed (blurred) to
 nearby cells.
--]]
local weights3d = {}
local weights2d = {}
local function smooth_charge(pa_charge, xi,yi,zi, rho, smooth_radius)
  -- Note: an improvement possibly could weight the blur matrix in a
  -- way to account more precisely for cylindrical symmetry (maybe
  -- useful).

  local smooth_size = smooth_radius * 3 -- three standard deviations.

  local M = smooth_radius*2+1
  local xm,ym,zm = xi-smooth_radius-1,yi-smooth_radius-1,zi-smooth_radius

  -- Rebuild matrix of weights.
  -- This is cached only updated if radius changed (faster).
  if #weights3d ~= M then
    -- print 'DEBUG:rebuild weights'
    local sum3d = 0
    local sum2d = 0
    weights3d = {}
    weights2d = {}
    local c = 1+smooth_radius
    for dz=1,M do  local wdz = {}; weights3d[dz] = wdz
    for dy=1,M do  local wdy = {};       wdz[dy] = wdy
    for dx=1,M do
      local val = math.exp(-((dx-c)^2+(dy-c)^2+(dz-c)^2) / (2 * smooth_size^2))
      wdy[dx] = val
      sum3d = sum3d + val
      if dz == c then
        sum2d = sum2d + val
      end
    end end end
    for dy=1,M do for dx=1,M do  -- normalize 2D
      weights2d[dy]         = weights2d[dy] or {}
      weights2d[dy][dx]     = weights3d[c ][dy][dx] / sum2d
    end end
    for dz=1,M do for dy=1,M do for dx=1,M do  -- normalize 3D
      weights3d[dz][dy][dx] = weights3d[dz][dy][dx] / sum3d
    end end end
  end

  -- Distribute charges in cell and nearby cells based on weights in
  -- matrix.
  if pa_charge.nz == 1 then -- 2D
    for dy=1,M do  local wdy = weights2d[dy]
    for dx=1,M do
      local x,y,z = xm+dx,ym+dy,zi
      if x < 0 then x = -x-1 end
      if y < 0 then y = -y-1 end
      pa_charge:potential_add(x,y,z, wdy[dx] * rho)
    end end
  else -- 3D
    for dz=1,M do  local wdz = weights3d[dz]
    for dy=1,M do  local wdy =       wdz[dy]
    for dx=1,M do
      local x,y,z = xm+dx,ym+dy,zm+dz
      if x < 0 then x = -x-1 end
      if y < 0 then y = -y-1 end
      if z < 0 then z = -z-1 end
      pa_charge:potential_add(x,y,z, wdy[dx] * rho)
    end end end
  end
end


--[[
 Deposits particle charge to space-charge array.
 Particle is assumed to be located at ion_px_mm,ion_py_mm,ion_pz_mm
 and represent physical charge ion_cwf.
 Units of ion_cwf are
   For 2D planar,      C mm^-1
   For 2D cylindrical, C
   For 3D,             C

 This is designed to be called within a SIMION other_actions segment.

 `pa_charge` - space charge array object to deposit charge into.
 `PIC_scale` - number of grid cells in each dimension of the PA per
   each grid cell in the charge density array
 `smooth_radius` - radius (Gaussian standard deviation) in charge array
   grid units to smooth over.  Set to 0 for no smoothing.
--]]
local TWO_PI = 2*math.pi
local function deposit_charge(pa_charge, PIC_scale, smooth_radius)
  -- Obtain space-charge array cell to store charge into.
  local xgu,ygu,zgu =
    pa_coords_to_array_coords(wb_coords_to_pa_coords(
      ion_px_mm,ion_py_mm,ion_pz_mm))
  local xi,yi,zi = floor(xgu/PIC_scale),floor(ygu/PIC_scale),
                   floor(zgu/PIC_scale)
  -- print('DEBUG', ion_number, xi,yi,zi, ion_mm_per_grid_unit)
  -- assert(pa_charge:inside(xi,yi,zi))

  -- Charge represented by this time-step.
  local q = ion_cwf
  if mode == 'current' then
    -- charge = current * time
    q = q * (ion_time_step * 1E-6) * M.current_ratio
      -- for 3D or 2D cylindrical: C; for 2D planar: (C/mm)
  end

  -- Charge density to add to array.
  local rho  -- C/mm^3
  if pa_charge.symmetry_type == '2dcylindrical' then
    local vol_mm3 = volumes[ion_instance] * (2*ygu+1)
      -- cell volume (mm^3).
      -- note: vol(ygu) = dx_mm*PI*(dy_mm*(ygu+1))^2 - dx_mm*PI*(dy_mm*ygu)^2
      --                = (dx_mm*PI*dy_mm^2)*(2*ygu+1)
      --       so vol(ygu) = vol(0)*(2*ygu+1)
    rho = q / vol_mm3  -- C / mm^3
  else  -- 2D planar or 3D
    rho = q / volumes[ion_instance]
      -- for 3D: C / mm^3; for 2D planar: (C/mm)/(mm^2)
  end

  -- Update charge array.
  if smooth_radius == 0 then
    pa_charge:potential_add(xi,yi,zi, rho)   -- to this cell only.
  else
    smooth_charge(pa_charge, xi,yi,zi, rho, smooth_radius)
  end
end
M.deposit_charge = deposit_charge

--[[
  Added by pst.
  Same as deposit_charge() but with instance of charge pa as argument.
  This allows to deposit charges on a charge_pa which does not correspond 
  to current ion_instance. For example, charges can be deposited to 
  coarse charge pa while ion is flying in fine pa.
--]]
function M.deposit_charge2(pa_charge_inst, PIC_scale)
  -- Obtain space-charge array cell to store charge into.
  local xgu,ygu,zgu =
    pa_charge_inst:pa_to_array_coords(pa_charge_inst:wb_to_pa_coords(
      ion_px_mm,ion_py_mm,ion_pz_mm))
  local xi,yi,zi = floor(xgu/PIC_scale),floor(ygu/PIC_scale),
                   floor(zgu/PIC_scale)
  -- print('DEBUG', ion_number, xi,yi,zi, ion_mm_per_grid_unit)
  -- assert(pa_charge:inside(xi,yi,zi))

  -- Charge represented by this time-step.
  local q = ion_cwf
  if mode == 'current' then
    -- charge = current * time
    q = q * (ion_time_step * 1E-6) * M.current_ratio
      -- for 3D or 2D cylindrical: C; for 2D planar: (C/mm)
  end

  -- Charge density to add to array.
  local rho  -- C/mm^3
   -- Compute volume (or area in 2D planar) of charge cell.
  local dx_mm = pa_charge_inst.pa.dx_mm
  local dy_mm = pa_charge_inst.pa.dy_mm
  local dz_mm = pa_charge_inst.pa.dz_mm
  local volumes2 =
    pa_charge_inst.pa.symmetry_type == '3dplanar'      and dx_mm*dy_mm*dz_mm or
    pa_charge_inst.pa.symmetry_type == '2dcylindrical' and dx_mm*dy_mm*dy_mm*math.pi or
    pa_charge_inst.pa.symmetry_type == '2dplanar'      and dx_mm*dy_mm or error'assert'

  if pa_charge_inst.pa.symmetry_type == '2dcylindrical' then
    local vol_mm3 = volumes2 * (2*ygu+1)
      -- cell volume (mm^3).
      -- note: vol(ygu) = dx_mm*PI*(dy_mm*(ygu+1))^2 - dx_mm*PI*(dy_mm*ygu)^2
      --                = (dx_mm*PI*dy_mm^2)*(2*ygu+1)
      --       so vol(ygu) = vol(0)*(2*ygu+1)
    rho = q / vol_mm3  -- C / mm^3
  else  -- 2D planar or 3D
    rho = q / volumes2
      -- for 3D: C / mm^3; for 2D planar: (C/mm)/(mm^2)
  end

  -- Update charge array.
  pa_charge_inst.pa:potential_add(xi,yi,zi, rho)   -- to this cell only.
end

-- `initialize_run` segment for 'charge' mode.
function M.charge_segment.initialize_run()
  -- Refine with no space charges (clears any previous refine that may
  -- have had charge data).  This is the initial refine.  It could
  -- take longer than subsequent refines since we are calculating
  -- from scratch here.
  for i, val in pairs(pa_info) do
    if type(val) == 'userdata' then
      local pa = simion.wb.instances[i].pa
      pa:refine{convergence=PIC_refine_convergence}
    end
  end
  
  -- initialize variables for next run
  last_tof = 0
  rebuild = {}
  ndisplay = 0
  nstep = 0
end


-- `other_actions` segment for 'charge' mode.
-- Called for each time step and each particle.
function M.charge_segment.other_actions()
  if PIC_enable == 0 then return end -- exit if Poisson effects disabled
  if ion_time_of_flight <= ion_time_of_birth then return end -- not yet born
  local charge_pa = pa_info[ion_instance]
  print(pa_info[ion_instance])
  if not charge_pa then return end
  if charge_pa == 'charge' then
    error('Particle must not see fields in charge density array instance number ' .. ion_instance .. '. ' ..
          'Make sure each charge density array is outside of the particle flying region ' ..
          'or fully overlapped by and lower in priority than its corresponding PA instance ' ..
          'to Poisson solve.')
  end

  -- When next time step...
  if ion_time_of_flight > last_tof then
    last_tof = ion_time_of_flight
    nstep = nstep + 1
    -- When charge data was rebuilt in previous time step...
    if rebuild[ion_instance] then
      rebuild[ion_instance] = false

      -- Refine PA given space-charge data.
      --
      -- Note: some performance may be gained by disabling
      -- skipped-point refining (skipped_point=false) when the
      -- space-charge is slowly and continuously changing.
      -- In practice, it doesn't seem to make much difference overall.
      -- SIMION does a non-skipped refine on the first step anyway.
      local pa = insts[ion_instance].pa
      --if not pcall(pa.refine, pa,
      --             {charge=pa_charge,
      --              convergence=PIC_refine_convergence,
      --              max_iterations=100,
      --              skipped_point=false})
      --then
      --  print 'Skipped point refine...'
      pa:refine{charge=charge_pa,
                convergence=PIC_refine_convergence,
                skipped_point=M.refine_skipped_point}
      --end

      -- Periodically update PE/contour display.
      if ndisplay > PIC_display_period then
        ndisplay = 1

        -- Note: space-charge array data is no longer needed, but
        -- we may re-scale it for display purposes so that
        -- charge density contours can be shown on the PE and
        -- contour views.
        if PIC_display_scaled_charge_1 ~= PIC_display_scaled_charge_2 then
          scale_charge_for_display(charge_pa,
            PIC_display_scaled_charge_1, PIC_display_scaled_charge_2)
        end

        -- Refresh sceeen for new potential values.
        -- note: sim_update_pe_surface is not sufficient since "pa"
        -- is not an active PA, and that doesn't update contour views.
        redraw_screen()
      else
        ndisplay = ndisplay + 1
      end
    end
    if nstep == PIC_refine_period then
      nstep = 0
      assert(M.is_init, "initialize segment not called")
      -- Start rebuilding space-charge array...
      rebuild[ion_instance] = true
      charge_pa:clear()
    end

    --print(ion_number, ion_time_of_flight)
  end

  -- Rebuild charge data.
  if rebuild[ion_instance] then
    if charge_pa then
      deposit_charge(charge_pa, charge_skip[ion_instance], PIC_smooth_radius)
    end
  end
end


-- `flym` segment for 'charge' mode.
function M.charge_segment.flym()
  -- Refine with no space charges (clears any previous refine that may
  -- have had charge data).  This is the initial refine.  It could
  -- take longer than subsequent refines since we are calculating
  -- from scratch here.
  if PIC_clear_refine ~= 0 then -- pst
	  for i, val in pairs(pa_info) do
		if type(val) == 'userdata' then
		  local pa = simion.wb.instances[i].pa
		  pa:refine{convergence=PIC_refine_convergence}
		end
	  end
	  redraw_screen()
  end -- pst
  run()
end

-- `flym` segment for 'current' mode.
function M.current_segment.flym()
  -- Refine with no space charges (clears any previous refine that may
  -- have had charge data).  This is the initial refine.  It could
  -- take longer than subsequent refines since we are calculating
  -- from scratch here.
  if PIC_clear_refine ~= 0 then -- pst
	  for i, val in pairs(pa_info) do
		if type(val) == 'userdata' then
		  local pa = simion.wb.instances[i].pa
		  pa:refine{convergence=PIC_refine_convergence}
		end
	  end
	  redraw_screen()
  end -- pst
  
  -- Standard single-run behavior if Poisson effects disabled.
  if PIC_enable == 0 then run() return end

  if sim_repulsion ~= 'none' then
    error('\n\nCharge repulsion effects (on Particles tab) should be disabled '..
          'when Poisson piclib is used since space-charge would be counted twice.\n')
  end
  
  for i=1, PIC_iterations do
    M.current_ratio = math.min(1, PIC_current_ratio_start + PIC_current_ratio_step * (i-1))
    sim_rerun_flym = (i>=PIC_iterations-1) and 0 or 1
       -- preserve trajectories in last two runs
    print('Iteration ', i, 'current_ratio=', M.current_ratio)
    run()
  end
end


-- `other_actions` segment for 'current' mode.
function M.current_segment.other_actions()
  if PIC_enable == 0 then return end -- exit if Poisson effects disabled
  if ion_time_of_flight <= ion_time_of_birth then return end -- not yet born
  local charge_pa = pa_info[ion_instance]
  if not charge_pa then return end
  if charge_pa == 'charge' then
    error('Particle must not see fields in charge density array instance number ' .. ion_instance .. '. ' ..
          'Make sure each charge density array is outside of the particle flying region ' ..
          'or fully overlapped by and lower in priority than its corresponding PA instance ' ..
          'to Poisson solve.')
  end

  deposit_charge(charge_pa, charge_skip[ion_instance], PIC_smooth_radius)
end


-- `initialize_run` segment for 'current' mode.
function M.current_segment.initialize_run()
  -- Clear space-charge density PA at start of rerun.
  for _, val in ipairs(pa_info) do
    if type(val) == 'userdata' then
      local pa = val
      pa:clear()
    end
  end
end


-- `initialize` segment for 'current' mode.
function M.current_segment.initialize()
  if PIC_enable == 0 then return end -- exit if Poisson effects disabled
  if ion_run == PIC_iterations then
    ion_color = 2  -- change color of trajectories in last run.
  end
end


-- `terminate_run` segment for 'current' mode.
function M.current_segment.terminate_run()
  if PIC_enable == 0 then return end -- exit if Poisson effects disabled

  -- Refine electric field PA at end of run, using space-charge
  -- determined during run.
  for i, val in ipairs(pa_info) do
    if type(val) == 'userdata' then
      local pa = simion.wb.instances[i].pa
      local charge_pa = val
      pa:refine{charge=charge_pa,
            convergence=PIC_refine_convergence,
            skipped_point=M.refine_skipped_point}
    end
  end
  
  -- Note: space-charge array data is no longer needed, but
  -- we may re-scale it for display purposes so that
  -- charge density contours can be shown on the PE and
  -- contour views.
  for i, val in ipairs(pa_info) do
    if type(val) == 'userdata' then
      local charge_pa = val
      if PIC_display_scaled_charge_1 ~= PIC_display_scaled_charge_2 then
        scale_charge_for_display(charge_pa,
        PIC_display_scaled_charge_1, PIC_display_scaled_charge_2)
      end
    end
  end

  -- Refresh sceeen for new potential values.
  -- note: sim_update_pe_surface is not sufficient since "pa"
  -- is not an active PA, and that doesn't update contour views.
  redraw_screen()
end


-- Incorporates segments in table `newsegment` into regular `segment` table.
function M.merge_segments(newsegment)
  for name, seg in pairs(newsegment) do
    local oldseg = segment[name]
    if oldseg then
      segment[name] = function() oldseg(); seg() end -- merge
    else
      segment[name] = seg
    end
  end
end


--[[
 Tells the PIC code that PA instance number `iinstance` should be
 Poisson solved with space-charge PA `chargepa`.  If `chargepa` is `nil`,
 then a space-charge PA is created in memory (reusing if it already exists).
 This need not be called explicitly if PA's are found with detect_charge_pas.
--]]
function M.add_charge_array(iinstance, chargepa)
  local pa = simion.wb.instances[iinstance].pa
  
  -- Find or create if necessary.
  if chargepa == nil then
    chargepa = M.find_or_create_charge_pa(pa)
  end
  
  -- Check.
  local mx = (pa.nx - 1) / (chargepa.nx) -- improve?
  local my = (pa.ny - 1) / (chargepa.ny)
  local mz = (pa.nz - 1) / (chargepa.nz)
  if math.floor(mx) ~= mx or math.floor(my) ~= my or
     not(pa.nz == 1 and chargepa.nz == 1) and math.floor(mz) ~= mz or
     mx ~= my or pa.nz ~= 1 and mx ~= mz
  then
        error('Array size for charge density pa '..tostring(chargepa)..' is inconsistent with array size for '..
          'potential array instance '..tostring(pa)..' to Poisson solve.  For potential array size '..
          '(nx,ny,nz) points, charge density array size should be '..
          '((nx-1)/2^N, (ny-1)/2^N, (nz == 1) and 1 or (nz-1)/2^N) for some integer N >= 0.')
  end
  
  charge_skip[iinstance] = mx
  pa_info[iinstance] = chargepa

  local jinstance
  for i=1,#simion.wb.instances do
    if simion.wb.instances[i].pa == chargepa then jinstance = i end
  end
  if jinstance then
    if jinstance > iinstance then
      error('Charge density PAs must preceed PAs to Poisson solve in the'..
            'PA Instances list on the PAs tab.  Note: Charge density instance number '..iinstance..
            ' corresponds to PA instance number '..jinstance..'.  Use L+/L- buttons to re-order.')
    end
    pa_info[jinstance] = 'charge'
  end
  
  -- Compute volume (or area in 2D planar) of charge cell.
  local dx_mm = chargepa.dx_mm
  local dy_mm = chargepa.dy_mm
  local dz_mm = chargepa.dz_mm
  volumes[iinstance] =
    chargepa.symmetry_type == '3dplanar'      and dx_mm*dy_mm*dz_mm or
    chargepa.symmetry_type == '2dcylindrical' and dx_mm*dy_mm*dy_mm*math.pi or
    chargepa.symmetry_type == '2dplanar'      and dx_mm*dy_mm or error'assert'
end


--[[
 If space-charge PA, returns filename without the terminal '-charge.pa'.
 Otherwise, returns `nil`.
--]]
local function find_charge_pa_basename(pa)
  if pa.filename then
    local base = pa.filename:match'^(.*)-[cC][hH][aA][rR][gG][eE]%.[pP][aA]$'
    return base or nil
  end
  return nil
end


--[[
 Determines electric field PA instance cooresponding to
 space-charge PA `chargepa`.
--]]
local function painstance_from_chargearray(chargepa)
  local basename = find_charge_pa_basename(chargepa)
  if not basename then
    local filename = chargepa.filename or ''
    error('charge PA file name "'..filename..'" should end with -charge.pa')
  end
  for i=1, #simion.wb.instances do
    local painstance = simion.wb.instances[i]
    if get_pa_basename(painstance.pa) == basename then
      return painstance, i
    end
  end
  return nil, nil, "'basename.pa*' not found in workbench PA instances list."
end


--[[
 Detects charge density arrays and corresponding PA's to Poisson solve.
 Only PA's added to the workbench PAs tab PA Instances list are auto-detected
 (others may be added explicitly via `add_chage_array`).
 
 It is critical that the particles never see the charge density array.
 You can ensure this by one of these:
   - don't add the space-charge PA to the workbench.  See `add_charge_array`.
   - move the space-charge array to some location where
     particles don't fly
   - overlap the space-charge array instance with its
     corresponding electric field PA instance.  As a safety, in the
     PA instances list, the charge density arrays are required to preceed
     the regular PA's.
--]]
function M.detect_charge_pas()
  -- Locate PA instances with name ending with '-charge.pa' and process...
  for i=1, #simion.wb.instances do
    local chargeinstance = simion.wb.instances[i]
    local base = find_charge_pa_basename(chargeinstance.pa)
    if base then -- file name ends with '-charge.pa'
      local painstance_, j, err = painstance_from_chargearray(chargeinstance.pa)
      if err then error(err) end
      M.add_charge_array(j, chargeinstance.pa)
    end
  end
end


--[[
 Finds charge PA in memory suitable for PA `pa` to Poisson solve.
 This will look first in the PA instances list on the workbench
 and then in the list of all PAs in memory.
 In the PA Instances list, any PA with name ending with '-charge.pa' is recognized
 as a charge-density PA.  It will be used when Poisson solving any PA with the
 same base name (e.g. 'test-charge.pa' will be used on 'test.pa').
--]]
function M.find_charge_pa(pa)
  if not pa.filename then return nil end -- must have name to match
  -- First look in PA instances in workbench.
  for i=1,#simion.wb.instances do
    local binst = simion.wb.instances[i]
    local basename = find_charge_pa_basename(binst.pa.filename)
    if basename and pa.filename:lower() == basename:lower() then
      if binst.scale ~= 1 then
        error("Space charge PA instance " .. binst.pa.filename ..
          " must have scale = 1.  Set mm/gu factor in PA file instead of PA instance.")
      end
      return binst.pa
    end
  end
  -- Next look in all memory.
  for i=1,#simion.pas do
    local bpa = simion.pas[i]
    local basename = find_charge_pa_basename(bpa)
    if basename and get_pa_basename(pa):lower() == basename:lower() then
      return bpa
    end
  end
  return nil -- not found
end


--[[
 Creates a space-charge PA compatible with given electric field PA `pa`.
 If `chargepa` is not `nil`, then the space-charge PA is created
 in that object.
--]]
function M.create_charge_pa(pa, chargepa)
  chargepa = chargepa or simion.pas:open()
  local nx,ny,nz = pa:size()
  nx = nx - 1
  ny = ny - 1
  if nz > 1 then nz = nz - 1 end
  chargepa:size(nx,ny,nz)
  chargepa.symmetry = pa.symmetry
  chargepa.dx_mm,chargepa.dy_mm,chargepa.dz_mm = pa.dx_mm,pa.dy_mm,pa.dz_mm
  chargepa.filename = get_pa_basename(pa)..'-charge.pa'
  chargepa.refinable = false
  return chargepa
end


--[[
 Returns a space-charge PA compatible with electric field PA `pa`.
 Reuses if possible; creates if necessary.
--]]
function M.find_or_create_charge_pa(pa)
  local chargepa = M.find_charge_pa(pa)
  chargepa = M.create_charge_pa(pa, chargepa)
  return chargepa
end


--[[
 Returns the charge density (C/mm^3) at position (x,y,z) in workbench coordinates.
 This consults charge density PA's associated with PA's on the workbench and
 honors PA instance priorities.
 On use for this is contour plotting:
   simion.import'../contour/contourlib81.lua'.plot(PIC.charge_density)
--]]
function M.charge_density(x,y,z)
  for i=#simion.wb.instances,1,-1 do -- highest to lowest priority
    local painst = simion.wb.instances[i]
    local pa = painst.pa
    local chargepa = pa_info[i] print(chargepa)
    if pa.potential_type == 'electric' and type(chargepa) == 'userdata' then
      local xp,yp,zp = painst:wb_to_pa_coords(x,y,z)
      if pa:inside_vc(xp,yp,zp) then
        local xa,ya,za = painst:pa_to_array_coords(xp,yp,zp)
        xa = math.floor(xa)
        ya = math.floor(ya)
        za = math.floor(za)
        if xa == pa.nx-1 then xa = xa-1 end
        if ya == pa.ny-1 then ya = ya-1 end
        if za == pa.nz-1 and pa.nz~=1 then za = za-1 end
        return chargepa:potential(xa,ya,za)
        -- IMPROVE: this could likely be improved to interpolate better
      end
    end
  end
end


--[[
 Initializes PIC code.  This includes auto-detecting PA's to use and applying
 segments.  This should be called at the top-level, before running segments.
 
 `mode_` - either 'charge' or 'current'.
   If 'current', particles are treated as rays of time-independent current.
     The Poisson solver is invoked between reruns.
   If 'charge', particles are reated as point charges (possibly time-dependent).
     Particles are flown grouped and the Poisson solver is invoked whever
     particles have moved significantly.
--]]
function M.load(mode_)
  -- Set mode.
  mode = mode_ or 'current'
  assert(mode == 'current' or mode == 'charge', 'invalid mode')
  M.segment = (mode == 'current') and M.current_segment or M.charge_segment

  -- Force grouped flying in 'charge' mode.  
  if mode == 'charge' then
    if sim_grouped == 0 then sim_trajectory_quality = 0 end
    sim_grouped = 1
  end

  -- Detect space-charge and Poisson solving PA's.
  M.detect_charge_pas()

  insts = simion.wb.instances

  -- Install segments.
  M.merge_segments(mode == 'charge' and M.charge_segment or M.current_segment)

end

_G.PIC = M

return M
