--[[
 Helper script to build and refine PAs.

Patrick Sturm
 (c) 2020 TOFWERK
--]]

local filename = "some-file"

local convergence = 1E-5 -- convergence objective to use.

simion.pas:close()  -- remove all PA's from RAM.

simion.command("gem2pa "..filename ..".gem "..filename..".pa#")

local pa = simion.pas:open(filename..".pa#")

pa:refine{convergence=convergence}

-- pas = {0,1,3,7}
-- for i=1,#pas do
-- 	local pa = simion.pas:open(filename..".pa#")
-- 	pa:refine{convergence=convergence, solutions = {pas[i]}}
-- end

simion.wb:load(filename..".iob")
