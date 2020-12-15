--[[
 Adds a grid plane to a PA.

Patrick Sturm
 (c) 2020 TOFWERK
--]]

local convergence = 1E-5 -- convergence objective to use.

simion.pas:close()  -- remove all PA's from RAM.
local pa = simion.pas:open("CTOF.pa#")

yi = 403
for zi = 72,132 do
	for xi = 46,138 do
		pa:point(xi,yi,zi, 0, true)
	end
end
pa:save("CTOF_withgrid.pa#")


pa:refine{convergence=convergence}