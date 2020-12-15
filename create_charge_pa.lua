--[[
SIMION Lua helper script to create a charge PA.

First remove all PAs from memory, then load the PA which you want to create a charge PA for into memory,
and then run this script.

Patrick Sturm
(c) 2019 TOFWERK
 --]]

simion.early_access(8.2)

-- load pa utility functions.
local PAL = simion.import('palib_pst.lua')

local pa = simion.pas[1]
local charge_pa = PAL.create_charge_pa(pa)
charge_pa:save()
