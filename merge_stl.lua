--[[
Function to merge multiple STL files into one.

STL files must be in binary format.

Example:
merge_binary_stl('final-1.stl', {"some-file.STL",  "some-other-file.STL"})
merge_binary_stl('final-2.stl', {"yet-another-file.STL"})

Adapted from SIMION SL Tools documentation.

Patrick Sturm
(c) 2020 TOFWERK
--]]

local function merge_binary_stl(dest, srcs)
  local destf = assert(io.open(dest, 'wb'))
  local function string2int(s)
    assert(#s == 4)
    local a,b,c,d = s:byte(1,4)
    return ((d*256 + c)*256 + b)*256 + a
  end
  local function int2string(i)
    assert(i >= 0 and i < 2^32)
    local a = i % 256; i = (i-a)/256
    local b = i % 256; i = (i-b)/256
    local c = i % 256; i = (i-c)/256
    local d = i % 256
    return string.char(a,b,c,d)
  end

  local srcfs = {}
  local ntriangles = {}
  local heads = {}
  local ntriangles_total = 0

  for i,src in ipairs(srcs) do
    srcfs[i] = assert(io.open(src, 'rb'))
    heads[i] = assert(srcfs[i]:read(80))
    ntriangles[i] = string2int(assert(srcfs[i]:read(4)))
    ntriangles_total = ntriangles_total + ntriangles[i]
    print(src, "ntriangles=" .. ntriangles[i])
  end

  destf:write(heads[1])
  destf:write(int2string(ntriangles_total))
  for i,src in ipairs(srcs) do
    local triangle_data = assert(srcfs[i]:read(ntriangles[i] * 50))
    srcfs[i]:close()
    destf:write(triangle_data)
  end

  destf:close()
end
