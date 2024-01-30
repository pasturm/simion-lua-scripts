--[[
 gnuplotlib.lua
 Library for simplified recording and plotting of data in gnuplot.
 This is intended to be simpler and more rubust than directly invoking
 the underlying gnuplot command.
 This interface of this library resembles that of examples\excellib.lua ,
 and the two libraries are often interchangeable.

 Example:
   local PLOT = simion.import 'gnuplotlib.lua'
   local t = {{10, 20}, {20,25}, {30,27}, title = 'testing'}
   local plot =PLOT.plot(t)

   simion.sleep(5) -- wait five seconds before updating again
   t[1][2] = 30    -- change data
   plot:update_data(t)

 D.Manura,
 (c) 2011-2017 Scientific Instrument Services, Inc. Licensed under SIMION 8.1.

 with enhancements by
 Patrick Sturm
 (c) 2020 TOFWERK
--]]

local GNUPLOT = {_VERSION='20170508'}
local Plot = {}

local CHART_TYPES = {}
CHART_TYPES.scatter = 'points'            -- scatterplot (xlXYScatter)
CHART_TYPES.scatter_lines = 'linespoints' -- scatterplot with connecting lines (xlXYScatterLines)
CHART_TYPES.bar_vertical = 'boxes'        -- vertical bar chart (xlColumnClustered)
CHART_TYPES.area = 'filledcurves'         -- area chart (xlArea)
CHART_TYPES.dots = 'dots'                 -- dots  -- added by pst
CHART_TYPES.lines = 'lines'               -- lines  -- added by pst
CHART_TYPES.none = 'none'                 -- no chart


-- List of commands to try to invoke gnuplot with.  You may add to this list.
-- Prefix with '[sh]' to invoke command via *nix-like shell ('/bin/sh -c').
GNUPLOT.paths = {
  (simion._internal.simion_root or '') .. '/gnuplot/bin/gnuplot.exe',
  (simion._internal.simion_root or '') .. '/gnuplot/bin/pgnuplot.exe',
     -- gnuplot version installed in a "gnuplot" subdirectory of SIMION program folder.
     -- note: piping to pgnuplot.exe isn't working under Wine, so prefer gnuplot.exe.
  'c:/cygwin/bin/gnuplot.exe',
     -- gnuplot in typical Cygwin location.
  'gnuplot',
     -- gnuplot in current PATH.
  'sh:gnuplot',
     -- gnuplot in *nix-like shell PATH (e.g. native gnuplot from Wine).
  'C:/Program Files/gnuplot/bin/gnuplot.exe'
     -- use gnuplot 5.2.8 instead of gnuplot-sis-20171205 (5.2rc) where axis labels are not rotated properly.  -- pst
}


-- Determines if file is readable (and exists).
local function file_readable(path)
  local fh = io.open(path)
  if fh then fh:close(); return true end
  return false
end


-- Executes command, returning all text outputted.
local function readpipe(cmd)
  local fh = io.popen(cmd)
  local out = fh and fh:read'*a' or ''
  if fh then fh:close() end
  -- note: In Lua 5.1 (unlike 5.2beta), fh:close() for pipe returns no error info
  return out
end


local gnuplot_cmd   -- command to run gnuplot [*1]
local osname        -- OS version reported by gnuplot (nil if unknown)


local function wrapcmd(cmd, args)  --[*1]
  return cmd:match'^sh:' and
         ('/bin/sh -c "%s"'):format((cmd:gsub('^sh:', '') .. args):gsub('"', "'")) or
         '""' .. cmd .. '"' .. args .. '"'  --[*2]
end

--[[
 Locates gnuplot program, setting internal variables
 and returning command name [*1].  Might also return OS name reported by
 gnuplot via second return value (else nil).
 Program names in GNUPLOT.paths will be searched in order.
 Values cached from previous call are used unless force=true.
--]]
function GNUPLOT.find(force)
  -- Detect *nix-like shell.
  local has_sh = file_readable'/bin/sh'

  -- Search for gnuplot in each path.
  local log = ''
  if not gnuplot_cmd or force then 
    for _, path in ipairs(GNUPLOT.paths) do
      -- As shortcut, skip non-existent non-relative EXE's
      -- and *nix-like shell paths without *nix-like shell.
      if not(path:match'[/\\].*%.exe$' and not file_readable(path)) and
         not(path:match'^sh:' and not has_sh)
      then
        -- Try to execute gnuplot, reading OS platform info.
        -- Note: pgnuplot.exe and wgnuplot.exe (unlike gnuplot) do not output stderr,
        --   and io.popen/pclose do not output a status code in Lua 5.1 (unlike 5.2), 
        --   so the "echo ok" ensures that something is outputted is on success.
        --   However, "&& echo ok" oddly always executes under Wine (1.3.15), even
        --   if the preceeding command fails.  However, Wine will print 'File not found'.
        -- Note: desipte the '2>', Wine can still display a "wine: cannot find" warning
        --   message to the terminal.  It's not obvious how to suppress this.
        local cmd = wrapcmd(path, [[ -e "show version" 2>&1 && echo ok]])
        local out = readpipe(cmd)
        log = log .. '$ ' .. cmd .. '\n' .. out .. '\n'
        if out:find'ok' and not out:find'File not found' then
          gnuplot_cmd = path
          osname = out:match'System:%s*([^\r\n]+)'
          break -- found
        end
      end
    end
  end
  if not gnuplot_cmd then
    error('Could not find gnuplot.\n' ..
          'Used search paths: ' .. table.concat(GNUPLOT.paths, ';') .. '\n' ..
          'Log: ' .. log
          , 2)
  end
  return gnuplot_cmd, osname
end


--[[
 Opens and returns file handle to gnuplot input stream.
--]]
function GNUPLOT.open()
  local prog, osname = GNUPLOT.find()

  local fh
  if GNUPLOT_TRACE then
    -- only write gnuplot commands to text file that user can subsequently use
    -- in gnuplot.
    print 'WARNING: GNUPLOT_TRACE is ON'
    fh = assert(io.open('gnuplot.txt', 'w'))
  elseif osname and not (osname:match'MS%-Windows' or osname:match'CYGWIN') then
    -- When the Windows version of SIMION calls a Linux version of gnuplot
    -- via Wine, command piping doesn't work well when if the pipe is kept
    -- open for a long time.  Command piping doesn't seem to work at all
    -- for the OS X version of gnuplot, even outside of CrossOver (e.g.
    -- "echo 'plot sin(x); pause -1' | gnuplot" doesn't work on OS X).
    -- However, SIMION can communicate with the Linux or OS X version of
    -- gnuplot by using a named pipe created on the non-Windows side.
    -- Yes, this code is ugly.  [*1]
    local tempdir = readpipe([[/bin/sh -c "rm -fr /tmp/simion.gnuplot.*; ]] ..
                             [[mktemp -d /tmp/simion.gnuplot.XXXXXXXXXXXXXXX"]])
                             :gsub('[\r\n].*', '')
        -- note: filesystem on current directory might not support pipes
        --   (even if directory is writable), but /tmp probably does.
        -- The rm helps cleanup temporary folders from previous invocations.
    assert(tempdir:match'simion%.gnuplot%................$', tempdir)
    local pipepath = tempdir .. '/pipe'
    local cmd = (osname and osname:match'Darwin' and 'DISPLAY=:1 ' or '') ..
                gnuplot_cmd:gsub('^sh:', '') .. ' ' .. pipepath
                -- CrossOver on OS X (Darwin) changes the DISPLAY environment
                -- variable, which prevents the x11 terminal type from working
                -- in the OS X version of gnuplot.  "DISPLAY=:1" typically
                -- should fix that.  If you make gnuplot use aquaterm (and
                -- your version of gnuplot was compiled to support that),
                -- then DISPLAY shouldn't matter.
    local cmd = ([[/bin/sh -c "rm -f %s ; mkfifo %s; (%s &)"]])
                :format(pipepath, pipepath, cmd)
    if os.execute(cmd) ~= 0 then error('command failed: ' .. cmd) end
    simion.sleep(1) -- gnuplot must be reading from pipe before we write to it.
                    -- Is there a more robust way to wait?
    fh = assert(io.open(pipepath, 'w'), 'failed writing to pipe')
  else
    fh = assert(io.popen(wrapcmd(prog, ''), 'w'))
  end

  if osname and osname:match'CYGWIN' then
    -- nothing displays on non-xterm Cygwin without this.
    fh:write('set terminal ggi\n')
  else
    fh:write('set terminal wxt\n')
  end
  
  -- Keep the handle around in global memory.  Garbage collection of
  -- the handle will close the graph.
  _G.gnuplot_fh = fh
  collectgarbage()  -- collect any previous handle
  return fh
end



--[[
  Plots data from table t in gnuplot.
  t has a very similar format to EXCEL.plot in excellib.lua.
  
  t is normally an array of rows of column values.  Alternately, t
  may be an array of arrays of rows of column values.
  t may optionally also contain these fields:

    header - array of column headers.
    title - title for plot.  Defaults to none if omitted.
    xlabel - x label for plot.  Defaults to t.header[1] if omitted.
    ylabel - y label for plot.  Defaults to t.header[2] if omitted.
    chart_type - any chart type name string in CHART_TYPES above
                   (defaults to 'scatter').  'none' for no chart.
    lines - if true and chart_type is nil then chart_type is set to 'scatter_lines'
    xmin,xmax,ymin,ymax - minimum and xmaximum values for x and y axes.
      Defaults to automatically defined if omitted.
    xtics,ytics - steps between major tick markers.
      Defaults to automatically defined if omitted.
    xlog,ylog - if true then x or y axis is logarithmically scaled
    gap - gap between bars (for chart_type='bar_vertical'). Ratio
      between 0 to 1 (0 is no gap, 1 is large gap).
    gnuplot_commands - any additional gnuplot command string
      to execute prior to plot (defaults to nil).  You may use this
      to further customize the plotting.

  Returns a "Plot" object, which may subsequently be used to call
  Plot:update_data.  It also has "fh", which is the pipe file handle
  used to communicate with gnuplot.,

  Examples:

    -- Plot points (1,2) and (4,5)
    GNUPLOT.plot {{1,2}, {4,5}}

    -- Plot with additional parameters
    GNUPLOT.plot {header={'time', 'speed'},
                title='my plot', xlabel='t', ylabel='s', lines=true,
                {1,2}, {4,5}}

    -- Plot two data series having same X values.
    -- First column is X.  Second and third columns are Y for each series.
    GNUPLOT.plot {{1,2,3}, {4,5,6}}

    -- Plot two data series having possibly different X values.
    -- Two independent sets of data:
    -- First column is X.  Second column is Y.
    GNUPLOT.plot {{{1,2},{4,5}}, {{2,6},{3,7}}}
--]]
function GNUPLOT.plot(t)
  -- Set chart type.
  local chart_type = t.chart_type
  if t.lines then -- Compatibility with older versions of this library.
    assert(chart_type == nil)
    chart_type = 'scatter_lines'
  end
  chart_type = chart_type or 'scatter'
  if not CHART_TYPES[chart_type] then
    error("undefined chart type " .. tostring(chart_type), 2)
  end
  if chart_type == 'none' then return end
  
  local fh = GNUPLOT.open()
  
  local plot = setmetatable({fh=fh, chart_type=chart_type, chart_title=t.title or ''}, {__index = Plot})  

  fh:write(('set xlabel "%s"\n'):format(t.xlabel or ''))
  fh:write(('set ylabel "%s"\n'):format(t.ylabel or ''))
  fh:write(('set xrange [%s:%s]\n'):format(t.xmin or '*', t.xmax or '*'))
  fh:write(('set yrange [%s:%s]\n'):format(t.ymin or '*', t.ymax or '*'))
  fh:write(('set xtics %s\n'):format(t.xtics or 'auto'))
  fh:write(('set ytics %s\n'):format(t.ytics or 'auto'))
  if t.gap then
    fh:write(('set boxwidth %g relative\n'):format(1/(1 + t.gap)))
  end
  if t.gridlines == nil or t.gridlines then
    fh:write('set grid xtics\n')
    fh:write('set grid ytics\n')
  end
  local logscale = (t.xlog and 'x' or '') .. (t.ylog and 'y' or '')
  fh:write(logscale == '' and 'unset logscale' or 'set logscale '..logscale, '\n')
  fh:write(t.gnuplot_commands or '', '\n')

  plot:update_data(t)
  
  return plot
end


--[[
  Plots new data.
  Data table `t` is in the same format as in `plot`, except that
  only row/column data values are used.
--]]

function Plot:update_data(t)
  -- Normalize table.
  local datasets = (t[1] and t[1][1] and type(t[1][1]) ~= 'table') and {t} or t

  local fh = self.fh
  fh:write(('set title "%s"\n'):format(self.chart_title or ''))
  fh:write('plot ')
  local iseries = 1
  for i,dataset in ipairs(datasets) do
    local header = dataset.header or {}
    local ncols = #(dataset.header or dataset[1])
    local color = dataset.color or {}  -- added by pst
    for c=2,ncols do
      -- fh:write(('%s"-" using 1:2 title "%s" '):format(
      --           (iseries==1 and '' or ', '), header[c] or ''))
      fh:write(('%s"-" using 1:2 title "%s" lc rgb "%s" '):format(
                (iseries==1 and '' or ', '), header[c] or '', 
                color[c-1] or 'dark-violet'))  -- added by pst
      fh:write('with ', CHART_TYPES[self.chart_type], ' ')  -- added by pst
      iseries = iseries + 1
    end
  end
  -- fh:write('with ', CHART_TYPES[self.chart_type], ' ')
  fh:write('\n')

  
  -- Write each data set.
  for i,dataset in ipairs(datasets) do
    -- number of rows and columns in data.
    local nrows = #dataset
    local ncols = #(dataset.header or dataset[1])

    for c=2,ncols do
      for r=1,nrows do
        fh:write(dataset[r][1], " ", dataset[r][c], '\n')
      end
      fh:write('e\n')
    end
  end
  
  fh:flush()
end

--[[ OLD VERSION using temporary file and replot.
-- Writes array datafiles used by gnuplot.
local function update_data_files(t)
  -- Normalize table.
  local datasets = (t[1] and t[1][1] and type(t[1][1]) ~= 'table') and {t} or t
  
  -- Write each data set.
  for i,dataset in ipairs(datasets) do
    local fh = assert(io.open('tmp' .. i .. '.dat', 'w'))
    -- number of rows and columns in data.
    local nrows = #dataset
    local ncols = #(dataset.header or dataset[1])
    for r=1,nrows do
      for c=1,ncols do
        fh:write(dataset[r][c], " ")
      end
      fh:write('\n')
    end
    fh:close()
  end
end
function Plot:update_data(t)
  update_data_files(t)
  
  local fh = self.fh
  fh:write('replot\n')
  fh:flush()
end
--]]


--[[
  Updates plot title.
--]]
function Plot:title(name)
  --local fh = self.fh
  --fh:write(('set title "%s"\n'):format(name or ''))
  --fh:write('replot\n')
  --fh:flush()
  self.chart_title = name
end

--[[
  Closes the plot (if open).
--]]
function Plot:close()
end


--[[
  Plots data as 2D surface plot.

  Example:
    local plot = GNUPLOT.plot_surface {
         xs={100,200,300},ys={500,550,600},
         xlabel='V1',ylabel='V2',zlabel='beam size',zmin=0,zmax=10}
    for i=1,3 do for j=1,3 do
      plot:update(i,j, i*j)
    end end
--]]
function GNUPLOT.plot_surface(t)
  local xs = assert(t.xs)
  local ys = assert(t.ys)
  local title = t.title or nil
  local xlabel = t.xlabel or 'x'
  local ylabel = t.ylabel or 'y'
  local zlabel = t.zlabel or 'z'
  local zmin = t.zmin or nil
  local zmax = t.zmax or nil
  --unused in gnuplot implementation:
  --  local wb_old = t.wb or nil  -- optionally reuse (excel) workbook object
  
  local plot = {}
  plot._queue = {}
  plot.queued_updates = false -- whether updates are queued
  plot.last_update_time = os.clock()

  local fh = GNUPLOT.open()
  
  
  local function isnum(o)
    return type(o) == 'number' and o == o and o ~= math.huge and o ~= -math.huge
  end

  local function doplot()
    local is2d = #xs > 1 and #ys > 1
    fh:write(('set title "%s"\n'):format(title or ''))
    if is2d then
      -- fh:write("splot '-' using 1:2:3 matrix title '' with image\n")
      fh:write("splot '-' using 1:2:3 title '' with image\n")  -- added by pst
    else
      fh:write("plot '-' using 1 title '' with boxes\n")
    end
    for yi=1,#ys do
      for xi=1,#xs do
        local v = (t[xi] or {})[yi] or 0
        -- fh:write( isnum(v) and ("%g "):format(v) or 'nan ',
        --   ((xi==#xs or #ys==1) and '\n' or ''))
        fh:write(("%g %g %g\n"):format(xs[xi], ys[yi], v)) -- added by pst
      end
    end
    fh:write(is2d and 'e\ne\n' or 'e\n')
    fh:flush()
  end

  function plot:force_update()
    for i=#plot._queue,1,-1 do
      local data = plot._queue[i]
      if not t[data.i] then t[data.i] = {} end
      t[data.i][data.j] = data.value
      table.remove(plot._queue)
    end
    doplot()
    plot.last_update_time = os.clock()
  end
  
  -- Sets value at index [i,j] to given value.
  -- Updates will be queued (e.g. for performance) if plot.queue_updates is true.
  function plot:update(i, j, value)
    table.insert(plot._queue, {i=i,j=j,value=value})
    if not plot.queued_updates or os.clock() - plot.last_update_time > 2 then
      plot:force_update()
    end
  end
  
  -- Clears all values in plot to zero.
  function plot:clear()
    for yi=1,#ys do
    for xi=1,#xs do
      t[xi][yi] = 0
    end end
    doplot()
  end
   
  -- Changes plot title.
  function plot:title(name)
    title = name
  end
  
  -- Brings plot to top if multiple plots exist.
  function plot:activate()
    -- do nothing in gnuplot implementation.
  end

  local is2d = #xs > 1 and #ys > 1

  local function settics(name, vals)
    fh:write('set '..name..'tics (')
    for i=1,#vals do
      fh:write(('"%g" %d'):format(vals[i], i-1), i==#vals and '' or ', ')
    end
    fh:write(')\n')
  end
  if is2d then
    settics('x', xs)
    settics('y', ys)
  else
    settics('x', #ys == 1 and xs or ys)
  end
  
  fh:write(('set xtics %s\n'):format(t.xtics or 'auto'))
  fh:write(('set ytics %s\n'):format(t.ytics or 'auto'))
  if is2d then
    fh:write(('set xlabel "%s"\n'):format(xlabel))
    fh:write(('set ylabel "%s"\n'):format(ylabel))
    fh:write(('set zlabel "%s"\n'):format(zlabel))
  else
    local label = (#ys == 1) and xlabel or ylabel
    fh:write(('set xlabel "%s"\n'):format(label))
    fh:write(('set ylabel "%s"\n'):format(zlabel))
  end
  if is2d and (zmin and zmax) then
    fh:write(('set cbrange [%g:%g]\n'):format(zmin, zmax))
  end
  if not is2d then
    fh:write('set style fill solid\n')
    fh:write('set boxwidth 0.5\n')
  end

  fh:write(t.gnuplot_commands or '', '\n')
  fh:write('set view map\n')
  doplot()
  
  return plot
end



return GNUPLOT

-- [*1] command names prefixed by 'sh:' are invoked via the *nix-like shell
--      ('/bin/sh -c').  Typically this is for invoking native Linux programs
--      from the Windows version of SIMION via Wine.
-- [*2] http://stackoverflow.com/questions/2642551/windows-c-system-call-with-spaces-in-command
