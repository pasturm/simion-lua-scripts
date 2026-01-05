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

 with modifications by Patrick Sturm
 (c) 2020 TOFWERK
--]]

local GNUPLOT = {_VERSION='20170508'}
local Plot = {}

local CHART_TYPES = {}
CHART_TYPES.scatter = 'points'            -- scatterplot (xlXYScatter)
CHART_TYPES.scatter_lines = 'linespoints' -- scatterplot with connecting lines (xlXYScatterLines)
CHART_TYPES.bar_vertical = 'boxes'        -- vertical bar chart (xlColumnClustered)
CHART_TYPES.area = 'filledcurves'         -- area chart (xlArea)
CHART_TYPES.dots = 'dots'                 -- dots
CHART_TYPES.lines = 'lines'               -- lines
CHART_TYPES.none = 'none'                 -- no chart


local path = 'C:/Program Files/gnuplot/bin/gnuplot.exe'


--[[
 Opens and returns file handle to gnuplot input stream.
--]]
function GNUPLOT.open(plot_index)

  local fh = assert(io.popen('""'..path..'""', 'w'))

  -- fh:write('set terminal wxt\n')
    fh:write(('set terminal wxt %d\n'):format(plot_index or 0))
    -- fh:write('set terminal windows\n')
    -- fh:write('set terminal qt\n')

  -- Keep the handle around in global memory.  Garbage collection of
  -- the handle will close the graph.
  _G["gnuplot_fh"..(plot_index or 0)] = fh
  collectgarbage()  -- collect any previous handle
  return fh
end


--[[
  Plots data from table t in gnuplot.
  
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
  
  local fh = GNUPLOT.open(t.plot_index)  -- plot_index added by pst
  
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
    local color = dataset.color or {}
    for c=2,ncols do
      fh:write(('%s"-" using 1:2 title "%s" lc rgb "%s" '):format(
                (iseries==1 and '' or ', '), header[c] or '', 
                color[c-1] or 'dark-violet'))
      fh:write('with ', CHART_TYPES[self.chart_type], ' ')
      iseries = iseries + 1
    end
  end
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

-- Scatter plot with variable point sizes. Added by pst.
function GNUPLOT.plot_points_ps_variable(t)
  local fh = GNUPLOT.open(t.plot_index)

  fh:write(('set xlabel "%s"\n'):format(t.xlabel or ''))
  fh:write(('set ylabel "%s"\n'):format(t.ylabel or ''))
  fh:write(('set xrange [%s:%s]\n'):format(t.xmin or '*', t.xmax or '*'))
  fh:write(('set yrange [%s:%s]\n'):format(t.ymin or '*', t.ymax or '*'))
  fh:write(('set xtics %s\n'):format(t.xtics or 'auto'))
  fh:write(('set ytics %s\n'):format(t.ytics or 'auto'))
  fh:write(('set title "%s"\n'):format(t.chart_title or ''))
  fh:write('set grid xtics\n')
  fh:write('set grid ytics\n')
  fh:write('unset key\n')
  fh:write(t.gnuplot_commands or '', '\n')

  fh:write('plot "-" using 1:2:3 with points ps variable pt 1\n')

  dataset = t[1]
  local nrows = #dataset
  for r=1,nrows do
    fh:write(dataset[r][1], " ", dataset[r][2], " ", dataset[r][3], '\n')
  end
  fh:write('e\n')
  
  fh:flush()
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
      fh:write("splot '-' using 1:2:3 title '' with image\n")
    else
      fh:write("plot '-' using 1 title '' with boxes\n")
    end
    for yi=1,#ys do
      for xi=1,#xs do
        local v = (t[xi] or {})[yi] or 0
        fh:write(("%g %g %g\n"):format(xs[xi], ys[yi], v))
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
