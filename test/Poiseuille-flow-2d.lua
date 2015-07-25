NX = 128
NY = 64

PARAMETERS = lb.LBGKParameters()

PARAMETERS.tau = 0.8

PY = 0

VT = 0.0 
VB = 0.0

START = 0 
END = 2048
FREQ = 128

cmap = lb.colormap()
cmap:append_color({1, 0, 0})
cmap:append_color({1, 1, 1})
cmap:append_color({0, 0, 1})

poiseuilleVelocity = function(iY, Ly, U) 
	y = iY/Ly
	return 4.*U*(y-y*y)
end

poiseuillePressure = function(iX, Lx, Ly, U, nu)
	Lx = Lx - 1
	Ly = Ly - 1
	return 8.*nu*U/(Ly*Ly)*(Lx/2.-iX);
end

density2rgb = function(rho, ux, uy)
	return cmap:map_value((rho - 2.7)/3)
end

initializer = function(x, y)
	return 2.5 
--	return 2.5 + math.random(), 0, 0
end

initialize = function(simulation, initializer, pinfo)
	--local pinfo = lb.LBPartitionInfo()
	--simulation:partition_info(pinfo)
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			local ux = poiseuilleVelocity(y, ny-1, 0.01)
			simulation:set_averages(x, y, 0,
						initializer(x + x0, y + y0), 
						ux, 
						0.0)
		end
	end

	simulation:set_equilibrium()
end

report_progress = function(start_wtime, t0, t1, t)
	local secs = lb.wtime() - start_wtime
	local steps = t - t0
	local sps = secs/steps

	print()
	io.write(string.format("<>> %d steps done in %f seconds (%f sec/step)\n",
			       steps, secs, sps))
	io.write(string.format("<>> ETA: %s\n", lb.wtime_string((t1 - t)*sps)))
end

make_filename = function(t, ext)
	local res = ""
	local digits = 1

	for n = 1, 8 do
		pp = math.pow(10, n)
		if math.floor(t / pp) == 0 then
			break
		else
			digits = digits + 1
		end
	end

	todo = 7 - digits
	for n = 0, todo do
		res = res .. "0"
	end

	return res .. t .. ext
end

render1 = function(callback, simulation, filename, pinfo)
--	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local rgb = lb.rgb(nx, ny)
	
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			local ux = poiseuilleVelocity(y, ny-1, 0.01)
			local rho = initializer(x + x0, y + y0)
			local uy = 0.0
			
			simulation:get_averages(rho, ux, uy)
			
			local c = callback(rho, ux, uy)
			
			rgb:set_pixel(x, ny - y - 1, c)
		end
	end

	rgb:save(filename)
	collectgarbage() -- this should on C side
end

filename_IJ = function(base, i, j)
	return string.format("%s-%d-%d.png", base, i, j)
end

renderN = function(callback, simulation, filename, pinfo)
--	local pinfo = simulation:partition_info()

	render1(callback, simulation,
		filename_IJ(filename, pinfo:processor_coords(0),
		    pinfo:processor_size(1) - pinfo:processor_coords(1) - 1))
	lb.mpi_barrier()

	if pinfo:processor_rank() == 0 then
		local cmd, filez
		cmd = string.format("montage +frame +shadow +label" ..
				    " -tile %dx%d -geometry %dx%d+0+0",
				    pinfo:processor_size(0),
				    pinfo:processor_size(1),
				    pinfo:size(0), pinfo:size(1))

		filez = ""
		for j = 0, pinfo:processor_size(1) - 1 do
			for i = 0, pinfo:processor_size(0) - 1 do
				filez = filez .. " "
				    .. filename_IJ(filename, i, j)
			end
		end

		cmd = cmd .. " " .. filez .. " -depth 8 " .. filename .. " "
		print(cmd)
		os.execute(cmd)

		cmd = "rm -f " .. filez
--		print(cmd)
--		os.execute(cmd)
	end
end

render = function(callback, simulation, t, pinfo)
	-- local filename = make_filename(t, ".h5")
	local filename = make_filename(t, ".png")

	if lb.is_parallel() == 1 then
		renderN(callback, simulation, filename, pinfo)
	else
		render1(callback, simulation, filename, pinfo)
	end
end

-- ##### Begin simulation code #####

simulation = lb.d2q9_BGK(1, NX, NY, 1, PY)
pinfo = lb.LBPartitionInfo()
simulation:partition_info(pinfo)
--math.randomseed(pinfo:processor_rank() + 1)

--pinfo.size[0] = 1

--print(pinfo.size[0])

simulation:set_parameters(PARAMETERS)

--[[if not PY then
	simulation:set_walls_speed(VT, VB)
end--]]

initialize(simulation, initializer, pinfo)

t0 = lb.wtime()
last_report = t0

--simulation:dump("init.h5")

for t = START, END do
	simulation:advance()

	if pinfo:processor_rank() == 0 then
		local now = lb.wtime()
		if now - last_report > 7 then
			report_progress(t0, START, END, t)
			last_report = now		
		end
	end

	if math.fmod(t, FREQ) == 0 then
		render(density2rgb, simulation, t, pinfo)
--		simulation:dump(make_filename(t))
		print(simulation:mass())
	end
end

par = lb.LBGKParameters()
simulation:get_parameters(par)

print("TAU = ", par.tau)
