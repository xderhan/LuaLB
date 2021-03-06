-- Helper routines

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
	--local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local rgb = lb.rgb(nx, ny)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			local rho, ux, uy = simulation:get_averages(x, y, 1)
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
	--local pinfo = simulation:partition_info()

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

render = function(callback, simulation, pinfo, t, ext)
	local filename = make_filename(t, ext)

	if lb.is_parallel() == 1 then
		renderN(callback, simulation, filename, pinfo)
	else
		render1(callback, simulation, filename, pinfo)
	end
end
