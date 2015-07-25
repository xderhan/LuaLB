--
--  parameters
--

NX = 8
NY = 32

PY = 0

START = 0
END = 2*2024

PARAMETERS = lb.LBGKParameters()

PARAMETERS.tau = 0.8

init_data = function(x, y)
	return 2.5705, 0.0, 0.0
end

--
--  the script
--

init = function(s)
	local pinfo = s:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
        local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
                        local global_x = x0 + x
                        local global_y = y0 + y

			s:set_averages(x, y, 0, init_data(global_x, global_y))
		end
	end

	s:set_equilibrium()
end

-----------------------------------------------------------------------------

dump = function(s, name)
	local pinfo = s:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(0)
	local x = 1

	file = io.open(name, "w");

	for y = 0, ny - 1 do
		rho, ux, uy = s:get_averages(x, y, 0)
		file:write(string.format("%d %f %f %f\n", y, rho, ux, uy))
	end
end

-----------------------------------------------------------------------------

simulation = lb.d2q9_BGK(1, NX, NY, 1, PY)

pinfo = lb.LBPartitionInfo()
simulation.partition_info(pinfo)

init(simulation)

simulation:set_parameters(PARAMETERS)
simulation:set_walls_speed(0.1, 0.0)

pinfo = simulation:partition_info()
rank = pinfo:processor_rank()

dump(simulation, "start.dat")

for i = START, END do
	simulation:advance()
	M = simulation:mass()
	if rank == 0 and math.mod(i, 16) == 0 then
		io.write(string.format("(%d) ==> mass = %f\n", i, M))
	end
end

dump(simulation, "stop.dat")
