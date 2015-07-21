require("utils.lua")

NX = 128
NY = 128

PARAMETERS = lb.LBGKParameters()

tau = 1.0
a = -0.1
b = 0.1
k = 0.05

PARAMETERS:set(tau)

PY = 1

VT = 0.0 
VB = 0.0

START = 0 
END = 2048
FREQ = 64

cmap = lb.colormap()
cmap:append_color(1, 0, 0)
cmap:append_color(1, 1, 1)
cmap:append_color(0, 0, 1)

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
	return cmap:map_value( (rho + 1)/2 )
end

initializer = function(x, y)
--	return 2.5 
	return 0.001*math.random(-1.0,1.0)
end

initialize = function(simulation, initializer)
	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			simulation:set_averages(x, y, 0,
					1.+initializer(x + x0, y + y0),0,0)
					
			simulation:set_averages(x, y, 1,
					initializer(x + x0, y + y0),0,0)
		end
	end

	simulation:set_equilibrium()
end

-- ##### Begin simulation code #####
my_collide = function(simulation)
	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do

			p,ux,uy = simulation:get_averages(x, y, 1)
			dxp,dyp,ddp = simulation:diff_opt(x, y, 1)
			
			r,ux,uy = simulation:get_averages(x, y, 0)
			
			local pI  = 0.5*a*p*p + 0.75*b*p*p*p*p - k*( p*ddp + 0.5*(dxp*dxp+dyp*dyp) ) + r/3.

			local pxx = pI + k*dxp*dxp
			local pyy = pI + k*dyp*dyp
			local pxy =      k*dxp*dyp
			
			local mu =  a*p + b*p*p*p - k*ddp
	
			-- 1st component
			local rux = r*ux
			local ruy = r*uy
	
			local ruxx = rux*ux + pxx
			local ruyy = ruy*uy + pyy
	
			local ruxy = 2.*rux*uy + 2.*pxy
			local rusq = ruxx + ruyy
			
			-- v0 - v3 
			simulation:collide_node(0,x,y,0, (1./9.)*( 3.*rux+4.5*ruxx-1.5*rusq))
			simulation:collide_node(1,x,y,0, (1./9.)*( 3.*ruy+4.5*ruyy-1.5*rusq))
			simulation:collide_node(2,x,y,0, (1./9.)*(-3.*rux+4.5*ruxx-1.5*rusq))
			simulation:collide_node(3,x,y,0, (1./9.)*(-3.*ruy+4.5*ruyy-1.5*rusq))

			-- v4 - v7 
			simulation:collide_node(4,x,y,0, (1./36.)*(3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
			simulation:collide_node(5,x,y,0, (1./36.)*(3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
			simulation:collide_node(6,x,y,0, (1./36.)*(3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
			simulation:collide_node(7,x,y,0, (1./36.)*(3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
			
			-- v8
			simulation:collide_node(8,x,y,0, r + (4./9.)*(-1.5*rusq));
			
			-- 2nd component
			rux = p*ux
			ruy = p*uy
	
			ruxx = rux*ux + mu
			ruyy = ruy*uy + mu
	
			ruxy = 2.*rux*uy
			rusq = ruxx + ruyy
			-- v0 - v3 
			simulation:collide_node(0,x,y,1, (1./9.)*( 3.*rux+4.5*ruxx-1.5*rusq))
			simulation:collide_node(1,x,y,1, (1./9.)*( 3.*ruy+4.5*ruyy-1.5*rusq))
			simulation:collide_node(2,x,y,1, (1./9.)*(-3.*rux+4.5*ruxx-1.5*rusq))
			simulation:collide_node(3,x,y,1, (1./9.)*(-3.*ruy+4.5*ruyy-1.5*rusq))

			-- v4 - v7 
			simulation:collide_node(4,x,y,1, (1./36.)*(3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
			simulation:collide_node(5,x,y,1, (1./36.)*(3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
			simulation:collide_node(6,x,y,1, (1./36.)*(3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
			simulation:collide_node(7,x,y,1, (1./36.)*(3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
			
			-- v8
			simulation:collide_node(8,x,y,1, p + (4./9.)*(-1.5*rusq));

		end
	end
end

simulation = lb.d2q9_BGK(2, NX, NY, 1, PY)
pinfo = simulation:partition_info()
--math.randomseed(pinfo:processor_rank() + 1)

simulation:set_parameters(PARAMETERS)

--[[if not PY then
	simulation:set_walls_speed(VT, VB)
end--]]

initialize(simulation, initializer)

t0 = lb.wtime()
last_report = t0

for t = START, END do
	--simulation:advance()
	simulation:stream()
	simulation:average()
	my_collide(simulation)
	
	if pinfo:processor_rank() == 0 then
		local now = lb.wtime()
		if now - last_report > 7 then
			report_progress(t0, START, END, t)
			last_report = now		
		end
	end

	if math.mod(t, FREQ) == 0 then
		render(density2rgb, simulation, t)
		--simulation:dump(make_filename(t,".h5"))
		--print(simulation:mass())
	end
end

simulation:destroy()
